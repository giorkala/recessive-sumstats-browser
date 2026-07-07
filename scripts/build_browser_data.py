#!/usr/bin/env python3
"""Build static browser data for BRaVa recessive summary statistics.

The script intentionally uses only the Python standard library so it can run on
machines without a scientific Python environment. It reads the supplementary
XLSX directly as XML, extracts phenotype metadata and gene-symbol mappings, then
streams large gzip TSV files once to produce QC summaries, top-hit indices, and
gene-bucketed gzip TSV chunks.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import heapq
import io
import json
import math
import os
import re
import shutil
import sys
import time
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
from xml.etree import ElementTree as ET


XLSX_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
XLSX_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS = {"main": XLSX_MAIN_NS, "rel": XLSX_REL_NS}

PRE_META_EXPECTED_COLUMNS = [
    "BIOBANK",
    "ANCESTRY",
    "SEX",
    "TRAIT",
    "ENCODING",
    "ANNOTATION",
    "ID",
    "BETA",
    "SE",
    "P",
    "N_CASE",
    "N_CTRL",
    "N_EFF",
    "AF",
    "AC",
    "N",
]

META_EXPECTED_COLUMNS = [
    "SEX",
    "TRAIT",
    "ID",
    "ANNOTATION",
    "ENCODING",
    "Z",
    "P",
    "number_of_pvals",
    "min_pvalue",
    "min_p_id",
    "CCT",
    "BIOBANK",
    "ANCESTRY",
]

COHORT_ALIASES = {
    "UK Biobank (UKB)": "UKB",
    "UK Biobank": "UKB",
    "UKB": "UKB",
    "All of Us (AOU)": "AOU",
    "All of Us": "AOU",
    "AOU": "AOU",
    "100k Genomes Project": "GEL",
    "100kGP": "GEL",
    "GEL": "GEL",
    "G&H": "GNH",
    "GNH": "GNH",
    "BBJ": "BBJ",
    "BioMe": "BioMe",
}

NULL_STRINGS = {"", "NA", "N/A", "NaN", "nan", "None", "null"}


def eprint(*parts: object) -> None:
    print(*parts, file=sys.stderr, flush=True)


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    cleaned = cleaned.strip("._")
    return cleaned or "unknown"


def parse_float(value: str) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text in NULL_STRINGS:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if math.isnan(parsed) or math.isinf(parsed):
        return None
    return parsed


def parse_int(value: str) -> Optional[int]:
    parsed = parse_float(value)
    if parsed is None:
        return None
    return int(parsed)


def json_number(value: str):
    parsed = parse_float(value)
    if parsed is None:
        return None
    if parsed.is_integer():
        return int(parsed)
    return parsed


def normalize_null(value: str):
    text = str(value).strip() if value is not None else ""
    return None if text in NULL_STRINGS else text


def fnv1a_32(text: str) -> int:
    value = 2166136261
    for byte in text.encode("utf-8"):
        value ^= byte
        value = (value * 16777619) & 0xFFFFFFFF
    return value


def bucket_for_id(gene_id: str, bucket_count: int) -> int:
    return fnv1a_32(gene_id) % bucket_count


def col_index(cell_ref: str) -> int:
    match = re.match(r"([A-Z]+)", cell_ref or "")
    if not match:
        return 0
    value = 0
    for char in match.group(1):
        value = value * 26 + ord(char) - 64
    return value - 1


class XlsxReader:
    def __init__(self, path: Path):
        self.path = path
        self.zip = zipfile.ZipFile(path)
        self.shared_strings = self._read_shared_strings()
        self.sheets = self._read_sheet_targets()

    def close(self) -> None:
        self.zip.close()

    def _read_shared_strings(self) -> List[str]:
        if "xl/sharedStrings.xml" not in self.zip.namelist():
            return []
        root = ET.fromstring(self.zip.read("xl/sharedStrings.xml"))
        return ["".join(si.itertext()) for si in root.findall("main:si", NS)]

    def _read_sheet_targets(self) -> Dict[str, str]:
        workbook = ET.fromstring(self.zip.read("xl/workbook.xml"))
        rels_root = ET.fromstring(self.zip.read("xl/_rels/workbook.xml.rels"))
        rels = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
        sheets = {}
        for sheet in workbook.findall("main:sheets/main:sheet", NS):
            rid = sheet.attrib[f"{{{XLSX_REL_NS}}}id"]
            target = rels[rid]
            if not target.startswith("/"):
                target = "xl/" + target
            else:
                target = target.lstrip("/")
            sheets[sheet.attrib["name"]] = target
        return sheets

    def _cell_value(self, cell: ET.Element) -> str:
        cell_type = cell.attrib.get("t")
        if cell_type == "inlineStr":
            return "".join(cell.itertext())
        value = cell.find("main:v", NS)
        if value is None:
            return ""
        text = value.text or ""
        if cell_type == "s":
            try:
                return self.shared_strings[int(text)]
            except (ValueError, IndexError):
                return text
        return text

    def rows(self, sheet_name: str) -> List[List[str]]:
        target = self.sheets[sheet_name]
        root = ET.fromstring(self.zip.read(target))
        rows = []
        for row in root.findall("main:sheetData/main:row", NS):
            values: List[str] = []
            for cell in row.findall("main:c", NS):
                idx = col_index(cell.attrib.get("r", ""))
                while len(values) < idx:
                    values.append("")
                values.append(self._cell_value(cell).strip())
            rows.append(values)
        width = max((len(row) for row in rows), default=0)
        for row in rows:
            row.extend([""] * (width - len(row)))
        return rows


def clean_cohort_name(value: str) -> Optional[str]:
    text = (value or "").strip()
    if not text:
        return None
    return COHORT_ALIASES.get(text, text)


def fill_cohort_headers(row: List[str]) -> List[Optional[str]]:
    filled: List[Optional[str]] = []
    current: Optional[str] = None
    for value in row:
        cohort = clean_cohort_name(value)
        if cohort:
            current = cohort
        filled.append(current)
    return filled


def extract_supplement_metadata(path: Path) -> Tuple[dict, dict]:
    reader = XlsxReader(path)
    try:
        phenotypes: Dict[str, dict] = {}

        if "ST1" in reader.sheets:
            rows = reader.rows("ST1")
            cohorts = fill_cohort_headers(rows[0]) if rows else []
            for row in rows[4:]:
                trait_id = row[1].strip() if len(row) > 1 else ""
                trait_name = row[0].strip() if row else ""
                if not trait_id or not trait_name:
                    continue
                record = {
                    "id": trait_id,
                    "label": trait_name,
                    "type": "binary",
                    "female_specific": normalize_null(row[2] if len(row) > 2 else ""),
                    "tested_lassen_2024": normalize_null(row[3] if len(row) > 3 else ""),
                    "total_cases": json_number(row[4] if len(row) > 4 else ""),
                    "total_controls": json_number(row[5] if len(row) > 5 else ""),
                    "case_control_ratio": json_number(row[6] if len(row) > 6 else ""),
                    "cohort_count": json_number(row[7] if len(row) > 7 else ""),
                    "cohort_ancestry": [],
                }
                ancestry_header = rows[3] if len(rows) > 3 else []
                for idx in range(8, len(row)):
                    ancestry = (ancestry_header[idx] if idx < len(ancestry_header) else "").strip()
                    cohort = cohorts[idx] if idx < len(cohorts) else None
                    value = parse_int(row[idx])
                    if cohort and ancestry and value is not None:
                        record["cohort_ancestry"].append(
                            {"cohort": cohort, "ancestry": ancestry, "n_cases": value}
                        )
                phenotypes[trait_id] = record

        if "ST2" in reader.sheets:
            rows = reader.rows("ST2")
            cohorts = fill_cohort_headers(rows[0]) if rows else []
            ancestry_header = rows[1] if len(rows) > 1 else []
            for row in rows[2:]:
                trait_id = row[1].strip() if len(row) > 1 else ""
                trait_name = row[0].strip() if row else ""
                if not trait_id or not trait_name:
                    continue
                record = {
                    "id": trait_id,
                    "label": trait_name,
                    "type": "quantitative",
                    "total_samples": json_number(row[2] if len(row) > 2 else ""),
                    "cohort_ancestry": [],
                }
                for idx in range(3, len(row)):
                    ancestry = (ancestry_header[idx] if idx < len(ancestry_header) else "").strip()
                    cohort = cohorts[idx] if idx < len(cohorts) else None
                    value = parse_int(row[idx])
                    if cohort and ancestry and value is not None:
                        record["cohort_ancestry"].append(
                            {"cohort": cohort, "ancestry": ancestry, "n_samples": value}
                        )
                phenotypes[trait_id] = record

        gene_map: Dict[str, dict] = {}
        if "ST4" in reader.sheets:
            rows = reader.rows("ST4")
            if rows:
                header = [cell.strip().lower() for cell in rows[0]]
                symbol_idx = header.index("gene") if "gene" in header else 0
                id_idx = None
                for candidate in ("ensembl_id", "ensmbl_id", "gene_id"):
                    if candidate in header:
                        id_idx = header.index(candidate)
                        break
                transcript_idx = header.index("transcript_id") if "transcript_id" in header else None
                if id_idx is not None:
                    for row in rows[1:]:
                        gene_id = row[id_idx].strip() if id_idx < len(row) else ""
                        if not gene_id:
                            continue
                        symbol = row[symbol_idx].strip() if symbol_idx < len(row) else ""
                        transcript = (
                            row[transcript_idx].strip()
                            if transcript_idx is not None and transcript_idx < len(row)
                            else ""
                        )
                        gene_map[gene_id] = {
                            "id": gene_id,
                            "symbol": symbol or None,
                            "transcript_id": transcript or None,
                            "source": "ST4",
                        }
        return phenotypes, gene_map
    finally:
        reader.close()


class ChunkWriters:
    def __init__(self, out_dir: Path, source_label: str, header: List[str], bucket_count: int):
        self.out_dir = out_dir
        self.source_label = source_label
        self.header = header
        self.bucket_count = bucket_count
        self.handles: Dict[int, io.TextIOWrapper] = {}
        self.paths: Dict[int, Path] = {}
        self.rows = Counter()
        self.out_dir.mkdir(parents=True, exist_ok=True)

    def write(self, bucket: int, row: List[str], symbol: Optional[str]) -> None:
        if bucket not in self.handles:
            path = self.out_dir / f"{bucket:03d}.tsv.gz"
            handle = gzip.open(path, "wt", encoding="utf-8", compresslevel=6, newline="")
            handle.write("\t".join(self.header) + "\n")
            self.handles[bucket] = handle
            self.paths[bucket] = path
        self.handles[bucket].write("\t".join(row) + "\n")
        self.rows[bucket] += 1

    def close(self) -> List[dict]:
        for handle in self.handles.values():
            handle.close()
        files = []
        for bucket, path in sorted(self.paths.items()):
            files.append(
                {
                    "bucket": bucket,
                    "path": str(path.name),
                    "rows": self.rows[bucket],
                    "bytes": path.stat().st_size,
                }
            )
        return files


def make_empty_qc(source_label: str, path: Path, header: List[str]) -> dict:
    return {
        "source": source_label,
        "input": path.name,
        "header": header,
        "rows": 0,
        "categories": defaultdict(Counter),
        "numeric": defaultdict(lambda: {"missing": 0, "invalid": 0, "min": None, "max": None}),
        "pvalue": {"missing": 0, "invalid": 0, "outside_0_1": 0, "min": None, "max": None},
        "se": {"missing": 0, "invalid": 0, "nonpositive": 0, "min": None, "max": None},
        "af": {"missing": 0, "invalid": 0, "outside_0_1": 0, "min": None, "max": None},
        "ac_n": {"checked": 0, "ac_gt_2n": 0, "negative": 0},
        "missing_id": 0,
        "missing_symbol": 0,
        "malformed_rows": 0,
        "adjacent_duplicate_keys": 0,
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "finished_at": None,
    }


def update_min_max(stats: dict, value: float) -> None:
    stats["min"] = value if stats["min"] is None else min(stats["min"], value)
    stats["max"] = value if stats["max"] is None else max(stats["max"], value)


def update_numeric(qc: dict, field: str, value: str) -> Optional[float]:
    stats = qc["numeric"][field]
    parsed = parse_float(value)
    if value.strip() in NULL_STRINGS:
        stats["missing"] += 1
        return None
    if parsed is None:
        stats["invalid"] += 1
        return None
    update_min_max(stats, parsed)
    return parsed


def compact_categories(categories: defaultdict) -> dict:
    compact = {}
    for key, counter in categories.items():
        compact[key] = dict(sorted(counter.items(), key=lambda item: (-item[1], item[0])))
    return compact


def compact_numeric(numeric: defaultdict) -> dict:
    return {key: value for key, value in numeric.items()}


def maybe_push_top(
    heap: List[Tuple[float, int, dict]],
    top_n: int,
    serial: int,
    p_value: Optional[float],
    row_dict: dict,
) -> None:
    if top_n <= 0 or p_value is None:
        return
    item = (-p_value, serial, row_dict)
    heapq.heappush(heap, item)
    if len(heap) > top_n:
        heapq.heappop(heap)


def process_summary_file(
    *,
    source_label: str,
    path: Path,
    expected_header: List[str],
    out_dir: Path,
    gene_map: Dict[str, dict],
    gene_ids_seen: set,
    bucket_count: int,
    top_n: int,
    max_rows: Optional[int],
    write_chunks: bool,
) -> Tuple[dict, List[dict], List[dict]]:
    source_dir = out_dir / "chunks" / source_label / "gene-buckets"
    top_heap: List[Tuple[float, int, dict]] = []
    chunk_files: List[dict] = []
    writers: Optional[ChunkWriters] = None
    serial = 0
    last_key = None

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"{path} is empty")
        header = header_line.rstrip("\n\r").split("\t")
        qc = make_empty_qc(source_label, path, header)
        qc["expected_header"] = expected_header
        qc["header_matches_expected"] = header == expected_header

        index = {name: idx for idx, name in enumerate(header)}
        if "ID" not in index or "P" not in index:
            raise ValueError(f"{path} must contain ID and P columns")

        category_fields = [
            field
            for field in ("BIOBANK", "ANCESTRY", "SEX", "TRAIT", "ENCODING", "ANNOTATION", "ID")
            if field in index
        ]
        numeric_fields = [
            field
            for field in (
                "BETA",
                "SE",
                "Z",
                "P",
                "N_CASE",
                "N_CTRL",
                "N_EFF",
                "AF",
                "AC",
                "N",
                "number_of_pvals",
                "min_pvalue",
                "CCT",
            )
            if field in index
        ]

        if write_chunks:
            writers = ChunkWriters(source_dir, source_label, header, bucket_count)

        for line_no, line in enumerate(handle, start=2):
            if max_rows is not None and qc["rows"] >= max_rows:
                break
            row = line.rstrip("\n\r").split("\t")
            if len(row) != len(header):
                qc["malformed_rows"] += 1
                continue

            qc["rows"] += 1
            serial += 1
            row_dict = {name: row[idx] for name, idx in index.items()}
            gene_id = row[index["ID"]].strip()
            if not gene_id:
                qc["missing_id"] += 1
                continue
            gene_ids_seen.add(gene_id)
            gene_record = gene_map.get(gene_id)
            symbol = gene_record.get("symbol") if gene_record else None
            if not symbol:
                qc["missing_symbol"] += 1

            key_fields = [
                row[index[field]]
                for field in ("BIOBANK", "ANCESTRY", "SEX", "TRAIT", "ENCODING", "ANNOTATION", "ID")
                if field in index
            ]
            key = tuple(key_fields)
            if key == last_key:
                qc["adjacent_duplicate_keys"] += 1
            last_key = key

            for field in category_fields:
                value = row[index[field]].strip()
                if value:
                    qc["categories"][field][value] += 1

            parsed_values = {}
            for field in numeric_fields:
                parsed_values[field] = update_numeric(qc, field, row[index[field]])

            p_value = parsed_values.get("P")
            if row[index["P"]].strip() in NULL_STRINGS:
                qc["pvalue"]["missing"] += 1
            elif p_value is None:
                qc["pvalue"]["invalid"] += 1
            else:
                update_min_max(qc["pvalue"], p_value)
                if p_value < 0 or p_value > 1:
                    qc["pvalue"]["outside_0_1"] += 1

            if "SE" in index:
                se = parsed_values.get("SE")
                if row[index["SE"]].strip() in NULL_STRINGS:
                    qc["se"]["missing"] += 1
                elif se is None:
                    qc["se"]["invalid"] += 1
                else:
                    update_min_max(qc["se"], se)
                    if se <= 0:
                        qc["se"]["nonpositive"] += 1

            if "AF" in index:
                af = parsed_values.get("AF")
                if row[index["AF"]].strip() in NULL_STRINGS:
                    qc["af"]["missing"] += 1
                elif af is None:
                    qc["af"]["invalid"] += 1
                else:
                    update_min_max(qc["af"], af)
                    if af < 0 or af > 1:
                        qc["af"]["outside_0_1"] += 1

            if "AC" in index and "N" in index:
                ac = parsed_values.get("AC")
                n = parsed_values.get("N")
                if ac is not None and n is not None:
                    qc["ac_n"]["checked"] += 1
                    if ac < 0 or n < 0:
                        qc["ac_n"]["negative"] += 1
                    if ac > 2 * n:
                        qc["ac_n"]["ac_gt_2n"] += 1

            bucket = bucket_for_id(gene_id, bucket_count)
            if writers:
                writers.write(bucket, row, symbol)

            top_record = {
                "source": source_label,
                "gene_id": gene_id,
                "symbol": symbol,
                "trait": row[index["TRAIT"]] if "TRAIT" in index else None,
                "annotation": row[index["ANNOTATION"]] if "ANNOTATION" in index else None,
                "encoding": row[index["ENCODING"]] if "ENCODING" in index else None,
                "sex": row[index["SEX"]] if "SEX" in index else None,
                "p": p_value,
                "biobank": row[index["BIOBANK"]] if "BIOBANK" in index else None,
                "ancestry": row[index["ANCESTRY"]] if "ANCESTRY" in index else None,
                "beta": parsed_values.get("BETA"),
                "se": parsed_values.get("SE"),
                "z": parsed_values.get("Z"),
                "n": parsed_values.get("N"),
                "n_eff": parsed_values.get("N_EFF"),
                "n_case": parsed_values.get("N_CASE"),
                "n_ctrl": parsed_values.get("N_CTRL"),
                "number_of_pvals": parsed_values.get("number_of_pvals"),
                "min_pvalue": parsed_values.get("min_pvalue"),
                "min_p_id": row[index["min_p_id"]] if "min_p_id" in index else None,
                "bucket": bucket,
            }
            maybe_push_top(top_heap, top_n, serial, p_value, top_record)

            if qc["rows"] % 1_000_000 == 0:
                eprint(f"{source_label}: processed {qc['rows']:,} rows")

    if writers:
        chunk_files = writers.close()

    qc["categories"] = compact_categories(qc["categories"])
    qc["numeric"] = compact_numeric(qc["numeric"])
    qc["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())

    top_hits = [item[2] for item in sorted(top_heap, key=lambda item: item[2]["p"] or 1)]
    return qc, top_hits, chunk_files


def write_json(path: Path, data) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_qc_markdown(path: Path, qc_records: List[dict]) -> None:
    lines = ["# QC Summary", ""]
    for qc in qc_records:
        lines.extend(
            [
                f"## {qc['source']}",
                "",
                f"- Rows: {qc['rows']:,}",
                f"- Header matches expected: {qc.get('header_matches_expected')}",
                f"- Malformed rows: {qc['malformed_rows']:,}",
                f"- Missing gene IDs: {qc['missing_id']:,}",
                f"- Rows without ST4 gene symbol: {qc['missing_symbol']:,}",
                f"- Adjacent duplicate keys: {qc['adjacent_duplicate_keys']:,}",
                f"- P-value min/max: {qc['pvalue']['min']} / {qc['pvalue']['max']}",
                f"- P-values outside [0, 1]: {qc['pvalue']['outside_0_1']:,}",
            ]
        )
        if qc.get("se"):
            lines.append(f"- Non-positive SE values: {qc['se']['nonpositive']:,}")
        if qc.get("af"):
            lines.append(f"- AF values outside [0, 1]: {qc['af']['outside_0_1']:,}")
        if qc.get("ac_n"):
            lines.append(f"- AC > 2N checks: {qc['ac_n']['ac_gt_2n']:,} / {qc['ac_n']['checked']:,}")
        lines.extend(["", "### Categories", ""])
        for field, values in qc["categories"].items():
            preview = ", ".join(f"{key} ({value:,})" for key, value in list(values.items())[:12])
            lines.append(f"- {field}: {preview}")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def rel_chunk_files(files: List[dict], prefix: str) -> List[dict]:
    return [
        {
            **item,
            "url": f"data/chunks/{prefix}/gene-buckets/{item['path']}",
        }
        for item in files
    ]


def build_gene_records(
    gene_map: Dict[str, dict],
    gene_ids_seen: set,
    bucket_count: int,
) -> List[dict]:
    records = []
    for gene_id in sorted(gene_ids_seen):
        mapped = gene_map.get(gene_id, {})
        records.append(
            {
                "id": gene_id,
                "symbol": mapped.get("symbol"),
                "transcript_id": mapped.get("transcript_id"),
                "bucket": bucket_for_id(gene_id, bucket_count),
                "mapping_source": mapped.get("source"),
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pre-meta", type=Path, help="Pre-meta cohort-level summary statistics .tsv.gz")
    parser.add_argument("--meta", type=Path, help="Meta-analysis summary statistics .tsv.gz")
    parser.add_argument("--supplement", type=Path, required=True, help="Supplementary tables .xlsx")
    parser.add_argument("--out", type=Path, required=True, help="Output directory, usually data")
    parser.add_argument("--bucket-count", type=int, default=128)
    parser.add_argument("--top-n", type=int, default=5000)
    parser.add_argument("--max-rows", type=int, help="Process at most this many rows per summary file")
    parser.add_argument(
        "--metadata-only",
        action="store_true",
        help="Build JSON metadata and QC only; do not write gene-bucket chunks",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove existing generated output directory before writing",
    )
    args = parser.parse_args()

    if args.bucket_count <= 0:
        raise SystemExit("--bucket-count must be positive")

    if args.clean and args.out.exists():
        shutil.rmtree(args.out)
    args.out.mkdir(parents=True, exist_ok=True)

    phenotypes, gene_map = extract_supplement_metadata(args.supplement)
    eprint(f"Loaded {len(phenotypes):,} phenotype records and {len(gene_map):,} ST4 gene mappings")

    qc_records = []
    chunk_manifest = {}
    top_by_source = {}
    gene_ids_seen = set()

    if args.pre_meta:
        qc, top_hits, chunks = process_summary_file(
            source_label="pre_meta",
            path=args.pre_meta,
            expected_header=PRE_META_EXPECTED_COLUMNS,
            out_dir=args.out,
            gene_map=gene_map,
            gene_ids_seen=gene_ids_seen,
            bucket_count=args.bucket_count,
            top_n=args.top_n,
            max_rows=args.max_rows,
            write_chunks=not args.metadata_only,
        )
        qc_records.append(qc)
        top_by_source["pre_meta"] = top_hits
        chunk_manifest["pre_meta"] = rel_chunk_files(chunks, "pre_meta")

    if args.meta:
        qc, top_hits, chunks = process_summary_file(
            source_label="meta",
            path=args.meta,
            expected_header=META_EXPECTED_COLUMNS,
            out_dir=args.out,
            gene_map=gene_map,
            gene_ids_seen=gene_ids_seen,
            bucket_count=args.bucket_count,
            top_n=args.top_n,
            max_rows=args.max_rows,
            write_chunks=not args.metadata_only,
        )
        qc_records.append(qc)
        top_by_source["meta"] = top_hits
        chunk_manifest["meta"] = rel_chunk_files(chunks, "meta")

    gene_records = build_gene_records(gene_map, gene_ids_seen, args.bucket_count)
    symbol_index: Dict[str, List[str]] = defaultdict(list)
    for record in gene_records:
        if record["symbol"]:
            symbol_index[record["symbol"].upper()].append(record["id"])

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "bucket_count": args.bucket_count,
        "max_rows": args.max_rows,
        "sources": {
            qc["source"]: {
                "input": qc["input"],
                "rows": qc["rows"],
                "header": qc["header"],
                "header_matches_expected": qc["header_matches_expected"],
                "categories": qc["categories"],
            }
            for qc in qc_records
        },
        "chunks": chunk_manifest,
        "counts": {
            "phenotypes": len(phenotypes),
            "genes": len(gene_records),
            "genes_with_symbols": sum(1 for record in gene_records if record["symbol"]),
            "symbols": len(symbol_index),
        },
        "files": {
            "phenotypes": "data/phenotypes.json",
            "genes": "data/genes.json",
            "symbol_index": "data/symbol_index.json",
            "top_hits": "data/top_hits.json",
            "qc_summary": "data/qc/summary.md",
        },
    }

    preferred_top = top_by_source.get("meta") or top_by_source.get("pre_meta") or []

    write_json(args.out / "phenotypes.json", sorted(phenotypes.values(), key=lambda item: item["id"]))
    write_json(args.out / "genes.json", gene_records)
    write_json(args.out / "symbol_index.json", dict(sorted(symbol_index.items())))
    write_json(args.out / "top_hits.json", preferred_top)
    for source, hits in top_by_source.items():
        write_json(args.out / f"top_hits.{source}.json", hits)
    write_json(args.out / "manifest.json", manifest)
    write_json(args.out / "qc" / "summary.json", qc_records)
    write_qc_markdown(args.out / "qc" / "summary.md", qc_records)

    eprint(f"Wrote browser data to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

