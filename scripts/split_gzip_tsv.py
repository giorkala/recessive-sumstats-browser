#!/usr/bin/env python3
"""Split a gzip TSV into line-preserving gzip TSV parts.

Each part repeats the original header and is independently decompressible.
The splitter checks the compressed output size periodically and starts a new
part when the current file reaches the requested target size.
"""

from __future__ import annotations

import argparse
import gzip
import io
import sys
from pathlib import Path


def open_part(out_dir: Path, prefix: str, index: int, header: str):
    path = out_dir / f"{prefix}.part_{index:03d}.tsv.gz"
    raw = path.open("wb")
    gz = gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=6, mtime=0)
    text = io.TextIOWrapper(gz, encoding="utf-8", newline="")
    text.write(header)
    if not header.endswith("\n"):
        text.write("\n")
    return path, raw, gz, text


def close_part(text, gz, raw) -> None:
    text.flush()
    text.detach()
    gz.close()
    raw.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--target-mib", type=float, default=90.0)
    parser.add_argument("--check-every", type=int, default=10000)
    args = parser.parse_args()

    target_bytes = int(args.target_mib * 1024 * 1024)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    part_index = 1
    part_rows = 0
    total_rows = 0

    with gzip.open(args.input, "rt", encoding="utf-8", newline="") as source:
        header = source.readline()
        if not header:
            raise SystemExit(f"{args.input} is empty")
        path, raw, gz, text = open_part(args.out_dir, args.prefix, part_index, header)

        for line in source:
            text.write(line)
            part_rows += 1
            total_rows += 1

            if part_rows % args.check_every == 0:
                text.flush()
                gz.flush()
                raw.flush()
                if raw.tell() >= target_bytes:
                    close_part(text, gz, raw)
                    print(f"Wrote {path} ({part_rows:,} rows)", file=sys.stderr)
                    part_index += 1
                    part_rows = 0
                    path, raw, gz, text = open_part(args.out_dir, args.prefix, part_index, header)

        close_part(text, gz, raw)
        print(f"Wrote {path} ({part_rows:,} rows)", file=sys.stderr)

    print(f"Total rows split: {total_rows:,}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

