# BRaVa Recessive Summary Statistics Browser

Static GitHub Pages browser for BRaVa recessive association summary statistics.

The public site is served directly from the repository root:

```text
recessive-browser/
  index.html
  assets/
  data/
  images/
  README.md
```

The browser is plain HTML, CSS, and JavaScript. Summary statistics are stored as small gzip-compressed TSV chunks plus lightweight JSON metadata, so no backend server is required.

## Local Preview

Because the app uses `fetch`, preview it with a local server from the repository root:

```bash
python3 -m http.server 8000
```

Then open `http://localhost:8000`.

## Editable Site Text

Landing-page wording, the paper link, contact email, and footer note live in:

```text
data/site_content.json
```

Edit that JSON file when you want to adjust public-facing copy without touching the browser code.

## GitHub Pages

This repository is ready for GitHub Pages. After pushing to GitHub:

1. Open the repository settings.
2. Go to `Pages`.
3. Set `Build and deployment` to `GitHub Actions`.
4. Push to the `main` branch, or run the `Deploy GitHub Pages` workflow manually.

The included workflow publishes the repository root, so `index.html`, `assets/`, `data/`, and `images/` are served directly.

## Data Size Notes

GitHub blocks regular Git files larger than 100 MiB. The generated browser chunks are intentionally split well below that threshold; this checkout is also just under GitHub's ideal repository-size guidance of 1 GB.

## Inputs For Rebuilding Data

Expected source files, kept outside this GitHub-ready static site folder:

- `BRaVa_recessive.sumstats.combined.UKB_AOU_GEL_GNH_BioMe_BBJ.090126.txt.gz`
  cohort-level, pre-meta summary statistics.
- `BRaVa_recessive.220525.meta_cauchy.txt.gz`
  meta-analysis summary statistics.
- `Brava_Recessive_Supplementary_Tables.xlsx`
  phenotype metadata and gene symbol mapping.

## Build Browser Data

From this directory:

```bash
python3 scripts/build_browser_data.py \
  --pre-meta ../BRaVa_recessive.sumstats.combined.UKB_AOU_GEL_GNH_BioMe_BBJ.090126.txt.gz \
  --meta ../BRaVa_recessive.220525.meta_cauchy.txt.gz \
  --supplement ../Brava_Recessive_Supplementary_Tables.xlsx \
  --out data \
  --bucket-count 128 \
  --top-n 5000
```

For a quick test run:

```bash
python3 scripts/build_browser_data.py \
  --pre-meta ../BRaVa_recessive.sumstats.combined.UKB_AOU_GEL_GNH_BioMe_BBJ.090126.txt.gz \
  --meta ../BRaVa_recessive.220525.meta_cauchy.txt.gz \
  --supplement ../Brava_Recessive_Supplementary_Tables.xlsx \
  --out data \
  --max-rows 200000
```

## Split Raw Downloads

If you want line-preserving raw downloadable parts below GitHub's 100 MiB file limit:

```bash
python3 scripts/split_gzip_tsv.py \
  --input ../BRaVa_recessive.sumstats.combined.UKB_AOU_GEL_GNH_BioMe_BBJ.090126.txt.gz \
  --out-dir data/downloads/pre_meta \
  --prefix BRaVa_recessive.pre_meta \
  --target-mib 90
```

Each output part is a valid `.tsv.gz` file with the original header repeated.

## Citation And Feedback

Please cite the associated BRaVa recessive association study when using this browser. Feedback links in the site send email to the contact listed in `data/site_content.json`.
