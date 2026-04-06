# Interactive Volcano Plot Tool

A lightweight R-based tool for generating interactive volcano plot webpages from differential expression (DEG) result tables.

## Features

- Supports flexible unique ID columns:
  - `trinity_id`
  - `gene_id`
  - `transcript_id`
  - `feature_id`
- Supports `gene_name` when available
- Generates interactive volcano plot HTML pages
- Generates an index page for navigation
- Supports webpage sharing through GitHub Pages

## Input format

### Required columns
- one unique ID column:
  - `trinity_id`, `gene_id`, `transcript_id`, or `feature_id`
- `log2FoldChange`
- `padj`

### Recommended columns
- `gene_name`
- `pvalue`
- `baseMean`
- `sampleA`
- `sampleB`
- `baseMeanA`
- `baseMeanB`

## Example input
`data/Example data.xlsx`

## How to run

```r
source("code/universal_deg_webpage.r")

make_deg_webpage(
  input_file = "data/Example data.xlsx",
  out_dir    = "docs"
)
```

## Output

The tool generates:
- `index.html`
- comparison-specific interactive volcano plot HTML
- static volcano plot figures

## Online demo

Add your GitHub Pages link here.

## License

MIT License
