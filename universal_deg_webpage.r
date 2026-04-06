

make_deg_webpage <- function(input_file,
                             out_dir,
                             lfc_cut = 2,
                             padj_cut = 0.05,
                             p_cap = 50) {
  
  ## =========================
  ## 1. Packages
  ## =========================
  pkgs <- c("readxl", "plotly", "htmlwidgets", "stringr")
  need <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(need)) {
    install.packages(need, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
  
  ## =========================
  ## 2. Output folder
  ## =========================
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## =========================
  ## 3. Color palette
  ## =========================
  color_palette <- c(
    "Up"   = "red",
    "Down" = "blue",
    "NS"   = "#D3D3D3"
  )
  
  ## =========================
  ## 4. Helper functions
  ## =========================
  safe_numeric <- function(x) {
    x <- as.character(x)
    x <- gsub(",", "", x)
    x <- gsub(" ", "", x)
    suppressWarnings(as.numeric(x))
  }
  
  clean_text <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    x[x %in% c("", "NA", "N/A", "na", "null", "NULL")] <- NA
    x
  }
  
  find_first_existing_col <- function(df, candidates, required = TRUE) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) {
      if (required) {
        stop(
          paste0(
            "Required column not found. Acceptable columns are: ",
            paste(candidates, collapse = ", ")
          ),
          call. = FALSE
        )
      } else {
        return(NA_character_)
      }
    }
    hit[1]
  }
  
  sanitize_filename <- function(x) {
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
    x <- stringr::str_replace_all(x, "_+", "_")
    x <- stringr::str_replace_all(x, "^_|_$", "")
    if (nchar(x) == 0) x <- "volcano_plot"
    x
  }
  
  ## =========================
  ## 5. Read input
  ## =========================
  ext <- tolower(tools::file_ext(input_file))
  
  if (ext == "xlsx") {
    dat <- readxl::read_xlsx(input_file)
  } else if (ext == "csv") {
    dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    dat <- read.delim(input_file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported input format. Please use .xlsx, .csv, .tsv, or .txt", call. = FALSE)
  }
  
  dat <- as.data.frame(dat, stringsAsFactors = FALSE)
  
  ## =========================
  ## 6. Detect columns
  ## =========================
  id_col <- find_first_existing_col(
    dat,
    candidates = c("trinity_id", "gene_id", "transcript_id", "feature_id"),
    required = TRUE
  )
  
  log2fc_col <- find_first_existing_col(
    dat,
    candidates = c("log2FoldChange"),
    required = TRUE
  )
  
  padj_col <- find_first_existing_col(
    dat,
    candidates = c("padj"),
    required = TRUE
  )
  
  gene_name_col <- find_first_existing_col(
    dat,
    candidates = c("gene_name"),
    required = FALSE
  )
  
  pvalue_col <- find_first_existing_col(
    dat,
    candidates = c("pvalue"),
    required = FALSE
  )
  
  baseMean_col <- find_first_existing_col(
    dat,
    candidates = c("baseMean"),
    required = FALSE
  )
  
  sampleA_col <- find_first_existing_col(
    dat,
    candidates = c("sampleA"),
    required = FALSE
  )
  
  sampleB_col <- find_first_existing_col(
    dat,
    candidates = c("sampleB"),
    required = FALSE
  )
  
  baseMeanA_col <- find_first_existing_col(
    dat,
    candidates = c("baseMeanA"),
    required = FALSE
  )
  
  baseMeanB_col <- find_first_existing_col(
    dat,
    candidates = c("baseMeanB"),
    required = FALSE
  )
  
  ## =========================
  ## 7. Standardize columns
  ## =========================
  dat$feature_id     <- clean_text(dat[[id_col]])
  dat$log2FoldChange <- safe_numeric(dat[[log2fc_col]])
  dat$padj           <- safe_numeric(dat[[padj_col]])
  
  dat$gene_name <- if (!is.na(gene_name_col)) clean_text(dat[[gene_name_col]]) else NA_character_
  dat$pvalue    <- if (!is.na(pvalue_col))    safe_numeric(dat[[pvalue_col]]) else NA_real_
  dat$baseMean  <- if (!is.na(baseMean_col))  safe_numeric(dat[[baseMean_col]]) else NA_real_
  dat$sampleA   <- if (!is.na(sampleA_col))   clean_text(dat[[sampleA_col]]) else NA_character_
  dat$sampleB   <- if (!is.na(sampleB_col))   clean_text(dat[[sampleB_col]]) else NA_character_
  dat$baseMeanA <- if (!is.na(baseMeanA_col)) safe_numeric(dat[[baseMeanA_col]]) else NA_real_
  dat$baseMeanB <- if (!is.na(baseMeanB_col)) safe_numeric(dat[[baseMeanB_col]]) else NA_real_
  
  ## gene display logic: gene_name first, then ID
  dat$display_gene <- ifelse(
    !is.na(dat$gene_name) & dat$gene_name != "",
    dat$gene_name,
    dat$feature_id
  )
  
  ## filter bad rows
  dat <- dat[!is.na(dat$feature_id) & dat$feature_id != "", ]
  dat <- dat[!is.na(dat$log2FoldChange) & !is.na(dat$padj), ]
  
  if (nrow(dat) == 0) {
    stop("No valid rows remain after filtering.", call. = FALSE)
  }
  
  ## =========================
  ## 8. Auto comparison name
  ## =========================
  comparison_name <- NA_character_
  
  sA <- unique(na.omit(dat$sampleA))
  sB <- unique(na.omit(dat$sampleB))
  
  if (length(sA) >= 1 && length(sB) >= 1) {
    comparison_name <- paste0(sA[1], " vs ", sB[1], " Comparison")
  }
  
  if (is.na(comparison_name) || comparison_name == "") {
    comparison_name <- tools::file_path_sans_ext(basename(input_file))
  }
  
  ## =========================
  ## 9. DEG classification
  ## =========================
  dat$status <- "NS"
  dat$status[dat$log2FoldChange >  lfc_cut & dat$padj < padj_cut] <- "Up"
  dat$status[dat$log2FoldChange < -lfc_cut & dat$padj < padj_cut] <- "Down"
  
  n_up   <- sum(dat$status == "Up", na.rm = TRUE)
  n_down <- sum(dat$status == "Down", na.rm = TRUE)
  n_ns   <- sum(dat$status == "NS", na.rm = TRUE)
  
  ## =========================
  ## 10. Y-axis processing
  ## =========================
  dat$padj_adj  <- ifelse(dat$padj <= 0, 1e-300, dat$padj)
  dat$neglog10p <- -log10(dat$padj_adj)
  dat$y_plot    <- ifelse(dat$neglog10p > p_cap, p_cap, dat$neglog10p)
  
  ## =========================
  ## 11. Hover text
  ## =========================
  gene_text <- ifelse(!is.na(dat$gene_name) & dat$gene_name != "", dat$gene_name, "NA")
  padj_text <- ifelse(!is.na(dat$padj), format(dat$padj, scientific = TRUE, digits = 3), "NA")
  pval_text <- ifelse(!is.na(dat$pvalue), format(dat$pvalue, scientific = TRUE, digits = 3), "NA")
  bm_text   <- ifelse(!is.na(dat$baseMean), round(dat$baseMean, 3), "NA")
  bmA_text  <- ifelse(!is.na(dat$baseMeanA), round(dat$baseMeanA, 3), "NA")
  bmB_text  <- ifelse(!is.na(dat$baseMeanB), round(dat$baseMeanB, 3), "NA")
  
  sampleA_name <- if (length(sA) >= 1) sA[1] else "SampleA"
  sampleB_name <- if (length(sB) >= 1) sB[1] else "SampleB"
  
  dat$hover <- paste0(
    "<b>Gene name:</b> ", gene_text, "<br>",
    "<b>Feature ID:</b> ", dat$feature_id, "<br>",
    "<b>Log2FC:</b> ", round(dat$log2FoldChange, 3), "<br>",
    "<b>Padj:</b> ", padj_text, "<br>",
    "<b>P-value:</b> ", pval_text, "<br>",
    "<b>baseMean:</b> ", bm_text, "<br>",
    "<b>", sampleA_name, " mean:</b> ", bmA_text, "<br>",
    "<b>", sampleB_name, " mean:</b> ", bmB_text
  )
  
  ## =========================
  ## 12. Title
  ## =========================
  title_html <- list(
    text = paste0(
      "<b style='font-size:30px; font-family:Arial;'>", comparison_name, "</b><br>",
      "<span style='font-size:18px; color:#555;'>",
      "Up: <b style='color:red;'>", n_up, "</b> | ",
      "Down: <b style='color:blue;'>", n_down, "</b> | ",
      "NS: ", n_ns, "</span>"
    ),
    x = 0.5,
    xanchor = "center",
    y = 0.96,
    yanchor = "top"
  )
  
  ## =========================
  ## 13. Volcano plot
  ## =========================
  p <- plotly::plot_ly(
    data = dat,
    x = ~log2FoldChange,
    y = ~y_plot,
    color = ~status,
    colors = color_palette,
    text = ~hover,
    hoverinfo = "text",
    type = "scatter",
    mode = "markers",
    marker = list(size = 7, opacity = 0.7)
  ) %>%
    plotly::layout(
      title = title_html,
      margin = list(t = 120, b = 60, l = 80, r = 80),
      xaxis = list(
        title = "<b>Log2 Fold Change</b>",
        showgrid = FALSE,
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 2,
        ticks = "outside"
      ),
      yaxis = list(
        title = "<b>-Log10 (Adjusted P-value)</b>",
        showgrid = FALSE,
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "black",
        linewidth = 2,
        ticks = "outside"
      ),
      shapes = list(
        list(
          type = "line",
          x0 = -lfc_cut, x1 = -lfc_cut,
          y0 = 0, y1 = 1,
          xref = "x", yref = "paper",
          line = list(color = "#CCC", dash = "dot", width = 1)
        ),
        list(
          type = "line",
          x0 = lfc_cut, x1 = lfc_cut,
          y0 = 0, y1 = 1,
          xref = "x", yref = "paper",
          line = list(color = "#CCC", dash = "dot", width = 1)
        ),
        list(
          type = "line",
          x0 = 0, x1 = 1,
          y0 = -log10(padj_cut), y1 = -log10(padj_cut),
          xref = "paper", yref = "y",
          line = list(color = "#CCC", dash = "dot", width = 1)
        )
      ),
      plot_bgcolor = "white",
      paper_bgcolor = "white",
      showlegend = FALSE
    ) %>%
    plotly::toWebGL()
  
  ## =========================
  ## 14. Save volcano html
  ## =========================
  volcano_file <- paste0(sanitize_filename(comparison_name), ".html")
  volcano_path <- file.path(out_dir, volcano_file)
  htmlwidgets::saveWidget(p, volcano_path, selfcontained = TRUE)
  
  ## =========================
  ## 15. Build index page
  ## =========================
  index_html <- paste0(
    "<!DOCTYPE html>
<html>
<head>
  <meta charset='UTF-8'>
  <meta name='viewport' content='width=device-width, initial-scale=1.0'>
  <title>RNA-Seq DEG Visualization Portal</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      background: #fdfdfd;
      margin: 0;
      padding: 60px 20px;
      display: flex;
      justify-content: center;
    }
    .card {
      max-width: 900px;
      width: 100%;
      background: white;
      padding: 40px;
      border-radius: 18px;
      box-shadow: 0 10px 30px rgba(0,0,0,0.05);
      border: 1px solid #eee;
    }
    h1 {
      color: #222;
      margin-top: 0;
      margin-bottom: 12px;
      font-size: 32px;
    }
    p {
      color: #555;
      line-height: 1.6;
      font-size: 16px;
    }
    .summary {
      margin: 25px 0 30px 0;
      padding: 20px;
      background: #fafafa;
      border: 1px solid #eee;
      border-radius: 12px;
    }
    .btn {
      display: inline-block;
      padding: 15px 28px;
      background: #333;
      color: white;
      text-decoration: none;
      border-radius: 50px;
      transition: 0.3s;
      font-weight: bold;
    }
    .btn:hover {
      background: red;
      transform: scale(1.03);
    }
    .meta {
      margin-top: 30px;
      font-size: 14px;
      color: #777;
    }
  </style>
</head>
<body>
  <div class='card'>
    <h1>RNA-Seq DEG Visualization Portal</h1>
    <p>
      This report provides an interactive volcano plot for differential expression analysis.
      Each point includes both a unique feature ID and gene name when available.
    </p>

    <div class='summary'>
      <p><b>Comparison:</b> ", comparison_name, "</p>
      <p><b>Up-regulated:</b> <span style='color:red;'>", n_up, "</span></p>
      <p><b>Down-regulated:</b> <span style='color:blue;'>", n_down, "</span></p>
      <p><b>Non-significant:</b> ", n_ns, "</p>
      <p><b>Thresholds:</b> |log2FC| > ", lfc_cut, " and adjusted P-value < ", padj_cut, "</p>
    </div>

    <a class='btn' href='", volcano_file, "'>Open interactive volcano plot</a>

    <div class='meta'>
      Generated automatically by make_deg_webpage().
    </div>
  </div>
</body>
</html>"
  )
  
  writeLines(index_html, file.path(out_dir, "index.html"))
  
  ## =========================
  ## 16. Summary
  ## =========================
  cat("\n--- All tasks completed successfully! ---\n")
  cat("Input file: ", input_file, "\n", sep = "")
  cat("Output directory: ", out_dir, "\n", sep = "")
  cat("Comparison: ", comparison_name, "\n", sep = "")
  cat("Volcano HTML: ", volcano_path, "\n", sep = "")
  cat("Index HTML: ", file.path(out_dir, "index.html"), "\n", sep = "")
  cat("Detected unique ID column: ", id_col, "\n", sep = "")
  cat("Detected gene_name column: ", ifelse(is.na(gene_name_col), "None", gene_name_col), "\n", sep = "")
  cat("Up: ", n_up, " | Down: ", n_down, " | NS: ", n_ns, "\n", sep = "")
  
  ## return useful objects
  invisible(list(
    data = dat,
    comparison_name = comparison_name,
    volcano_file = volcano_path,
    index_file = file.path(out_dir, "index.html"),
    detected_id_col = id_col,
    detected_gene_name_col = ifelse(is.na(gene_name_col), NA, gene_name_col),
    counts = c(Up = n_up, Down = n_down, NS = n_ns)
  ))
}