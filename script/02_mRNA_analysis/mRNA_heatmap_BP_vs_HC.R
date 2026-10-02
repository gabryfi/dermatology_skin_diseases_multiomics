#!/usr/bin/env Rscript

# =============================================================================
# BAMBU GENE COUNTS -> DESeq2 -> mRNA HEATMAP (3 FIGURE FILES)
# =============================================================================
#
# Purpose
#   Reproduce the BP-versus-HC heatmap starting from bambu
#   counts_gene.txt. The script creates exactly three figure files in a new
#   output directory: PNG, TIFF and SVG. It writes no tables or log files.
#   The internal groups PB/SANO are displayed as BP/HC in the figure.
#
# Identifier policy
#   Ensembl gene IDs remain the primary, unique identifiers throughout import,
#   DESeq2, filtering and ranking. Gene symbols are obtained
#   from the same GTF used by bambu and are used only as graphical row labels.
#   Unmapped IDs remain Ensembl IDs. Duplicate symbols are displayed as
#   "SYMBOL [ENSEMBL]" so that no identifiers are collapsed.
#
# Statistical/graphical choices
#   - bambu count estimates are rounded before DESeq2;
#   - DESeq2 model: ~ condition (no batch term);
#   - Wald test, PB versus SANO by default;
#   - DEG definition: BH FDR < 0.05 and |log2FC| >= log2(1.5);
#   - top 50 ranked by FDR, nominal p value, |log2FC|, Ensembl ID;
#   - complete DESeq2 VST with blind = FALSE for visualization only;
#   - pheatmap scale = "row" (per-gene z scores);
#   - rows: Pearson correlation; columns: Euclidean distance;
#     complete linkage for both, matching the miRNA heatmaps;
#   - the color scale is centered at zero without clipping;
#   - internal treatment/control labels PB/SANO remain unchanged for DESeq2.
#
# Primary references
#   Love MI, Huber W, Anders S. Genome Biol. 2014;15:550.
#   doi:10.1186/s13059-014-0550-8
#   Chen Y et al. Nat Methods. 2023;20:1187-1195.
#   doi:10.1038/s41592-023-01908-w
#   Official bambu workflow:
#   https://bioconductor.org/packages/bambu
#
# Command for the PB/SANO counts file used in this project
#   Rscript bambu_DESeq2_heatmap_BP_vs_HC_publication_vst.R \
#     --counts counts_gene.txt \
#     --gtf BULK_HUMAN_G38.gtf \
#     --outdir BPvsHC_mRNA_heatmap_publication_vst
#
# The current bambu column names contain "healthy" and "pem"; therefore the
# script can infer SANO/PB automatically. For any other naming convention,
# provide --metadata with these tab-separated columns:
#   count_column    sample_id    condition
#
# Run with --help for all options.
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

SCRIPT_VERSION <- "2.0.0"

usage <- function() {
  cat(
    paste0(
      "bambu_DESeq2_heatmap_BP_vs_HC_publication_vst.R v", SCRIPT_VERSION, "\n\n",
      "Required:\n",
      "  --counts FILE       bambu counts_gene.txt (GENEID + sample columns)\n\n",
      "  --gtf FILE          exact GTF used by bambu; symbols are used only\n",
      "                      as heatmap labels\n\n",
      "Optional:\n",
      "  --metadata FILE     TSV: count_column, sample_id, condition\n",
      "  --outdir DIR        output directory [BPvsHC_mRNA_heatmap_publication_vst]\n",
      "  --reference LABEL   reference condition [SANO]\n",
      "  --treatment LABEL   treatment condition [PB]\n",
      "  --fdr NUMBER        BH FDR cutoff [0.05]\n",
      "  --fc NUMBER         absolute linear fold-change cutoff [1.5]\n",
      "  --top INTEGER       maximum number of genes displayed [50]\n",
      "  --dpi INTEGER       PNG and TIFF resolution [600]\n",
      "  --help              print this help\n\n",
      "Example:\n",
      "  Rscript bambu_DESeq2_heatmap_BP_vs_HC_publication_vst.R \\\n",
      "    --counts counts_gene.txt \\\n",
      "    --gtf BULK_HUMAN_G38.gtf \\\n",
      "    --outdir BPvsHC_mRNA_heatmap_publication_vst\n\n",
      "Outputs:\n",
      "  Heatmap_top50_DEG_BP_vs_HC_vst_rowZ.png\n",
      "  Heatmap_top50_DEG_BP_vs_HC_vst_rowZ.tiff\n",
      "  Heatmap_top50_DEG_BP_vs_HC_vst_rowZ.svg\n"
    )
  )
}

parse_arguments <- function(args) {
  defaults <- list(
    counts = NULL,
    gtf = NULL,
    metadata = NULL,
    outdir = "BPvsHC_mRNA_heatmap_publication_vst",
    reference = "SANO",
    treatment = "PB",
    fdr = "0.05",
    fc = "1.5",
    top = "50",
    dpi = "600"
  )

  if (length(args) == 0L || any(args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0L)
  }

  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!grepl("^--", token)) {
      stop("Unexpected positional argument: ", token, "\nUse --help.")
    }

    if (grepl("=", token, fixed = TRUE)) {
      pieces <- strsplit(sub("^--", "", token), "=", fixed = TRUE)[[1L]]
      key <- pieces[[1L]]
      value <- paste(pieces[-1L], collapse = "=")
      i <- i + 1L
    } else {
      key <- sub("^--", "", token)
      if (i == length(args)) stop("Missing value for --", key)
      value <- args[[i + 1L]]
      i <- i + 2L
    }

    if (!key %in% names(defaults)) {
      stop("Unknown option --", key, "\nUse --help.")
    }
    defaults[[key]] <- value
  }

  defaults
}

as_single_number <- function(x, name, lower = -Inf, upper = Inf,
                             lower_open = FALSE, integer = FALSE) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || !is.finite(value)) {
    stop(name, " must be one finite number")
  }
  lower_ok <- if (lower_open) value > lower else value >= lower
  if (!lower_ok || value > upper) {
    stop(name, " is outside the permitted range")
  }
  if (integer && value != round(value)) stop(name, " must be an integer")
  value
}

args <- parse_arguments(commandArgs(trailingOnly = TRUE))

if (is.null(args$counts) || !nzchar(args$counts)) {
  stop("--counts is required. Use --help.")
}
if (is.null(args$gtf) || !nzchar(args$gtf) ||
    tolower(args$gtf) %in% c("none", "na")) {
  stop("--gtf is required so that heatmap labels use the exact bambu annotation.")
}
if (identical(args$reference, args$treatment)) {
  stop("Reference and treatment labels must differ")
}

fdr_cutoff <- as_single_number(
  args$fdr, "--fdr", lower = 0, upper = 1, lower_open = TRUE
)
linear_fc_cutoff <- as_single_number(
  args$fc, "--fc", lower = 1, lower_open = TRUE
)
top_n <- as.integer(as_single_number(
  args$top, "--top", lower = 2, integer = TRUE
))
tiff_dpi <- as.integer(as_single_number(
  args$dpi, "--dpi", lower = 300, integer = TRUE
))
abs_log2fc_cutoff <- log2(linear_fc_cutoff)

counts_file <- normalizePath(args$counts, mustWork = TRUE)
gtf_file <- normalizePath(args$gtf, mustWork = TRUE)
metadata_file <- NULL
if (!is.null(args$metadata) && nzchar(args$metadata) &&
    !tolower(args$metadata) %in% c("none", "na")) {
  metadata_file <- normalizePath(args$metadata, mustWork = TRUE)
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(args$outdir, mustWork = TRUE)

required_packages <- c("DESeq2", "pheatmap", "rtracklayer")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing R packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall DESeq2/rtracklayer with BiocManager and pheatmap from CRAN."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(pheatmap)
})

message("[1/6] Reading and validating bambu counts")
count_table <- read.delim(
  counts_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "NaN", "")
)

if (!"GENEID" %in% names(count_table)) {
  stop("The counts file must contain a column named GENEID")
}
if (nrow(count_table) < 2L) stop("The counts file contains fewer than two genes")
if (ncol(count_table) < 3L) stop("The counts file contains fewer than two samples")

ensembl_ids <- trimws(as.character(count_table$GENEID))
if (anyNA(ensembl_ids) || any(!nzchar(ensembl_ids))) {
  stop("GENEID contains missing or empty identifiers")
}
if (anyDuplicated(ensembl_ids)) {
  duplicate_ids <- unique(ensembl_ids[duplicated(ensembl_ids)])
  stop(
    "Duplicated GENEID values (first examples): ",
    paste(head(duplicate_ids, 10L), collapse = ", ")
  )
}

sample_columns <- setdiff(names(count_table), "GENEID")
count_frame <- count_table[, sample_columns, drop = FALSE]
numeric_columns <- lapply(count_frame, function(x) suppressWarnings(as.numeric(x)))
count_matrix_double <- as.matrix(
  data.frame(numeric_columns, check.names = FALSE)
)
rownames(count_matrix_double) <- ensembl_ids
colnames(count_matrix_double) <- sample_columns

if (anyNA(count_matrix_double)) {
  stop("The sample columns contain missing or non-numeric values")
}
if (any(!is.finite(count_matrix_double))) {
  stop("The count matrix contains non-finite values")
}
if (any(count_matrix_double < 0)) stop("The count matrix contains negative values")

max_rounding_difference <- max(abs(count_matrix_double - round(count_matrix_double)))
if (max_rounding_difference > 1e-8) {
  message("bambu counts include decimal estimates; values were rounded for DESeq2 ",
          "(maximum absolute change: ",
          signif(max_rounding_difference, 4L), ").")
}
count_matrix <- round(count_matrix_double)
if (max(count_matrix) > .Machine$integer.max) {
  stop("At least one rounded count exceeds R's integer limit")
}
storage.mode(count_matrix) <- "integer"
if (any(colSums(count_matrix) == 0L)) {
  stop(
    "Samples with zero total counts after rounding: ",
    paste(colnames(count_matrix)[colSums(count_matrix) == 0L], collapse = ", ")
  )
}
if (sum(rowSums(count_matrix) > 0L) < 2L) {
  stop("Fewer than two genes have non-zero counts")
}

extract_barcode_number <- function(x) {
  patterns <- c(
    "barcode_barcode([0-9]+)",
    "barcode([0-9]+)",
    "(?:^|[_-])([0-9]{1,3})(?:[_-]|$)"
  )
  output <- rep(NA_character_, length(x))
  for (pattern in patterns) {
    unresolved <- which(is.na(output))
    if (length(unresolved) == 0L) break
    matches <- regexec(pattern, x[unresolved], perl = TRUE, ignore.case = TRUE)
    values <- regmatches(x[unresolved], matches)
    has_value <- lengths(values) >= 2L
    if (any(has_value)) {
      output[unresolved[has_value]] <- vapply(
        values[has_value], `[[`, character(1L), 2L
      )
    }
  }
  output
}

infer_metadata <- function(columns, reference, treatment) {
  lower_names <- tolower(columns)
  reference_match <- grepl(
    "healthy|sano|control|(^|[_-])hc([_-]|$)",
    lower_names,
    perl = TRUE
  )
  treatment_match <- grepl(
    "pem|bullous|(^|[_-])pb([_-]|$)",
    lower_names,
    perl = TRUE
  )

  if (any(reference_match & treatment_match)) {
    stop("Automatic metadata inference assigned both groups to at least one sample")
  }
  if (any(!reference_match & !treatment_match)) {
    stop(
      "Automatic condition inference failed for: ",
      paste(columns[!reference_match & !treatment_match], collapse = ", "),
      "\nProvide a metadata TSV with --metadata."
    )
  }

  condition <- ifelse(reference_match, reference, treatment)
  barcode <- extract_barcode_number(columns)
  sample_id <- columns
  if (all(!is.na(barcode))) {
    width <- max(2L, max(nchar(barcode)))
    padded <- sprintf(paste0("%0", width, "d"), as.integer(barcode))
    sample_id <- paste0(padded, "_", condition)
  }
  if (anyDuplicated(sample_id)) sample_id <- columns

  data.frame(
    count_column = columns,
    sample_id = sample_id,
    condition = condition,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

message("[2/6] Preparing sample metadata")
if (is.null(metadata_file)) {
  samples <- infer_metadata(
    sample_columns,
    reference = args$reference,
    treatment = args$treatment
  )
} else {
  samples <- read.delim(
    metadata_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  required_metadata <- c("count_column", "sample_id", "condition")
  missing_metadata <- setdiff(required_metadata, names(samples))
  if (length(missing_metadata) > 0L) {
    stop(
      "Metadata file is missing columns: ",
      paste(missing_metadata, collapse = ", ")
    )
  }
  samples <- samples[, required_metadata, drop = FALSE]
  for (column in names(samples)) samples[[column]] <- as.character(samples[[column]])
}

required_metadata <- c("count_column", "sample_id", "condition")
if (anyNA(samples[, required_metadata, drop = FALSE]) ||
    any(!nzchar(as.matrix(samples[, required_metadata, drop = FALSE])))) {
  stop("Sample metadata contains missing or empty required values")
}
if (anyDuplicated(samples$count_column)) stop("Duplicated count_column in metadata")
if (anyDuplicated(samples$sample_id)) stop("Duplicated sample_id in metadata")
if (!setequal(samples$count_column, sample_columns)) {
  stop("Metadata count_column values do not exactly match the counts file")
}

samples <- samples[match(sample_columns, samples$count_column), , drop = FALSE]
if (any(!samples$condition %in% c(args$reference, args$treatment))) {
  stop("This script requires exactly the reference and treatment conditions")
}
if (!all(c(args$reference, args$treatment) %in% samples$condition)) {
  stop("Reference or treatment condition is absent from the metadata")
}

samples$condition <- factor(
  samples$condition,
  levels = c(args$reference, args$treatment)
)
if (any(table(samples$condition) < 2L)) {
  stop("Reference and treatment must each contain at least two replicates")
}
if (identical(args$reference, "SANO") && identical(args$treatment, "PB") &&
    any(table(samples$condition) != 6L)) {
  stop("This manuscript heatmap expects exactly 6 SANO and 6 PB samples")
}

colnames(count_matrix) <- samples$sample_id
rownames(samples) <- samples$sample_id

design_formula <- ~ condition
design_matrix <- model.matrix(design_formula, data = samples)
if (qr(design_matrix)$rank < ncol(design_matrix)) {
  stop("The condition-only design is not full rank")
}

message("[3/6] Running DESeq2")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = samples,
  design = design_formula
)
dds <- DESeq(dds, quiet = TRUE)
res <- results(
  dds,
  contrast = c("condition", args$treatment, args$reference),
  alpha = fdr_cutoff,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

map_gene_symbols <- function(ids, gtf_path = NULL) {
  output <- rep(NA_character_, length(ids))
  names(output) <- ids
  ambiguous_count <- 0L

  if (is.null(gtf_path)) {
    warning(
      "No GTF supplied: Ensembl IDs will be retained as heatmap labels. ",
      "For the paper, rerun with the exact GTF used by bambu."
    )
    return(list(symbol = output, ambiguous = ambiguous_count))
  }

  message("[4/6] Mapping Ensembl IDs to display-only gene symbols")
  gtf <- rtracklayer::import(gtf_path)
  gtf_metadata <- as.data.frame(S4Vectors::mcols(gtf))
  if (!all(c("type", "gene_id", "gene_name") %in% names(gtf_metadata))) {
    stop("The GTF must provide type, gene_id and gene_name attributes")
  }

  gene_rows <- !is.na(gtf_metadata$type) & gtf_metadata$type == "gene"
  if (!any(gene_rows)) stop("The supplied GTF contains no gene features")
  annotation <- data.frame(
    gene_id = sub("\\.[0-9]+$", "", as.character(gtf_metadata$gene_id[gene_rows])),
    gene_symbol = trimws(as.character(gtf_metadata$gene_name[gene_rows])),
    stringsAsFactors = FALSE
  )
  annotation <- annotation[
    !is.na(annotation$gene_id) & nzchar(annotation$gene_id) &
      !is.na(annotation$gene_symbol) & nzchar(annotation$gene_symbol),
    , drop = FALSE
  ]
  annotation <- unique(annotation)

  symbols_by_id <- split(annotation$gene_symbol, annotation$gene_id)
  unique_symbol <- vapply(
    symbols_by_id,
    function(x) {
      values <- unique(x[!is.na(x) & nzchar(x)])
      if (length(values) == 1L) values[[1L]] else NA_character_
    },
    character(1L)
  )
  ambiguous_count <- sum(is.na(unique_symbol))
  clean_ids <- sub("\\.[0-9]+$", "", ids)
  output <- unname(unique_symbol[clean_ids])
  names(output) <- ids

  list(symbol = output, ambiguous = ambiguous_count)
}

symbol_mapping <- map_gene_symbols(rownames(res), gtf_file)

result <- data.frame(
  ENSEMBL = rownames(res),
  gene_symbol = unname(symbol_mapping$symbol[rownames(res)]),
  as.data.frame(res),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  row.names = NULL
)
result$foldChangeRatio <- ifelse(
  is.na(result$log2FoldChange), NA_real_, 2^result$log2FoldChange
)
result$signedLinearFoldChange <- ifelse(
  is.na(result$log2FoldChange),
  NA_real_,
  ifelse(
    result$log2FoldChange >= 0,
    2^result$log2FoldChange,
    -(2^(-result$log2FoldChange))
  )
)
result <- result[
  order(result$padj, result$pvalue, -abs(result$log2FoldChange),
        result$ENSEMBL, na.last = TRUE),
  , drop = FALSE
]

eligible <- result[
  is.finite(result$padj) & result$padj < fdr_cutoff &
    is.finite(result$log2FoldChange) &
    abs(result$log2FoldChange) >= abs_log2fc_cutoff,
  , drop = FALSE
]
eligible <- eligible[
  order(eligible$padj, eligible$pvalue, -abs(eligible$log2FoldChange),
        eligible$ENSEMBL),
  , drop = FALSE
]
eligible$DEG_rank <- seq_len(nrow(eligible))

if (nrow(eligible) < 2L) {
  stop("Fewer than two genes satisfy the declared FDR/fold-change criteria")
}
message("[5/6] Applying VST and selecting heatmap genes")
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
eligible_vst <- assay(vsd)[eligible$ENSEMBL, , drop = FALSE]
eligible_sd <- apply(eligible_vst, 1L, stats::sd)
usable <- is.finite(eligible_sd) & eligible_sd > 0
if (sum(usable) < 2L) {
  stop("Fewer than two eligible genes have non-zero finite VST variance")
}

selected_ids <- rownames(eligible_vst)[usable]
selected_ids <- head(selected_ids, min(top_n, length(selected_ids)))
selected_vst <- assay(vsd)[selected_ids, , drop = FALSE]

selected_table <- eligible[match(selected_ids, eligible$ENSEMBL), , drop = FALSE]
selected_table$heatmap_rank <- seq_len(nrow(selected_table))

display_labels <- selected_table$gene_symbol
use_symbol <- !is.na(display_labels) & nzchar(display_labels)
display_labels[!use_symbol] <- selected_table$ENSEMBL[!use_symbol]
duplicate_labels <- duplicated(display_labels) |
  duplicated(display_labels, fromLast = TRUE)
display_labels[duplicate_labels & use_symbol] <- paste0(
  display_labels[duplicate_labels & use_symbol],
  " [", selected_table$ENSEMBL[duplicate_labels & use_symbol], "]"
)
selected_table$heatmap_label <- display_labels
rownames(selected_vst) <- display_labels

# This matrix fixes the color limits only. pheatmap applies scale='row' once.
# The color range contains every observed z score: no clipping.
row_z <- t(scale(t(selected_vst)))
if (any(!is.finite(row_z))) {
  stop("Non-finite values were produced during row-wise z scoring")
}
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100L)
z_limit <- max(abs(row_z)) + 1e-8
heatmap_breaks <- seq(-z_limit, z_limit,
                      length.out = length(heatmap_colors) + 1L)

reference_display <- if (identical(args$reference, "SANO")) "HC" else args$reference
treatment_display <- if (identical(args$treatment, "PB")) "BP" else args$treatment
if (identical(reference_display, treatment_display)) {
  stop("Reference and treatment cannot have the same display label")
}
comparison <- paste0(treatment_display, "_vs_", reference_display)
display_group <- ifelse(
  as.character(samples$condition) == args$treatment,
  treatment_display,
  reference_display
)

# Default inferred IDs (such as 01_PB and 02_SANO) become BP_01 and HC_02.
# Metadata-supplied IDs without this pattern remain unchanged.
display_sample_ids <- as.character(samples$sample_id)
if (identical(args$reference, "SANO")) {
  display_sample_ids <- sub("^([0-9]+)_SANO$", "HC_\\1", display_sample_ids)
}
if (identical(args$treatment, "PB")) {
  display_sample_ids <- sub("^([0-9]+)_PB$", "BP_\\1", display_sample_ids)
}
if (anyDuplicated(display_sample_ids)) {
  stop("Displayed sample names are not unique; check --metadata")
}
colnames(selected_vst) <- display_sample_ids
message("Sample mapping (count column -> figure label [group]):")
for (i in seq_len(nrow(samples))) {
  message("  ", samples$count_column[[i]], " -> ",
          display_sample_ids[[i]], " [", display_group[[i]], "]")
}

annotation_col <- data.frame(
  condition = factor(
    display_group,
    levels = c(reference_display, treatment_display)
  ),
  row.names = display_sample_ids,
  check.names = FALSE
)
annotation_colors <- list(
  condition = stats::setNames(
    c("skyblue", "pink"),
    c(reference_display, treatment_display)
  )
)

message("[6/6] Building and exporting the heatmap")
heatmap_object <- pheatmap(
  selected_vst,
  color = heatmap_colors,
  breaks = heatmap_breaks,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = TRUE,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 9,
  fontsize_row = 7.5,
  fontsize_col = 8,
  angle_col = "45",
  treeheight_row = 45,
  treeheight_col = 45,
  legend = TRUE,
  annotation_legend = TRUE,
  main = "mRNA expression (VST; row z-score)",
  silent = TRUE
)

figure_width_in <- 10
figure_height_in <- 11
figure_stem <- file.path(
  output_dir,
  paste0("Heatmap_top", nrow(selected_table), "_DEG_", comparison,
         "_vst_rowZ")
)

draw_heatmap <- function(filename, device = c("png", "tiff", "svg")) {
  device <- match.arg(device)
  if (device == "png") {
    grDevices::png(
      filename, width = figure_width_in, height = figure_height_in,
      units = "in", res = tiff_dpi, bg = "white"
    )
  } else if (device == "tiff") {
    grDevices::tiff(
      filename, width = figure_width_in, height = figure_height_in,
      units = "in", res = tiff_dpi, compression = "lzw", bg = "white"
    )
  } else {
    grDevices::svg(
      filename, width = figure_width_in, height = figure_height_in,
      onefile = TRUE, bg = "white", family = "sans"
    )
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(heatmap_object$gtable)
}

draw_heatmap(paste0(figure_stem, ".png"), "png")
draw_heatmap(paste0(figure_stem, ".tiff"), "tiff")
draw_heatmap(paste0(figure_stem, ".svg"), "svg")

message("Completed")
message(
  "DEGs passing declared criteria: ", nrow(eligible), "\n",
  "Genes displayed: ", nrow(selected_table), "\n",
  "PNG:  ", paste0(figure_stem, ".png"), "\n",
  "TIFF: ", paste0(figure_stem, ".tiff"), "\n",
  "SVG:  ", paste0(figure_stem, ".svg")
)
