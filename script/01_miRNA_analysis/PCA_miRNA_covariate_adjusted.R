# PCA of miRNA raw counts with covariate adjustment for visualization
# Input files must be in the current working directory.
# The script creates exactly six figure files: PNG, TIFF and SVG for each PCA.

required_packages <- c("DESeq2", "limma", "readxl", "ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

counts_file <- "Dermatologia_miRNA_PAPER mirna pirna expression values(1).xlsx"
metadata_file <- "metadati.csv"

if (!file.exists(counts_file)) {
  stop("Input file not found: ", counts_file, call. = FALSE)
}
if (!file.exists(metadata_file)) {
  stop("Input file not found: ", metadata_file, call. = FALSE)
}

# Read raw miRNA counts. The first two columns contain miRNA annotations.
count_table <- as.data.frame(
  readxl::read_excel(counts_file, sheet = "miRNA Count"),
  check.names = FALSE
)

if (ncol(count_table) < 3L) {
  stop("The 'miRNA Count' sheet does not contain sample columns.", call. = FALSE)
}

feature_ids <- trimws(as.character(count_table[[1L]]))
if (anyNA(feature_ids) || any(feature_ids == "") || anyDuplicated(feature_ids)) {
  stop("The miRNA names in the first column must be present and unique.", call. = FALSE)
}

count_matrix <- data.matrix(count_table[, -(1:2), drop = FALSE])
rownames(count_matrix) <- feature_ids

if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
  stop("The count matrix contains missing or non-finite values.", call. = FALSE)
}
if (any(count_matrix < 0) || any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("The 'miRNA Count' sheet must contain non-negative integer counts.", call. = FALSE)
}

count_matrix <- round(count_matrix)
storage.mode(count_matrix) <- "integer"

# Read and validate metadata. The prurito column is intentionally not used.
metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  strip.white = TRUE,
  na.strings = c("", "NA", "N/A")
)

required_columns <- c("ID", "condition", "eta", "sesso", "ospedale")
missing_columns <- setdiff(required_columns, colnames(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata$ID <- trimws(as.character(metadata$ID))
metadata$condition <- trimws(as.character(metadata$condition))
metadata$sesso <- toupper(trimws(as.character(metadata$sesso)))
metadata$ospedale <- trimws(as.character(metadata$ospedale))
metadata$eta <- suppressWarnings(as.numeric(as.character(metadata$eta)))

if (anyDuplicated(metadata$ID)) {
  stop("Metadata IDs must be unique.", call. = FALSE)
}
if (!setequal(colnames(count_matrix), metadata$ID)) {
  missing_in_metadata <- setdiff(colnames(count_matrix), metadata$ID)
  missing_in_counts <- setdiff(metadata$ID, colnames(count_matrix))
  stop(
    paste0(
      "Sample IDs do not match between counts and metadata. ",
      "Missing in metadata: ", paste(missing_in_metadata, collapse = ", "),
      "; missing in counts: ", paste(missing_in_counts, collapse = ", ")
    ),
    call. = FALSE
  )
}

# Put metadata in exactly the same order as the count-matrix columns.
metadata <- metadata[match(colnames(count_matrix), metadata$ID), , drop = FALSE]
rownames(metadata) <- metadata$ID

condition_order <- c("HC", "BP", "AD", "Pso")
unexpected_conditions <- setdiff(unique(metadata$condition), condition_order)
if (length(unexpected_conditions) > 0L) {
  stop(
    "Unexpected condition labels: ",
    paste(unexpected_conditions, collapse = ", "),
    call. = FALSE
  )
}

metadata$condition <- factor(metadata$condition, levels = condition_order)

# Binary coding requested for sex: F = 0, M = 1.
if (anyNA(metadata$sesso) || !all(metadata$sesso %in% c("F", "M"))) {
  stop("The sesso column must contain only F or M, with no missing values.", call. = FALSE)
}
metadata$sex_M <- ifelse(metadata$sesso == "M", 1, 0)

# Hospital is handled as a two-level factor, not as a continuous number.
if (anyNA(metadata$ospedale) || !all(metadata$ospedale %in% c("1", "2"))) {
  stop("The ospedale column must contain only 1 or 2, with no missing values.", call. = FALSE)
}
metadata$ospedale <- factor(metadata$ospedale, levels = c("1", "2"))

# No imputation is performed: all 31 samples require a verified age.
if (anyNA(metadata$eta)) {
  stop(
    "Missing eta for: ",
    paste(metadata$ID[is.na(metadata$eta)], collapse = ", "),
    ". Enter the verified age before running the PCA.",
    call. = FALSE
  )
}

if (ncol(count_matrix) != 31L || nrow(metadata) != 31L) {
  stop("This analysis expects exactly 31 samples in both input files.", call. = FALSE)
}

condition_colors <- c(
  "HC" = "#7F7F7F",
  "BP" = "#D55E00",
  "AD" = "#0072B2",
  "Pso" = "#009E73"
)

make_adjusted_pca <- function(counts, meta, conditions_to_keep, plot_title) {
  use_samples <- meta$condition %in% conditions_to_keep
  meta_sub <- meta[use_samples, , drop = FALSE]
  meta_sub$condition <- droplevels(meta_sub$condition)
  meta_sub$ospedale <- droplevels(meta_sub$ospedale)
  counts_sub <- counts[, rownames(meta_sub), drop = FALSE]

  if (nlevels(meta_sub$condition) < 2L) {
    stop("At least two conditions are required for PCA.", call. = FALSE)
  }
  if (nlevels(meta_sub$ospedale) < 2L) {
    stop("Both hospital levels must be represented in each PCA subset.", call. = FALSE)
  }
  if (stats::sd(meta_sub$eta) == 0) {
    stop("Age has zero variance in a PCA subset.", call. = FALSE)
  }

  # Remove only completely uninformative miRNAs; selection of the 500 most
  # variable features is performed after transformation and adjustment,
  # matching the default DESeq2 plotPCA strategy.
  counts_sub <- counts_sub[rowSums(counts_sub) > 0L, , drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = ~ condition
  )

  vst_object <- DESeq2::varianceStabilizingTransformation(
    dds,
    blind = TRUE
  )
  vst_matrix <- SummarizedExperiment::assay(vst_object)

  # Preserve the biological condition while removing the additive effects of
  # hospital, standardized age and binary sex from the transformed matrix.
  condition_design <- stats::model.matrix(~ condition, data = meta_sub)
  numeric_covariates <- cbind(
    age_z = as.numeric(scale(meta_sub$eta)),
    sex_M = meta_sub$sex_M
  )
  hospital_design <- stats::model.matrix(~ ospedale, data = meta_sub)[, -1L, drop = FALSE]
  full_design <- cbind(condition_design, hospital_design, numeric_covariates)

  if (qr(full_design)$rank < ncol(full_design)) {
    stop(
      "Condition and covariates are not jointly estimable in the selected samples.",
      call. = FALSE
    )
  }

  adjusted_matrix <- limma::removeBatchEffect(
    vst_matrix,
    batch = meta_sub$ospedale,
    covariates = numeric_covariates,
    design = condition_design
  )

  if (anyNA(adjusted_matrix) || any(!is.finite(adjusted_matrix))) {
    stop("Covariate adjustment produced non-finite values.", call. = FALSE)
  }

  feature_variance <- apply(adjusted_matrix, 1L, stats::var)
  n_top <- min(500L, length(feature_variance))
  selected_features <- order(feature_variance, decreasing = TRUE)[seq_len(n_top)]

  pca <- stats::prcomp(
    t(adjusted_matrix[selected_features, , drop = FALSE]),
    center = TRUE,
    scale. = FALSE
  )
  percent_variance <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

  pca_data <- data.frame(
    PC1 = pca$x[, 1L],
    PC2 = pca$x[, 2L],
    condition = meta_sub$condition,
    sample = rownames(meta_sub),
    check.names = FALSE
  )

  ggplot2::ggplot(
    pca_data,
    ggplot2::aes(x = PC1, y = PC2, color = condition)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.3,
      color = "grey85"
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.3,
      color = "grey85"
    ) +
    ggplot2::geom_point(size = 3.2, alpha = 0.95) +
    ggplot2::scale_color_manual(
      values = condition_colors,
      breaks = levels(meta_sub$condition),
      drop = FALSE,
      name = "Condition"
    ) +
    ggplot2::labs(
      title = plot_title,
      x = sprintf("PC1 (%.1f%%)", percent_variance[1L]),
      y = sprintf("PC2 (%.1f%%)", percent_variance[2L])
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

save_plot_triplet <- function(plot_object, filename_stem) {
  width_mm <- 180
  height_mm <- 140

  ggplot2::ggsave(
    filename = paste0(filename_stem, ".png"),
    plot = plot_object,
    device = "png",
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )

  grDevices::tiff(
    filename = paste0(filename_stem, ".tiff"),
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 600,
    compression = "lzw",
    bg = "white"
  )
  print(plot_object)
  grDevices::dev.off()

  grDevices::svg(
    filename = paste0(filename_stem, ".svg"),
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    bg = "white",
    onefile = FALSE
  )
  print(plot_object)
  grDevices::dev.off()
}

# Both plots are fully built before any output file is written.
pca_all <- make_adjusted_pca(
  counts = count_matrix,
  meta = metadata,
  conditions_to_keep = c("HC", "BP", "AD", "Pso"),
  plot_title = "All samples"
)

pca_hc_bp <- make_adjusted_pca(
  counts = count_matrix,
  meta = metadata,
  conditions_to_keep = c("HC", "BP"),
  plot_title = "HC vs BP"
)

save_plot_triplet(pca_all, "PCA_all_samples_covariate_adjusted")
save_plot_triplet(pca_hc_bp, "PCA_HC_vs_BP_covariate_adjusted")
