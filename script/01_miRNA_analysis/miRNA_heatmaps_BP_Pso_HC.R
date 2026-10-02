#!/usr/bin/env Rscript

# Heatmap miRNA per il manoscritto: HC, BP e Pso (AD esclusa).
# Input: gli stessi due file Excel dell'analisi precedente, nella stessa cartella.
# Output: due heatmap, ciascuna in PNG e SVG, nella cartella indicata sotto.
# Il test DESeq2 usa i conteggi interi originali. Solo per la visualizzazione,
# i conteggi vengono trasformati con VST (blind=FALSE) e pheatmap applica
# lo z-score separatamente a ciascun miRNA (scale='row').

suppressPackageStartupMessages({
  library(DESeq2)
  library(readxl)
  library(pheatmap)
})

counts_file <- "Dermatologia_miRNA_PAPER mirna pirna expression values(8).xlsx"
gene_globe_file <- "Pso vs. HC, BP vs. HC, AD vs. HC_selection(5).xlsx"
output_dir <- "Heatmaps_miRNA_BP_Pso_HC_publication_vst"

for (file in c(counts_file, gene_globe_file)) {
  if (!file.exists(file)) stop("File di input non trovato: ", file)
}

# 1. Conteggi grezzi; i nove campioni AD non vengono inclusi nel modello.
tbl <- as.data.frame(
  read_excel(counts_file, sheet = "miRNA Count", .name_repair = "minimal"),
  check.names = FALSE
)
if (!"Name" %in% names(tbl)) stop("Manca la colonna Name nei conteggi.")
feature_ids <- trimws(as.character(tbl$Name))
if (anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
    anyDuplicated(feature_ids)) {
  stop("Gli identificativi miRNA devono essere presenti e univoci.")
}

sample_headers <- grep(
  "^(SANO|PB|PSO)_[^ ]+ count$", names(tbl), value = TRUE
)
sample_ids <- sub(" count$", "", sample_headers)
prefix <- sub("_.*$", "", sample_ids)
condition_map <- c(SANO = "HC", PB = "BP", PSO = "Pso")
condition <- factor(unname(condition_map[prefix]),
                    levels = c("HC", "BP", "Pso"))

observed <- table(condition)
expected <- c(HC = 6L, BP = 7L, Pso = 9L)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("Numero inatteso di campioni: ",
       paste(names(observed), observed, collapse = ", "))
}
sample_numbers <- sub("^[^_]+_([^_]+)_.*$", "\\1", sample_ids)
sample_labels <- paste0(as.character(condition), "_", sample_numbers)
if (anyDuplicated(sample_ids) || anyDuplicated(sample_labels)) {
  stop("Identificativi o etichette di campione duplicati.")
}

count_matrix <- as.matrix(tbl[, sample_headers, drop = FALSE])
if (!is.numeric(count_matrix) || anyNA(count_matrix) ||
    any(!is.finite(count_matrix)) || any(count_matrix < 0) ||
    any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
  stop("La matrice deve contenere conteggi grezzi interi e non negativi.")
}
if (max(count_matrix) > .Machine$integer.max) {
  stop("Almeno un conteggio supera il limite per gli interi R.")
}
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- feature_ids
colnames(count_matrix) <- sample_ids

metadata <- data.frame(
  condition = condition,
  row.names = sample_ids
)
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ condition
)

# 2. Come nello script precedente, i primi 50 derivano dal test omnibus
#    LRT (HC, BP, Pso). Un test a tre gruppi non ha un unico log2 fold change.
dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)
lrt <- as.data.frame(results(dds, alpha = 0.05))
passing <- rownames(lrt)[!is.na(lrt$padj) & lrt$padj < 0.05]
if (length(passing) < 2L) {
  stop("Meno di due miRNA superano FDR < 0.05 nel test omnibus.")
}
passing <- passing[order(lrt[passing, "padj"],
                         lrt[passing, "pvalue"], passing)]
top_miRNAs <- head(passing, 50L)

# 3. I 30 miRNA condivisi vengono selezionati dalla tabella GeneGlobe:
#    FC lineare con modulo > 1.5 e FDR <= 0.05 in entrambi i confronti.
#    L'Excel fornito e' gia' una lista prefiltrata di 30 righe: per verificare
#    che non manchi alcun miRNA serve il report GeneGlobe completo.
gg <- as.data.frame(
  read_excel(gene_globe_file, sheet = "Sheet1", .name_repair = "minimal"),
  check.names = FALSE
)
needed <- c(
  "Name",
  "BP vs. HC Fold change", "BP vs. HC FDR p-value",
  "Pso vs. HC Fold change", "Pso vs. HC FDR p-value"
)
if (!all(needed %in% names(gg))) {
  stop("Mancano colonne nella tabella GeneGlobe: ",
       paste(setdiff(needed, names(gg)), collapse = ", "))
}
for (column in needed[-1L]) {
  if (!is.numeric(gg[[column]]) || anyNA(gg[[column]]) ||
      any(!is.finite(gg[[column]]))) {
    stop("Colonna GeneGlobe non numerica o incompleta: ", column)
  }
}
keep <- abs(gg[["BP vs. HC Fold change"]]) > 1.5 &
        gg[["BP vs. HC FDR p-value"]] <= 0.05 &
        abs(gg[["Pso vs. HC Fold change"]]) > 1.5 &
        gg[["Pso vs. HC FDR p-value"]] <= 0.05
shared_miRNAs <- trimws(as.character(gg$Name[keep]))
if (length(shared_miRNAs) != 30L || anyNA(shared_miRNAs) ||
    any(!nzchar(shared_miRNAs)) || anyDuplicated(shared_miRNAs)) {
  stop("La tabella GeneGlobe non produce esattamente 30 miRNA univoci.")
}
if (!all(c(top_miRNAs, shared_miRNAs) %in% rownames(dds))) {
  stop("Alcuni miRNA selezionati non sono nella matrice dei conteggi.")
}

# 4. Figura: VST sui conteggi dei 22 campioni, poi z-score per riga.
#    La VST e lo z-score NON cambiano i p-value o le liste dei miRNA.
#    Correlazione per le righe, distanza euclidea per i campioni,
#    complete linkage. Scala blu-bianco-rosso centrata su zero; la scala
#    dei colori e' identica nelle due figure e non taglia alcun valore.
# La trasformazione completa non richiede il minimo di 1000 righe espresse
# imposto dalla funzione veloce DESeq2::vst() con nsub predefinito.
vst_counts <- SummarizedExperiment::assay(
  varianceStabilizingTransformation(dds, blind = FALSE)
)
annotation_col <- data.frame(
  condition = factor(condition, levels = c("HC", "BP", "Pso")),
  row.names = sample_labels
)
annotation_colors <- list(
  condition = c(HC = "skyblue", BP = "pink", Pso = "orange")
)
hm_colors <- colorRampPalette(c("blue", "white", "red"))(100)

all_plotted <- unique(c(top_miRNAs, shared_miRNAs))
sd_plotted <- apply(vst_counts[all_plotted, , drop = FALSE], 1L, sd)
if (any(!is.finite(sd_plotted) | sd_plotted == 0)) {
  stop("MiRNA selezionati con varianza nulla dopo VST: ",
       paste(all_plotted[!is.finite(sd_plotted) | sd_plotted == 0],
             collapse = ", "))
}
# Il calcolo qui serve solo a fissare gli stessi limiti cromatici per le due
# figure; pheatmap(scale='row') esegue la scalatura una sola volta nel grafico.
all_row_z <- t(scale(t(vst_counts[all_plotted, , drop = FALSE])))
z_limit <- max(abs(all_row_z)) + 1e-8
hm_breaks <- seq(-z_limit, z_limit, length.out = length(hm_colors) + 1L)

dir.create(output_dir, showWarnings = FALSE)

save_figure <- function(gtable, filename, format, height) {
  if (format == "png") {
    png(filename, width = 10, height = height, units = "in",
        res = 600, bg = "white")
  } else {
    svg(filename, width = 10, height = height, onefile = TRUE,
        bg = "white", family = "sans")
  }
  on.exit(dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(gtable)
}

plot_heatmap <- function(miRNAs, stem, height, font_size) {
  mat <- vst_counts[miRNAs, , drop = FALSE]
  colnames(mat) <- sample_labels

  hm <- pheatmap(
    mat,
    scale = "row",
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = hm_colors,
    breaks = hm_breaks,
    main = "miRNA expression (VST; row z-score)",
    show_rownames = TRUE,
    show_colnames = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    fontsize = 9,
    fontsize_row = font_size,
    fontsize_col = 8,
    angle_col = "45",
    border_color = NA,
    silent = TRUE
  )
  for (format in c("png", "svg")) {
    save_figure(hm$gtable, file.path(output_dir, paste0(stem, ".", format)),
                format, height)
  }
}

plot_heatmap(
  top_miRNAs,
  paste0("miRNA_BP_Pso_HC_Heatmap_top", length(top_miRNAs),
         "_overall_LRT_vst_rowZ"),
  height = 11,
  font_size = 7.5
)
plot_heatmap(
  shared_miRNAs,
  "miRNA_BP_Pso_HC_Heatmap_30_shared_BP_Pso_vst_rowZ",
  height = 9.4,
  font_size = 9
)

message(
  "Completato: ", length(top_miRNAs), " miRNA omnibus, ",
  length(shared_miRNAs), " miRNA condivisi; ",
  paste(names(observed), observed, collapse = ", "),
  "; quattro figure PNG/SVG (VST, z-score per riga) in ", output_dir
)
