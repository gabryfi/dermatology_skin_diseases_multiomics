#!/bin/bash

set -euo pipefail
shopt -s nullglob

# ==================================================
# =============== PARAMETRI UTENTE =================
# ==================================================

# I 4 percorsi delle corse/librerie distinte
BATCH1_LIB1="/media/user/8Tb/endometriosi/RNAseq_endometriosi/lib1_endo/20260616_1431_P2S-02262-A_PAY70792_d633faf2/fastq_pass"
BATCH1_LIB2="/media/user/8Tb/endometriosi/RNAseq_endometriosi/lib2_endo/20260616_1431_P2S-02262-B_PBC10666_3aee24a9/fastq_pass"
BATCH2_LIB1="/media/user/8Tb/endometriosi/endometriosis_RNAseq2/lib1_endo/20260622_1633_P2S-02262-A_PAY70792_aad78502/fastq_pass"
BATCH2_LIB2="/media/user/8Tb/endometriosi/endometriosis_RNAseq2/lib2_endo/20260622_1633_P2S-02262-B_PBC10666_e854015b/fastq_pass"

CONTROL_LABEL="NO_PAIN"
TREATMENT_LABEL="PAIN"

CONTROL_BARCODES=(07 08 13 14 15 16)
TREATMENT_BARCODES=(01 02 03 04 05 06 09 11 12)

BATCH_NAME="Endometriosi"

# Nuova cartella di output dedicata ai 4 batch completi
OUTPUT_DIR="/media/user/8Tb/endometriosi_analisi_4batches"

THREADS=30

# Reference genome e GTF
REF="/home/user/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/home/user/reference/gtf_human/BULK_HUMAN_G38.gtf"

CIBERSORTX_GENE_ID_TYPE="gene_name"

# Filtri fastcat
MIN_READ_LENGTH=500
MIN_QSCORE=10

# Parametri featureCounts
FEATURECOUNTS_STRAND=0

# Pychopper
RUN_PYCHOPPER="yes"

# ==================================================
# ============ PARAMETRI AUTOMATICI ================
# ==================================================

# Definiamo i nomi e i percorsi dei 4 batch tecnici
BATCH_NAMES=("Batch1_Lib1" "Batch1_Lib2" "Batch2_Lib1" "Batch2_Lib2")
BATCH_PATHS=("$BATCH1_LIB1" "$BATCH1_LIB2" "$BATCH2_LIB1" "$BATCH2_LIB2")

RUN_DIR_B1_L1="$(dirname "$BATCH1_LIB1")"
RUN_DIR_B1_L2="$(dirname "$BATCH1_LIB2")"
RUN_DIR_B2_L1="$(dirname "$BATCH2_LIB1")"
RUN_DIR_B2_L2="$(dirname "$BATCH2_LIB2")"

RUN_USER="${SUDO_USER:-$(id -un)}"
RUN_GROUP="$(id -gn "$RUN_USER")"

LOG_DIR="$OUTPUT_DIR/00_logs"
METADATA_DIR="$OUTPUT_DIR/00_run_metadata"
MERGED_FASTQ_DIR="$OUTPUT_DIR/01_merged_fastq"
NANOPLOT_DIR="$OUTPUT_DIR/02_nanoplot"
FILTERED_FASTQ_DIR="$OUTPUT_DIR/03_fastcat_filtered_fastq"
PYCHOPPER_DIR="$OUTPUT_DIR/04_pychopper"
ALIGNMENT_DIR="$OUTPUT_DIR/05_minimap2_sam"
BAM_DIR="$OUTPUT_DIR/06_samtools_bam"
FEATURECOUNTS_DIR="$OUTPUT_DIR/07_featureCounts"
CIBERSORTX_DIR="$OUTPUT_DIR/08_CIBERSORTx"
DESEQ2_DIR="$OUTPUT_DIR/09_DESeq2"

SAMPLE_MANIFEST="$METADATA_DIR/sample_manifest.tsv"

ALL_BARCODES=("${CONTROL_BARCODES[@]}" "${TREATMENT_BARCODES[@]}")
# Il totale dei campioni ora tiene conto di tutte e 4 le corse indipendenti (15 barcode x 4 = 60)
TOTAL_SAMPLES=$((${#ALL_BARCODES[@]} * 4))

sample_name() {
    local barcode="$1"
    local condition="$2"
    local batch="$3"
    printf "%s_%s_barcode_%s_%s" "$BATCH_NAME" "$batch" "$barcode" "$condition"
}

progress_bar() {
    local step="$1"
    local current="$2"
    local total="$3"
    local label="$4"
    local width=30
    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))
    local percent=$(( current * 100 / total ))
    local bar_done
    local bar_empty

    bar_done="$(printf "%${filled}s" "" | tr " " "#")"
    bar_empty="$(printf "%${empty}s" "" | tr " " "-")"

    printf "\n[%s] [%s%s] %d/%d (%d%%) %s\n" \
        "$step" "$bar_done" "$bar_empty" "$current" "$total" "$percent" "$label"
}

check_command() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "Errore: comando non trovato nel PATH: $cmd"
        exit 1
    fi
}

check_fastq_inputs() {
    local batch_name="$1"
    local batch_path="$2"
    local barcode="$3"
    
    local barcode_dir="$batch_path/barcode${barcode}"
    if [ ! -d "$barcode_dir" ]; then
        echo "Errore: cartella barcode non trovata: $barcode_dir ($batch_name)"
        exit 1
    fi

    local files=("$barcode_dir"/*.fastq.gz)
    if [ "${#files[@]}" -eq 0 ]; then
        echo "Errore: nessun file .fastq.gz trovato in: $barcode_dir ($batch_name)"
        exit 1
    fi
}

write_manifest_row() {
    local barcode="$1"
    local condition="$2"
    local batch="$3"
    local sample
    sample="$(sample_name "$barcode" "$condition" "$batch")"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$barcode" \
        "$condition" \
        "$batch" \
        "$MERGED_FASTQ_DIR/${sample}_combined.fastq" \
        "$FILTERED_FASTQ_DIR/${sample}_filtered.fastq" \
        "$PYCHOPPER_DIR/${sample}_full_length.fastq" \
        "$ALIGNMENT_DIR/${sample}_aln.sam" \
        "$BAM_DIR/${sample}.bam" \
        "$BAM_DIR/${sample}_sorted.bam" \
        "$FEATURECOUNTS_DIR/${sample}_counts.txt" \
        >> "$SAMPLE_MANIFEST"
}

process_group() {
    local batch="$1"
    local condition="$2"
    shift 2
    local barcodes=("$@")
    local barcode

    for barcode in "${barcodes[@]}"; do
        write_manifest_row "$barcode" "$condition" "$batch"
    done
}

echo "======================================="
echo "STEP 0 - CONTROLLI E PREPARAZIONE"
echo "======================================="

for cmd in sudo zcat NanoPlot fastcat pychopper minimap2 samtools featureCounts Rscript; do
    check_command "$cmd"
done

# Controlli di integrità su tutti i 4 batch di input
for i in "${!BATCH_NAMES[@]}"; do
    check_fastq_inputs "${BATCH_NAMES[$i]}" "${BATCH_PATHS[$i]}" "$barcode"
done

for reference_file in "$REF" "$GTF"; do
    if [ ! -f "$reference_file" ]; then
        echo "Errore: file reference non trovato: $reference_file"
        exit 1
    fi
done

sudo mkdir -p "$OUTPUT_DIR"
sudo chown -R "$RUN_USER":"$RUN_GROUP" "$RUN_DIR_B1_L1" "$RUN_DIR_B1_L2" "$RUN_DIR_B2_L1" "$RUN_DIR_B2_L2" "$OUTPUT_DIR"
sudo chmod -R u+rwX "$RUN_DIR_B1_L1" "$RUN_DIR_B1_L2" "$RUN_DIR_B2_L1" "$RUN_DIR_B2_L2" "$OUTPUT_DIR"

mkdir -p \
    "$LOG_DIR" \
    "$METADATA_DIR" \
    "$MERGED_FASTQ_DIR" \
    "$NANOPLOT_DIR" \
    "$FILTERED_FASTQ_DIR" \
    "$PYCHOPPER_DIR" \
    "$ALIGNMENT_DIR" \
    "$BAM_DIR" \
    "$FEATURECOUNTS_DIR" \
    "$CIBERSORTX_DIR" \
    "$DESEQ2_DIR"

Rscript -e 'pkgs <- c("DESeq2", "pheatmap", "ggplot2", "matrixStats"); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing) > 0) stop("Pacchetti R mancanti: ", paste(missing, collapse = ", "))' \
    > "$LOG_DIR/R_package_check.log" 2>&1 || {
        echo "Errore: controllo pacchetti R fallito. Vedi: $LOG_DIR/R_package_check.log"
        exit 1
    }

printf "sample\tbarcode\tcondition\tbatch\tmerged_fastq\tfiltered_fastq\tfull_length_fastq\tsam_file\tbam_file\tsorted_bam_file\tcount_file\n" \
    > "$SAMPLE_MANIFEST"

for batch in "${BATCH_NAMES[@]}"; do
    process_group "$batch" "$CONTROL_LABEL" "${CONTROL_BARCODES[@]}"
    process_group "$batch" "$TREATMENT_LABEL" "${TREATMENT_BARCODES[@]}"
done

echo "Cartella output: $OUTPUT_DIR"
echo "Manifest campioni: $SAMPLE_MANIFEST"

#########################################
############### QC ######################
#########################################

echo "======================================="
echo "STEP 1 - CONCATENAZIONE CHUNKS + NANOPLOT + FASTCAT"
echo "======================================="

current=0
for i in "${!BATCH_NAMES[@]}"; do
    batch="${BATCH_NAMES[$i]}"
    b_path="${BATCH_PATHS[$i]}"

    for condition in "$CONTROL_LABEL" "$TREATMENT_LABEL"; do
        if [ "$condition" = "$CONTROL_LABEL" ]; then
            barcodes=("${CONTROL_BARCODES[@]}")
        else
            barcodes=("${TREATMENT_BARCODES[@]}")
        fi

        for barcode in "${barcodes[@]}"; do
            current=$((current + 1))
            sample="$(sample_name "$barcode" "$condition" "$batch")"
            files=("$b_path/barcode${barcode}"/*.fastq.gz)

            progress_bar "QC/FILTER" "$current" "$TOTAL_SAMPLES" "$sample"

            # zcat qui unisce i singoli chunk fastq.gz generati da ONT per questa specifica corsa
            zcat "${files[@]}" \
                > "$MERGED_FASTQ_DIR/${sample}_combined.fastq" \
                2> "$LOG_DIR/${sample}_zcat.log"

            NanoPlot \
                -t "$THREADS" \
                --fastq "$MERGED_FASTQ_DIR/${sample}_combined.fastq" \
                --outdir "$NANOPLOT_DIR/Nanoplot_${sample}" \
                --plots hex \
                > "$LOG_DIR/${sample}_nanoplot.log" \
                2>&1

            fastcat \
                --min_length "$MIN_READ_LENGTH" \
                --min_qscore "$MIN_QSCORE" \
                "$MERGED_FASTQ_DIR/${sample}_combined.fastq" \
                > "$FILTERED_FASTQ_DIR/${sample}_filtered.fastq" \
                2> "$LOG_DIR/${sample}_fastcat.log"
        done
    done
done

#########################################
############ PYCHOPPER ##################
#########################################

echo "======================================="
echo "STEP 2 - PYCHOPPER"
echo "======================================="

current=0
for batch in "${BATCH_NAMES[@]}"; do
    for condition in "$CONTROL_LABEL" "$TREATMENT_LABEL"; do
        if [ "$condition" = "$CONTROL_LABEL" ]; then
            barcodes=("${CONTROL_BARCODES[@]}")
        else
            barcodes=("${TREATMENT_BARCODES[@]}")
        fi

        for barcode in "${barcodes[@]}"; do
            current=$((current + 1))
            sample="$(sample_name "$barcode" "$condition" "$batch")"

            progress_bar "PYCHOPPER" "$current" "$TOTAL_SAMPLES" "$sample"

            if [ "$RUN_PYCHOPPER" = "yes" ]; then
                pychopper \
                    -t "$THREADS" \
                    -r "$PYCHOPPER_DIR/${sample}_pychopper_report.txt" \
                    -u "$PYCHOPPER_DIR/${sample}_unclassified.fastq" \
                    "$FILTERED_FASTQ_DIR/${sample}_filtered.fastq" \
                    "$PYCHOPPER_DIR/${sample}_full_length.fastq" \
                    > "$LOG_DIR/${sample}_pychopper.log" \
                    2>&1
            else
                cp "$FILTERED_FASTQ_DIR/${sample}_filtered.fastq" \
                    "$PYCHOPPER_DIR/${sample}_full_length.fastq"
                printf "Pychopper skipped. Filtered FASTQ copied to full_length FASTQ.\n" \
                    > "$LOG_DIR/${sample}_pychopper.log"
            fi
        done
    done
done

#########################################
############ ALIGNMENT ##################
#########################################

echo "======================================="
echo "STEP 3 - MINIMAP2 ALIGNMENT"
echo "======================================="

current=0
for batch in "${BATCH_NAMES[@]}"; do
    for condition in "$CONTROL_LABEL" "$TREATMENT_LABEL"; do
        if [ "$condition" = "$CONTROL_LABEL" ]; then
            barcodes=("${CONTROL_BARCODES[@]}")
        else
            barcodes=("${TREATMENT_BARCODES[@]}")
        fi

        for barcode in "${barcodes[@]}"; do
            current=$((current + 1))
            sample="$(sample_name "$barcode" "$condition" "$batch")"

            progress_bar "MINIMAP2" "$current" "$TOTAL_SAMPLES" "$sample"

            minimap2 \
                -t "$THREADS" \
                -a \
                -x splice \
                -uf \
                --secondary=no \
                "$REF" \
                "$PYCHOPPER_DIR/${sample}_full_length.fastq" \
                > "$ALIGNMENT_DIR/${sample}_aln.sam" \
                2> "$LOG_DIR/${sample}_minimap2.log"
        done
    done
done

#########################################
######### SAM -> SORTED BAM #############
#########################################

echo "======================================="
echo "STEP 4 - SAMTOOLS BAM PROCESSING"
echo "======================================="

current=0
for batch in "${BATCH_NAMES[@]}"; do
    for condition in "$CONTROL_LABEL" "$TREATMENT_LABEL"; do
        if [ "$condition" = "$CONTROL_LABEL" ]; then
            barcodes=("${CONTROL_BARCODES[@]}")
        else
            barcodes=("${TREATMENT_BARCODES[@]}")
        fi

        for barcode in "${barcodes[@]}"; do
            current=$((current + 1))
            sample="$(sample_name "$barcode" "$condition" "$batch")"

            progress_bar "SAMTOOLS" "$current" "$TOTAL_SAMPLES" "$sample"

            samtools view \
                -S -b \
                "$ALIGNMENT_DIR/${sample}_aln.sam" \
                > "$BAM_DIR/${sample}.bam" \
                2> "$LOG_DIR/${sample}_samtools_view.log"

            samtools sort \
                --threads "$THREADS" \
                -o "$BAM_DIR/${sample}_sorted.bam" \
                "$BAM_DIR/${sample}.bam" \
                > "$LOG_DIR/${sample}_samtools_sort.log" \
                2>&1

            samtools index \
                "$BAM_DIR/${sample}_sorted.bam" \
                > "$LOG_DIR/${sample}_samtools_index.log" \
                2>&1

            samtools quickcheck \
                "$BAM_DIR/${sample}_sorted.bam" \
                > "$LOG_DIR/${sample}_samtools_quickcheck.log" \
                2>&1

            rm "$ALIGNMENT_DIR/${sample}_aln.sam"
        done
    done
done

#########################################
########## FEATURECOUNTS ################
#########################################

echo "======================================="
echo "STEP 5 - FEATURECOUNTS"
echo "======================================="

current=0
for batch in "${BATCH_NAMES[@]}"; do
    for condition in "$CONTROL_LABEL" "$TREATMENT_LABEL"; do
        if [ "$condition" = "$CONTROL_LABEL" ]; then
            barcodes=("${CONTROL_BARCODES[@]}")
        else
            barcodes=("${TREATMENT_BARCODES[@]}")
        fi

        for barcode in "${barcodes[@]}"; do
            current=$((current + 1))
            sample="$(sample_name "$barcode" "$condition" "$batch")"

            progress_bar "FEATURECOUNTS" "$current" "$TOTAL_SAMPLES" "$sample"

            featureCounts \
                -L \
                -T "$THREADS" \
                --maxMOp 100 \
                -s "$FEATURECOUNTS_STRAND" \
                -a "$GTF" \
                -o "$FEATURECOUNTS_DIR/${sample}_counts.txt" \
                -t exon \
                -g gene_id \
                -B \
                "$BAM_DIR/${sample}_sorted.bam" \
                > "$LOG_DIR/${sample}_featureCounts.log" \
                2>&1
        done
    done
done

#########################################
######## INPUT PER CIBERSORTX ###########
#########################################

echo "======================================="
echo "STEP 6 - CIBERSORTx INPUT MATRIX"
echo "======================================="

Rscript - \
    "$SAMPLE_MANIFEST" \
    "$CIBERSORTX_DIR" \
    "$GTF" \
    "$CIBERSORTX_GENE_ID_TYPE" \
    > "$LOG_DIR/cibersortx_matrix.log" 2>&1 <<'RSCRIPT'

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- args[1]
cibersortx_out_dir <- args[2]
gtf_file <- args[3]
cibersortx_id_type <- args[4]

dir.create(cibersortx_out_dir, showWarnings = FALSE, recursive = TRUE)

samples <- read.delim(manifest_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

missing_files <- samples$count_file[!file.exists(samples$count_file)]
if (length(missing_files) > 0) {
    stop("File featureCounts mancanti: ", paste(missing_files, collapse = ", "))
}

read_featurecounts <- function(file, sample_name) {
    x <- read.delim(file, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
    count_col <- colnames(x)[ncol(x)]
    out <- data.frame(Geneid = x$Geneid, Length = x$Length, count = as.numeric(x[[count_col]]), stringsAsFactors = FALSE)
    colnames(out)[3] <- sample_name
    out
}

fc_list <- lapply(seq_len(nrow(samples)), function(i) {
    read_featurecounts(samples$count_file[i], samples$sample[i])
})

length_table <- fc_list[[1]][, c("Geneid", "Length")]
for (i in seq_along(fc_list)) {
    length_table <- merge(length_table, fc_list[[i]][, c("Geneid", "Length")], by = "Geneid", all = TRUE, suffixes = c("", paste0(".", i)))
}

length_cols <- grep("^Length", colnames(length_table), value = TRUE)
length_matrix <- as.matrix(length_table[, length_cols, drop = FALSE])
length_final <- apply(length_matrix, 1, function(z) z[!is.na(z)][1])
length_final <- data.frame(Geneid = length_table$Geneid, Length = as.numeric(length_final))

count_list <- lapply(fc_list, function(x) x[, c("Geneid", colnames(x)[3])])
counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)
counts[is.na(counts)] <- 0
counts <- merge(length_final, counts, by = "Geneid", all.y = TRUE)
counts <- counts[!is.na(counts$Length) & counts$Length > 0, ]

count_matrix <- as.matrix(counts[, samples$sample, drop = FALSE])
rownames(count_matrix) <- counts$Geneid

length_kb <- counts$Length / 1000
rpk <- sweep(count_matrix, 1, length_kb, "/")
scaling_factor <- colSums(rpk, na.rm = TRUE) / 1e6
tpm <- sweep(rpk, 2, scaling_factor, "/")
tpm[!is.finite(tpm)] <- 0

extract_gtf_attr <- function(attr, key) {
    pattern <- paste0(key, " \\\"([^\\\"]+)\\\"")
    m <- regexec(pattern, attr, perl = TRUE)
    reg <- regmatches(attr, m)
    vapply(reg, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1))
}

if (cibersortx_id_type == "gene_name") {
    gtf_lines <- readLines(gtf_file, warn = FALSE)
    gtf_lines <- gtf_lines[!grepl("^#", gtf_lines)]
    gtf <- read.delim(text = gtf_lines, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
    gtf <- gtf[gtf$V3 == "exon", ]
    gene_map <- data.frame(Geneid = extract_gtf_attr(gtf$V9, "gene_id"), Gene = extract_gtf_attr(gtf$V9, "gene_name"), stringsAsFactors = FALSE)
    gene_map <- gene_map[!is.na(gene_map$Geneid) & !is.na(gene_map$Gene) & gene_map$Gene != "", ]
    gene_map <- gene_map[!duplicated(gene_map$Geneid), ]

    tpm_df <- data.frame(Geneid = rownames(tpm), tpm, check.names = FALSE)
    tpm_df <- merge(gene_map, tpm_df, by = "Geneid", all.y = TRUE)
    tpm_df <- tpm_df[!is.na(tpm_df$Gene) & tpm_df$Gene != "", ]
    expr_tpm <- aggregate(. ~ Gene, data = tpm_df[, c("Gene", samples$sample)], FUN = sum)

    counts_df <- data.frame(Geneid = rownames(count_matrix), count_matrix, check.names = FALSE)
    counts_df <- merge(gene_map, counts_df, by = "Geneid", all.y = TRUE)
    counts_df <- counts_df[!is.na(counts_df$Gene) & counts_df$Gene != "", ]
    expr_counts <- aggregate(. ~ Gene, data = counts_df[, c("Gene", samples$sample)], FUN = sum)
} else {
    expr_tpm <- data.frame(Gene = rownames(tpm), tpm, check.names = FALSE)
    expr_counts <- data.frame(Gene = rownames(count_matrix), count_matrix, check.names = FALSE)
}

write.table(expr_tpm, file = file.path(cibersortx_out_dir, "cibersortx_mixture_TPM.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.csv(samples[, c("sample", "barcode", "condition", "batch")], file = file.path(cibersortx_out_dir, "cibersortx_samples_metadata.csv"), row.names = FALSE)

RSCRIPT

#########################################
############### DESEQ2 ##################
#########################################

echo "======================================="
echo "STEP 7 - DESEQ2 (CON SUPPORTO COLLAPSING ANTI-PSEUDOREPLICAZIONE)"
echo "======================================="

Rscript - \
    "$SAMPLE_MANIFEST" \
    "$DESEQ2_DIR" \
    "$TREATMENT_LABEL" \
    "$CONTROL_LABEL" \
    > "$LOG_DIR/deseq2.log" 2>&1 <<'RSCRIPT'

suppressMessages({
  library(DESeq2)
  library(pheatmap)
  library(ggplot2)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- args[1]
deseq2_out_dir <- args[2]
treatment <- args[3]
control <- args[4]

# ----------------------------------------------------------------------
# PARAMETRO CRITICO PER REPLICHE TECNICHE
# ----------------------------------------------------------------------
# TRUE (Consigliato): Genera grafici QC (PCA/Heatmap) con tutte e 60 le corse 
# separate, ma unisce (somma) i conteggi prima del test dei DEGs per evitare
# la pseudo-replicazione. Il design finale passerà automaticamente a '~ condition'.
# FALSE: Mantiene i 60 campioni indipendenti usando il design '~ batch + condition'.
COLLAPSE_REPLICATES <- TRUE
# ----------------------------------------------------------------------

dir.create(deseq2_out_dir, showWarnings = FALSE, recursive = TRUE)
samples <- read.delim(manifest_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

read_featurecounts_for_deseq <- function(file, sample_name) {
    df <- read.delim(file, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
    counts <- df[, ncol(df), drop = FALSE]
    rownames(counts) <- df$Geneid
    colnames(counts) <- sample_name
    counts
}

count_list <- lapply(seq_len(nrow(samples)), function(i) read_featurecounts_for_deseq(samples$count_file[i], samples$sample[i]))
count_matrix <- do.call(cbind, count_list)

coldata <- data.frame(
    condition = factor(samples$condition, levels = c(control, treatment)),
    batch = factor(samples$batch),
    barcode = factor(samples$barcode), # Mappa l'identità del paziente originale
    row.names = samples$sample,
    stringsAsFactors = FALSE
)

count_matrix <- round(as.matrix(count_matrix))
storage.mode(count_matrix) <- "integer"

# Costruiamo l'oggetto base con tutte e 60 le corse separate
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ batch + condition)
dds <- dds[rowSums(counts(dds)) > 0, ]

# 1. GRAFICI QC COMPLETI (Sulle 60 corse per valutare visivamente i 4 batch)
message("Generazione grafici QC su tutte e 60 le corse...")
rld <- rlog(dds, blind = FALSE)

pcaData <- plotPCA(rld, intgroup = c("condition", "batch"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = batch, label = barcode)) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_text(vjust = -0.8, size = 2.5, show.legend = FALSE) +
    xlab(paste0("PC1: ", percentVar[1], "%")) +
    ylab(paste0("PC2: ", percentVar[2], "%")) +
    theme_minimal() +
    ggtitle("PCA - 60 Corse distinte (4 Batch Tecnici)")

ggsave(filename = file.path(deseq2_out_dir, "PCA_plot_60_runs.png"), plot = pca_plot, width = 9, height = 6, dpi = 300)

rld_mat <- assay(rld)
top_genes <- head(order(rowVars(rld_mat), decreasing = TRUE), 30)
mat <- rld_mat[top_genes, , drop = FALSE]
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("condition", "batch"), drop = FALSE])
pheatmap(mat, annotation_col = anno, fontsize_row = 8, fontsize_col = 6, cluster_cols = TRUE, filename = file.path(deseq2_out_dir, "Heatmap_top30_genes_60_runs.png"))


# 2. STRATEGIA DI ANALISI DELLE REPLICHE
if (COLLAPSE_REPLICATES) {
    message("COLLAPSE_REPLICATES = TRUE: Aggrego le repliche tecniche per barcode...")
    # Somma i conteggi delle corse dello stesso barcode
    dds_analyzed <- collapseReplicates(dds, groupby = dds$barcode)
    # Poiché il disegno è perfettamente bilanciato (tutti i campioni sono in tutti i batch),
    # la somma pulisce l'effetto batch. Il design passa a ~ condition.
    design(dds_analyzed) <- ~ condition
    dds_analyzed <- DESeq(dds_analyzed)
} else {
    message("COLLAPSE_REPLICATES = FALSE: Analizzo i 60 campioni come indipendenti con design ~ batch + condition...")
    dds_analyzed <- DESeq(dds)
}

# Estrazione risultati statistici finali
res <- results(dds_analyzed, contrast = c("condition", treatment, control))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, res_df$pvalue, na.last = TRUE), ]

write.csv(res_df, file = file.path(deseq2_out_dir, "deseq_results_ordered_by_padj.csv"), row.names = FALSE)
write.csv(as.data.frame(counts(dds_analyzed, normalized = TRUE)), file = file.path(deseq2_out_dir, "normalized_counts_DESeq2.csv"))

# Grafico Volcano finale dei DEGs
volcano_df <- res_df
volcano_df$threshold <- "Not significant"
volcano_df$threshold[!is.na(volcano_df$padj) & volcano_df$padj < 0.05 & abs(volcano_df$log2FoldChange) > 1] <- "padj<0.05 & |log2FC|>1"
volcano_df$neglog10padj <- -log10(volcano_df$padj)

volcano_plot <- ggplot(volcano_df, aes(log2FoldChange, neglog10padj, color = threshold)) +
    geom_point(alpha = 0.6, na.rm = TRUE) +
    scale_color_manual(values = c("Not significant" = "grey", "padj<0.05 & |log2FC|>1" = "red")) +
    theme_minimal() +
    xlab(paste0("log2 Fold Change (", treatment, " vs ", control, ")")) +
    ylab("-log10 adjusted p-value") +
    ggtitle("Volcano plot dei Geni Differenzialmente Espressi") +
    theme(legend.position = "bottom")

ggsave(filename = file.path(deseq2_out_dir, "Volcano_plot.png"), plot = volcano_plot, width = 7, height = 5, dpi = 300)

sink(file.path(deseq2_out_dir, "sessionInfo.txt"))
sessionInfo()
sink()

cat("Analisi DESeq2 completata. Repliche collassate:", COLLAPSE_REPLICATES, "\n")
RSCRIPT

echo "======================================="
echo "PIPELINE COMPLETATA CON SUCCESSO!"
echo "======================================="