#####################################################################
##### cDNA SEQ P2-solo
#### mkdir PNI_FGF9
## PNI
zcat /mnt/c/data/fastq_pass_PNI_part1/barcode01/*.fastq.gz > part1_barcode01_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part1/barcode02/*.fastq.gz > part1_barcode02_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part1/barcode03/*.fastq.gz > part1_barcode03_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part1/barcode04/*.fastq.gz > part1_barcode04_pni_combined.fastq

zcat /mnt/c/data/fastq_pass_PNI_part2/barcode01/*.fastq.gz > part2_barcode01_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part2/barcode02/*.fastq.gz > part2_barcode02_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part2/barcode03/*.fastq.gz > part2_barcode03_pni_combined.fastq
zcat /mnt/c/data/fastq_pass_PNI_part2/barcode04/*.fastq.gz > part2_barcode04_pni_combined.fastq

## Fgf9
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode05/*fastq.gz > part1_barcode05_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode06/*fastq.gz > part1_barcode06_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode07/*fastq.gz > part1_barcode07_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode08/*fastq.gz > part1_barcode08_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode09/*fastq.gz > part1_barcode09_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part1/barcode10/*fastq.gz > part1_barcode10_fgf9_combined.fastq

zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode05/*.fastq.gz > part2_barcode05_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode06/*.fastq.gz > part2_barcode06_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode07/*.fastq.gz > part2_barcode07_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode08/*.fastq.gz > part2_barcode08_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode09/*.fastq.gz > part2_barcode09_fgf9_combined.fastq
zcat /mnt/c/data/fastq_pass_FGF9_part2/barcode10/*.fastq.gz > part2_barcode10_fgf9_combined.fastq

### merge part 1 and 2
## mkdir fastqs_combined
cat mnt/c/data/PNI_FGF9_analysis/fastqs_combined/*_barcode01_pni_combined.fastq > barcode01_pni_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode02_pni_combined.fastq > barcode02_pni_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode03_pni_combined.fastq > barcode03_pni_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode04_pni_combined.fastq > barcode04_pni_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode05_fgf9_combined.fastq > barcode05_fgf9_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode06_fgf9_combined.fastq > barcode06_fgf9_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode07_fgf9_combined.fastq > barcode07_fgf9_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode08_fgf9_combined.fastq > barcode08_fgf9_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode09_fgf9_combined.fastq > barcode09_fgf9_combined_full.fastq
cat /mnt/c/data/drg_pni_analysis/fastqs_combined/*_barcode10_fgf9_combined.fastq > barcode10_fgf9_combined_full.fastq



NanoPlot -t 60 --fastq combined_barcode01_pni_combined.fastq --outdir Nanoplot_01_pni --plots hex dot
NanoPlot -t 60 --fastq combined_barcode02_pni_combined.fastq --outdir Nanoplot_02_pni --plots hex dot
NanoPlot -t 60 --fastq combined_barcode03_pni_combined.fastq --outdir Nanoplot_03_pni --plots hex dot
NanoPlot -t 60 --fastq combined_barcode04_pni_combined.fastq --outdir Nanoplot_04_pni --plots hex dot
NanoPlot -t 60 --fastq combined_barcode05_fgf9_combined.fastq --outdir Nanoplot_05_fgf9 --plots hex dot
NanoPlot -t 60 --fastq combined_barcode06_fgf9_combined.fastq --outdir Nanoplot_06_fgf9 --plots hex dot
NanoPlot -t 60 --fastq combined_barcode07_fgf9_combined.fastq --outdir Nanoplot_07_fgf9 --plots hex dot
NanoPlot -t 60 --fastq combined_barcode08_fgf9_combined.fastq --outdir Nanoplot_08_fgf9 --plots hex dot
NanoPlot -t 60 --fastq combined_barcode09_fgf9_combined.fastq --outdir Nanoplot_09_fgf9 --plots hex dot
NanoPlot -t 60 --fastq combined_barcode10_fgf9_combined.fastq --outdir Nanoplot_10_fgf9 --plots hex dot

####################################################################
### Fastcat and pychopper
fastcat --min_length 500 --min_qscore 10 barcode01_pni_combined_full.fastq > barcode01_pni_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode02_pni_combined_full.fastq > barcode02_pni_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode03_pni_combined_full.fastq > barcode03_pni_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode04_pni_combined_full.fastq > barcode04_pni_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode05_fgf9_combined_full.fastq > barcode05_fgf9_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode06_fgf9_combined_full.fastq > barcode06_fgf9_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode07_fgf9_combined_full.fastq > barcode07_fgf9_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode08_fgf9_combined_full.fastq > barcode08_fgf9_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode09_fgf9_combined_full.fastq > barcode09_fgf9_filtered.fastq
fastcat --min_length 500 --min_qscore 10 barcode10_fgf9_combined_full.fastq > barcode10_fgf9_filtered.fastq

pychopper -t 60 -r barcode01_pni_pychopper_report.txt -u barcode01_pni_unclassified.fastq barcode01_pni_filtered.fastq barcode01_pni_full_length.fastq
pychopper -t 60 -r barcode02_pni_pychopper_report.txt -u barcode02_pni_unclassified.fastq barcode02_pni_filtered.fastq barcode02_pni_full_length.fastq
pychopper -t 60 -r barcode03_pni_pychopper_report.txt -u barcode03_pni_unclassified.fastq barcode03_pni_filtered.fastq barcode03_pni_full_length.fastq
pychopper -t 60 -r barcode04_pni_pychopper_report.txt -u barcode04_pni_unclassified.fastq barcode04_pni_filtered.fastq barcode04_pni_full_length.fastq

pychopper -t 60 -r barcode05_fgf9_pychopper_report.txt -u barcode05_fgf9_unclassified.fastq barcode05_fgf9_filtered.fastq barcode05_fgf9_full_length.fastq
pychopper -t 60 -r barcode06_fgf9_pychopper_report.txt -u barcode06_fgf9_unclassified.fastq barcode06_fgf9_filtered.fastq barcode06_fgf9_full_length.fastq
pychopper -t 60 -r barcode07_fgf9_pychopper_report.txt -u barcode07_fgf9_unclassified.fastq barcode07_fgf9_filtered.fastq barcode07_fgf9_full_length.fastq
pychopper -t 60 -r barcode08_fgf9_pychopper_report.txt -u barcode08_fgf9_unclassified.fastq barcode08_fgf9_filtered.fastq barcode08_fgf9_full_length.fastq
pychopper -t 60 -r barcode09_fgf9_pychopper_report.txt -u barcode09_fgf9_unclassified.fastq barcode09_fgf9_filtered.fastq barcode09_fgf9_full_length.fastq
pychopper -t 60 -r barcode10_fgf9_pychopper_report.txt -u barcode10_fgf9_unclassified.fastq barcode10_fgf9_filtered.fastq barcode10_fgf9_full_length.fastq

####################################################################
### Minimap2 alignment
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode01_pni_full_length.fastq  > barcode01_pni_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode02_pni_full_length.fastq  > barcode02_pni_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode03_pni_full_length.fastq  > barcode03_pni_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode04_pni_full_length.fastq  > barcode04_pni_aln.sam

minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode05_fgf9_full_length.fastq  > barcode05_fgf9_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode06_fgf9_full_length.fastq  > barcode06_fgf9_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode07_fgf9_full_length.fastq  > barcode07_fgf9_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode08_fgf9_full_length.fastq  > barcode08_fgf9_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode09_fgf9_full_length.fastq  > barcode09_fgf9_aln.sam
minimap2 -t 60 -a -x splice /home/nasslab/ref/Mus_musculus.GRCm39.dna.primary_assembly.fa barcode10_fgf9_full_length.fastq  > barcode10_fgf9_aln.sam


####################################################################
#### Samtools
samtools view -S -b barcode01_pni_aln.sam > barcode01_pni.bam 
samtools view -S -b barcode02_pni_aln.sam > barcode02_pni.bam 
samtools view -S -b barcode03_pni_aln.sam > barcode03_pni.bam 
samtools view -S -b barcode04_pni_aln.sam > barcode04_pni.bam 
samtools view -S -b barcode05_fgf9_aln.sam > barcode05_fgf9.bam 
samtools view -S -b barcode06_fgf9_aln.sam > barcode06_fgf9.bam 
samtools view -S -b barcode07_fgf9_aln.sam > barcode07_fgf9.bam  
samtools view -S -b barcode08_fgf9_aln.sam > barcode08_fgf9.bam 
samtools view -S -b barcode09_fgf9_aln.sam > barcode09_fgf9.bam 
samtools view -S -b barcode10_fgf9_aln.sam > barcode10_fgf9.bam  

samtools sort -o barcode01_pni_sorted.bam barcode01_pni.bam 
samtools sort -o barcode02_pni_sorted.bam barcode02_pni.bam 
samtools sort -o barcode03_pni_sorted.bam barcode03_pni.bam 
samtools sort -o barcode04_pni_sorted.bam barcode04_pni.bam 
samtools sort -o barcode05_fgf9_sorted.bam barcode05_fgf9.bam 
samtools sort -o barcode06_fgf9_sorted.bam barcode06_fgf9.bam 
samtools sort -o barcode07_fgf9_sorted.bam barcode07_fgf9.bam 
samtools sort -o barcode08_fgf9_sorted.bam barcode08_fgf9.bam 
samtools sort -o barcode09_fgf9_sorted.bam barcode09_fgf9.bam 
samtools sort -o barcode10_fgf9_sorted.bam barcode10_fgf9.bam 

samtools index barcode01_pni_sorted.bam
samtools index barcode02_pni_sorted.bam
samtools index barcode03_pni_sorted.bam
samtools index barcode04_pni_sorted.bam
samtools index barcode05_fgf9_sorted.bam
samtools index barcode06_fgf9_sorted.bam
samtools index barcode07_fgf9_sorted.bam
samtools index barcode08_fgf9_sorted.bam
samtools index barcode09_fgf9_sorted.bam
samtools index barcode10_fgf9_sorted.bam

#####################################################################
#### FeatureCounts quantification
#mkdir featurecounts
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode01pni_counts.txt -t exon -g gene_id -B ../barcode01_pni_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode02pni_counts.txt -t exon -g gene_id -B ../barcode02_pni_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode03pni_counts.txt -t exon -g gene_id -B ../barcode03_pni_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode04pni_counts.txt -t exon -g gene_id -B ../barcode04_pni_sorted.bam

featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode05fgf9_counts.txt -t exon -g gene_id -B ../barcode05_fgf9_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode06fgf9_counts.txt -t exon -g gene_id -B ../barcode06_fgf9_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode07fgf9_counts.txt -t exon -g gene_id -B ../barcode07_fgf9_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode08fgf9_counts.txt -t exon -g gene_id -B ../barcode08_fgf9_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode09fgf9_counts.txt -t exon -g gene_id -B ../barcode09_fgf9_sorted.bam
featureCounts -L -T 60 --maxMOp 100 -a /home/nasslab/ref/Mus_musculus.GRCm39.113.gtf  -o barcode10fgf9_counts.txt -t exon -g gene_id -B ../barcode10_fgf9_sorted.bam


####################################################################
####################################################################
####################################################################
##################### into R #######################################
#### Libraries needed
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(RColorBrewer) 
library(enrichR)
library(pheatmap)
library(biomaRt)
library(ggrepel)
library(limma)


####################################################################
### Loading data
b2_old = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_17122024/barcode02veh_counts.txt",header=TRUE)
b3_old = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_17122024/barcode03veh_counts.txt",header=TRUE)
b11 = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_10012025/barcode11veh_counts.txt",header=TRUE)     
b12 = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_10012025/barcode12veh_counts.txt",header=TRUE)    
b13 = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_10012025/barcode13veh_counts.txt",header=TRUE)  
b14 = read.table("/home/nasslab/PNI_FGF9_analysis/featuresCounts_fgf9_veh_10012025/barcode14veh_counts.txt",header=TRUE)   

b1 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode01pni_counts.txt",header=TRUE)
b2 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode02pni_counts.txt",header=TRUE)
b3 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode03pni_counts.txt",header=TRUE)
b4 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode04pni_counts.txt",header=TRUE)

b5 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode05fgf9_counts.txt",header=TRUE)
b6 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode06fgf9_counts.txt",header=TRUE)     
b7 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode07fgf9_counts.txt",header=TRUE)    
b8 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode08fgf9_counts.txt",header=TRUE)  
b9 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode09fgf9_counts.txt",header=TRUE)   
b10 = read.table("/home/nasslab/PNI_FGF9_analysis/fastqs_combined/featurecounts/barcode10fgf9_counts.txt",header=TRUE)   


####################################################################
### Prep DESeq2 object PNI PNI PNI PNI PNI PNI v1

counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b13[,7],b14[,7],b1[,7],b2[,7],b3[,7],b4[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode13","barcode14","barcode01","barcode02","barcode03","barcode04")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","veh","veh","pni","pni","pni","pni")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_pni_fgf9_21_1_25_labelled.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 

####################################################################
### Prep DESeq2 object PNI PNI PNI PNI PNI PNI v2

counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b13[,7],b1[,7],b2[,7],b3[,7],b4[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode13","barcode01","barcode02","barcode03","barcode04")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","veh","pni","pni","pni","pni")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_pni_fgf9_21_1_25_labelled.v2.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 

####################################################################
### Differential expression analysis and Gene symbol convertion
dds <- DESeq(dds)  # Run the analysis
results <- results(dds)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)  # Mouse dataset
results <- as.data.frame(results)
ensembl_ids <- rownames(results)
gene_info <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl)
results <- merge(results, gene_info, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(results)[1] <- "Ensembl_ID"
filtered_df <- results[!is.na(results$padj) & results$padj < 0.1, ]
filtered_df <- filtered_df[filtered_df$log2FoldChange > 0.58 | filtered_df$log2FoldChange < -0.58, ]

###########################################################################################################
####################################################################
### Prep DESeq2 object FGF9 3.5 H 

counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b13[,7],b14[,7],b5[,7],b6[,7],b7[,7],b8[,7],b9[,7],b10[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode13","barcode14","barcode05","barcode06","barcode07","barcode08","barcode09","barcode10")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","veh","veh","fgf9","fgf9","fgf9","fgf9","fgf9","fgf9")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_pni_fgf9_21_1_25_labelled.v3.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 



###########################################################################################################
####################################################################
### Prep DESeq2 object FGF9 3.5 H V4

counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b13[,7],b5[,7],b7[,7],b8[,7],b9[,7],b10[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode13","barcode05","barcode07","barcode08","barcode09","barcode10")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","veh","fgf9","fgf9","fgf9","fgf9","fgf9")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_pni_fgf9_21_1_25_labelled.v4.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 

###########################################################################################################
####################################################################
### Prep DESeq2 object FGF9 3.5 H V5

counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b7[,7],b8[,7],b9[,7],b10[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode07","barcode08","barcode09","barcode10")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","fgf9","fgf9","fgf9","fgf9")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_pni_fgf9_21_1_25_labelled.v5.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 

dds <- DESeq(dds)  # Run the analysis
results <- results(dds)
ensembl <- useMart("ensembl")
ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)  # Mouse dataset
results <- as.data.frame(results)
ensembl_ids <- rownames(results)
gene_info <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                   filters = "ensembl_gene_id",
                   values = ensembl_ids,
                   mart = ensembl)
results <- merge(results, gene_info, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(results)[1] <- "Ensembl_ID"
filtered_df <- results[!is.na(results$padj) & results$padj < 0.01, ]
filtered_df <- filtered_df[filtered_df$log2FoldChange > 1 | filtered_df$log2FoldChange < -1, ]


library(dplyr)

# Conta il numero di righe duplicate (uguali in tutte le colonne)
filtered_dfu = filtered_df
filtered_dfu$mgi_symbol= NULL
duplicates <- filtered_dfu %>%
  group_by(across(everything())) %>%
  summarize(count = n()) %>%
  filter(count > 1)

# Visualizza le righe duplicate
duplicates



counts = cbind(b2_old[,7],b3_old[,7],b11[,7],b12[,7],b13[,7],b14[,7],b1[,7],b2[,7],b3[,7],b4[,7],b5[,7],b6[,7],b7[,7],b8[,7],b9[,7],b10[,7])
rownames(counts) = b2[,1]
colnames(counts) = c("barcode02veh","barcode03veh","barcode11","barcode12","barcode13","barcode14","barcode01","barcode02","barcode03","barcode04","barcode05","barcode06","barcode07","barcode08","barcode09","barcode10")         
counts = as.data.frame(counts)
sample_names <- colnames(counts) # Create sample names
conditions <- c("veh","veh","veh","veh","veh","veh","pni","pni","pni","pni","fgf9","fgf9","fgf9","fgf9","fgf9","fgf9")
batch = c("Lib1","Lib1","Lib2","Lib2","Lib2","Lib2","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3","Lib3")
metadata <- data.frame(sample = sample_names, condition = conditions,batch= batch)
rownames(metadata) <- sample_names
colnames(counts) <- row.names(metadata) # Make sure the row names of metadata match the column names of the count matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ condition)
dds # Check the DESeqDataSet object
dds$condition <- relevel(dds$condition, ref = "veh")

#### PCA
vst_data <- varianceStabilizingTransformation(dds)
mat <- assay(vst_data)
mat_corrected <- limma::removeBatchEffect(mat, batch = colData(dds)$batch)
pca <- prcomp(t(mat_corrected))
pca_df <- as.data.frame(pca$x)
pca_df$condition <- colData(dds)$condition

svg("PCA_all.svg")
ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  # Size of points
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, max.overlaps = 10) +  # Add labels with rownames
  labs(title = "Fgf9 vs veh",
       x = paste0("PC1: ", round(100 * summary(pca)$importance[2, 1], 1), "% variance"),
       y = paste0("PC2: ", round(100 * summary(pca)$importance[2, 2], 1), "% variance")) +
  theme_minimal()
dev.off() 