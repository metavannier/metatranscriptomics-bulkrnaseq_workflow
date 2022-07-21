# #################################################################
# 
#                               DESeq2
#                          
# #################################################################

library(DESeq2)
library(readr)

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Loading the parameters
project=snakemake@params[["project"]]
samples=snakemake@params[["samples"]]
condition=snakemake@params[["condition"]]
ref_level=snakemake@params[["ref_level"]]
normalized_counts_file=snakemake@output[["normalized_counts_file"]]

# Rename column name of the count matrix as coldata
cts <- as.matrix(read.table(snakemake@input[["cts"]], header=T, row.names = 1))
cts_ordered <- matrix(, nrow = length(rownames(cts)), ncol = length(project))
colnames(cts_ordered) <- condition
for (i in 1:length(project)) {
  colnames(cts) <- lapply(colnames(cts), sub, pattern = paste(project[[i]],"_",sep=""), replacement = "")
  colnames(cts) <- lapply(colnames(cts), sub, pattern = paste("_",samples[[i]],"\\.bam",sep=""), replacement = "")
}
colnames(cts) <- lapply(colnames(cts), sub, pattern = paste("\\.",sep=""), replacement = "-")
# Reorder the column as sample.tsv 
genomic_idx <-  match(colnames(cts_ordered),colnames(cts))
cts  <- cts[,genomic_idx]

# Metadata for deseq2
# coldata.tsv and sample.tsv need the same order of sample
coldata_read <- read.delim(snakemake@input[["coldata"]], header=TRUE, comment.char="#", quote="")
coldata <- coldata_read[,-1]
rownames(coldata) <- coldata_read[,1]
coldata$condition <- factor(coldata_read$condition)
coldata$type <- factor(coldata_read$type)

# Remove samples outlier
rmproj_list = as.list(strsplit(snakemake@params[["rmproj_list"]], ",")[[1]])

if(length(rmproj_list)!=0){
  for (i in 1:length(rmproj_list)) {
      name <- rmproj_list[[i]]
      coldata <- coldata[-match((name), table = rownames(coldata)), ]
  }
}

# Check that sample names match in both files
if (all(colnames(cts) %in% rownames(coldata)) & all(colnames(cts) == rownames(coldata))){
  # Create the DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
} else {
  print("sample names doesn't match in both files")
}

# Remove uninformative columns (to do when filter not already done with the CPM threshold)
#dds <- dds[ rowSums(counts(dds)) > 10, ]

# Specifying the reference level
dds$condition <- relevel(dds$condition, ref = ref_level)

# DESeq: Normalization and preprocessing (counts divided by sample-specific size factors
# determined by median ratio of gene counts relative to geometric mean per gene)
dds <- DESeq(dds, parallel=parallel)
# To save the object in a file for later use
saveRDS(dds, file=snakemake@output[["rds"]])

# Already done in the DESeq function
#dds <- estimateSizeFactors( dds)
print(sizeFactors(dds))
# Save the normalized data matrix
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=normalized_counts_file, sep="\t", quote=F, col.names=NA)
