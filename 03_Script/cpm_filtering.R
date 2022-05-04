# #################################################################
# 
#        Calcul the cpm (Counts Per Million) and
#            filter the low expressed genes
#              
# #################################################################


library(MASS)
library(edgeR)

# Browse the list of file count to produce a matrix
file_list=list.files(snakemake@params[["path"]],pattern = "_count.txt$", full.names=TRUE)
output_count=snakemake@output[["count_df"]]
output_cpm=snakemake@output[["cpm"]]
output_filter_count=snakemake@output[["output_filter_count"]]

# Dataframe with the gene count for selected sample
dataframe_total_count <- data.frame()
# Remove the duplicate sequenced library with error in the fastqc report

rmrun_list = as.list(strsplit(snakemake@params[["rmrun_list"]], ",")[[1]])
if(length(rmrun_list)!=0){
  for (i in 1:length(rmrun_list)) {
      name <- rmrun_list[[i]]
      rmrun_file <- as.numeric(grep(pattern = name, file_list))
      file_list <- file_list[-rmrun_file]
  }
}

gene_name <- read.delim(file_list[[1]], header=FALSE, comment.char="#", quote="")
dataframe_total_count <- cbind(gene_name[,1])
for (i in 1:length(file_list)) {
  count <- read.delim(file_list[[i]], header=FALSE, comment.char="#", quote="")
  count[1,] <- lapply(count[1,], sub, pattern = ".*\\/.*\\/", replacement = "")
  dataframe_total_count <- cbind(dataframe_total_count, count[,7])
}

# write the dataframe with count value of each samples
write.table(dataframe_total_count, file = output_count, sep="\t", quote = FALSE,row.names = FALSE, col.names = FALSE)

# Matrix transformation for cpm calculation
mtx_total_count <- as.matrix(dataframe_total_count[-1,-1])
class(mtx_total_count) <- "numeric"
cpm <- cpm(mtx_total_count)
df_cpm <- as.data.frame(cpm)
df_cpm <- cbind(Gene_id=dataframe_total_count[-1,1], df_cpm)

## Remove row with n value < cpm
df_cpm[] <- lapply(df_cpm, function(x) {
    if(is.factor(x)) as.numeric(as.character(x)) else x
})
# Dataframe with the cpm filtrated
df_cpm_filter <- data.frame()
for (j in 1:nrow(df_cpm)) {
    flag = 0
    for (k in 2:ncol(df_cpm)) {
      if(df_cpm[j,k] >= snakemake@params[["thresh_cpm"]]){
        flag <- flag+1
      }
    }
    
    if(flag >= snakemake@params[["thresh_sample"]]){
        df_cpm_filter <- rbind(df_cpm_filter,df_cpm[j,])
    }
}

df_cpm_filter <- rbind(dataframe_total_count[1,],df_cpm_filter)
write.table(df_cpm_filter, file = output_cpm, sep=" ", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Dataframe with the count value of the genes without low expressed genes
colnames_genes = colnames(dataframe_total_count[1,])
dataframe_filtered_count <- data.frame()
for (i in 1:nrow(df_cpm_filter)) {
  core_genes=df_cpm_filter[i,1]
  coregenes=as.character(core_genes)
  ligne=subset(dataframe_total_count,dataframe_total_count[,1]==coregenes)
  dataframe_filtered_count <- rbind(dataframe_filtered_count,ligne)
}
# Matrix transformation for cpm calculation
write.table(dataframe_filtered_count, file = output_filter_count, quote = FALSE, row.names = FALSE, col.names = FALSE)


