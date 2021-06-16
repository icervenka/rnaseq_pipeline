#!/usr/bin/Rscript

## install packages
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("gridExtra")
#install.packages("pheatmap")
#install.packages("rjson")
#source("https://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("Biobase")
#biocLite("DESeq2")


# parse commandline arguments
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("usage: rnaseq-htseq.R directory

      Arguments:
      directory   - character, 'directory of source files for RNASeq'

      Example:
      ./rnaseq-htseq.R /home/rnaseq")

  q(save="no")
}


directory = args[1]


#plug in libraries
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(Biobase))
suppressMessages(library(biomaRt))
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(rjson))
suppressMessages(library(readr))
suppressMessages(library(apeglm))
suppressMessages(library(tibble))

# function declaration ----------------------------------------------------------------
#convert summary from DESeq into array of values (gene numbers with significant fold change)
summary_match = function(summary) {
  arr = unname(sapply(summary, function(x) {str_trim(str_match(x,":(.*)%" )[2])}))
  arr = arr[!is.na(arr)]
  arr = unname(sapply(arr, strsplit, ", "))
  return(arr)
}

scaleFUN <- function(x) sprintf("%.2f", x)

update_json = function(arr, key, value) {
  arr[key] = value
  return(arr)
}

# parameters --------------------------------------------------------------------------

#set working directory
#setwd(directory)
setwd("/home/igocer/rnaseq/nk_cells/")
directory = getwd()

# import script variables from json -------------------------------------------------------------
global_vars <- fromJSON(file="global_vars.json")

# log files
df_hisat = read.csv("hisat_log.csv")
df_htseq = read.csv("htseq_log.csv")

# load text file with experimental design and add appropriate header
col_data <- read.csv(global_vars$pheno_file, header = FALSE)
#names(col_data) <- c("sample", "file", "condition")
names(col_data) <- c("sample", "cytokine", "donor", "rna", "treatment")


#col_data_2factor = data.frame(animal = c(rep("pymt", 9), rep("wt", 6)), treatment = c(rep("cxhigh", 3), rep("cxlow", 3), rep("ctrl", 3), rep("cxlow", 3), rep("ctrl", 3)))
#col_data_2factor = data.frame(col_data, col_data_2factor)
#col_data = cbind.data.frame(col_data, subject = factor(c(rep(1,5), rep(2, 6), rep(3, 6))))

# graph palette
volPallete <- c("#B0B0B0", "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33")

# contrast workflow ------------------------------------------------

#count_data <- read.csv("nk_cells_hg38_counts.txt", sep="\t", skip = 1, row.names = "Geneid", 
#                 colClasses = c("character", rep("NULL", 5), rep("integer", length(col_data$sample))))
count_data <- read.csv("GSE77808_dataCount.csv", row.names = "gene_id")

#names(count_data) <- c(as.character(col_data$sample))
col_data$sample_id = factor(names(count_data))
col_data = col_data %>% dplyr::select(sample_id, everything())
col_data = col_data %>% dplyr::filter(rna == "Polysomal") %>% dplyr::mutate(condition = factor(paste0(cytokine, "_", treatment)))
count_data = count_data %>% select(col_data$sample)

dds = DESeqDataSetFromMatrix(count_data, col_data, design = ~ donor + condition)
#dds <- DESeqDataSetFromHTSeqCount(col_data, design = ~ subject + condition)

# remove genes with zero counts
dds <- dds[rowSums(counts(dds)) > 10,]
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

# set up reference condition
#dds$condition <- relevel(dds$condition, ref=global_vars$ref_level)

#levels = c("wt","wt.cx", "pymt", "pymt.cxlo", "pymt.cxhi")
dds$condition <- relevel(dds$condition, ref="IL2_A")
#dds$condition <- factor(dds$condition, levels = levels)

# add Genomic Ranges to calulate fpkm
txdb = makeTxDbFromGFF("/shared/genome/hg38/annotation/hg38.gtf", format="gtf")
txdb_exons = exonsBy(txdb, by="gene")
txdb_exons = txdb_exons[names(txdb_exons) %in% rownames(dds)]
rowRanges(dds) = txdb_exons

dds  = dds[ all(!seqnames(dds) %in% c("Y")), ]


# data transformations
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

#sample_expression = counts(dds2, normalized=T)


# DESeq
dds2 <- DESeq(dds)



# create results
group_levels = nlevels(col_data$condition)
level_comb = combn(levels(dds$condition), 2)

if(group_levels == 2) {
  res_array = results(dds2)
} else if(group_levels > 2) {
  res_array = apply(level_comb, 2, function(x) { results(dds2, contrast=c("condition", x[2], x[1]))})
} else {
  print("Error, only one group found.")
}

# 
# if(group_levels == 2) {
#   res_array = list(lfcShrink(dds2, coef=length(resultsNames(dds2)), type="apeglm"))
# } else if(group_levels > 2) {
#   #res_array = apply(level_comb, 2, function(x) { lfcShrink(dds2, contrast=c("condition", x[2], x[1]))})
#   res_array = apply(level_comb, 2, function(x) { lfcShrink(dds2, coef=length(resultsNames(dds2)), type="apeglm")})
# } else {
#   print("Error, only one group found.")
# }



res_array_names = apply(level_comb, 2, function(x) { paste(x[2], x[1], sep="_")})
names(res_array) = res_array_names

res_summary <- lapply(res_array, function(x) {capture.output(summary((x), global_vars$fdr_value))})

res_array = lapply(res_array, function(x) {
  rownames_to_column(data.frame(x))
})

# res_array = lapply(names(res_array), function(x) {
#   rownames_to_column(data.frame(res_array[[x]]))
# })

# order results by adjusted p-value
res_array_ordered <- lapply(res_array, function(x) {x[order(x$padj),]})





# summary log for numbers of upregulated and downregulated genes
logs = lapply(res_summary, summary_match)
summary_log = lapply(logs, function(x) {data.frame(matrix(unlist(x), nrow=1, byrow=T))})
summary_log = do.call("rbind", summary_log)
summary_header = c("logFC_up", "%logFC", "logFC_down", "%logFC_down", "outliers", "%outliers", "low_counts", "%lowcounts")
summary_index = row.names(summary_log)
if(dim(summary_log)[1] == 1) {
  summary_log = apply(summary_log, 2, as.character) %>% as.numeric()
  summary_log = data.frame(t(summary_log), stringsAsFactors = F)
} else {
  #TODO test on more data sets
  summary_log = apply(summary_log, 2, as.character) %>% apply(2, as.numeric)
}

row.names(summary_log) = summary_index
colnames(summary_log) = summary_header


# load biomart database
# assign gene names
# retrieve mouse gene list
write("getting gene names from to ensembl", stderr())
if(global_vars$animal == 'm') {
  mart <- useMart('ensembl', 'mmusculus_gene_ensembl')
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"), values=rownames(dds), mart= mart)
  global_vars$symbol = "mgi_symbol"
} else if(global_vars$animal == 'h') {
  mart <- useMart('ensembl', 'hsapiens_gene_ensembl')
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=rownames(dds), mart= mart)
  global_vars$symbol = "hgnc_symbol"
} else {
  print("Species not defined")
}

# add gene names to ensembl IDs
names(G_list) = c("ensembl_gene_id", "symbol")
nms <- data.frame(ensembl_gene_id = rownames(dds), stringsAsFactors = FALSE)
nms_join <- left_join(nms, G_list)
nms_join = nms_join[!duplicated(nms_join$ensembl_gene_id),]
na_gene_id <- nms_join$ensembl_gene_id[is.na(nms_join$symbol)]
nms_join$symbol[is.na(nms_join$symbol)] <- na_gene_id
empty_gene_id <- nms_join$ensembl_gene_id[nms_join$symbol == ""]
nms_join$symbol[nms_join$symbol == ""] <- empty_gene_id

# join ensembl ids to gene names and sort by adj p value
res_joined <- lapply(res_array_ordered, function(x) { data.frame(ensembl_gene_id = row.names(x), data.frame(x), stringsAsFactors = FALSE) })
res_joined <- lapply(res_joined, function(x) {inner_join(nms_join, x)})
res_joined <- lapply(res_joined, function(x) { x[!duplicated(x$ensembl_gene_id),]})
res_joined_ordered <- lapply(res_joined, function(x) {x[order(x$padj),]})

res_joined = res_array
res_joined_ordered = res_array_ordered
# filter results by adjusted p value
res_joined_ordered_filtered = lapply(res_joined_ordered, function(x) {x[!is.na(x$padj) & x$padj < global_vars$fdr_value, ]})

#seperate dir creation to different script
unique_dir_id = format(Sys.time(), "%y%m%d%H%M")
output_dir = paste0(directory, "/",unique_dir_id, "_", global_vars$file_prefix)
global_vars = update_json(global_vars, 'unique_dir', paste0(unique_dir_id, "_", global_vars$file_prefix))
global_vars = update_json(global_vars, 'output_dir', output_dir)
dir.create(output_dir)
setwd(output_dir)

#result_filenames = names(res_joined)
result_filenames = names(res_array)

write("exporting files", stderr())
# export files with all hits
sapply(result_filenames, function (x) write.table(res_joined_ordered[[x]], file=paste0(x, "_all.txt"), row.names = F, quote = F, sep = "\t" ))

# export filtered results
sapply(result_filenames, function (x) write.table(res_joined_ordered_filtered[[x]], file=paste0(x, "_filtered.txt"), row.names = F, quote = F, sep = "\t" ))

# export fpkm values
fpkm = cbind.data.frame(nms_join, fpkm(dds, robust = F))
write.table(fpkm, "rpkm.csv", quote = F, row.names = F, sep = "\t")

# export deseq corrected sample expression values
sample_expression = cbind.data.frame(nms_join, assays(dds2)[["mu"]])
write.table(sample_expression, "sample_expression.csv", quote = F, row.names = F, sep = "\t")

setwd(directory)
save.image(paste0(global_vars$file_prefix, ".RData"))
toJSON(global_vars) %>% writeLines("global_vars.json")

write("generating report", stderr())
#source(paste0(global_vars$exec_path,"/render_deseq.R"), local = TRUE)
source(paste0("~/programming/scripts/python/r_rnaseq","/render_deseq.R"), local = TRUE)
