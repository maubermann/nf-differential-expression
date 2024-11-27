#!/usr/bin/env Rscript

require("tximport")
require("DESeq2")
library("rtracklayer")

args <- commandArgs(trailingOnly = TRUE)

#quant_parent_dir <- "/Users/mo/code/assignments/BIOF/results/salmon"
#experiment_info_file <- "/Users/mo/code/assignments/BIOF/experiment_info.txt"
#gtf_file <- "/Users/mo/code/assignments/BIOF/Mus_musculus_c57bl6nj.C57BL_6NJ_v1.113.gtf.gz"

quant_parent_dir <- args[1]
experiment_info_file <- args[2]
gtf_file <- args[3]


#mostly following tutorial on bioconductor
experiment_info <- read.table(experiment_info_file, header = TRUE, stringsAsFactors = FALSE)
samples <- experiment_info$sample
conditions <- experiment_info$condition

quant_files <- file.path(quant_parent_dir, paste0(samples, "_quant"), "quant.sf")

#get transcript to gene mapping from gtx file 
gtf <- import(gtf_file)
tx2gene <- as.data.frame(gtf)[, c("transcript_id", "gene_id")]
tx2gene <- tx2gene[(!is.na(tx2gene$transcript_id) )& (!is.na(tx2gene$gene_id)),]
tx2gene <- tx2gene[!duplicated(tx2gene),]



txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

colData <- data.frame(run = samples, condition = factor(conditions),row.names = samples)
dds <- DESeqDataSetFromTximport(txi, colData, ~ condition)

dds <- DESeq(dds)
res <- results(dds)



res_with_pvalue <- res[!is.na(res$pvalue) ,]
res_with_padj <- res[!is.na(res$padj) ,]

#transcriptome wide significant:
res_padj_sig <-res_with_padj[res_with_padj$padj < 0.05,]


write.csv(data.frame(genes = rownames(res_padj_sig), padj = res_padj_sig$padj), file = "genes_padj_significant.csv")

pdf("MA_plot.pdf", width = 7, height = 5) 
plotMA(res, ylim = c(-2, 2))             
dev.off()   

