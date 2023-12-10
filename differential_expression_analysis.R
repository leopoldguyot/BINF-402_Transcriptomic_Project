############### IMPORTS ##################
library(Rsubread)
library(DESeq2)
library(tidyverse)
############### BODY ####################

files <- list.files(path = "data_output/mapped_reads",
                    pattern = "^v2_filtered.*\\.bam$",
                    full.names = TRUE)
fc <- featureCounts(files,
                    annot.inbuilt = "hg38")

sample_table <- read.csv(file = "data/sample_table_links.csv")
rownames(sample_table) <- sample_table[[3]]
sample_table <- sample_table[order(row.names(sample_table)), ]
colnames(fc$counts) <- sub("^v2_filtered_(.*?)\\.fastq\\.gz_reads_mapped\\.bam$", "\\1", colnames(fc$counts))

desd <- DESeqDataSetFromMatrix(countData = fc$counts,
                              colData = sample_table,
                              design = ~ tissue_name)

as_tibble(assay(desd), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

desd <- desd[rowSums(assay(desd)) > 0, ] #removing rows containing only zeros
desd <- DESeq(desd)
rldesd <- rlogTransformation(desd)
plotPCA(rldesd, intgroup = "tissue_name")
# kmeans unsupervised clustering
sizeFactors(desd) #extreme value for ENCFF486DIS other like ENCFF689CZI are quite low

plotDispEsts(desd)
desd$tissue_name # gastrocnemius_medialis is the reference level
res_pancreas_vs_gastro <- results(desd,
               name = "tissue_name_pancreas_body_vs_gastrocnemius_medialis")
res_prostate_vs_gastro <- results(desd,
                                  name = "tissue_name_prostate_gland_vs_gastrocnemius_medialis")
res_tbl_pancreas_vs_gastro <- as_tibble(res_pancreas_vs_gastro,rownames = "gene")
res_tbl_prostate_vs_gastro <- as_tibble(res_prostate_vs_gastro,rownames = "gene")

metadata(res_pancreas_vs_gastro)$filterThreshold # really high
as_tibble(metadata(res_pancreas_vs_gastro)$filterNumRej) %>%
  ggplot(aes(x = theta, y = numRej)) +
  geom_point() +
  geom_vline(xintercept = 0.872449,
             color = 'red')

hist(res_tbl_pancreas_vs_gastro$pvalue)
plotMA(res_pancreas_vs_gastro)

res_tbl_pancreas_vs_gastro %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  theme(legend.position = "bottom")


# Maybe plot count plot of higher expressed miRNA
