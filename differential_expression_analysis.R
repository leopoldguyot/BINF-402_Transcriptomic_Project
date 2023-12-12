############### IMPORTS ##################
library(Rsubread)
library(DESeq2)
library(tidyverse)
library(patchwork)
############### BODY ####################

files <- list.files(path = "data_output/mapped_reads",
                    pattern = "^v2_filtered.*\\.bam$",
                    full.names = TRUE)
fc <- featureCounts(files,
                    annot.inbuilt = "hg38") 

# assembling the count matrix with the data sets names 
# (to combine the count matrix with the metadata into a DESeq object)
sample_table <- read.csv(file = "data/sample_table_links.csv")
rownames(sample_table) <- sample_table[[3]]
sample_table <- sample_table[order(row.names(sample_table)), ]
colnames(fc$counts) <- sub("^v2_filtered_(.*?)\\.fastq\\.gz_reads_mapped\\.bam$", "\\1", colnames(fc$counts))

#creation of the DESeq object
desd <- DESeqDataSetFromMatrix(countData = fc$counts, 
                               colData = sample_table,
                               design = ~ tissue_name)

# plot count distribution (did not used in the final report)
as_tibble(assay(desd), rownames = "gene") %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:13) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ sample)

# plot of the sequence depth for each data set
SD <- tibble(sample = colnames(desd), SeqDepth = colSums(assay(desd))) %>%
  ggplot(aes(x = sample, y = SeqDepth)) +
  geom_col()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(plot = SD,"Figures/differential_analysis/depth.pdf", width = 4, height = 2.5)

desd <- desd[rowSums(assay(desd)) > 0, ] #removing rows containing only zeros

desd$tissue_name # gastrocnemius_medialis is the reference level
#change to have prostate
desd$tissue_name <- relevel(desd$tissue_name, ref = "pancreas_body")
desd <- DESeq(desd) #Running the differential analysis (normalization + stats test)

saveRDS(desd, file = "data_output/differential_results.rds")
rldesd <- rlogTransformation(desd) # optimise for PCA
#retrieve PC1 and PC2 coordinates 
pca_res <- plotPCA(rldesd, intgroup = "tissue_name", returnData = TRUE) 
# unsupervised clustering
kmeans_res <- kmeans(pca_res[,1:2], 3, 100) # between_SS / total_SS =  99.0 %

ggsave("Figures/differential_analysis/pca.pdf",
       plotPCA(rldesd, intgroup = "tissue_name"),
       height = 3, width = 6) #PCA plot

sizeFactors(desd) 
#extreme value for ENCFF486DIS other like ENCFF689CZI are quite low (information equivalent to seq depth)

# Ploting the estimate dispersion plot
pdf("Figures/differential_analysis/dispersion.pdf", width = 7, height = 4.5)
plotDispEsts(desd)
dev.off() 


#retrieving the results for each pair of comparison
res_pancreas_vs_gastro <- results(desd,
               name = "tissue_name_gastrocnemius_medialis_vs_pancreas_body")
res_pancreas_vs_prostate <- results(desd,
                                  name = "tissue_name_prostate_gland_vs_pancreas_body")
res_tbl_pancreas_vs_gastro <- as_tibble(res_pancreas_vs_gastro,rownames = "gene")
res_tbl_pancreas_vs_prostate <- as_tibble(res_pancreas_vs_prostate,rownames = "gene")


# Volcano plots
volcano_pan_vs_gas <- res_tbl_pancreas_vs_gastro %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  ggtitle("Pancreas compared to Gastrocnemius medialis")+
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))

volcano_pan_vs_pro <- res_tbl_pancreas_vs_prostate %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj),
             color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_colour_manual(values = c("gray", "red")) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = 1) +
  geom_vline(xintercept = -1) +
  ggtitle("Pancreas compared to Prostate")+
  theme(legend.position = "bottom",plot.title = element_text(hjust = 0.5))

ggsave(file = "Figures/differential_analysis/volcano.pdf",
       plot = volcano_pan_vs_pro + volcano_pan_vs_gas,
       height = 5.5, width = 10)
# Count number of features with padj under 0.05 and abs(log2FoldChange) above 1

n_vs_gastro <- res_tbl_pancreas_vs_gastro %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  nrow()

n_vs_prost <- res_tbl_pancreas_vs_prostate %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  nrow()
tot_vs_gastro <- nrow(res_pancreas_vs_gastro)
tot_vs_prost <- nrow(res_pancreas_vs_prostate)

# create a contengent table to use a chi test
table <- as.table(rbind(c(n_vs_gastro, tot_vs_gastro),c(n_vs_prost, tot_vs_prost)))
dimnames(table) <- list(pair = c("vs_gastro", "vs_prost"), dif_miRNA = c(TRUE, FALSE))

chisq.test(table) #p-value = 0.2705
