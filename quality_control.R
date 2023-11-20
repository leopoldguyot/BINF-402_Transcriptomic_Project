library(Rqc)
library(tidyverse)

qcRes <- rqc(path = "data/mirnaseq_data/",
             pattern = "^ENC.*fastq.gz$",
             openBrowser = FALSE,
             outdir = "data_output/QC/")


statistics_per_cycle <- rqcCycleQualityBoxCalc(qcRes)
statistics_per_reads <- rqcReadQualityCalc(qcRes)
test <- statistics_per_cycle %>%
  filter(filename == "ENCFF327IXJ.fastq.gz")

ggplot(test, aes(x = cycle, ymin = ymin, ymax = ymax, lower = lower, middle = middle, upper = upper))+
  geom_boxplot(stat = "identity")

rqcCycleQualityPlot(qcRes) # Use this representation instead of the classic boxplot ?


############# EXPORTS ###############

ggsave("Figures/QC_plots/global_cycle_quality_boxplot.pdf",rqcCycleQualityBoxPlot(qcRes))
ggsave("Figures/QC_plots/global_reads_quality_plot.pdf",rqcReadQualityPlot(qcRes))
