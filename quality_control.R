library(Rqc)
library(tidyverse)
library(patchwork)

stats_by_cycle_plot <- function(qcRes){
  groups <- unique(perFileInformation(qcRes)$group)
  for (group in groups) {
    qcRes.sub <- subsetByGroup(qcRes, group)
    print(rqcCycleAverageQualityPlot(qcRes.sub))
  }
}



qcRes <- rqc(path = "data/mirnaseq_data/",
             pattern = "^ENC.*fastq.gz$",
             openBrowser = FALSE,
             outdir = "data_output/QC/",
             file = "rqc_report_pre_filtered")

qcRes_filtered <- rqc(path = "data_output/filtered_reads/",
                      pattern = "^filtered.*fastq.gz$",
                      openBrowser = FALSE,
                      outdir = "data_output/QC/",
                      file = "rqc_report_post_filtered")

qcRes_filtered_v2 <- rqc(path = "data_output/filtered_reads/",
                      pattern = "^v2_filtered.*fastq.gz$",
                      openBrowser = FALSE,
                      outdir = "data_output/QC/",
                      file = "rqc_report_post_filtered_v2")
statistics_per_cycle <- rqcCycleQualityBoxCalc(qcRes)
statistics_per_reads <- rqcReadQualityCalc(qcRes)
statistics_per_cycle_filtered <- rqcCycleQualityBoxCalc(qcRes_filtered)
statistics_per_reads_filtered <- rqcReadQualityCalc(qcRes_filtered)

stats_by_cycle_plot(qcRes)
stats_by_cycle_plot(qcRes_filtered)
stats_by_cycle_plot(qcRes_filtered_v2)


ggplot(test, aes(x = cycle, ymin = ymin, ymax = ymax, lower = lower, middle = middle, upper = upper))+
  geom_boxplot(stat = "identity")

rqcCycleQualityPlot(qcRes) # Use this representation instead of the classic boxplot ?


############# EXPORTS ###############

ggsave("Figures/QC_plots/global_cycle_quality_boxplot.pdf",rqcCycleQualityBoxPlot(qcRes))
ggsave("Figures/QC_plots/global_reads_quality_plot.pdf",rqcReadQualityPlot(qcRes))
ggsave("Figures/QC_plots/global_cycle_quality_boxplot_filtered.pdf",rqcCycleQualityBoxPlot(qcRes_filtered))
ggsave("Figures/QC_plots/global_reads_quality_plot_filtered.pdf",rqcReadQualityPlot(qcRes_filtered))
