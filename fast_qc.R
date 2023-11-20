library(fastqcr)

fastqc(fq.dir = "data/mirnaseq_data/", qc.dir = "data_output/fast_qc/")
?qc_report
qc = qc_read("data_output/fast_qc/ENCFF327IXJ_fastqc.zip")
test <- qc_plot(qc, "Per base sequence quality")
