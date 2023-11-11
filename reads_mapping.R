library(Rsubread)

buildindex("data_output/index/hg38_index", "data/genome/hg38.fa.gz",
           indexSplit = TRUE)

mapping <- function(filename){
  basename <- basename(filename)
  align(index = "data_output/index/hg38_index", readfile1 = filename, type = "rna",
        output_file = paste("data_output/mapped_reads/",basename,"_reads_mapped.bam", sep = ""), nthreads = 4)
}

files <- list.files(path = "data/mirnaseq_data", full.names = TRUE)
lapply(files, mapping)
