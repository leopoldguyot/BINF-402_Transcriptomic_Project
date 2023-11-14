############ IMPORTS ###############
library(Rsubread)
library(ShortRead)


############ FUNCTIONS ############

quality_control <- function(filename){
  print(paste("Quality control processing :", filename))
  fq <- readFastq(filename)
  t_fq = trimEnds(fq, "4")
  output_file <- file.path("data_output/filtered_reads",
                           paste("filtered_",basename(filename), sep = ""))
  writeFastq(t_fq, output_file)
  return(NULL)# will keep full result in memory (between iteration) if not added
}

mapping <- function(filename){
  print(paste("mapping :", filename))
  basename <- basename(filename)
  align(index = "data_output/index/hg38_index", readfile1 = filename, type = "rna",
        output_file = paste("data_output/mapped_reads/",basename,"_reads_mapped.bam", sep = ""), nthreads = 4)
  return(NULL)# will keep full result in memory (between iteration) if not added
}


########### MAIN CODE ############

buildindex("data_output/index/hg38_index", "data/genome/hg38.fa.gz",
           indexSplit = TRUE)

reads_files <- list.files(path = "data/mirnaseq_data", full.names = TRUE)
lapply(reads_files, quality_control)

filtered_reads_files <- list.files(path = "data_output/filtered_reads", full.names = TRUE)
lapply(filtered_reads_files, mapping)
