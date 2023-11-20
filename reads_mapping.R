############ IMPORTS ###############
library(Rsubread)
library(ShortRead)
library(QuasR)

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

destination_filtered <- sapply(reads_files,
                               function(name) paste("data_output/filtered_reads/filtered_",
                                                    basename(name), sep = ""
                                                    )
                               )

preprocessReads(reads_files,
                outputFilename = destination_filtered,
                truncateStartBases = 5,
                nBases = 0)

reads_files_to_be_mapped <- c(list.files(path = "data_output/filtered_reads",
                                         full.names = TRUE),
                              reads_files)
lapply(reads_files_to_be_mapped, mapping)
