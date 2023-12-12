############ IMPORTS ###############
library(Rsubread)
library(ShortRead)
library(QuasR)
library(tidyverse)
############ FUNCTIONS ############

# LEGACY
# quality_control <- function(filename){
#   print(paste("Quality control processing :", filename))
#   fq <- readFastq(filename)
#   t_fq = trimEnds(fq, "4")
#   output_file <- file.path("data_output/filtered_reads",
#                            paste("filtered_",basename(filename), sep = ""))
#   writeFastq(t_fq, output_file)
#   return(NULL)# will keep full result in memory (between iteration) if not added
# }

mapping <- function(filename){
  print(paste("mapping :", filename))
  basename <- basename(filename)
  align(index = "data_output/index/hg38_index", readfile1 = filename, type = "rna",
        output_file = paste("data_output/mapped_reads/",basename,"_reads_mapped.bam", sep = ""), nthreads = 4)
  return(NULL)# will keep full result in memory (between iteration) if not added
}


########### BODY ############

buildindex("data_output/index/hg38_index", "data/genome/hg38.fa.gz",
           indexSplit = TRUE)


# FILTERING
reads_files <- list.files(path = "data/mirnaseq_data", full.names = TRUE)

destination_filtered <- sapply(reads_files,
                               function(name) paste("data_output/filtered_reads/filtered_",
                                                    basename(name), sep = ""
                                                    )
                               )
# Cut the 5 first start bases (adapter)
preprocessReads(reads_files,
                outputFilename = destination_filtered,
                truncateStartBases = 5,
                truncateEndBases = 7,
                nBases = 0)

# Remove read with bad quality

for (fastqFile in destination_filtered) {
  print(paste("processing file : ",fastqFile ))
  # Set up FastqStreamer for the current file
  f <- FastqStreamer(fastqFile, readerBlockSize = 1000)

  while (length(fq <- yield(f))) {
    qPerBase = as(quality(fq), "matrix")
    qcount = rowMeans(qPerBase)

    writeFastq(fq[qcount > 20],
               paste("data_output/filtered_reads/","v2_",basename(fastqFile), sep = ""),
               mode = "a")
  }
}

reads_files_to_be_mapped <- c(list.files(path = "data_output/filtered_reads",
                                         full.names = TRUE),
                              reads_files)
lapply(reads_files_to_be_mapped, mapping) # applying the mapping on all files (36)

mapping_files <- list.files(path = "data_output/mapped_reads",
                    pattern = "*\\.bam$",
                    full.names = TRUE) # selection of .bam files (36)

# create a dataframe with information of read mapped for each data set
mapping_proportion <- propmapped(mapping_files) 
mapping_proportion["processed"] <- c(rep("No", 12),
                                       rep("Trimmed",12),
                                       rep("Trimmed+ReadsFiltered",12))
saveRDS(mapping_proportion, file = "data_output/mapped_reads/mapping_proportion.rds")

#plot mapping props
mapping_prop_plot <- ggplot(mapping_proportion, aes(x = processed, y = PropMapped, color = processed))+
  geom_boxplot() +
  theme(legend.position = "none")+
  ylab("Proportion mapped")+
  xlab("Type of processing")
ggsave(filename = "Figures/mapping_props.pdf",plot = mapping_prop_plot, height = 3, width = 4)
