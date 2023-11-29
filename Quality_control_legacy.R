############## IMPORTS #################
library(ShortRead)
library(ggplot2)
library(stringr)
############# FUNCTIONS ###############

load_quality <- function(filename){
  fq <- readFastq(filename)
  return(quality(fq))
}

extract_read_means <- function(quality){
  alphabetScore(quality)/width(quality)
}

extract_cycle_means <- function(quality){
  tabel <- alphabetByCycle(quality)
  rownames(tabel) <- sapply(row.names(tabel),
                            function(rowname) utf8ToInt(rowname)-33)
  tabel <- as.data.frame(tabel)
  colnames(tabel) <- 1:ncol(tabel)
  Ncol <- colSums(tabel)
  for (row in 1:nrow(tabel)) {
    tabel[row,] <- tabel[row,] * as.integer(rownames(tabel[row,]))
  }
  i = 0
  means <- apply(tabel, 2, function(col) {
    i = i+1
    sum(col)/Ncol[i]})
  return(means)
}


read_mean_plot <- function(quality,name){
  means <- extract_read_means(quality)
  #Name
  ggsave(filename = paste("Figures/legacy_plots/","read_mean_", name,".pdf", sep = ""),
         ggplot(mapping = aes(x = means, after_stat(scaled))) +
            geom_density())
}

cycle_mean_plot <- function(quality, name){
  means <- extract_cycle_means(quality)
  means <- as.data.frame(means)
  ggsave(filename = paste("Figures/legacy_plots/","cycle_mean_", name,".pdf", sep = ""),
         ggplot(means, aes(x = 1:length(means), y = means)) + geom_line())
}
run_qc_plots <- function(fastqfiles){
  for (file in fastqfiles){
    base <- str_remove(basename(file),".fastq.gz")
    print(base)
    quality <- load_quality(file)
    read_mean_plot(quality,base)
    cycle_mean_plot(quality,base)
  }

}
############### BODY #######################

reads_files <- list.files(path = "data/mirnaseq_data", full.names = TRUE)

run_qc_plots(reads_files)
