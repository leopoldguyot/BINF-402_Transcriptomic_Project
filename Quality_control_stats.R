############## IMPORTS #################
library(ShortRead)
library(tidyverse)
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


run_group_means <- function(fastqfiles, groupname, df) {
  result_df <- df

  for (file in fastqfiles) {
    base <- str_remove(basename(file), ".fastq.gz")
    print(base)
    quality <- load_quality(file)
    means <- extract_cycle_means(quality)

    if (length(means) < 50) {
      means <- c(rep(NA, 50 - length(means)), means)
    }

    result_df <- bind_rows(result_df,
                           data.frame(file_name = base,
                                      group = groupname,
                                      as.data.frame(t(means))))
  }

  return(result_df)
}

############### BODY #######################

reads_files <- list.files(path = "data/mirnaseq_data", full.names = TRUE)

run_qc_plots(reads_files)

groups_cycle_means_df <- data.frame("file_name"= character(),"group" = character())

vanilla <- list.files("data/mirnaseq_data/", full.names = TRUE)
groups_cycle_means_df <-  run_group_means(vanilla, "vanilla", groups_cycle_means_df)
trimmed <- list.files("data_output/filtered_reads", full.names = TRUE, pattern = "^filtered")
groups_cycle_means_df <- run_group_means(trimmed, "trimmed", groups_cycle_means_df)
trimmed_filtered <- list.files("data_output/filtered_reads", full.names = TRUE, pattern = "^v2")
groups_cycle_means_df <- run_group_means(trimmed_filtered, "trimmed + filtered", groups_cycle_means_df)

groups_cycle_means_df_pl <- groups_cycle_means_df %>% 
  pivot_longer(cols = 3:52,
               names_to = "cycle",
               values_to = "mean")

ggplot(groups_cycle_means_df_pl, aes(x = cycle, y = mean, color = group))+
  geom_line()