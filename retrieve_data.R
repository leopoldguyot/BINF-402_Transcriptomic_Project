start_time <- Sys.time()
options(timeout = 5000)
# Human genome download :
download.file(
  'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz',
  'data/genome/hg38.fa.gz')

# mRNA data download (16 different tissues):

mrna_metadata <- read.csv(file = "data/mrnaseq_links.csv")


apply(mrna_metadata, 1, function(row) {
  download.file(row["link"], file.path("data", "mrnaseq_data", row["code"]))
})

end_time <- Sys.time()
end_time - start_time
