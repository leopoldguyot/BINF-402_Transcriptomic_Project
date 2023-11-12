options(timeout = 5000)
# Human genome download :
download.file(
  'https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz',
  'data/genome/hg38.fa.gz')

# miRNA data download (4 replicates for 3 tissues):

sample_table <- read.csv(file = "data/sample_table_links.csv")

apply(sample_table, 1, function(row) {
  base <- basename(row["link"])
  download.file(row["link"], file.path("data", "mirnaseq_data", base))
})
