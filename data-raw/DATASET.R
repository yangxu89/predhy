# data-raw/process.R
# Data import and processing pipeline


input_geno <- read.table(file="data-raw/geno.hmp.txt",sep="\t",header=TRUE)
input_geno1 <- read.csv(file="data-raw/geno.num.csv",row.names = 1)

hybrid_phe <- read.csv(file="data-raw/hybrid.phe.csv")


use_data(input_geno, input_geno1, hybrid_phe, overwrite = T)
