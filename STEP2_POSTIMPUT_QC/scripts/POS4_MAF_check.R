print("R script started")

args <- commandArgs(trailingOnly = TRUE)

maf_file <- args[1]

if (is.na(maf_file) || maf_file == "") {
  stop("ERROR: No input file provided to R script")
}

maf_freq <- read.table(maf_file, header = TRUE, as.is = TRUE)

pdf("MAF_distribution.pdf")
hist(maf_freq[,5],
     main = "MAF distribution",
     xlab = "MAF")
dev.off()

print("R script finished")