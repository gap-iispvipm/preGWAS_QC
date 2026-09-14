print("R script started")

args <- commandArgs(trailingOnly = TRUE)

imiss_file <- args[1]
lmiss_file <- args[2]

if (is.na(imiss_file) || imiss_file == "") {
  stop("ERROR: No .imiss file provided to R script")
}
if (is.na(lmiss_file) || lmiss_file == "") {
  stop("ERROR: No .lmiss file provided to R script")
}

indmiss <- read.table(file = imiss_file, header = TRUE)
snpmiss <- read.table(file = lmiss_file, header = TRUE)

pdf("histimiss.pdf")
hist(indmiss[, 6], main = "Histogram individual missingness", xlab = "Missing rate")
dev.off()

pdf("histlmiss.pdf")
hist(snpmiss[, 5], main = "Histogram SNP missingness", xlab = "Missing rate")
dev.off()

print("R script finished")




