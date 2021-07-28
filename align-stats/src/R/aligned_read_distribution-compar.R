#
# Simple comparison of transcriptome aligned read length between two libraries
#
# TODO: Remove the hard-coded parts and make it read args
#

library("rio")

fileIn1 <- "results/lengths/hsa.dRNASeq.HeLa.polyA.1.alignedReadLen.transcriptome.txt"
fileIn2 <- "results/lengths/hsa.dRNASeq.HeLa.polyA.REL5.long.1.alignedReadLen.transcriptome.txt"
compar<-c("CTRL", "5P-Poly(A)")

a <- rio::import(fileIn1, format = "tsv") # CTRL
b <- rio::import(fileIn2, format = "tsv") # 5P-Poly(A)

plotName <- paste(compar[1], "vs", compar[2])

ofile <- paste0(strsplit(x = basename(fileIn1), split = ".aligned")[[1]][1], "vs", strsplit(x = basename(fileIn2), split = ".aligned")[[1]][1],
                ".alignedReadLen.transcriptome.pdf")
ofile<-paste(dirname(fileIn1), ofile, sep="/")
print(ofile)

# Get median length
med1 <- median(rep(a$V2, a$V1)) # Expand the frequency table
med2 <- median(rep(b$V2, b$V1)) # Expand the frequency table

a$V3<-a$V1/sum(a$V1)
b$V3<-b$V1/sum(b$V1)

pdf(ofile)
  plot(
    x = a$V2, y = a$V3, cex = 0.2, col = "red", type = "l",
    xlim = c(0, 5000),
    xlab = "Read length (capped at 5k)", ylab = "Fraction of reads",
    main = paste0("Aligned read length distribution\n", "Median length:\n", compar[1], ": ", med1, "; ", compar[2], ": ", med2)
  )
  lines(x = b$V2, y = b$V3, cex = 0.2, col = "blue")
  legend("topright",
    legend = c(compar[1], compar[2]),
    col = c("red", "blue"), lty = 1, cex = 1.2
  )
dev.off()
