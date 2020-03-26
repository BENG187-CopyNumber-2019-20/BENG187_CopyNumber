#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(DNAcopy)

cn <- read.table(args[1], header=T)
CNA.object <- CNA(genomdat = cn$adjusted_log_ratio, chrom = cn$chrom, maploc = cn$chr_start, data.type = 'logratio')
CNA.smoothed <- smooth.CNA(CNA.object)
segment <- segment(CNA.smoothed, verbose=0, min.width=2, undo.SD=3)
p.segment <- segments.p(segment)
write.table(p.segment[,2:10], file = paste(args[1],"segmentation", sep = "."), append = FALSE, quote = FALSE, row.names = FALSE)
q(save="no")
