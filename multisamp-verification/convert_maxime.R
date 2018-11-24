args <- commandArgs(trailingOnly=TRUE)
infn <- args[1]
outfn <- args[2]

indata <- get(load(infn))[[1]]
M <- nrow(indata$omega)

outdata <- data.frame(
  "id" = paste('s', seq(0,M-1), sep=''),
  "name" = paste(indata$chromosome[,1], indata$position[,1], sep='_'),
  "var_reads" = apply(indata$mutCount, 1, paste, collapse=','),
  "total_reads" = apply(indata$mutCount + indata$WTCount, 1, paste, collapse=','),
  "var_read_prob" = apply(indata$omega, 1, paste, collapse=',')
)
write.table(outdata, file=outfn, quote=FALSE, sep='\t', row.names = FALSE)
