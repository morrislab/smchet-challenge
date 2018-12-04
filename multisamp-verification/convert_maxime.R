args <- commandArgs(trailingOnly=TRUE)
infn <- args[1]
outfn <- args[2]

indata <- get(load(infn))[[1]]
M <- nrow(indata$dataset$chromosome)

outdata <- data.frame(
  "id" = paste('s', seq(0,M-1), sep=''),
  "name" = paste(indata$dataset$chromosome[,1], indata$dataset$position[,1], sep='_'),
  "var_reads" = apply(indata$dataset$mutCount, 1, paste, collapse=','),
  "total_reads" = apply(indata$dataset$mutCount + indata$dataset$WTCount, 1, paste, collapse=','),
  "copy_number_adjustment" = apply(indata$dataset$copyNumberAdjustment, 1, paste, collapse=','),
  "total_copy_number" = apply(indata$dataset$totalCopyNumber, 1, paste, collapse=','),
  "var_read_prob" = apply(round(indata$dataset$copyNumberAdjustment) / round(indata$dataset$totalCopyNumber), 1, paste, collapse=','),
  "cluster" = indata$clustering$best.node.assignments
)

write.table(outdata, file=outfn, quote=FALSE, sep='\t', row.names = FALSE)
