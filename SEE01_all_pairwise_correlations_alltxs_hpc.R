args <- commandArgs(TRUE)
hapDT <- as.character(args[1])
#chem <- as.character(args[2])
gpositions <- as.character(args[2])
chemPos <- as.character(args[3])
outdir <- as.character(args[4])

library(tictoc)
library(data.table)
library(pcaPP)

hapFreqs <- fread(hapDT, header = T)
uniquePos <- fread(chemPos, header = F)

hapFreqsFilt <- hapFreqs[, grep(paste0("\\.", gpositions), colnames(hapFreqs)), with = FALSE]
chemPOS <- data.table(gpIdx = uniquePos[V1 %like% paste0("\\.", gpositions)])
chemPOS[, gp := gsub(".*-", "", gpIdx.V1)][, gp := gsub(paste0("\\.", gpositions), "", gp)]
gpSplit <- split(chemPOS, chemPOS$gp)

tic("Total time")
print(gpositions)
flush.console()
listofCors <- lapply(gpSplit, function(x) {
    #print(x)
    if(length(x) > 1) {print(x)}
    newCols <- gsub("-R[0-9][0-9]$", "", colnames(hapFreqsFilt))
    findPos <- which(newCols %in% x$gpIdx.V1)
    corrDT <- hapFreqsFilt[, c(findPos), with = FALSE]
    names(corrDT) <- gsub(paste0("\\.", gpositions), "", names(corrDT))
    corrDT <- na.omit(corrDT)
    cors <- cor(corrDT, use = "everything", method = "spearman")
    #cors <- cor.fk(corrDF)
    corsDF <- as.data.frame(cors, stringsAsFactors = FALSE)
    names(corsDF) <- unlist(lapply(names(corsDF), function(x) {paste(strsplit(x, "-")[[1]][c(1,3)], collapse = "-")}))
    corsDF
})
allCors <- as.data.frame(do.call(rbind, listofCors), stringsAsFactors = FALSE)
toc()
write.table(allCors, file=paste0(getwd(), "/", outdir,"/","Idx_", gpositions, "_spearman_cors.txt"), col.names = T, quote = F, sep = "\t")