args <- commandArgs(TRUE)
hapDT <- as.character(args[1])
chem <- as.character(args[2])
chemPos <- as.character(args[3])
outdir <- as.character(args[4])

library(tictoc)
library(data.table)
library(pcaPP)

hapFreqs <- fread(hapDT, header = T)
uniquePos <- fread(chemPos, header = F)

hapFreqsFilt <- hapFreqs[, grep(chem, colnames(hapFreqs)), with = FALSE]

hapFreqsFilt2 <- melt(hapFreqsFilt)

chemPOS <- uniquePos[V1 %like% chem]
tic("Total time")
print(chem)
flush.console()
listofCors <- lapply(chemPOS$V1, function(x) {
    #print(x)
    if(length(x) > 1) {print(x)}
    newCols <- gsub("-R[0-9][0-9]$", "", colnames(hapFreqsFilt))
    findPos <- which(newCols == x)
    corrDF <- hapFreqsFilt[, c(findPos), with = FALSE]
    corrDF <- na.omit(corrDF)
    cors <- cor(corrDF, use = "everything", method = "spearman")
    #cors <- cor.fk(corrDF)
    corsDF <- as.data.frame(cors, stringsAsFactors = FALSE)
    names(corsDF) <- unlist(lapply(names(corsDF), function(x) {paste(strsplit(x, "-")[[1]][c(1,3)], collapse = "-")}))
    corsDF
})
allCors <- as.data.frame(do.call(rbind, listofCors), stringsAsFactors = FALSE)
toc()
write.table(allCors, file=paste0(getwd(), "/", outdir,"/",chem, "_spearman_cors.txt"), col.names = T, quote = F, sep = "\t")