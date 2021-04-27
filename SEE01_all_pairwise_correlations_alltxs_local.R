projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))


library(tictoc)
library(data.table)

hapFreqs <- fread("Reformatted_SEE01_hap_freqs_v3.txt", header = T)
uniquePos <- fread("Unique_chemPos_v2.txt", header = F)

chemPOS <- data.table(gpIdx = uniquePos)
chemPOS[, gp := gsub(".*-", "", gpIdx.V1)][, gp := gsub(paste0("\\.", gpositions), "", gp)]
gpSplit <- split(chemPOS, chemPOS$gp)
newCols <- gsub("-R[0-9][0-9]$", "", colnames(hapFreqs))

job::job(corsTest = {
counter <- 0
listofCors <- lapply(gpSplit, function(x) {
    counter <<- counter + 1
    if(counter %% 500 == 0) {
        print(counter)
        flush.console()
    }
    findPos <- which(newCols %in% x$gpIdx.V1)
    corrDT <- hapFreqs[, c(findPos), with = FALSE]
    names(corrDT) <- gsub("\\.*", "", names(corrDT))
    corrDT <- na.omit(corrDT)
    cors <- cor(corrDT, use = "everything", method = "pearson")
    corsDF <- as.data.frame(cors, stringsAsFactors = FALSE)
    names(corsDF) <- unlist(lapply(names(corsDF), function(x) {paste(strsplit(x, "-")[[1]][c(1,3)], collapse = "-")}))
    corsDF
})
}, import = c(hapFreqs, gpSplit, newCols), packages = c("data.table")) 
allCors <- as.data.frame(do.call(rbind, corsTest$listofCors), stringsAsFactors = FALSE)

job::job({
write.table(allCors, file=paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/all_chems_all_reps_pearson_cors.txt"), col.names = T, quote = F, sep = "\t")
}, import = c(allCors, projectDir), packages = c("data.table"))