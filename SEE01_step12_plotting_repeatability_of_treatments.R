##############################################################################

# Project:     SEE01
# Author:      Robert Linder
# Date:        2021
# Title:       SEE01 paper 1 repeatability plots
# Description: Also include scatterplot of LOD scores vs correlation across all positions for all chemicals, fit a lm, then try doing each chemical separately to see if there are some chemicals where there is a correlation b/w LOD score and per-site correlation. May want to try to weight average correlations across replicates so that haplotypes that moved more are given more weight. Maybe can compute mean correlation across replicates for the top two or three moving haplotypes.

##############################################################################

# ============================================================================
# Load packages and sourced files
# Sourced files are kept in the default working directory	

library(tictoc)
library(data.table)
library(tidyverse)
library(scales)
library(ggplot2)
library(GGally)
library(ggpubr)
library(ggbeeswarm)
library(DescTools)
library(abind)
library(job)
library(ggdendro)
source('formatting/Haplotype_file_splitter.R')
source('utils/Writing_files_helper.R')
source('utils/Reading_files_helper.R')
source('calculating/Calculating_test_statistics.R')
source('formatting/Haplotype_file_reformatter.R')
source('seq/Position_offsetter.R')



# ============================================================================
# Set global options

defDir <- getwd()
#projectDir <- "/Users/robertlinder/Dropbox/Long_lab/DXQTL03/Primary_experiments/"
projectDir <- "/Users/robertlinder/Dropbox/Long_lab/SEE01/Primary_experiments/"
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

# ============================================================================
# Custom functions

rep.list <- function(object, repObject) {
    rep(list(object), length(repObject))
}

inflect <- function(avgDiffDT, colNames, threshold = 1){
    avgDiffDF <- as.data.frame(avgDiffDT, stringsAsFactors = FALSE)
    allCols <- lapply(colNames, function(col) {
        up   <- sapply(1:threshold, function(n) c(avgDiffDF[,col][-(seq(n))], rep(NA, n)))
        down <-  sapply(-1:-threshold, function(n) c(rep(NA,abs(n)), avgDiffDF[,col][-seq(length(avgDiffDF[,col]), length(avgDiffDF[,col]) - abs(n) + 1)]))
        a <- cbind(avgDiffDF[,col],up,down)
        minMaxList <- list(minima = which(apply(a, 1, min) == a[,1]), maxima = which(apply(a, 1, max) == a[,1]))
        names(minMaxList) <- c(paste0(col, "_minima"), paste0(col, "_maxima"))
        minMaxList
    } ) 
    allCols
}

add.sig.column <- function(avgDiffDT, filteredMinMaxList, LODcutoff) {
    avgDiffDT$significant <- 0
    minOnly <- unlist(filteredMinMaxList[[1]][1])
    maxOnly <- unlist(filteredMinMaxList[[1]][2])
    minRows <- unlist(minOnly)
    maxRows <- unlist(maxOnly)
    avgDiffDT$significant[minRows] <- 2
    avgDiffDT$significant[maxRows] <- 1
    avgDiffDT[LOD < LODcutoff, significant := 0]
    avgDiffDT
}

lod.support.intervals <- function(LODdrop) {
    function(DT) {
        chromSplit <- split(DT, DT$chr)
        LODintervals <- lapply(chromSplit, function(chromDT) {
            #print(chromDT$chr[1])
            #flush.console()
            chromDT[, row := .I]
            chromDT[, LOD2 := LOD - LODdrop]
            rightLod <- chromDT[chromDT, x.row-i.row, on = .(row > row, LOD <= LOD2), mult="first"]
            leftLod <- chromDT[chromDT, x.row-i.row, on = .(row < row, LOD <= LOD2), mult="last"]
            chromDT$rightBound <- rightLod
            chromDT$leftBound <- leftLod
            if(any(is.na(chromDT$rightBound))) {
                chromDT$rightBound[is.na(chromDT$rightBound)] <- nrow(chromDT) - chromDT$row[is.na(chromDT$rightBound)]}
            if(any(is.na(chromDT$leftBound))) {
                chromDT$leftBound[is.na(chromDT$leftBound)] <- 1 - chromDT$row[is.na(chromDT$leftBound)]}
            chromDT
        } )
        chromDTs <- do.call(rbind, LODintervals)
        chromDTs
    }
}

find.lod.dts <- function(DT, peakGP) {
    DT[, row := .I]
    counter <- 0
    findFlanks <- lapply(peakGP, function(peak) {
        counter <<- counter + 1
        peakRow <- DT[gp == peak, row]
        rowsDown <- DT[gp == peak, rightBound]
        rowsUp <- DT[gp == peak, leftBound]
        lodDT <- DT[row <= (peakRow + rowsDown) & row >= (peakRow + rowsUp)][, Idx := counter]
        lodDT
    } )
    flanksDT <- do.call(rbind, findFlanks)
    flanksDT
}

# ============================================================================
# Load data needed for downstream analyses

founderNames <- fread("Founder_names.txt", header = F)
founderNames <- founderNames[,V1]
nFounders <- length(founderNames)
indHapDiffsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_tx_ind_haps_sum_of_squared_diffs_tables/", analysisType = "haps_sq_diffs", samplePattern = "^SEE12B02")
indLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Hap_adjusted_LOD_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^SEE12B02")
indHapDevsDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Ind_rep_hap_dev_tables/", analysisType = "avg_hap_devs", samplePattern = "^SEE12B02")
repLODDTs <- reading_read.in.dts.as.list(projectRootDir = projectDir, folderPath = "Data_Analysis/Sequencing_analysis/Tables/Repeatability_tables/", analysisType = "hap_adjusted_LOD", samplePattern = "^18way")

treatKeyDT <- fread("treatment_key.txt", header = T)
offsets <- read.table("newoffsets.txt", header = T)
offsets$lines <- ceiling(offsets[,2]/50)
offsets$chrbytelengths <- rowSums(offsets[,c(2,3,5)])
offsets$chr <- 1:17
g_l <- c(0, cumsum(offsets$len))
ch.bounds <- c(0, g_l[1:17] + offsets[,2])
#correcting <- offsets$totaloffset[match(data$CHROM,offsets$chr)]
max.pos <- g_l[length(g_l)]
mid.ch <- diff(ch.bounds)/2
midpt.ch <- ch.bounds[2:18] - mid.ch
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
colorCodes = c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
hx = sapply(strsplit(colorCodes, ","), function(x) rgb(x[1], x[2], x[3], maxColorValue=255))
hx2 = hx[-c(5,21,24)]   
mycols = c(hx2,"#000000")

# ============================================================================
# Plot the average per haplotype deviation for each replicate genomewide

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Ind_reps_hap_dev_plots/")
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
popDT <- do.call(rbind, indHapDevsDTs)
maxDev <- max(popDT$sqrtAvgHapDifs)
txDTs <- split(popDT, popDT$chemWeek)


txLooping <- lapply(txDTs, function(DT) {
    tic('Total time')
    ylabel <- "Mean per haplotype deviation"
    plotCounter <- 0
    fileName <- DT$chemWeek[1]
    mainTitle <- fileName
    print(fileName)
    flush.console()
    DF <- as.data.frame(DT)
    png(paste0(filePath, "SEE01_", fileName, "_avg_hap_devs_", DT$Week[1], ".png"), width = 10, height = 10, units = 'in', res = 150)
    layout(matrix(1:length(unique(DF$id))))
    par(mar = c(0.7,2.5,0,6), mgp = c(1.5,0.7,0), oma = c(2,2.5,2,2))
    plotting <- lapply(unique(DF$id), function(x) {
        plotCounter <<- plotCounter + 1
        ymax <- 0.5
        ymin <- 0
        DT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(gp),max(gp)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
        DT[, box()]
        DT[id == x, points(gp, sqrtAvgHapDifs, col = "black", bty = "n", pch = 16, cex = 0.5)]
        chromSplit <- split(DT[id == x], DT[id == x, chr])
        smoothing <- lapply(chromSplit, function(ch) {
            ch[, lines(ksmooth(gp, sqrtAvgHapDifs, kernel = "normal", bandwidth = 100000), col = "red", lwd = 1)]
        } )
        lapply(seq(1, length(ch.bounds), 2), function(y) {
            rect(ch.bounds[y], ymin, ch.bounds[y+1], ymax, col = chromBarCol, border = NA, xpd = TRUE) 
        } )
        DT[, axis(2, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), las = 1, cex.axis = 0.75)]
        #DT[, mtext(text = yTickLab, side = 2, line = 0.5, cex = 0.6, las = 1)]
        print(x)
        flush.console()
        DT[, text(x = 11500000, y = ymax - ymax/8, labels = DT[id == x, unique(Replicate)], cex = 1)]
        if(plotCounter == 1) {
            DT[, title(main = mainTitle, cex.main = 1, line = 0.5, xpd = NA)]}
        if(length(unique(DF$id)) %% 2 == 1 & plotCounter == length(unique(DF$id))/2 + 0.5) {
            DT[, mtext(text = ylabel, side = 2, line = 3.5, cex = 0.75)]} else if(length(unique(DF$id)) %% 2 == 0 & plotCounter == round(length(unique(DF$id))/2, 0)) {DT[, text(grconvertX(-0.07, "npc", "user"), grconvertY(0.02, "npc", "user"), labels = ylabel, xpd = NA, cex = 1.25, srt = 90)]}
        if(plotCounter == length(unique(DF$id))) {DT[, axis(1, at=g_l[1:16] + offsets[[2]][1:16]/2, labels= as.roman(1:16), las = 1, cex.axis = 1)]}
    } )
    dev.off()
    toc()
} )

# ============================================================================
# Plot the average per haplotype deviation for each treatment genomewide
## Needs work

filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Ind_reps_hap_dev_plots/")
gray <- col2rgb('grey50')
chromBarCol <- rgb(gray[1],gray[2],gray[3], maxColorValue = 255, alpha = 100)
popDT <- do.call(rbind, indHapDevsDTs)
reps <- popDT[, .(reps = uniqueN(Replicate)), by = chemWeek]
avgDT <- popDT[, .(sqrtAvgHapDifs = mean(sqrtAvgHapDifs), sdHapDifs = sd(sqrtAvgHapDifs)), by = c("gp", "chemWeek")]
avgDTreps <- merge(avgDT, reps, by = "chemWeek")
avgDTreps[, c("upperCI", "lowerCI") := .(sqrtAvgHapDifs + 1.96*(sdHapDifs/sqrt(reps)), sqrtAvgHapDifs - 1.96*(sdHapDifs/sqrt(reps)))]
avgDTreps[, chr := rep(popDT$chr[c(1:11574)], uniqueN(chemWeek))]
maxDev <- max(avgDTreps$sqrtAvgHapDifs)
txDTs <- split(avgDTreps, avgDTreps$chemWeek)


ylabel <- "Mean per haplotype deviation"
plotCounter <- 0
DF <- as.data.frame(DT)
png(paste0(filePath, "SEE01_all_chems_avg_hap_devs_", DT$Week[1], ".png"), width = 10, height = 10, units = 'in', res = 150)
layout(matrix(1:length(txDTs)))
par(mar = c(0.7,2.5,0,6), mgp = c(1.5,0.7,0), oma = c(2,2.5,2,2))
plotting <- lapply(txDTs, function(DT) {
    plotCounter <<- plotCounter + 1
    tic()
    chemUsing <- gsub("18way_|_12", "", DT$chemWeek[1])
    chemUsing <- gsub("_", " ", chemUsing)
    print(chemUsing)
    flush.console()
    ymax <- 0.2
    ymin <- 0
    DT[, plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(min(gp),max(gp)), ylim = c(ymin, ymax), axes="FALSE", yaxs="i", xaxs="i")]
    DT[, box()]
    DT[, polygon(x = c(gp, rev(gp)), y = c(lowerCI, rev(upperCI)), col =  adjustcolor("dodgerblue", alpha.f = 0.25), border = NA)]
    DT[, points(gp, sqrtAvgHapDifs, col = "black", bty = "n", pch = 16, cex = 0.5)]
    chromSplit <- split(DT, DT$chr)
    smoothing <- lapply(chromSplit, function(ch) {
        ch[, lines(ksmooth(gp, sqrtAvgHapDifs, kernel = "normal", bandwidth = 100000), col = "red", lwd = 1)]
    } )
    lapply(seq(1, length(ch.bounds)-2, 2), function(y) {
        rect(ch.bounds[y], ymin, ch.bounds[y+1], ymax, col = chromBarCol, border = NA, xpd = TRUE) 
    } )
    DT[, text(x = 11300000, y = 0.19, labels = chemUsing, cex = 1)]
    DT[, axis(2, at=seq(0,0.2,0.05), labels=seq(0,0.2,0.05), las = 1, cex.axis = 0.75)]
    #DT[, mtext(text = yTickLab, side = 2, line = 0.5, cex = 0.6, las = 1)]
    #if(plotCounter == 1) {
        #DT[, title(main = mainTitle, cex.main = 1, line = 0.5, xpd = NA)]}
    if(length(txDTs) %% 2 == 1 & plotCounter == length(txDTs)/2 + 0.5) {
        DT[, mtext(text = ylabel, side = 2, line = 3.5, cex = 0.75)]} else if(length(txDTs) %% 2 == 0 & plotCounter == round(length(txDTs)/2, 0)) {DT[, text(grconvertX(-0.07, "npc", "user"), grconvertY(0.02, "npc", "user"), labels = ylabel, xpd = NA, cex = 1.25, srt = 90)]}
    if(plotCounter == length(txDTs)) {DT[, axis(1, at=g_l[1:16] + offsets[[2]][1:16]/2, labels= as.roman(1:16), las = 1, cex.axis = 1)]}
    toc()
} )
dev.off()

# ============================================================================
# Plot the average per-site correlations genome-wide for each chemical using Spearman's rho

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
files <- dir(path = directory, pattern = "_spearman_cors.txt$")
setwd(directory)
allCorDFs <- lapply(files, function(read) {
    read.table(read, header = T, sep = "\t")
} )
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

perSiteCor <- lapply(allCorDFs, function(corr) {
    identifier <- gsub("18way_|_12.*", "", rownames(corr)[1])
    print(identifier)
    flush.console()
    tic('total time')
    corr$gp <- substr(rownames(corr), regexpr("_12-", rownames(corr)) + 4, regexpr("-R", rownames(corr)) - 1)
    corrLst <- split(corr, corr$gp)
    rmveGP <- lapply(corrLst, function(x) {subset(x, select=-c(gp))})
    corrLstMats <- lapply(rmveGP, as.matrix)
    matMeans <- lapply(corrLstMats, function(x) mean(x[upper.tri(x)]))
    corrMeans <- unlist(matMeans)
    meanCorDT <- data.table(gp = as.numeric(names(corrMeans)), avgRho = corrMeans)[order(gp)][, Chemical := identifier]
    meanCorDT
} )


# ============================================================================
# Plot a pairwise correlation matrix of all replicates for each sample using Spearman correlation

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
files <- dir(path = directory, pattern = "_spearman_cors.txt$")
setwd(directory)
#allCorDFs <- lapply(files, function(read) {
    #read.table(read, header = T, sep = "\t")
#} )
allCorDF <- read.table("all_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

allCorDF <- allCorDF[-grep("fluconazole|nic", rownames(allCorDF)),]
allCorDF <- allCorDF[-grep("cadmium_chloride_12.*-R0[1-3]|YPD_12.*-R10|glacial_acetic_acid_12.*-R08|chlorpromazine_12.*-R12|chlorpromazine_12.*-R13", rownames(allCorDF)),]
allCorDF <- allCorDF[, -grep("fluconazole|nic|cadmium_chloride_12.R0[1-3]|YPD_12.R10|glacial_acetic_acid_12.R08|chlorpromazine_12.R12|chlorpromazine_12.R13", names(allCorDF))]

job::job(corJob = {
    chemReps <- unique(names(allCorDF))
    chemNames <- gsub("X18way_|_12.*", "", chemReps)
    names(chemReps) <- chemNames
    chemSplit <- split(chemReps, names(chemReps))
    chemLoop <- lapply(chemSplit, function(chem) {
        avgLoop <- lapply(chem, function(x) {
            tic()
            corrGrp <- allCorDF[,grep(x, names(allCorDF))]
            corrGrp <- as.data.frame(corrGrp)
            names(corrGrp) <- x
            rownames(corrGrp) <- rownames(allCorDF)
            corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
            corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
            #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
            corrGrpDT <- as.data.table(corrGrp)
            chemName <- gsub("X18way_|_12.*", "", names(corrGrp)[1])
            corrGrpDT <- corrGrpDT[grepl(chemName, id)]
            corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
            setnames(corrAvgs, c("id", x))
            toc()
            corrAvgs
        } )
        chemAvgs <- Reduce(function(x,y) merge(x = x, y = y, by = "id"), avgLoop)
        chemAvgs
    } )
}, import = c(allCorDF), packages = c("data.table", "tictoc") )

plotLooper <- lapply(corJob$chemLoop, function(x) {
    CC3 <- as.data.frame(x)
    rownames(CC3) <- x$id
    CC4 <- CC3[, -grep("^id", names(CC3))]
    names(CC4) <- gsub("X18way_", "", names(CC4))
    names(CC4) <- gsub("_12.", "-", names(CC4))
    chemTitle <- gsub("-.*", "", names(CC4))[1]
    names(CC4) <- gsub(".*-", "", names(CC4))
    rownames(CC4) <- gsub(".*-", "", rownames(CC4))
    data <- as.matrix(CC4)
    df <- reshape2::melt(data)
    gplot <- ggplot(data = df) + geom_tile(aes(x = factor(Var1), y = factor(Var2), fill = value)) + theme_bw(base_size = 8) + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkred") +  theme(text = element_text(size=10)) + theme(axis.title=element_blank()) + theme(legend.title=element_blank()) + ggtitle(chemTitle) + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank()) + theme(axis.text.x = element_text(angle = 45)) + coord_equal()
    gplot
} )
arrangePlot <- do.call(ggarrange, c(plotLooper, list(common.legend = TRUE, legend = "right")))
savePlot <- annotate_figure(arrangePlot, top = text_grob("SEE01_wk12_replicate_correlations", color = "black", size = 10, hjust = 0.5, vjust = -0.05))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/SEE01_wk12_haplotype_spearman_correlations.pdf"), savePlot, width = 8.5, height = 8.5, units = "in")

# ============================================================================
# Plot a dendrogram of all replicates for each sample using Spearman correlation

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
setwd(directory)
allCorDF <- read.table("all_chems_all_reps_spearman_cors.txt", header = T, sep = "\t")
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

allCorDF <- allCorDF[-grep("fluconazole|nic", rownames(allCorDF)),]
allCorDF <- allCorDF[-grep("cadmium_chloride_12.*-R0[1-3]|YPD_12.*-R10|glacial_acetic_acid_12.*-R08|chlorpromazine_12.*-R12|chlorpromazine_12.*-R13", rownames(allCorDF)),]
allCorDF <- allCorDF[, -grep("fluconazole|nic|cadmium_chloride_12.R0[1-3]|YPD_12.R10|glacial_acetic_acid_12.R08|chlorpromazine_12.R12|chlorpromazine_12.R13", names(allCorDF))]

job::job(corJob = {
    chemReps <- unique(names(allCorDF))
    corAvgs <- lapply(chemReps, function(x) {
        tic()
        corrGrp <- allCorDF[,grep(x, names(allCorDF))]
        corrGrp <- as.data.frame(corrGrp)
        names(corrGrp) <- x
        rownames(corrGrp) <- rownames(allCorDF)
        corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
        corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
        #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
        corrGrpDT <- as.data.table(corrGrp)
        corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
        setnames(corrAvgs, c("id", x))
        toc()
        corrAvgs
    } )
}, import = c(allCorDF), packages = c("data.table", "tictoc") )


filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/")
CC2 <- do.call(cbind, corJob$corAvgs)
CC3 <- as.data.frame(CC2)
rownames(CC3) <- CC2$id
CC4 <- CC3[, -grep("^id", names(CC3))]
names(CC4) <- gsub("X18way_", "", names(CC4))
names(CC4) <- gsub("_12.", "-", names(CC4))
out <- hclust(as.dist(1-CC4^2))

fileName <- paste0(filePath, "all_reps_Spearman_cors_v2.pdf")
mainTitle <- paste0("Spearman_correlations_v2")
pdf(fileName, width = 10, height = 10)
par(mar = c(7.5,3,1.5,0.5), mgp = c(2,0.7,0), oma = c(0,0,0,0))
par(cex = 0.75)
plot(as.dendrogram(out), xlab = "", ylab = "", main = "", sub = "", axes = FALSE, ylim = c(0, 1))
rect.hclust(out, h = 0.85, border = "red")
par(cex = 1)
title(main = mainTitle, ylab = "height")
axis(2)
dev.off()

# ============================================================================
# Plotting avg heterozygosity vs avg correlation per treatment (not significant)

hapDT <- fread(paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Summary_tables/SEE01_Reps_Using_Hets.txt"), header = T)
hapDT$treatment <- gsub("^18way_|_([0-1]?[0-9])$", "", hapDT$chemical)
avgHets <- hapDT[, mean(het), by = chemical]
setnames(avgHets, c("chemWeek", "avgHet"))
avgHets$id <- unlist(lapply(avgHets$chemWeek, function(x) {
    splitting <- strsplit(x,"_")[[1]]
    if(length(splitting) == 3) {return(splitting[2])} else{
        return(paste(splitting[2:3], collapse = "_"))}
} ) )


avgChemCorDFs <- lapply(corJob$chemLoop, function(x) {
    x <- as.data.frame(x)
    rownames(x) <- x$id
    names(x)[2:length(names(x))] <- x$id
    chemName <- gsub("-.*", "", names(x)[2])
    x$id <- NULL
    avgCor <- mean(x[upper.tri(x)])
    avgCorDF <- data.frame(chemical = chemName, avgCor = avgCor, stringsAsFactors = FALSE)
    avgCorDF
} )

avgCorsDF <- do.call(rbind, avgChemCorDFs)
avgCorDT <- as.data.table(avgCorsDF)
avgHets[, chemical := gsub("18way_|_12", "", chemWeek)][, chemWeek := NULL]
avgCorHetDF <- merge(avgHets, avgCorDT, by = "chemical")
avgCorHetDF[chemical %like% "YPD", c("chemical", "id") := .('ypd', 'ypd')]
avgCorHetDF <- avgCorHetDF[order(avgCorHetDF$chemical),]
avgCorHetDF <- na.omit(avgCorHetDF)

#avgCorHetDF$relHet <- avgCorHetDF$avgHet/hapDT$Heterozygosity[hapDT$Treatment == "BAS01"]

corrlm <- lm(avgCor ~ avgHet, data = avgCorHetDF)
modsum <- summary(corrlm)
r2 <- modsum$adj.r.squared
mylabel <- bquote(italic(R)^2 == .(format(r2, digits = 2)))
filename <- paste0("SEE01_wk12_avgHetvsCor.png")
filePath <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/")
png(paste0(filePath, filename), width = 10, height = 10, units = 'in', res = 150)
par(mar = c(1,1.5,0,0), mgp = c(1.5,0.7,0), oma = c(2,2,1,1))
plot(0,0,type="n", ann=FALSE, xaxt='n', yaxt='n', xlim=c(0, 1), ylim = c(0, 1), axes="FALSE", yaxs="i", xaxs="i")
counter <- 0
lapply(avgCorHetDF$id, function(p) {
    counter <<- counter + 1
    points(x = avgCorHetDF$avgHet[avgCorHetDF$id == p], y = avgCorHetDF$avgCor[avgCorHetDF$id == p], col = mycols[counter], pch = 20, bty = "n", cex = 1)
} )
abline(corrlm)
mtext(text = "average within-treatment correlation", side = 2, line = 2, cex = 0.85)
mtext(text = "average within-treatment heterozygosity", side = 1, line = 2, cex = 0.85)
text(x = 0.55, y = 0.55, labels = mylabel)
axis(2, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), las = 1, cex.axis = 0.75)
axis(1, at=c(0, 0.25, 0.5, 0.75, 1), labels=c(0, 0.25, 0.5, 0.75, 1), las = 1, cex.axis = 0.75)
legend(0.05, 0.95, legend = c(avgCorHetDF$id), col = mycols[1:nrow(avgCorHetDF)], pch = 20, xpd = NA, cex = 1, bty = "n")
dev.off()

# ============================================================================
# Find the most significant peaks, the middle peak, and a smaller peak for each drug. Find the 2.5 LOD support interval surrounding the peak LOD scores to subset these positions for plotting.

popDT <- do.call(rbind, indHapDiffsDTs)
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
repShapes <- data.frame(reps = c(paste0("R0", 1:9), paste0("R1", 0:6)), shapes = c(0:6, 8:12, 15:18))
    
colNames <- rep(list('LOD'), length(indLODDTs))
#threshold <- rep(list(25), length(indLODDTs))
threshold <- rep(list(50), length(indLODDTs))
allPeaks <- Map(inflect, indLODDTs, colNames, threshold)
peakCol <- Map(add.sig.column, indLODDTs, allPeaks, rep.list(5, indLODDTs))
top <- lapply(peakCol, function(x) {
    peaks <- x[significant == 1]
    quants <- quantile(peaks$LOD, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1))
    quantsIdx <- sapply(quants, function(y) which.min(abs(peaks$LOD - y)))
    quantsDT <- peaks[quantsIdx][order(LOD, decreasing = TRUE)][, Idx := 1:.N]
    quantsDT
} )
peakPos <- lapply(top, function(x) {x$gp})
lodInts <- lod.support.intervals(2.5)
addLodInts <- Map(lodInts, peakCol)
lodWins <- Map(find.lod.dts, addLodInts, peakPos)
allLodWins <- do.call(rbind, lodWins)
orderedWins <- split(allLodWins, allLodWins$Idx)

chemIndexer <- lapply(orderedWins, function(idx) {
    popDTsub <- merge(popDT, idx[, c("gp", "Chemical", "Idx")], by = c("gp", "Chemical"))[order(Chemical, Replicate)]
    popDTsub
} )

# ============================================================================
# Find the two most increased haplotypes within each window.

idxLooper <- lapply(chemIndexer, function(idx)  {
    freqDifs <- idx[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- as.data.frame(idx[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    maxIdx <- apply(freqDifs, 1, which.max)
    idxDF <- data.frame(row = 1:length(maxIdx), col = maxIdx)
    maxHap <- collHaps[as.matrix(idxDF)]
    idx[, c("maxHap", "maxChange") := .(maxHap, maxVals)]
    nxtMax <- apply(freqDifs, 1, function(x) rev(sort(x))[2])
    nxtMaxIdx <- apply(freqDifs, 1, function(x) which(x == rev(sort(x))[2])[1])
    idxDF2 <- data.frame(row = 1:length(nxtMaxIdx), col = unlist(nxtMaxIdx))
    maxHap2 <- collHaps[as.matrix(idxDF2)]
    idx[, c("maxHap2", "maxChange2") := .(maxHap2, nxtMax)]
    idx <- idx[, c("chr", "pos", "gp", "Chemical", "Replicate", "maxHap", "maxChange", "maxHap2", "maxChange2", "Idx")]
    idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    idxMelt1 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical", "Replicate", "Idx"), measure.vars = c("maxChange", "maxChange2"))
    idxMelt2 <- melt(idx, id.vars = c("chr", "pos", "gp", "Chemical", "Replicate", "Idx"), measure.vars = c("maxHap", "maxHap2"))
    setnames(idxMelt2, old = c("variable", "value"), new = c("hapType", "haplotype"))
    idxMrge <- idxMelt1[, c("hapType", "haplotype") := .(idxMelt2$hapType, idxMelt2$haplotype)]
    setnames(idxMrge, old = c("variable", "value"), new = c("changeType", "frequencyChange"))
    hapColors <- match(idxMrge$haplotype, topHaps)
    hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
    Colors <- mycols[hapColors]
    reps <- match(idxMrge$Replicate, repShapes$reps)
    Shapes <- repShapes$shapes[reps]
    idxMrge[, c("color.codes", "Shape") := .(Colors, Shapes)][, "pos(kb)" := pos/1000][, "pos" := NULL]
    idxMrge
})

# ============================================================================
# Find the sd of the most increased haplotypes at each position genomewide.

popDT <- do.call(rbind, indHapDiffsDTs)
freqDifs <- popDT[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
collHaps <- as.data.frame(popDT[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
maxIdx <- apply(freqDifs, 1, which.max)
idxDF <- data.frame(row = 1:length(maxIdx), col = maxIdx)
maxHap <- collHaps[as.matrix(idxDF)]
popDT[, c("maxHap", "maxChange") := .(maxHap, maxVals)]
popDT2 <- popDT[, c("chr", "pos", "gp", "Chemical", "Replicate", "maxHap", "maxChange")]
popDT2[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
maxHapsDT <- popDT2[, .(maxHaps = uniqueN(maxHap), avgFreqChange = mean(maxChange)), by = c("Chemical", "gp")] ### can try seeing if higher LOD translates to same haplotype increasing most in frequency genomewide- score is number of unique maxHaps over number of replicates

### now calculate the sd for positions where only a single haplotype is the most frequent

oneHap <- maxHapsDT[maxHaps == 1, gp]
popDT3 <- popDT2[gp %in% oneHap]
sdDT <- popDT3[, .(sdMax = sd(maxChange)), by = c("Chemical", "gp")]

# ============================================================================
# Find the Spearman correlation coefficient for the most significant peaks, the middle peak, and a smaller peak for each drug. Make a dataframe that has the Chemical and averaged correlation coefficient over that interval. 

directory <- paste0(projectDir, "Data_Analysis/Sequencing_analysis/Tables/Correlation_tables/")
files <- dir(path = directory, pattern = "_spearman_cors.txt$")
setwd(directory)
allCorDFs <- lapply(files, function(read) {
    read.table(read, header = T, sep = "\t")
} )
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

allIdx <- do.call(rbind, idxLooper)
job::job(allChems = {
    byChems <- lapply(allCorDFs, function(corr) {
    identifier <- strsplit(rownames(corr)[1], "_")[[1]][2]
    print(identifier)
    flush.console()
    tic('total time')
    corr$gp <- substr(rownames(corr), regexpr("_12-", rownames(corr)) + 4, regexpr("-R", rownames(corr)) - 1)
    corrLst <- split(corr, corr$gp)
    rmveGP <- lapply(corrLst, function(x) {subset(x, select=-c(gp))})
    corrLstMats <- lapply(rmveGP, as.matrix)
    matMeans <- unlist(lapply(corrLstMats, function(x) mean(x[upper.tri(x)])))
    meanCorDT <- data.table(gp = as.numeric(names(matMeans)), avgRho = matMeans)[order(gp)]
    chemIdx <- allIdx[grepl(identifier, Chemical)]
    idxSplit <- split(chemIdx, chemIdx$Idx)
    idxLoop <- lapply(idxSplit, function(idx) {
            idxMeanCorDT <- meanCorDT[grepl(paste(as.character(unique(idx$gp)), collapse = "|"), gp)]
            meanCor <- idxMeanCorDT[, mean(avgRho)]
            meanCorDT2 <- data.table(Chemical = idx$Chemical[1], avgRho = meanCor, idx = idx$Idx[1])
            meanCorDT2
    } )
    toc()
    idxBind <- do.call(rbind, idxLoop)
    idxBind
} )
}, import = c(allIdx, allCorDFs), packages = c("data.table", "tictoc") )

chemCorsDT <- do.call(rbind, allChems$byChems)     

# ============================================================================
# Plot the frequency of the most changed haplotype (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots. Have these plotted as a separate sub-panel for each chemical using facet. Show the most significant peak, a middle peak, and a smaller peak that is likely still real to show how repeatability changes as a function of the LOD score, with each type of peak a separate panel (A-C). Include as text the standard deviation.

counter <- 0
plotLooper <- lapply(idxLooper, function(idce) {
    counter <<- counter + 1
    chrDT <- idce[, .(chr = chr[1]), by = "Chemical"]
    idce <- idce[, hapGrps := haplotype][!hapGrps %in% topHaps, hapGrps := "other"]
    gplot <- ggplot(data=idce, aes(`pos(kb)`, frequencyChange, group = interaction(Replicate, changeType), colour = hapGrps, shape = Replicate)) + facet_wrap(~as.factor(Chemical), scales = "free_x") + coord_cartesian(ylim = c(-1, 1)) + xlab("position (kb)") + ylab("haplotype frequency change") + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(idce$color.codes, idce$hapGrps)) + scale_shape_manual(values = setNames(idce$Shape, idce$Replicate)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + geom_text(data = chrDT, aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.25, label = paste0('chr', as.roman(chr))), inherit.aes = FALSE) + labs(tag = LETTERS[counter])
    ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Hap_change_repeatability_by_peak_size_", counter, ".pdf"), gplot, width = 8.5, height = 9, units = "in")
} )
#savePlot <- do.call(ggarrange, c(plotLooper, list(common.legend = TRUE, legend = "right")))

# ============================================================================
# Data wrangle frequency change of all haplotypes (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots.

idxLooperAll <- lapply(chemIndexer, function(idx)  {
    chemSplit <- split(idx, idx$Chemical)
    chemLooper <- lapply(chemSplit, function(chem) {
        freqDifs <- chem[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
        collHaps <- chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)]
        allHaps <- cbind(chem$chr, chem$pos, chem$gp, chem$Chemical, chem$Replicate, chem$Idx, collHaps, freqDifs)
        freqsDF <- as.data.frame(freqDifs)
        newFreqs <- lapply(freqsDF, function(y) gsub("^-\\d.*|^\\d.*", 1, y))
        newFreqsDF <- as.data.frame(do.call(cbind, newFreqs))
        newFreqsDF2 <- sapply(newFreqsDF, as.numeric)
        allHaps$reps <- rowSums(newFreqsDF2, na.rm = T)
        positions <- rep(unlist(allHaps[,2]), allHaps$reps) 
        gpositions <- rep(unlist(allHaps[,3]), allHaps$reps)
        replicates <- sort(rep(unlist(allHaps[,5]), allHaps$reps))
        allfreqs <- as.numeric(unlist(t(freqDifs)))
        allfreqs2 <- na.omit(allfreqs)
        hapVec <- as.character(unlist(t(collHaps)))
        hapVec2 <- na.omit(hapVec)
        idxDT <- data.table(chr = allHaps$V1[1], pos = positions, gp = gpositions, Chemical = allHaps$V4[1], Replicate = replicates, Idx = allHaps$V6[1], haps = hapVec2, freqDifs = allfreqs2)
        idxDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
        hapColors <- match(idxDT$haps, topHaps)
        hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
        Colors <- mycols[hapColors]
        reps <- match(idxDT$Replicate, repShapes$reps)
        Shapes <- repShapes$shapes[reps]
        idxDT[, c("color.codes", "Shape") := .(Colors, Shapes)][, "pos(kb)" := pos/1000][, "pos" := NULL]
        idxDT
    } )
    idxDTs <- do.call(rbind, chemLooper)
    idxDTs
} )


# ============================================================================
# Create a data table of only adjacent positions with different combos to plot as line segments below.

 diffDTs <- lapply(idxLooperAll, function(idx) {
     chemSplit <- split(idx, idx$Chemical)
     chemLooper <- lapply(chemSplit, function(chem) {
         repSplit <- split(chem, chem$Replicate)
         repLooper <- lapply(repSplit, function(repl) {
             posLooper <- lapply(unique(repl$gp)[-length(unique(repl$gp))], function(position) {
                 hapCombos1 <- repl$haps[repl$gp == position]
                 hapCombos2 <- repl$haps[repl$gp == (position + 1000)]
                 if(any(!hapCombos2 %in% hapCombos1)) {
                     diffHapDT <- repl[gp %in% c(position, position + 1000)]
                     diffHaps1 <- hapCombos1[which(!hapCombos1 %in% hapCombos2)]
                     diffHaps2 <- hapCombos2[which(!hapCombos2 %in% hapCombos1)]
                     diffHapsQuery <- gsub("(?<=\\d)([A-B])", "|\\1\\2", diffHaps1, perl=T)
                     matching <- lapply(diffHapsQuery, function(x) {grep(x, diffHaps2)})
                     hapDTs <- lapply(1:length(matching), function(x) {
                        DT <-  data.table(haps1 = diffHaps1[[x]], haps2 = diffHaps2[matching[[x]]])
                        DT
                     } )
                     hapDT <- do.call(rbind, hapDTs)
                     hapDT[, c("chr", "Replicate", "Chemical", "Idx") := .(diffHapDT$chr[1], diffHapDT$Replicate[1], diffHapDT$Chemical[1], diffHapDT$Idx[1])]
                     pos1 <- unlist(lapply(hapDT$haps1, function(x) {diffHapDT$`pos(kb)`[x == diffHapDT$haps]}))
                     pos2 <- unlist(lapply(hapDT$haps2, function(x) {diffHapDT$`pos(kb)`[x == diffHapDT$haps]}))
                     hapDT[, c("pos(kb)1", "pos(kb)2") := .(pos1, pos2)]
                     freq1 <- unlist(lapply(hapDT$haps1, function(x) {diffHapDT$freqDifs[x == diffHapDT$haps]}))
                     freq2 <- unlist(lapply(hapDT$haps2, function(x) {diffHapDT$freqDifs[x == diffHapDT$haps]}))
                     hapDT[, c("freq1", "freq2") := .(freq1, freq2)]
                     hapDT
                 }
            } )
            allPos <- do.call(rbind, posLooper)
            allPos
         } )
         allReps <- do.call(rbind, repLooper)
         allReps
    } )
    allChems <- do.call(rbind, chemLooper)    
    allChems
} )


# ============================================================================
# Plot the frequency change of all haplotypes (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots.

### Still a bug when grouping (look at sodium chloride, index 1)
counter <- 0
plotLooper <- lapply(idxLooperAll, function(indexing) {
    counter <<- counter + 1
    diffDT <- diffDTs[[counter]]
    print(counter)
    flush.console()
    chrDT <- indexing[, .(chr = chr[1]), by = "Chemical"]
    indexing <- indexing[, hapGrps := haps][!hapGrps %in% topHaps, hapGrps := "other"][order(Replicate, haps, gp)][, c("difCol") := .(c(1000, diff(gp))), by = .(Replicate, haps)][, row := .I]
    indexing[, grpCol := paste0(haps, "_", difCol)]
    indexing[, grp := with(rle(as.vector(grpCol)), rep(seq_along(lengths), lengths)), by = .(Replicate)]
    reGrp <- indexing[difCol != 1000, row]
    nextGrp <- indexing[row %in% (reGrp + 1), difCol]
    nextRep <- indexing[row %in% (reGrp + 1), Replicate]
    grpRight <- which(nextGrp == 1000) 
    grping <- indexing[row %in% reGrp[grpRight], grp := indexing[row %in% (reGrp[grpRight] + 1), grp]]
    repSplit <- split(indexing, indexing$Replicate)
    repLoop <- lapply(repSplit, function(repl) {
        rleGrps <- rle(repl$grp)
        rleGrps$values <- 1:length(rleGrps$values)
        repl$grp <- rep(rleGrps$values, rleGrps$lengths)
        repl
    } )
    allRepsDT <- do.call(rbind, repLoop)
    correctGrps <- rle(allRepsDT$grp)
    correctGrps$values <- 1:length(correctGrps$values)
    allRepsDT$grp <- rep(correctGrps$values, correctGrps$lengths)
    corDT <- chemCorsDT[idx %in% indexing$Idx[1]]
    corDT[, rhoLabel := sprintf("italic(R) == %.2f", avgRho)]
    gplot <- ggplot(data=allRepsDT, aes(`pos(kb)`, freqDifs, colour = hapGrps, shape = Replicate, group = interaction(grp, Replicate))) + facet_wrap(~as.factor(Chemical), scales = "free_x") + coord_cartesian(ylim = c(-1, 1.1)) + xlab("position (kb)") + ylab("haplotype frequency change") + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(allRepsDT$color.codes, allRepsDT$hapGrps)) + scale_shape_manual(values = setNames(allRepsDT$Shape, allRepsDT$Replicate)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + geom_text(data = chrDT, aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.25, label = paste0('chr', as.roman(chr))), inherit.aes = FALSE) + geom_text(data = corDT, aes(x = -Inf, y = Inf, hjust = -0.15, vjust = 1.25, label = rhoLabel), parse = TRUE, inherit.aes = FALSE)
    ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Hap_change_repeatability_by_peak_size_Idx_", counter, ".pdf"), gplot, width = 8.5, height = 9, units = "in")
} )

# ============================================================================
# Plot the frequency change of all haplotypes (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots, each sample plotted separately.

counter <- 0
plotLooper <- lapply(idxLooperAll, function(indexing) {
    counter <<- counter + 1
    diffDT <- diffDTs[[counter]]
    print(counter)
    flush.console()
    chemSplit <- split(indexing, indexing$Chemical)
    chemLooper <- lapply(chemSplit, function(chem){
        print(chem$Chemical[1])
        flush.console()
        spliceDT <- diffDT[Chemical == chem$Chemical[1]]
        chem <- chem[, hapGrps := haps][!hapGrps %in% topHaps, hapGrps := "other"][order(Replicate, haps, gp)][, c("difCol") := .(c(1000, diff(gp))), by = .(Replicate, haps)][, row := .I]
        chem[, grpCol := paste0(haps, "_", difCol)]
        chem[, grp := with(rle(as.vector(grpCol)), rep(seq_along(lengths), lengths)), by = .(Replicate)]
        reGrp <- chem[difCol != 1000, row]
        nextGrp <- chem[row %in% (reGrp + 1), difCol]
        nextRep <- chem[row %in% (reGrp + 1), Replicate]
        grpRight <- which(nextGrp == 1000) 
        grping <- chem[row %in% reGrp[grpRight], grp := chem[row %in% (reGrp[grpRight] + 1), grp]]
        repSplit <- split(chem, chem$Replicate)
        repLoop <- lapply(repSplit, function(repl) {
            rleGrps <- rle(repl$grp)
            rleGrps$values <- 1:length(rleGrps$values)
            repl$grp <- rep(rleGrps$values, rleGrps$lengths)
            repl
        } )
        allRepsDT <- do.call(rbind, repLoop)
        correctGrps <- rle(allRepsDT$grp)
        correctGrps$values <- 1:length(correctGrps$values)
        allRepsDT$grp <- rep(correctGrps$values, correctGrps$lengths)
        corDT <- chemCorsDT[idx %in% chem$Idx[1] & Chemical %in% chem$Chemical[1]]
        corDT[, rhoLabel := sprintf("italic(R) == %.2f", avgRho)]
        gplot <- ggplot(data=allRepsDT, aes(`pos(kb)`, freqDifs, colour = hapGrps, shape = Replicate, group = interaction(grp,Replicate))) + coord_cartesian(ylim = c(-1, 1.1)) + xlab(paste0("chr", as.roman(allRepsDT$chr[1]), " position (kb)")) + ylab("haplotype frequency change") + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(allRepsDT$color.codes, allRepsDT$hapGrps)) + scale_shape_manual(values = setNames(allRepsDT$Shape, allRepsDT$Replicate)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + geom_text(data = corDT, aes(x = -Inf, y = Inf, hjust = -0.15, vjust = 1.25, label = rhoLabel), parse = TRUE, inherit.aes = FALSE) + geom_segment(data = spliceDT, aes(x = `pos(kb)1`, xend = `pos(kb)2`, y = freq1, yend = freq2), colour = "darkgray", linetype = 2, show.legend = FALSE, inherit.aes = FALSE)
        ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/Ind_Plots/Hap_change_repeatability_by_peak_size_spearman", chem$Chemical[1], "_", counter, ".pdf"), gplot, width = 8.5, height = 9, units = "in")
    } )
} )

# ============================================================================
# Calculating z-scores for log transformed LOD data and comparing the 3 vs 3 replicates genomewide via Spearmans' rho

popDT <- do.call(rbind, repLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
lodDF <- as.data.frame(lodDT)
lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
lodDT <- as.data.table(lodDF)

# ============================================================================
# Plot LOD scores for all pairs of treatments with trend line to look for correlations within each chemical; doesn't change the fact that there is less correlation

popDT <- do.call(rbind, repLODDTs)
popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
chemSplit <- split(popDT, popDT$chemWeek)
chemLoop <- lapply(chemSplit, function(chem) {
    print(chem$chemWeek[1])
    flush.console()
    replodDT <- dcast(chem, gp + chr + pos ~ Chemical, value.var = "LOD")
    names(replodDT) <- gsub("18way_", "", names(replodDT))
    names(replodDT) <- gsub("_", " ", names(replodDT))
    replodDF <- as.data.frame(replodDT)
    #replodDF[,4:ncol(replodDF)] <- scale(replodDF[,4:ncol(replodDF)])
    #gg1 = makePairs(replodDF[,-c(1:3)])
    #mega_iris = data.frame(gg1$all)
    panelPlot <- ggplot(replodDF, aes(replodDF[,4], replodDF[,5])) + xlab("") + ylab("") + facet_grid(names(replodDF)[4] ~ names(replodDF[5]), scales = "free", labeller = labeller(yvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) 
    panelPlot
} )
savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_zscore_correlations_spearman.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the Spearman correlation between LOD scores considering 1:n top peaks for each chemical separately


popDT <- do.call(rbind, repLODDTs)
topHapsDF <- read.table("topHapCombosAllChems.txt", header = T)
topHaps <- as.character(unlist(strsplit(topHapsDF$topHapCombos, ";")))
repShapes <- data.frame(reps = c(paste0("R0", 1:9), paste0("R1", 0:6)), shapes = c(0:6, 8:12, 15:18))

firstRep <- repLODDTs[seq(1, length(repLODDTs), 2)]
colNames <- rep(list('LOD'), length(firstRep))
#threshold <- rep(list(25), length(indLODDTs))
threshold <- rep(list(50), length(firstRep))
allPeaks <- Map(inflect, firstRep, colNames, threshold)
peakCol <- Map(add.sig.column, firstRep, allPeaks, rep.list(5, firstRep))
top <- lapply(peakCol, function(x) {
    peaks <- x[significant == 1]
    topPeaks <- peaks[order(-LOD)][1:20][, Idx := 1:20] ### looking at the top 20 peaks
    topPeaks
} )
peakPos <- lapply(top, function(x) {x$gp})
lodInts <- lod.support.intervals(2.5)
addLodInts <- Map(lodInts, peakCol)
lodWins <- Map(find.lod.dts, addLodInts, peakPos)
allLodWins <- do.call(rbind, lodWins)
popDTsub <- merge(popDT, allLodWins[, c("gp", "Idx", "chemWeek")], by = c("gp", "chemWeek"))[order(Chemical, gp)]

popDT[, logLOD := log(LOD), by = seq_len(nrow(popDT))]
#lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "logLOD")
lodDT <- dcast(popDT, gp + chr + pos ~ Chemical, value.var = "LOD")
names(lodDT) <- gsub("18way_", "", names(lodDT))
names(lodDT) <- gsub("_", " ", names(lodDT))
#lodDF <- as.data.frame(lodDT)
#lodDF[,4:ncol(lodDF)] <- scale(lodDF[,4:ncol(lodDF)])
#lodDT <- as.data.table(lodDF)
#lodDT2 <- melt(lodDT, id.vars = c("gp", "chr", "pos"), measure.vars = c(names(lodDT)[4:ncol(lodDT)]), variable.name = "Chemical", value.name = "zscore")
lodDT2 <- melt(lodDT, id.vars = c("gp", "chr", "pos"), measure.vars = c(names(lodDT)[4:ncol(lodDT)]), variable.name = "Chemical", value.name = "lod")
popDTsub[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
mrgeZs <- merge(popDTsub, lodDT2, by = c("gp", "Chemical"))
grp1 <- mrgeZs[Idx == 1][, grp := 1]
grp2 <- mrgeZs[Idx %like% 1:5][, grp := 2]
grp3 <- mrgeZs[Idx %like% 1:10][, grp := 3]
grp4 <- mrgeZs[Idx %like% 1:20][, grp := 4]
grp5 <- mrgeZs[Idx %like% 1:20][, grp := 5]
allGrps <- list(grp1, grp2, grp3, grp4, grp5)

idxLoop <- lapply(allGrps, function(idx) {
    print(idx$grp[1])
    flush.console()
    fileName <- paste0("within_chem_grp_", idx$grp[1], "_v2.pdf")
    #replodDT <- dcast(idx, gp + chr.x + pos.x + Idx ~ Chemical, value.var = "zscore")
    #names(replodDT) <- gsub("18way_", "", names(replodDT))
    #names(replodDT) <- gsub("_", " ", names(replodDT))
    #replodDF <- as.data.frame(replodDT)
    chemSplit <- split(idx, idx$chemWeek)
    chemLoop <- lapply(chemSplit, function(chem){
        replodDT <- dcast(chem, gp + chr.x + pos.x + grp ~ Chemical, value.var = "LOD")
        replodDF <- as.data.frame(replodDT)
        panelPlot <- ggplot(replodDF, aes(replodDF[,5], replodDF[,6])) + xlab("") + ylab("") + facet_grid(names(replodDF)[5] ~ names(replodDF[6]), scales = "free", labeller = labeller(yvar = label_wrap_gen(16))) + geom_point() + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + geom_smooth(span = 1) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1))
        panelPlot
    } )
    savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
    savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
    ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/", fileName), savingPlot, width = 8.5, height = 9, units = "in")
} )

# ============================================================================
# Find the average Spearman correlation per site genomewide for each chemical separately.

job::job(gpCorJob = {
    chemReps <- unique(names(allCorDF))
    chemNames <- gsub("X18way_|_12.*", "", chemReps)
    names(chemReps) <- chemNames
    chemSplit <- split(chemReps, names(chemReps))
    byChems <- lapply(chemSplit, function(chem) {
        ids <- paste(chem, collapse = "|")
        tic()
        corrGrp <- allCorDF[,grep(ids, names(allCorDF))]
        corrGrp <- as.data.frame(corrGrp)
        #rownames(corrGrp) <- rownames(allCorDF)
        gpos <- substr(rownames(corrGrp), 1, unlist(gregexpr(pattern ='\\.',rownames(corrGrp)))[1] - 1)
        gpos <- gsub("\\..*", "", gpos)
        corrGrp$id <- gsub(".*18way_", "", rownames(allCorDF))
        corrGrp$id <- gsub("_12.*-", "-", corrGrp$id)
        corrGrp$gp <- gpos
        #corrGrp$gp <- as.numeric(gsub(".*_12-|-R.*", "", rownames(allCorDF)))
        corrGrpDT <- as.data.table(corrGrp)
        chemName <- gsub("X18way_|_12.*", "", names(corrGrp)[1])
        print(chemName)
        flush.console()
        corrGrpDT <- corrGrpDT[grepl(chemName, id)]
        corrGrpDT2 <- corrGrpDT[, id := NULL]
        corrLst <- split(corrGrpDT2, corrGrpDT2$gp)
        rmveGP <- lapply(corrLst, function(x) {subset(x, select=-c(gp))})
        corrLstMats <- lapply(rmveGP, as.matrix)
        matMeans <- unlist(lapply(corrLstMats, function(x) mean(x[upper.tri(x)])))
        meanCorDT <- data.table(Chemical = chemName, gp = as.numeric(names(matMeans)), avgRho = matMeans)[order(gp)]
        toc()
        meanCorDT
    } )
}, import = c(allCorDF), packages = c("data.table", "tictoc") )

gpCorDT <- do.call(rbind, gpCorJob$byChems)     

# ============================================================================
# Plot the correlation between LOD score and spearman correlation per site genomewide. No correlation!

popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := gsub("18way_", "", Chemical)]
chem2 <- merge(popDT, gpCorDT, by = c("gp", "Chemical"))
chem2[, Chemical := gsub("_", " ", Chemical)]
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, avgRho)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/LOD_vs_avgRho_correlations.pdf"), panelPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation between LOD score and mean per site haplotype deviation genomewide.Again no correlation!

popDT <- do.call(rbind, indLODDTs)
devDT <- do.call(rbind, indHapDevsDTs)
avgDevDT <- devDT[, .(sqrtAvgHapDifs = mean(sqrtAvgHapDifs)), by = c("chemWeek", "gp")]
chem2 <- merge(popDT, avgDevDT, by = c("gp", "chemWeek"))
chem2[, Chemical := gsub("18way_|_", " ", Chemical)]
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, sqrtAvgHapDifs)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_correlations.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation between mean per site correlation and mean per site haplotype deviation genomewide.

devDT <- do.call(rbind, indHapDevsDTs)
avgDevDT <- devDT[, .(sqrtAvgHapDifs = mean(sqrtAvgHapDifs)), by = c("chemWeek", "gp")]
chem2 <- merge(gpCorDT, avgDevDT, by = c("gp", "chemWeek"))
chem2[, Chemical := gsub("18way_|_", " ", chemWeek)]
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(avgRho, sqrtAvgHapDifs)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_correlations.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation between LOD score and mean standard deviation of the most increased haplotype at each position genomewide.

popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
chem2 <- merge(popDT, sdDT, by = c("gp", "Chemical"))
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, sdMax)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_correlations.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation between LOD score and number of unique max haplotypes genomewide.
popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
chem2 <- merge(popDT, maxHapsDT, by = c("gp", "Chemical"))
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, maxHaps)) + facet_wrap(~as.factor(Chemical), scales = "free_x", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 1)) + stat_cor(aes(label = ..r.label..), label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) 
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/LOD_vs_maxHaps.pdf"), panelPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation between LOD score and avg max haplotypes frequency change per site genomewide.

popDT <- do.call(rbind, indLODDTs)
popDT[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
chem2 <- merge(popDT, maxHapsDT, by = c("gp", "Chemical"))
chemDF <- as.data.frame(chem2)
panelPlot <- ggplot(chemDF, aes(LOD, avgFreqChange)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + stat_cor(aes(label = ..r.label..), method = "spearman", label.x.npc = "left", label.y.npc = "top", size = 3) + scale_fill_continuous(type = "viridis") + geom_smooth(span = 0.3) + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_correlations.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

# ============================================================================
# Plot the correlation the number of max haplotypes and the avg max haplotypes frequency change per site genomewide.

chemDF <- as.data.frame(maxHapsDT)
panelPlot <- ggplot(chemDF, aes(maxHaps, avgFreqChange)) + facet_wrap(~as.factor(Chemical), scales = "free", labeller = labeller(xvar = label_wrap_gen(16))) + geom_hex(bins = 50) + scale_fill_continuous(type = "viridis") + stat_cor(aes(label = ..r.label..)) + geom_smooth(method = "lm") + theme_bw(base_size = 8.5) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust = 1)) + scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

savePlot <- do.call(ggarrange, c(chemLoop, list(common.legend = TRUE, legend = "right")))
savingPlot <- annotate_figure(savePlot, bottom = text_grob("LOD1", color = "black", size = 14), left = text_grob("LOD2", color = "black", rot = 90, size = 14))
ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Repeatability_plots/within_chem_correlations.pdf"), savingPlot, width = 8.5, height = 9, units = "in")

maxHapsDT <- popDT2[, .(maxHaps = uniqueN(maxHap)), by = c("Chemical", "gp")]

# ============================================================================
# Trouble-shooting

chemSplit <- split(popDT, popDT$chemWeek)
chemLoop <- lapply(chemSplit, function(chem) {
print(chem$chemWeek[1])
flush.console()
    

library(ggfortify)
