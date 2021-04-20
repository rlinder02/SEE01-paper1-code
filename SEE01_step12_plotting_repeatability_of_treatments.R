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
avgDT <- popDT[, .(sqrtAvgHapDifs = mean(sqrtAvgHapDifs)), by = chemWeek]
maxDev <- max(avgDT$sqrtAvgHapDifs)
txDTs <- split(avgDT, avgDT$chemWeek)


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
        DT[, axis(2, at=c(0,0.5), labels=c(0,0.5), las = 1, cex.axis = 0.75)]
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
allCorDFs <- lapply(files, function(read) {
    read.table(read, header = T, sep = "\t")
} )
setwd(paste0(projectDir, "Data_Storage/Sequencing_Data_Processed/"))

byChems <- lapply(allCorDFs, function(corr) {
    tic('total time')
    chemReps <- unique(names(corr))
    print(chemReps[1])
    flush.console()
    corAvgs <- lapply(chemReps, function(x) {
        corrGrp <- corr[,grep(x, names(corr))]
        corrGrp <- as.data.frame(corrGrp)
        names(corrGrp) <- x
        rownames(corrGrp) <- rownames(corr)
        corrGrp$id <- unlist(lapply(rownames(corrGrp), function(y) {strsplit(y, "-")[[1]][3]}))
        corrGrpDT <- as.data.table(corrGrp)
        corrAvgs <- corrGrpDT[, mean(get(x), na.rm = T), by = id]
        setnames(corrAvgs, c("id", x))
        corrAvgs
    } )
} )

chemCorDTs <- lapply(byChems, function(x) {
    bindReps <- do.call(cbind, x)
    bindRepsDF <- as.data.frame(bindReps, stringsAsFactors = FALSE)
    chemName <- strsplit(names(bindRepsDF)[2], "\\.")[[1]][1]
    idCols <- which(names(bindReps) == "id")
    bindReps2 <- bindRepsDF[,!duplicated(colnames(bindRepsDF))]
    rownames(bindReps2) <- bindReps2$id
    names(bindReps2)[2:length(names(bindReps2))] <- bindReps2$id
    bindReps2$id <- NULL
    chemDT <- list(bindReps2)
    names(chemDT) <- chemName
    chemDT
} )

ways <- sapply(chemCorDTs, function(x) {strsplit(names(x), "_")[[1]][1]})
waySplit <- split(chemCorDTs, ways)
wayLooper <- lapply(waySplit, function(wayList) {
    weeks <- sapply(wayList, function(x) {strsplit(names(x), "_")[[1]][length(strsplit(names(x), "_")[[1]])]})
    print(weeks)
    flush.console()
    weekSplit <- split(wayList, weeks)
    weekLooper <- lapply(weekSplit, function(weekList) {
        fileName1 <- gsub("^X", "", names(weekList[[1]]))
        fileName2 <- paste(strsplit(fileName1, "_")[[1]][c(1, length(strsplit(names(weekList[[1]]), "_")[[1]]))], collapse = "_wk")
        fileName3 <- paste0(fileName2, "_haplotype_spearman_correlations.pdf")
        print(fileName3)
        flush.console()
        plotLooper <- lapply(1:length(weekList), function(x) {
            data <- as.matrix(weekList[[x]][[1]])
            df <- reshape2::melt(data)
            mainTitle <- paste(strsplit(names(weekList[[x]]), "_")[[1]][2:(length(strsplit(names(weekList[[x]]), "_")[[1]])-1)], collapse = "_")
            gplot <- ggplot(data = df) + geom_tile(aes(x = factor(Var1), y = factor(Var2), fill = value)) + theme_bw(base_size = 8) + scale_fill_continuous(limits=c(0, 1), low = "white", high = "darkred") +  theme(text = element_text(size=10)) + theme(axis.title=element_blank()) + theme(legend.title=element_blank()) + ggtitle(mainTitle) + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.border = element_blank()) + theme(axis.text.x = element_text(angle = 45)) + coord_equal()
            gplot
        } )
        arrangePlot <- do.call(ggarrange, c(plotLooper, list(common.legend = TRUE, legend = "right")))
        savePlot <- annotate_figure(arrangePlot, top = text_grob(paste0(fileName2, "_replicate_correlations"), color = "black", size = 10, hjust = 0.5, vjust = -0.05))
        ggsave(file = paste0(projectDir, "Data_Analysis/Sequencing_analysis/Plots/Summary_plots/SEE01_", fileName3), savePlot, width = 8.5, height = 8.5, units = "in")
    } )
} )

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


avgChemCorDFs <- lapply(byChems, function(x) {
    bindReps <- do.call(cbind, x)
    bindRepsDF <- as.data.frame(bindReps, stringsAsFactors = FALSE)
    chemName <- strsplit(names(bindRepsDF)[2], "\\.")[[1]][1]
    idCols <- which(names(bindReps) == "id")
    bindReps2 <- bindRepsDF[,!duplicated(colnames(bindRepsDF))]
    rownames(bindReps2) <- bindReps2$id
    names(bindReps2)[2:length(names(bindReps2))] <- bindReps2$id
    bindReps2$id <- NULL
    avgCor <- mean(bindReps2[upper.tri(bindReps2)])
    avgCorDF <- data.frame(chemWeek = chemName, avgCor = avgCor, stringsAsFactors = FALSE)
    avgCorDF
} )

avgCorsDF <- do.call(rbind, avgChemCorDFs)
avgCorDT <- as.data.table(avgCorsDF)
avgCorDT$chemWeek <- gsub("^X", "", avgCorDT$chemWeek)
avgCorHetDF <- merge(avgHets, avgCorDT, by = "chemWeek")
avgCorHetDF[chemWeek %like% "YPD", c("chemWeek", "id") := .(gsub("YPD", "ypd", chemWeek), 'ypd')]
avgCorHetDF <- avgCorHetDF[order(avgCorHetDF$chemWeek),]
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
# Find the most increased haplotype within each window .

idxLooper <- lapply(chemIndexer, function(idx)  {
    freqDifs <- idx[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
    collHaps <- as.data.frame(idx[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
    maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
    #minVals <- apply(freqDifs, 1, min, na.rm = T)
    #maxMinDF <- data.frame(max = maxVals, min = abs(minVals))
    #maxRows <- apply(maxMinDF, 1, which.max)
    maxIdx <- apply(freqDifs, 1, which.max)
    #minIdx <- apply(freqDifs, 1, which.min)
    idxDF <- data.frame(row = 1:length(maxIdx), col = maxIdx)
    maxHap <- collHaps[as.matrix(idxDF)]
    idx[, c("maxHap", "maxChange") := .(maxHap, maxVals)]
    idx <- idx[, c("chr", "pos", "gp", "Chemical", "Replicate", "maxHap", "maxChange", "Idx")]
    idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
    hapColors <- match(idx$maxHap, topHaps)
    hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
    Colors <- mycols[hapColors]
    reps <- match(idx$Replicate, repShapes$reps)
    Shapes <- repShapes$shapes[reps]
    idx[, c("color.codes", "Shape") := .(Colors, Shapes)][, "pos(kb)" := pos/1000][, "pos" := NULL]
    idx
})

topCorLooper <- lapply(chemIndexer, function(idx)  {
    chemSplit <- split(idx, idx$Chemical)
        chemLooper <- lapply(chemSplit, function(chem) {
        freqDifs <- chem[, tstrsplit(freqDifs, split = ";", type.convert = TRUE, fixed = TRUE)]
        freqDifsDF <- as.data.frame(freqDifs)
        collHaps <- as.data.frame(chem[, tstrsplit(collapsedFounders, split = ";", type.convert = TRUE, fixed = TRUE)])
        collHaps$gp <- chem$gp
        maxVals <- apply(freqDifs, 1, max, na.rm = T) ## only looking at increasing values
        #minVals <- apply(freqDifs, 1, min, na.rm = T)
        #maxMinDF <- data.frame(max = maxVals, min = abs(minVals))
        #maxRows <- apply(maxMinDF, 1, which.max)
        maxIdx <- apply(freqDifs, 1, which.max)
        #minIdx <- apply(freqDifs, 1, which.min)
        idxDF <- data.frame(row = 1:length(maxIdx), col = maxIdx)
        maxHap <- collHaps[as.matrix(idxDF)]
        chem[, c("maxHap", "maxChange") := .(maxHap, maxVals)]
        maxHaps <- chem[, unique(maxHap)]
        gpSplit <- split(collHaps, collHaps$gp)
        gpLooper <- lapply(gpSplit, function(gpos) {
            maxHapsIdx <- data.frame(apply(gpos, 1, function(x) which(x %in% maxHaps)))
            names(maxHapsIdx) <- gsub("X", "", names(maxHapsIdx))
            columns <- nrow(maxHapsIdx)
            freqIdxDF <- data.frame(row = sort(rep(as.numeric(names(maxHapsIdx)), columns)), col = unlist(maxHapsIdx))
            topFreqs <- freqDifsDF[as.matrix(freqIdxDF)]
            hapVecs <- split(topFreqs, cut(seq_along(topFreqs), length(topFreqs)/columns, labels = FALSE))
            hapVecDF <- as.data.frame(hapVecs)
            cors <- cor(hapVecDF, use = "everything", method = "pearson")
            corsDF <- as.data.frame(cors)
            
            slimChem <- chem[, c("chr", "pos", "gp", "Chemical", "Replicate", "maxHap", "maxChange", "Idx")]

        
        eachPos <- slimChem[, maxHap, by = gp]
        idx[, Chemical := gsub("18way_", "", Chemical)][, Chemical := gsub("_", " ", Chemical)]
        hapColors <- match(idx$maxHap, topHaps)
        hapColors[is.na(hapColors)] <- (length(topHaps) + 1)
        Colors <- mycols[hapColors]
        reps <- match(idx$Replicate, repShapes$reps)
        Shapes <- repShapes$shapes[reps]
        idx[, c("color.codes", "Shape") := .(Colors, Shapes)][, "pos(kb)" := pos/1000][, "pos" := NULL]
        idx
})

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

chemCorsDT <- do.call(rbind, byChems)     

# ============================================================================
# Find the Pearson correlation coefficient for the most significant peaks, the middle peak, and a smaller peak for each drug that only looks at the most increased haplotype (if different replicates have different most changed haplotypes, use the mean correlation of these haplotypes). Make a dataframe that has the Chemical and averaged correlation coefficient over that interval. Need to convert Pearson coefficient to z-score, average the z-scores, then back-transform.

hapFreqs <- fread("Reformatted_SEE01_hap_freqs_v2.txt")




# ============================================================================
# Plot the frequency of the most changed haplotype (for each replicate), with each replicate a different shape and the haplotypes colored consistently with previous plots. Have these plotted as a separate sub-panel for each chemical using facet. Show the most significant peak, a middle peak, and a smaller peak that is likely still real to show how repeatability changes as a function of the LOD score, with each type of peak a separate panel (A-C). Include as text the Pearson correlation coefficient.

counter <- 0
plotLooper <- lapply(idxLooper, function(idce) {
    counter <<- counter + 1
    chrDT <- idce[, .(chr = chr[1]), by = "Chemical"]
    corDT <- chemCorsDT[idx %in% idce$Idx[1]]
    ggplot2::theme_update(plot.tag = element_text(face = "bold", colour = "black"))
    gplot <- ggplot(data=idce, aes(`pos(kb)`, maxChange, group = Replicate, colour = maxHap, shape = Replicate)) + facet_wrap(~as.factor(Chemical), scales = "free_x") + coord_cartesian(ylim = c(0, 1)) + xlab("position (kb)") + ylab("haplotype frequency change") + geom_point(size = 2) + geom_line() + scale_colour_manual(values=setNames(idce$color.codes, idce$maxHap)) + scale_shape_manual(values = setNames(idce$Shape, idce$Replicate)) + theme_bw(base_size = 12) + theme(panel.grid = element_blank()) + geom_text(data = chrDT, aes(x = Inf, y = Inf, hjust = 1.1, vjust = 1.25, label = paste0('chr', as.roman(chr))), inherit.aes = FALSE) + geom_text(parse = TRUE, data = corDT, aes(x = -Inf, y = Inf, hjust = -0.15, vjust = 1.25, label = round(avgRho, 2)), inherit.aes = FALSE) + labs(tag = LETTERS[counter])
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
# Trouble-shooting


library(ggfortify)
