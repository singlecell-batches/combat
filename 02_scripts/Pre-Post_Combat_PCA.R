setwd("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data")
getwd()
#g1 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group1.matrix.txt", header = F) #Read in counts matrices
#g2 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group2.matrix.txt", header = F)
#g3 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group3.matrix.txt", header = F)
#g4 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group4.matrix.txt", header = F)
#g5 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group5.matrix.txt", header = F)
#g6 <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\group6.matrix.txt", header = F)

#colnames(g1) <- c("TranscriptID", "GeneID", paste("c",1:1000, sep = "")) #add some column names
#colnames(g2) <- c("TranscriptID", "GeneID", paste("c",1001:2000, sep = ""))
#colnames(g3) <- c("TranscriptID", "GeneID", paste("c",2001:7000, sep = ""))
#colnames(g4) <- c("TranscriptID", "GeneID", paste("c",7001:11000, sep = ""))
#colnames(g5) <- c("TranscriptID", "GeneID", paste("c",11001:12500, sep = ""))
#colnames(g6) <- c("TranscriptID", "GeneID", paste("c",12501:13300, sep = ""))

#mg <- merge(g1, g2, by = c("TranscriptID", "GeneID")) #merging the data
#mg1 <- merge(mg, g3, by = c("TranscriptID", "GeneID"))
#mg2 <- merge(mg1, g4, by = c("TranscriptID", "GeneID"))
#mg3 <- merge(mg2, g5, by = c("TranscriptID", "GeneID"))
#mm <- merge(mg3, g6, by = c("TranscriptID", "GeneID"))

#merged_matrix <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\merged.txt", header = F)

Counts <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\10percentallgenes.txt", header = TRUE, row.names = NULL)
#Counts1 <- sapply(Counts[, -1], as.numeric)

head(Counts)
dim(Counts)

meta <- read.table("C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\10percentmeta.txt", header = T)
head(meta)
meta$Batch <- factor(meta$Batch)
meta$Batch

props <- apply(Counts[, -1],1,function(x) length(which(x==0))/length(x)) #finds gene expression proportion
w <- which(props > 0.8) #which gene did not express much?
CountsSubset <- Counts[-w,] #remove the low expression genes
CountsSubset1 <- sapply(CountsSubset[, -1], as.numeric)
dim(CountsSubset1)

write.table(CountsSubset1, file = "C:\\Users\\anbra\\Desktop\\Single_Cell_Batches\\Data\\TenpercentSubsetgenes.txt")

CountsSubsetScaled <- scale(CountsSubset1) #or cb if run through ComBat
head(CountsSubsetScaled[,1])

pcmm <- prcomp(t(CountsSubsetScaled))
attributes(pcmm)
head(pcmm$x)
?prcomp

library(ggplot2) #can use ggplot
library(ggfortify) #ggfortify only needed for autoplot

#autoplot(pcmm,data = meta[c(s1,s2,s3,s4,s5,s6),], colour = "Batch", main = "100 from each batch no scale") #takes too long

plot(pcmm$x[,1],pcmm$x[,2], col = meta$Batch, pch = 16, xlab = "PC1", ylab = "PC2",main = "Scaled Data")
legend(locator(1), legend = c("1","2","3","4","5","6"), col = levels(meta$Batch), pch = 16) #col is color

source("https://bioconductor.org/biocLite.R")
biocLite("sva")
library(sva)
?ComBat

cb <- ComBat(dat = CountsSubset1, batch=meta$Batch, mod=NULL, par.prior=TRUE, mean.only=FALSE)
head(cb[, 1])
dim(cb)

CbScaled <- scale(cb) #or cb if run through ComBat
pcmmCB <- prcomp(t(CbScaled))

plot(pcmmCB$x[,1],pcmmCB$x[,2], col = meta$Batch, pch = 16, xlab = "PC1", ylab = "PC2",main = "Combat Scaled Data")
legend(locator(1), legend = c("1","2","3","4","5","6"), col = levels(meta$Batch), pch = 16) #col is color

pcmmCBNoScale <- prcomp(t(cb))
plot(pcmmCBNoScale$x[,1],pcmmCBNoScale$x[,2], col = meta$Batch, pch = 16, xlab = "PC1", ylab = "PC2",main = "Combat Data")
legend(locator(1), legend = c("1","2","3","4","5","6"), col = levels(meta$Batch), pch = 16) #col is color
