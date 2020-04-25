y <- DGEList(count = counts.keep, group = group) # DGE list object in edgeR is used to store count data
ynorm <- calcNormFactors(y)

group
design <- model.matrix(~ 0 + group)
design
colnames(design) <- levels(group)
design

y$samples$lib.size
par(mfrow=c(1,1))
par(mar=c(7,7,4,2)+0.1,mgp=c(5,1,0)) # read about the description of this command
barplot(y$samples$lib.size,names=colnames(y),las=2) #barplot of the library size
abline(h=1.50e+07,col="red") # library size of all of the samples were found to be more than 2.75e+07
title(main = "Barplot of library sizes", xlab = "Sample number", ylab = "Number of reads")
dev.off()
logcounts <- cpm(y,log=TRUE) # log transformation of the library size has been performed so as to shrink the data so as make it easier for visualization
boxplot(logcounts, xlab="Sample number ", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue") # and the blue line is the median of all the library sizes
title("Boxplots of logCPMs (unnormalised)")
plotMDS(y)
title("MDS plot") # Multi dimentional scaling was performed on the unnormalized log transformes cpm's
sampleinfo <- read.delim("Bioneer_sampleInfo.txt")
sampleinfo
levels(sampleinfo$CellType)
col.cell <- c("purple","orange","dark green","blue","red","black")[sampleinfo$CellType]
data.frame(sampleinfo$CellType,col.cell)
par(mfrow=c(1,1))
plotMDS(y,col=col.cell, top = 1000)
legend("topright",fill=c("purple","orange","dark green"),legend=levels(sampleinfo$CellType),cex=0.5)
title("Cell type")

v <- voom(ynorm,design,plot = TRUE) #This plot can also tell us if there are any genes that look really variable in our data, and if we’ve filtered the low counts adequately.
#voom: transforms the readcounts into logCPM while taking into consideration the mean variance of the data
v
names(v)
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
fit <- lmFit(v,design) # lmFit function was used to find out the differentially expresssed genes. It uses voom transformed function along with the design matix that we have already specified in the command.
#fit <- glmQLFit(v)
#lmFit calculates group's mean based on the design matirix as well as gene wise variance.
names(fit)

#this contrast is basically on the basis of growth phase varoation
cont.matrix1 <- makeContrasts(Clone1 = ((ProDay3-WT2)+(ProDay1-WT1)), levels=design)
#NR is indicating noise reduction
#This is only for vertical contrast


fit.cont1 <- contrasts.fit(fit, cont.matrix1) #Then we fit the contrast matrix with the help of fit function
fit.cont1 <- eBayes(fit.cont1)
dim(fit.cont1)
summa.fit1 <- decideTests(fit.cont1)
summary(summa.fit1)

Contrast_Clone1 <- topTable(fit.cont1,coef="Clone1",sort.by="p", number = "Inf")

Contrast_Clone1 <- cbind(rownames(Contrast_Clone1), Contrast_Clone1)
rownames(Contrast_Clone1) <- NULL
colnames(Contrast_Clone1) <- c("SYMBOL", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

par(mfrow=c(1,4))
volcanoplot(fit.cont1,coef=1,highlight=0,names=rownames(fit.cont1), main = "Pro3")
write.table(topTable(fit.cont1,coef="Clone1",sort.by="p", number = "Inf"), "Clone1.txt")

#biocLite("piano") # Only once
#library(piano)
#biocLite("snowfall") # Only once
#library(snowfall); library(snow)

gsc_all <- loadGSC("Glycerolipid.gmt", type="gmt")
gsc_all <- loadGSC("Stress_related_pathways.gmt", type="gmt")
gsc_all <- loadGSC("c2.all.v6.2.symbols.gmt", type="gmt")
gsc_all <- loadGSC("c5.all.v6.2.symbols.gmt", type="gmt")

padj1 <- Contrast_Clone1$adj.P.Val # All these are CHo genes
logfc1 <- Contrast_Clone1$logFC # All these are CHo genes
names(padj1) <- names(logfc1) <- Contrast_Clone1$SYMBOL # All these are CHo genes
gsaRes_all1 <- runGSA(padj1, logfc1, gsc=gsc_all, gsSizeLim = c(1,Inf))
tvalue1 <- Contrast_Clone1$t # t-statistics has been calculated with the help significant CHO gene sets
names(tvalue1) <- Contrast_Clone1$SYMBOL # column names has been provided with the help of human ENTREZID
gsaRes1 <- runGSA(tvalue1, geneSetStat="mean", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes2 <- runGSA(tvalue1, geneSetStat="median", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes3 <- runGSA(tvalue1, geneSetStat="sum", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes4 <- runGSA(tvalue1, geneSetStat="maxmean", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes5 <- runGSA(padj1, logfc1, geneSetStat="fisher", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes6 <- runGSA(padj1, logfc1,geneSetStat="stouffer", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")
gsaRes7 <- runGSA(padj1, logfc1, geneSetStat="tailStrength", gsc=gsc_all, nPerm = 1000, gsSizeLim = c(10,800), adjMethod = "none")

resList1 <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
names(resList1) <- c("mean","median","sum","maxmean","fisher","stouffer","tailStrength")
ch1 <- consensusHeatmap(resList1, cutoff=30, method="mean", ncharLabel = 50)
#write.csv(ch1, "~/Desktop/PhD/RNA-seq1/Bioneer_data/PacBio/Lipid/ch1.csv")
GSAsummaryTable(gsaRes_all1, save=T, file="Contrast_Clone1_C2.txt")

cs1 <- consensusScores(resList1,class="non")
cs1_mu <- consensusScores(resList1,class="mixed",direction="up")
cs1_md <- consensusScores(resList1,class="mixed",direction="down")
cs1_dd <- consensusScores(resList1,class="distinct",direction="down")
cs1_du <- consensusScores(resList1,class="distinct",direction="up")
