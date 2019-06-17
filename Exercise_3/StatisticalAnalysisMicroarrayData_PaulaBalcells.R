## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----include=FALSE-------------------------------------------------------
if (!require(BiocManager)) install.packages("BiocManager")

installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
    BiocManager::install(pkg)
}else{
  require(pkg, character.only=T)
  }
}


## ----echo = T, results = 'hide', message=FALSE---------------------------
installifnot("pd.mogene.1.0.st.v1")
installifnot("mogene10sttranscriptcluster.db")
installifnot("oligo")
installifnot("limma")
installifnot("Biobase")
installifnot("arrayQualityMetrics")
installifnot("genefilter")
installifnot("multtest")
installifnot("annotate")
installifnot("xtable")
installifnot("gplots")
installifnot("scatterplot3d")
installifnot("affyPLM")
installifnot("hgug4112a.db")
installifnot("ReactomePA")

install.packages("survival")
install.packages("foreign")
install.packages("Hmisc")
install.packages("gdtools")
install.packages("gridSVG")
install.packages("ggrepel")


## ----echo = T, results = 'hide', message=FALSE---------------------------
workingDir <-getwd()
dataDir <- file.path(workingDir, "dades")
resultsDir <- file.path(workingDir, "results")


## ----echo = T, results = 'hide', message=FALSE---------------------------
require(GEOquery)
#getGEOSuppFiles(GEO = "GSE66597", makeDirectory = TRUE)
rawDir <- file.path(workingDir, "GSE66597/GSE66597_RAW")


## ------------------------------------------------------------------------
targets <- readTargets("targets.txt",row.names="FileName",sep="")
fileNames <- as.character(targets$FileName)
sampleNames <- as.character(targets$sampleName)
sampleGroup <- as.character(targets$Group)
sampleReplicate <- as.character(targets$Replicate)
sampleColor <- as.character(targets$Color)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  targets, booktabs = TRUE, caption = "Targets file"
)


## ----echo = T, results = 'hide', message=FALSE---------------------------
rawData <- read.maimages(targets,path=rawDir,source="agilent", green.only=TRUE)


## ----echo = T, results = 'hide'------------------------------------------
# require(arrayQualityMetrics)
# arrayQualityMetrics(rawData, outdir = file.path(resultsDir, "QCDir.Raw"), force=TRUE)


## ------------------------------------------------------------------------
rawDataLog <- log2(rawData$E)


## ------------------------------------------------------------------------
boxplot(rawDataLog, which="all",las=2, main="Intensity distribution of RAW data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(log2(rawDataLog))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of RawData", 
     cex=0.7,  hang=-1)


## ----message=FALSE-------------------------------------------------------
require(ggplot2)
require(ggrepel)

plotPCA3 <- function(X, labels, factor, title, scale, colors, size=1.5, glineas=0.25){
  data <- prcomp(t(X),scale=scale)
  dataDf <- data.frame(data$x)
  Group <- factor
  loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
  p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
    theme_classic() +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(aes(color = Group), alpha = 0.55, size = 3) +
    coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
    scale_fill_discrete(name = "Group")
  p1 + geom_text_repel(aes(y = PC2 + 0.25, label = labels), segment.size = 0.25, size = size) + 
    labs(x = c(paste("PC1",loads[1],"%")), y = c(paste("PC2", loads[2], "%"))) +
    ggtitle(paste("Principal Component Analysis for: ", title, sep=" ")) +
    scale_color_manual(values=colors)
}


## ------------------------------------------------------------------------
plotPCA3(rawDataLog, labels = targets$sampleName, factor = targets$Group, title = "Raw data", scale = FALSE, size = 3, colors = targets$Color)


## ------------------------------------------------------------------------
control <- rawData$genes$ControlType==1
filteredData <- rawData[!control,]
dup <- duplicated(filteredData$genes$ProbeName)
filteredData <- filteredData[!dup,]


## ----message=FALSE-------------------------------------------------------
bg.set <- limma::backgroundCorrect(filteredData, method = "normexp")
nset <- normalizeBetweenArrays(bg.set, method = "quantile")


## ------------------------------------------------------------------------
expressions <- data.frame(nset$E, row.names = nset$genes$ProbeName)
genes <- t(expressions)[0,]

#AssayData
require(Biobase)
expressions <- as.matrix(expressions)
eset <- ExpressionSet(expressions)

#Information about covariates
columnDesc <-  data.frame(labelDescription= c("File Name", "Sample Name", "Group", "Replicate", "Color"))
annot <- new("AnnotatedDataFrame", data=targets, varMetadata= columnDesc)

rownames(pData(annot))<-pData(annot)$sampleName
eset <- ExpressionSet(assayData=expressions, phenoData=annot)

#Information about features
eset <- ExpressionSet(assayData=expressions, phenoData=annot, featureNames =genes)
show(eset)


## ----echo = T, results = 'hide', message=FALSE---------------------------
require(arrayQualityMetrics)
arrayQualityMetrics(eset, outdir = file.path(resultsDir, "QCDir.Normalized"), force=TRUE)


## ------------------------------------------------------------------------
boxplot(exprs(eset), which="all",las=2, main="Intensity distribution of normalized data", 
        cex.axis=0.6, col=sampleColor, names=sampleNames)


## ------------------------------------------------------------------------
clust.euclid.average <- hclust(dist(t(exprs(eset))),method="average")
plot(clust.euclid.average, labels=sampleNames, main="Hierarchical clustering of normalized data",
     cex=0.7,  hang=-1)


## ------------------------------------------------------------------------
plotPCA3(exprs(eset), labels = targets$sampleName, factor = targets$Group, title = "Normalized data", scale = T, size = 3, colors = targets$Color)


## ------------------------------------------------------------------------
sds <- apply(eset, 1, sd)
sdsO<- sort(sds)
plot(1:length(sdsO), sdsO, main="Distribution of variability for all genes", 
     sub="Vertical lines represent 90% and 95% percentiles", xlab="Gene index (from least to most variable)",
     ylab="Standard deviation") +
abline(v=length(sds)*c(0.9,0.95))


## ------------------------------------------------------------------------
require(genefilter)
require(hgug4112a.db)

annotation(eset) <- "hgug4112a.db"
filtered <- nsFilter(eset, require.entrez = TRUE, remove.dupEntrez = TRUE, var.filter = TRUE, var.func = IQR, var.cutoff = 0.75, filterByQuantile = TRUE, feature.exclude = "^AFFX")

print(filtered$filter.log)


## ------------------------------------------------------------------------
eset_filtered <- filtered$eset


## ----include=FALSE-------------------------------------------------------
write.csv(exprs(eset), file="./results/normalized.Data.csv")
write.csv(exprs(eset_filtered), file="./results/normalized.Filtered.Data.csv")
save(eset, eset_filtered, file="./results/normalized.Data.Rda")


## ------------------------------------------------------------------------
treat <- targets$Group
lev <- factor(treat, levels = unique(treat))
design <-model.matrix(~0+lev)
colnames(design) <- levels(lev)
rownames(design) <- targets$sampleName
print(design)


## ------------------------------------------------------------------------
require(limma)
cont.matrix <- makeContrasts (
  hm6vscon6 = hm6h-con6h,
  hm12vscon12 = hm12h-con12h,
  hm24vscon24 = hm24h-con24h,
  hm24vshm6 = hm24h-hm6h,
  levels=design)
print(cont.matrix)


## ------------------------------------------------------------------------
fit<-lmFit(eset, design)
fit.main<-contrasts.fit(fit, cont.matrix)
fit.main<-eBayes(fit.main)


## ------------------------------------------------------------------------
topTab_hm6vscon6 <- topTable(fit.main, number=nrow(fit.main), coef="hm6vscon6", adjust="fdr")
tt1 <- head(topTab_hm6vscon6)

topTab_hm12vscon12 <- topTable(fit.main, number=nrow(fit.main), coef="hm12vscon12", adjust="fdr")
tt2 <-head(topTab_hm12vscon12)

topTab_hm24vscon24  <- topTable(fit.main, number=nrow(fit.main), coef="hm24vscon24", adjust="fdr")
tt3 <- head(topTab_hm24vscon24)

topTab_hm24vshm6  <- topTable(fit.main, number=nrow(fit.main), coef="hm24vshm6", adjust="fdr")
tt4 <- head(topTab_hm24vshm6)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tt1, booktabs = TRUE, caption = "Top table results for hm6vscon6 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tt2, booktabs = TRUE, caption = "Top table results for hm12vscon12 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tt3, booktabs = TRUE, caption = "Top table results for hm24vscon24 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tt4, booktabs = TRUE, caption = "Top table results for hm24vshm6 comparison"
)


## ------------------------------------------------------------------------
require(hgug4112a.db)
annotatedTopTable <- function(topTab, anotPackage){
  topTab <- cbind(PROBEID=rownames(topTab), topTab)
  myProbes <- rownames(topTab)
  thePackage <- eval(parse(text = anotPackage))
  geneAnots <- select(thePackage, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
  annotatedTopTab <- merge(x=geneAnots, y=topTab, by.x="PROBEID", by.y="PROBEID")
  return(annotatedTopTab)
}


## ------------------------------------------------------------------------
topAnnotated_hm6vscon6 <- annotatedTopTable(topTab_hm6vscon6, anotPackage = "hgug4112a.db")
tta1 <- head(topAnnotated_hm6vscon6[1:6,1:4])

topAnnotated_hm12vscon12 <- annotatedTopTable(topTab_hm12vscon12, anotPackage = "hgug4112a.db")
tta2 <- head(topAnnotated_hm12vscon12[1:6,1:4])

topAnnotated_hm24vscon24 <- annotatedTopTable(topTab_hm24vscon24, anotPackage = "hgug4112a.db")
tta3 <- head(topAnnotated_hm24vscon24[1:6,1:4])

topAnnotated_hm24vshm6 <- annotatedTopTable(topTab_hm24vshm6, anotPackage = "hgug4112a.db")
tta4 <- head(topAnnotated_hm24vshm6[1:6,1:4])


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tta1, booktabs = TRUE, caption = "Annotated top table results for hm6vscon6 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tta2, booktabs = TRUE, caption = "Annotated top table results for hm12vscon12 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tta3, booktabs = TRUE, caption = "Annotated top table results for hm24vscon24 comparison"
)


## ----echo=FALSE----------------------------------------------------------
knitr::kable(
  tta4, booktabs = TRUE, caption = "Annotated top table results for hm24vshm6 comparison"
)


## ----echo = T, results = 'hide', message=FALSE---------------------------
geneSymbols <- select(hgug4112a.db, rownames(fit.main), c("SYMBOL"))
SYMBOLS <- geneSymbols$SYMBOL


## ----showresults---------------------------------------------------------
volcanoplot(fit.main, coef=1, highlight = 10, names=SYMBOLS, main=paste("Differentially expressed genes", colnames(cont.matrix)[1], sep="\n")) +
abline(v=c(-1,1))

volcanoplot(fit.main, coef=2, highlight = 10, names=SYMBOLS, main=paste("Differentially expressed genes", colnames(cont.matrix)[2], sep="\n")) +
abline(v=c(-1,1))

volcanoplot(fit.main, coef=3, highlight = 10, names=SYMBOLS, main=paste("Differentially expressed genes", colnames(cont.matrix)[3], sep="\n")) +
abline(v=c(-1,1))

volcanoplot(fit.main, coef=4, highlight = 10, names=SYMBOLS, main=paste("Differentially expressed genes", colnames(cont.matrix)[4], sep="\n")) +
abline(v=c(-1,1))


## ------------------------------------------------------------------------
require(limma)
res <- decideTests(fit.main, method="separate", adjust.method="fdr", p.value=0.1, lfc=1)
sum.res.rows <- apply(abs(res),1,sum)
res.selected <- res[sum.res.rows!=0,]
print(summary(res))


## ------------------------------------------------------------------------
vennDiagram(res.selected[,1:3], cex=0.9) +
title("Genes in common between the three comparisons\n Genes selected with FDR < 0.1 and logFC > 1")


## ------------------------------------------------------------------------
listOfTables <- list(hm6vscon6 = topTab_hm6vscon6, 
                     hm12vscon12 = topTab_hm12vscon12, 
                     hm24vscon24 = topTab_hm24vscon24,
                     hm24vshm6 = topTab_hm24vshm6)
listOfSelected <- list()
i = 1
for (i in 1:length(listOfTables)){
  # select the toptable
  topTab <- listOfTables[[i]]
  # select the genes to be included in the analysis
  whichGenes<-topTab["adj.P.Val"]<0.15
  selectedIDs <- rownames(topTab)[whichGenes]
  # convert the ID to Entrez
  EntrezIDs <- select(hgug4112a.db, selectedIDs, c("ENTREZID"))
  EntrezIDs <- EntrezIDs$ENTREZID
  listOfSelected[[i]] <- EntrezIDs
  names(listOfSelected)[i] <- names(listOfTables)[i]
}
sapply(listOfSelected, length)


## ------------------------------------------------------------------------
EntrezUni <- topAnnotated_hm24vscon24$ENTREZID


## ----echo = T, results = 'hide', message=FALSE---------------------------
require(ReactomePA)

listOfData <- listOfSelected[1:4]
comparisonsNames <- names(listOfData)

organism_ <- "human"
universe <- as.character(EntrezUni)

for (i in 1:length(listOfData)){
  data <- listOfData[[i]]
  genesIn <- listOfSelected[[i]]
  comparison <- comparisonsNames[i]
  enrich.result <- enrichPathway(gene = genesIn,
                                 pvalueCutoff = 0.05,
                                 readable = T,
                                 organism =  organism_,
                                 universe = universe,
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 pAdjustMethod = "BH")
  
  if (length(rownames(enrich.result@result)) != 0) {
  write.csv(as.data.frame(enrich.result), 
             file =paste0("./results/","ReactomePA.Results.",comparison,".csv"), 
             row.names = FALSE)
  
  pdf(file=paste0("./results/","ReactomePABarplot.",comparison,".pdf"))
    print(barplot(enrich.result, showCategory = 15, font.size = 4, 
            title = paste0("Reactome Pathway Analysis for ", comparison,". Barplot")))
  dev.off()
  
  pdf(file = paste0("./results/","ReactomePAcnetplot.",comparison,".pdf"))
    print(cnetplot(enrich.result, categorySize = "geneNum", showCategory = 15, 
         vertex.label.cex = 0.75))
  dev.off()
  }
}



## ------------------------------------------------------------------------
Tab.react.6 <- read.csv2(file.path("./results/ReactomePA.Results.hm6vscon6.csv"), sep = ",", header = TRUE, row.names = 1)
Tab.react.6.short <- Tab.react.6[1:4, 1:5]

Tab.react.12 <- read.csv2(file.path("./results/ReactomePA.Results.hm12vscon12.csv"), sep = ",", header = TRUE, row.names = 1)
Tab.react.12.short <- Tab.react.12[1:4, 1:5]

Tab.react <- read.csv2(file.path("./results/ReactomePA.Results.hm24vscon24.csv"), sep = ",", header = TRUE, row.names = 1)
Tab.react.short <- Tab.react[1:4, 1:5]

Tab.react.hm <- read.csv2(file.path("./results/ReactomePA.Results.hm24vshm6.csv"), sep = ",", header = TRUE, row.names = 1)
Tab.react.hm.short <- Tab.react.hm[1:4, 1:5]


## ----echo=FALSE----------------------------------------------------------
knitr::kable(Tab.react.6.short, booktabs=TRUE,caption="First rows and columns for Reactome results on hm6vscon6 comparison")


## ----echo=FALSE----------------------------------------------------------
knitr::kable(Tab.react.12.short, booktabs=TRUE,caption="First rows and columns for Reactome results on hm12vscon12 comparison")


## ----echo=FALSE----------------------------------------------------------
knitr::kable(Tab.react.short, booktabs=TRUE,caption="First rows and columns for Reactome results on hm24vscon24 comparison")


## ----echo=FALSE----------------------------------------------------------
knitr::kable(Tab.react.hm.short, booktabs=TRUE,caption="First rows and columns for Reactome results on hm24vshm6 comparison")

