date()
getwd()
library(WGCNA);
library(RColorBrewer)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

# load expression data
data = read.table("All ADT RPKM.txt",header=TRUE, sep="\t");
# Take a quick look at what is in the data set:
dim(data);  #37223   157
names(data);
names(data)<-gsub("assemblies.|.bam","",names(data))
#  Make each row corresponds to a gene and column to a sample or auxiliary information.
datExprT = as.data.frame(t(data[, -1]));
names(datExprT) = data[,1];
dim(datExprT<-datExprT[-grep(".A$|.D$",rownames(datExprT),perl=TRUE ), ] )
# make subsets for each genome
A2Data<-datExprT[grep("A2",rownames(datExprT) ),]
D5Data<-datExprT[grep("D5",rownames(datExprT) ),]
TM1Data<-datExprT[grep("TM1",rownames(datExprT) ),]
YucData<-datExprT[grep("Yuc",rownames(datExprT) ),]
AD3Data<-datExprT[grep("AD3",rownames(datExprT) ),]

# We work with five sets:
nSets = 5
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = A2Data);  names(multiExpr[[1]]$data) = "A2";  names(multiExpr[[1]]$data)= data[,1];
multiExpr[[2]] = list(data = D5Data);  names(multiExpr[[1]]$data) = "D5";  names(multiExpr[[1]]$data)= data[,1];
multiExpr[[3]] = list(data = TM1Data); names(multiExpr[[1]]$data) = "TM1"; names(multiExpr[[1]]$data)= data[,1];
multiExpr[[4]] = list(data = YucData); names(multiExpr[[1]]$data) = "Yuc"; names(multiExpr[[1]]$data)= data[,1];
multiExpr[[5]] = list(data = AD3Data); names(multiExpr[[1]]$data) = "AD3"; names(multiExpr[[1]]$data)= data[,1];
# Check that the data has the correct format for many functions operating on multiple sets:
checkSets(multiExpr)
# $nSets 5
# $nGenes 37223
# $nSamples 12 12 12 12 12
# $structureOK  TRUE

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK
# Excluding 14670 genes from the calculation due to too many missing samples or zero variance.
# If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data:
if (!gsg$allOK)
{
    # Print information about the removed genes:
    if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], collapse = ", ")))
    for (set in 1:nSets)
    {
        if (sum(!gsg$goodSamples[[set]]))
        printFlush(paste("In set", setLabels[set], "removing samples",paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
        # Remove the offending genes and samples
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
    }
}
checkSets(multiExpr)
# $nSets 5
# $nGenes 22553
# $nSamples 12 12 12 12 12
# $structureOK  TRUE

# Still too many genes to make networks, reduce the size by applying anova test
dpa=as.factor(rep(seq(10,40,by=10),each=3))
panova<-function(x){anova(aov(x~dpa))$"Pr(>F)"[1]}  # p-value of anova tests
lanova<-data.frame(gene=names(multiExpr[[1]]$data))
for (set in 1:nSets){ lanova[,shortLabels[set]] <- apply(multiExpr[[set]]$data, 2, panova ) }
lanova$minP<-apply(lanova[,-1],1,min )
table(lanova$A2<0.05)   # FALSE  TRUE 9946 12607
table(lanova$D5<0.05)   # FALSE  TRUE 11900 10653
table(lanova$TM1<0.05)  # FALSE  TRUE 9564 12989
table(lanova$Yuc<0.05)  # FALSE  TRUE 13430  9123
table(lanova$AD3<0.05)  # FALSE  TRUE 7502 15051
table(lanova$minP<0.05) # FALSE  TRUE 2499 20054
multiExpr0 <- multiExpr  #save previous as 0
for (set in 1:nSets)
{
    multiExpr[[set]]$data = multiExpr[[set]]$data[, lanova$minP<0.05];
}
# Update exprSize
checkSets(multiExpr)
# $nSets 5
# $nGenes 20054
# $nSamples 12 12 12 12 12
# $structureOK  TRUE

# We now cluster the samples on their Euclidean distance, separately in each set.
pdf(file = "SampleClustering0.pdf", width = 12, height = 12);
par(mfrow=c(5,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr0[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
pdf(file = "SampleClusteringS.pdf", width = 12, height = 12);
par(mfrow=c(5,1))
par(mar = c(0, 4, 2, 0))
sampleTrees = list()
for (set in 1:nSets)
{
    sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
for (set in 1:nSets)
plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
xlab="", sub="", cex = 0.7);
dev.off()
# looking good, at least no changes of topology due to filtering;

save(multiExpr, lanova, file = "concensus-01-dataInput.RData")
