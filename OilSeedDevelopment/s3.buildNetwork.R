date()
library(WGCNA)
library(RColorBrewer)
library(ggplot2);
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
lnames = load(file = "concensus-01-dataInput.RData")
lnames
nSets<-checkSets(multiExpr)$nSets   #5
nGenes<-checkSets(multiExpr)$nGenes #20054
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3")

# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error    
    for (j in c(10, 20, 30) )
    {
        softPower = j
        print(paste("Start building network for ",shortLabels[set]," using soft threshold ",j,"......",sep=""))
        # Network construction
        
        net = blockwiseModules(
             # Input data
             subDat,
             # Data checking options
             checkMissingData = TRUE,
             
             # Options for splitting data into blocks
             blocks = NULL,
             randomSeed = 12345,
             maxBlockSize = nGenes,  # 5000 for 4G memory, 20000 for 16G, 300000 for 32 G
             
             # Network construction arguments: correlation options, use default
             # Adjacency and topology overlap function options
             power = j, networkType = "signed", TOMType = "signed",
             
             # Saving or returning TOM
             saveTOMs = TRUE,
             saveTOMFileBase = paste(shortLabels[set],"_power",j,"_TOM",sep=""),
             
             # Basic tree cut options
             minModuleSize = 30,
             
             # Thredhold to merge modules: a height cut of 0.25 corresponding to correlation of 0.75
             mergeCutHeight = 0.25,
             
             # others
             reassignThreshold = 0,
             numericLabels = TRUE,
             verbose = 3)
             
        assign(paste(shortLabels[set],"net",j,sep=""), net)
        }
}
save(list=grep(".+net.+",ls(), value=TRUE), file = "concensus-03-buildNetwork.RData")

