brary(WGCNA);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# import prefiltered datasets, use datExpr3
lnames = load(file = "oilseed-02-dataFilter.RData");
lnames
     
# work with individual genome, then work with different soft threshold
for (i in c("A2","D5","AD3","TM1","Yuc") )
{
     i="A2"
     # Extract total read counts for each genome
     subTraits <- datTraits[datTraits$genome == i & datTraits$sub=="T",]
     subDat    <-   datExpr3[datTraits$genome == i & datTraits$sub=="T",]

     for j in c(20, 25, 30, 40) )
     {
          j= 20  # test different softpower
          softPower = j
          
          # Network construction
          print(paste("Start building network for genome ",i," using soft thres
hold ",j,"......",sep=""))
          net = blockwiseModules(
                         # Input data
                         subDat,                         
                         # Data checking options
                         checkMissingData = TRUE,
                         
                         # Options for splitting data into blocks
                         blocks = NULL,
                         maxBlockSize = 20000,  # 5000 for 4G memory, 20000 for
 16G, 300000 for 32 G
                         
                         # Network construction arguments: correlation options,
 use default  
                         # Adjacency and topology overlap function options
                         power = j, networkType = "unsigned", TOMType = "unsign
ed", 
                         
                         # Thredhold to merge modules: a height cut of 0.25 cor
responding to correlation of 0.75
                         mergeCutHeight = 0.25,
                         
                         # others
                         minModuleSize = 30,
                         reassignThreshold = 0, 
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = paste(i,"_power",j,"_TOM",sep=""),
                         verbose = 3)
          save(net, file = paste(i,"_power",j,".RData",sep="") )
          table(net$colors)
          
          # open a graphics window
          sizeGrWindow(12, 9)
          # Convert labels to colors for plotting
          mergedColors = labels2colors(net$colors)
          # Plot the dendrogram and the module colors underneath
          pdf(paste(i,"_power",j,".pdf",sep=""))
          plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes
[[1]]], "Module colors",
                         dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, gu
ideHang = 0.05)
          dev.off()
     }
}
  
