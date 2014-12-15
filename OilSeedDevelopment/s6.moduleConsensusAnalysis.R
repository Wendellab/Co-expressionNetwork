############### Step 6.  Consensus Modules analysis  ###############
# nohup R CMD BATCH s6.moduleConsensusAnalysis.R &
################


date()
getwd()
library(WGCNA);
library(flashClust);
library(RColorBrewer);
library(ggplot2);
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

remove(list=ls())
load("concensus-01-dataInput.RData")          # lanova, multiExpr
load("concensus-04-networkTopology.RData")    # netA2, netD5, netTM1, netYuc, netAD3, intrakA2, intrakD5, intrakTM1, intrakYuc, intrakAD3,
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('multiplot.r', chdir = TRUE)


checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #5
nGenes<-checkSets(multiExpr)$nGenes #20054
nSamples<-checkSets(multiExpr)$nSamples[1] #12
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3")
softPower = 20

# put all net in one list in correct order
nets<-list(A2=netA2,D5=netD5,TM1=netTM1,Yuc=netYuc,AD3=netAD3)
names(nets)==shortLabels  # TRUE make sure correct order
all.colors<-cbind(netA2$colors, netD5$colors, netTM1$colors, netYuc$colors, netAD3$colors)
# put all anova p together
dpa=as.factor(rep(seq(10,40,by=10),each=3))
anovaP<-list()
for(i in 1:nSets)
{
    MEs<-nets[[i]]$MEs
    pval<-apply(MEs,2,function(x){round(anova(aov(x~dpa) )$"Pr(>F)"[1],4)})
    pval<-as.data.frame(pval)
    pval$symbol<-ifelse(pval$pval<0.05,"*"," ")
    pval$numeric<-as.numeric(substring(rownames(pval),3) )
    pval<-pval[order(pval$numeric),]
    pval$symbol[1]<-" "  # ME0 always meaningless
    anovaP[[i]]<-pval
}
names(anovaP)<-shortLabels


# for any consensus module analysis, define new CmultiExpr and TOMs
comparisons <- list( all =1:5,          ## all five consensus
                     diploid = c(1,2),  ## diploid consensus, A2 and D5
                     trio.TM1 =c(1,2,3),
                     trio.Yuc =c(1,2,4),
                     trio.AD3 =c(1,2,5),
                     AD1 =c(3,4),       ## G. hirsutum
                     polyploids =c(3,4,5) ) ## all three polyploids

# RUN MY HUGE LOOP FOR ALL THESE COMPARISONS
for(li in 1:7){
comp<-comparisons[[li]]
print(paste("Consensus module analysis for ",names(comparisons)[li], ": ",sep=""))
print(shortLabels[comp])


######### Start of my consensus analysis unit ###########
nSets<-length(comp)
CmultiExpr = vector(mode = "list", length = nSets)
for(i in 1:nSets){ CmultiExpr[[i]] = multiExpr[[comp[i]]]}
fileName<-paste(shortLabels[comp],sep="",collapse="_")

print("Start to read in TOMs.")
# Initialize an appropriate array to hold the TOMs
TOMs = array(0, dim = c(nSets, nGenes, nGenes));
# load and store TOMs from each individual data set
for (set in 1:nSets)
{
 #  load(paste("/Volumes/jfw-lab/home/jing/oilseedNetwork/concensusTotal/",shortLabels[comp[set]], "_power20_TOM-block.1.RData", sep="") )   #local use
    load(paste("~/jfw-lab/home/jing/oilseedNetwork/concensusTotal/",shortLabels[comp[set]], "_power20_TOM-block.1.RData", sep="") )   #biocrunch
    TOMs[set, , ]<-as.matrix(TOM)
}


# Step-by-step consensus network construction
# consensus is defined as the component-wise minimum of the multiple TOMs.

# first, TOMs need to be scaled to be comparable, using the 95th percentile
# BUT BE CAREFUL, for each individual network, original TOMs and scaled TOMs might generate different numbers of modules!!!!
print("Scale TOMs using 95th percentile: values, scaling powers")
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# prepare the TOM scaled
TOMscaled <-TOMs
# Loop over sets
for (set in 1:nSets)
{
    # Select the sampled TOM entries
    TOMScalingSamples[[set]] = as.dist(TOMs[set, , ])[scaleSample]
    # Calculate the 95th percentile
    scaleQuant[set] = quantile(TOMScalingSamples[[set]], probs = scaleP, type = 8)
    # Scale the other TOMs
    if (set>1)
    {
        scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
        TOMscaled[set, ,] = TOMs[set, ,]^scalePowers[set];
    }
}
# check the scaling achieved using a quantile-quantile plot
scaledTOMSamples = list();
for (set in 1:nSets)
scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window
pdf(paste("checkTOMscaling-",fileName,".pdf",sep=""))
pwset<-combn(nSets,2)
for(i in 1:choose(nSets,2))
{
    # qq plot of the unscaled samples
    qqUnscaled = qqplot(TOMScalingSamples[[pwset[1,i]]], TOMScalingSamples[[pwset[2,i]]], plot.it = TRUE, cex = 0.6, xlab = paste("TOM in", setLabels[comp[pwset[1,i]]]), ylab = paste("TOM in", setLabels[comp[pwset[2,i]]]), main = "Q-Q plot of TOM", pch = 20)
    # qq plot of the scaled samples
    qqScaled = qqplot(scaledTOMSamples[[pwset[1,i]]], scaledTOMSamples[[pwset[2,i]]], plot.it = FALSE)
    points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
    abline(a=0, b=1, col = "blue")
    legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
}
dev.off()
# the red scaled quantile plot should be closer than black unscaled plot to the theoretic blue line
print(scaleQuant)   # 0.2002059 0.2124997 0.2477379 0.2039054 0.2607832
print(scalePowers)  # 1.000000 1.038477 1.152664 1.011515 1.196674

# Second, calculate the consensus Topological Overlap by taking the component-wise ("parallel") minimum of the TOMs in individual sets.
# Between diploids - diploid consensus network vs A2 and D5 each
if(nSets==2) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ])
if(nSets==3) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ])
if(nSets==4) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ])
if(nSets==5) consensusTOM = pmin(TOMscaled[1, , ], TOMscaled[2, , ], TOMscaled[3, , ], TOMscaled[4, , ], TOMscaled[5, , ])
# Clustering
consTree = flashClust(as.dist(1-consensusTOM), method = "average")
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = as.matrix(1-consensusTOM), deepSplit = 2, cutHeight = 0.995, minClusterSize = minModuleSize, pamRespectsDendro = TRUE )
unmergedColors = labels2colors(unmergedLabels)
# a quick summary of the module detection,
print(paste("Before merging, ", length(unique(unmergedLabels)), " modules were resulted.",sep=""))
# Calculate module eigengenes
unmergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs);
# Cluster consensus modules
consMETree = flashClust(as.dist(consMEDiss), method = "average");
# merge modules
merge = mergeCloseModules(CmultiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
# Numeric module labels
moduleLabels = merge$colors;
# the corresponding colors for large module ID need to be adjusted, otherwise they cannot be plotted
sortedModules = sort(unique(moduleLabels))
moduleLabels.adjusted = match(moduleLabels,sortedModules)-1
# Convert labels to colors
moduleColors = labels2colors(moduleLabels.adjusted)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs;
# Calculate new module eigengenes
mergedMEs = multiSetMEs(CmultiExpr, colors = NULL, universalColors = moduleLabels)
# Calculate consensus dissimilarity of consensus module eigengenes, Cluster consensus modules
consMETree.merged = flashClust(as.dist(consensusMEDissimilarity(mergedMEs) ), method = "average");

# Save useful info for the consensus modules
save(consMEs, unmergedLabels, moduleColors, moduleLabels, moduleLabels.adjusted, consTree, file = paste("s6.moduleConsensus.", fileName, ".RData", sep="") )
print("Save consensus modules information and start to draw figures.")


# Plot the consensus module results
pdf(paste("consensus-",fileName,".pdf",sep=""))
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes before merging",xlab = "", sub = "")
abline(h=0.25, col = "red")
plot(consMETree.merged, main = "Consensus clustering of consensus module eigengenes after merging",xlab = "", sub = "")
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors), c("Unmerged", "Merged"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(consTree, cbind(moduleColors, labels2colors(all.colors[,comp]) ), c("Consensus",shortLabels[comp]), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Consensus gene dendrogram and module colors")
# Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
# Displaying module heatmap and the eigengene
# sizeGrWindow(8,7);

for(i in 1:nSets)
{   MEs<-consMEs[[i]]$data
    Nmodules<-dim(MEs)[2]
    module.names<-names(MEs)
    plots <- list()  # new empty list
    for(me in 1:Nmodules)
    {
        which.module=module.names[me]
        module.color=labels2colors(match(as.numeric(substring(which.module,3)),sortedModules)-1)
        #heatmap
        par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
        plotMat(t(scale(CmultiExpr[[i]]$data[,moduleColors==module.color ]) ), nrgcols=30,rlabels=T,rcols=module.color,
        main=paste(shortLabels[comp[i]],which.module, module.color, sep=": "), cex.main=2)
        #barplot
        par(mar=c(5, 4.2, 0, 0.7))
        barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
        ylab="eigengene expression",xlab="seed development (dpa)", names.arg=rep(c(10,20,30,40), each=3) )
        #line, anova
        df<-data.frame(ME=MEs[,me], dpa, module = which.module )
        fit<-aov(ME~dpa,df)
        dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "module"))
        plots[[me]]<- ggplot(dfc, aes(x=dpa, y=ME, group=module)) +
        geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.1) +
        geom_line(colour=module.color) + geom_point( ) +
        ggtitle(paste(which.module," ",module.color,", anova P=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
        theme(plot.title=element_text( size=8))
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }
}
dev.off()


# Plot consensus module preservation results
pdf(paste("modulePreservasion-",fileName,".pdf",sep=""),width=8, height=10)
# Recalculate consMEs to give them color names, or
# consMEs = multiSetMEs(multiExpr, universalColors = merge$colors);
#sizeGrWindow(8,10);
par(cex = 0.4)
plotEigengeneNetworks(consMEs, setLabels[comp], marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1), zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
# Characterizing consensus modules by differential expression of their corresponding eigengenes in the various time points. Red means over-expression, green under-expression; numbers in each cell give the corresponding t-test p-value. Each column corresponds to an eigengene and each row corresponds to a time point.
pwset<-combn(nSets,2)
par(mfrow=c(4,1) )
for(i in 1:choose(nSets,2) )
{
    # the pair for comparison
    MEs1<-consMEs[[pwset[1,i]]]$data
    MEs2<-consMEs[[pwset[2,i]]]$data
    diffTable <- apply(MEs1-MEs2,2, function(x){unlist(lapply(split(x,dpa),mean ))} )
    pvalTable <- matrix(0, nrow = 4, ncol = dim(MEs1)[2]);
    for(cc in 1:dim(MEs1)[2])
    {
        for(rr in 1:4)
        { pvalTable[rr,cc] = t.test(split(MEs1[,cc],dpa)[[rr]],split(MEs2[,cc],dpa)[[rr]])$p.value }
    }
    labeledHeatmap( Matrix = diffTable, xLabels = names(MEs1), yLabels = c(10,20,30,40),
    colorLabels = TRUE, colors = blueWhiteRed(50),
    textMatrix = round(as.matrix(pvalTable),2),
    cex.text = 0.7,  zlim = c(-0.8,0.8),setStdMargins = FALSE,
    main = paste("Consensus MEs differential expression: ", paste(shortLabels[comp[pwset[,i] ]],collapse=" vs " ), sep=""  ) )
}
dev.off()


# Plot marginal analysis between consensus modules and each individual dataset
pdf(paste("Consensus&marginal-",fileName,".pdf",sep=""),width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
consMEs.no<-length(unique(moduleLabels))
consMEs.module<-labels2colors(match(as.numeric(substring(names(consMEs[[1]]$data ),3) ), sortedModules)-1 )  #colors in order
# loop pairwise comparison
for(i in 1:5)  # compare to all 5 , no matter which consensus dealt with
{
    coln<-consMEs.no                  #put consensus modules in columns
    rown<-ncol(nets[[i]]$MEs )  # put each individual set modules in rows
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # color list of MEs in the color of decreasing numbers of memebers
    colModule  <- consMEs.module
    rowModule  <- labels2colors(as.numeric(names(table(nets[[i]]$colors)) ))
    # colors for each gene
    colColors  <- moduleColors
    rowColors  <- labels2colors(nets[[i]]$colors )
    # anova significance sybol
    rowP  <- anovaP[[i]]$symbol
    # Initialize tables of p-values and of the corresponding counts
    pTable = matrix(0, nrow = rown, ncol = coln);
    CountTbl = matrix(0, nrow = rown, ncol = coln);
    # Execute all pairwaise comparisons
    for (rmod in 1:rown)
    for (cmod in 1:coln)
    {
        rMembers = (rowColors == rowModule[rmod] );
        cMembers = (colColors == colModule[cmod] );
        pTable[rmod, cmod] = -log10(fisher.test(rMembers, cMembers, alternative = "greater")$p.value);
        CountTbl[rmod, cmod] = sum(rowColors == rowModule[rmod]  & colColors == colModule[cmod]  )
    }
    
    # display the p-value and counts in a color-coded table. The colors will indicate the p-value signicance
    # Truncate p values smaller than 10^{-50} to 10^{-50}
    pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)]);
    pTable[pTable>50 ] = 50 ;
    # Marginal counts (really module sizes)
    rModTotals = apply(CountTbl, 1, sum)
    cModTotals = apply(CountTbl, 2, sum)
    # Use function labeledHeatmap to produce the color-coded table with all the trimmings
    labeledHeatmap( Matrix = pTable, colorLabels = TRUE,
    xLabels = paste(" ", colModule), yLabels = paste(" ", rowModule),
    xSymbols = paste(names(consMEs[[1]]$data ),"-", colModule, ": ", cModTotals, " ", sep=""),
    ySymbols = paste(names(nets)[i], rowModule, ": ", rModTotals, rowP, sep=""),
    textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
    main = "Correspondence of dataset-specific to consensus modules ",
    cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
}
dev.off()

######### End of my consensus analysis unit ###########
}


