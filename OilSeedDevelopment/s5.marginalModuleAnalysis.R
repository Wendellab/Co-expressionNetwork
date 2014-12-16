############### Step 5.  Marginal Modules analysis  ###############
# nohup R CMD BATCH s5.marginalModuleAnalysis.R &
################

date()
getwd()
library(WGCNA);
library(flashClust);
library(RColorBrewer);
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

remove(list=ls())
load("concensus-01-dataInput.RData")          # lanova, multiExpr
load("concensus-04-networkTopology.RData")    # netA2, netD5, netTM1, netYuc, netAD3, intrakA2, intrakD5, intrakTM1, intrakYuc, intrakAD3,

checkSets(multiExpr)
nSets<-checkSets(multiExpr)$nSets   #5
nGenes<-checkSets(multiExpr)$nGenes #20054
nSamples<-checkSets(multiExpr)$nSamples[1] #12
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3")
softPower = 20

##### let us first conduct pairwise marginal module analysis
# Standard marginal module detection analysis studies data-set specific modules

pwset<-combn(nSets,2)
# put all net in one list in correct order
nets<-list(A2=netA2,D5=netD5,TM1=netTM1,Yuc=netYuc,AD3=netAD3)
names(nets)==shortLabels  # TRUE make sure correct order
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

# prepare plotting
pdf("marginalDetection.pdf",width=10,height=7)
par(mfrow=c(1,1));
par(cex = 1.0);
par(mar=c(8, 10.4, 2.7, 1)+0.3);
# loop pairwise comparison
for(i in 1:choose(nSets,2))
{
    coln<-ncol(nets[[pwset[1,i]]]$MEs )  # number of MEs
    rown<-ncol(nets[[pwset[2,i]]]$MEs )
    # color list of MEs in the color of decreasing numbers of memebers
    colModule  <- labels2colors(as.numeric(names(table(nets[[pwset[1,i]]]$colors)) ))
    rowModule  <- labels2colors(as.numeric(names(table(nets[[pwset[2,i]]]$colors)) ))
    # colors for each gene
    colColors  <- labels2colors(nets[[pwset[1,i]]]$colors )
    rowColors  <- labels2colors(nets[[pwset[2,i]]]$colors )   # colors for each gene
    # anova significance sybol
    colP  <- anovaP[[pwset[1,i]]]$symbol
    rowP  <- anovaP[[pwset[2,i]]]$symbol
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
         xSymbols = paste(names(nets)[pwset[1,i]], colModule, ": ", cModTotals, colP, sep=""),
         ySymbols = paste(names(nets)[pwset[2,i]], rowModule, ": ", rModTotals, rowP, sep=""),
         textMatrix = CountTbl, colors = blueWhiteRed(100)[50:100],
         main = "Correspondence of modules",
         cex.text = 1, cex.lab = 1, setStdMargins = FALSE      )
}
dev.off()



