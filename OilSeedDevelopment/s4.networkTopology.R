date()
getwd()
library(WGCNA);
library(RColorBrewer);
library(scatterplot3d);
library(ggplot2);
options(stringsAsFactors = FALSE);
enableWGCNAThreads();

remove(list=ls())
load("concensus-01-dataInput.RData")
load("concensus-03-buildNetwork.RData")
source('multiscalefreeplot.r', chdir = TRUE)
source('summarySE.r', chdir = TRUE)
source('multiplot.r', chdir = TRUE)

nSets<-checkSets(multiExpr)$nSets   #5
nGenes<-checkSets(multiExpr)$nGenes #20054
nSamples<-checkSets(multiExpr)$nSamples[1] #12
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("A2 seed development", "D5 seed development" ,"TM1 seed development", "Yuc seed development", "AD3 seed development")
shortLabels = c("A2", "D5" ,"TM1", "Yuc", "AD3")
powers<-c(10,20,30)

# work with individual genome, then work with different soft threshold
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    genome <-  shortLabels[set]
    net10<-get(paste(genome,"net10",sep="") )
    net20<-get(paste(genome,"net20",sep="") )
    net30<-get(paste(genome,"net30",sep="") )
    gg<-net20$goodGenes    #get the good genes descision
    
    pdf(paste(genome, "_connectivity.pdf", sep=""))

    # Network topology and other concepts need to be explored more!!!!!!!!!!!!!!!!!!
    # The following computes the network connectivity (Connectivity)
    # pdf("A2_network_connectivity.pdf")
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=3))
    names(degree)= c("power=10","power=20","power=30")
    for(j in 1:3 ){
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        # Creating Scale Free Topology Plots (SUPPLEMENTARY FIGURE S1 in our article)
        scaleFreePlot(k, truncated=T,main= paste("beta=",powers[j]))
        degree[,j]<-k
    }
    multiScaleFreePlot(degree, paste("Connectivity distribution - ", genome, sep=""))
    
    # Plot dendrograms of power=20 with all three clustering results
    plotDendroAndColors( net10$dendrograms[[1]], main = "Power=10 dendrogram", cbind(labels2colors(net10$colors[gg]),labels2colors(net20$colors[gg]),labels2colors(net30$colors[gg])),c("power=10","power=20","power=30"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
    plotDendroAndColors( net20$dendrograms[[1]], main = "Power=20 dendrogram",  cbind(labels2colors(net10$colors[gg]),labels2colors(net20$colors[gg]),labels2colors(net30$colors[gg])),c("power=10","power=20","power=30"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
    plotDendroAndColors( net30$dendrograms[[1]], main = "Power=30 dendrogram",  cbind(labels2colors(net10$colors[gg]),labels2colors(net20$colors[gg]),labels2colors(net30$colors[gg])),c("power=10","power=20","power=30"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)

    dev.off()
}

# plot power=20 for all genomes
pdf("All_Connectivity.pdf")
for(j in 1:3)
{
    degree<-as.data.frame(matrix(NA,nrow=nGenes,ncol=5))
    names(degree)<- shortLabels
    for(set in 1:nSets)
    {
        subDat <-  multiExpr[[set]]$data
        k= softConnectivity(subDat,power=powers[j], type = "signed")
        degree[,set]<-k
    }
    multiScaleFreePlot(degree, paste("Connectivity distribution, Power=",powers[j],sep=""))
}
dev.off()


# According to above results, power=20 looks resonable for following analysis
# work with individual genome, conduct module level analysis
for (set in 1:nSets )
{
    # Extract total read counts for each genome
    subDat    <-  multiExpr[[set]]$data
    subDat   <-  apply(subDat,2,as.numeric)  # important, otherwise report error
    genome <-  shortLabels[set]
    net<-get(paste(genome,"net20",sep="") )
    gg<-net20$goodGenes    #get the good genes descision
    adjacency = adjacency(subDat, power = 20, type = "signed")
    Nmodules= dim(net$MEs)[2]
    
    load(net$TOMFiles)
    print(net$TOMFiles)
    print(paste("Number of modules in ",genome," network is ",Nmodules,sep=""))
    colorsa<- labels2colors(net$colors)
    
    # Local network property
    # The function intramodularConnectivity computes the whole network connectivity kTotal, within module connectivity (kWithin). kOut=kTotal-kWithin and kDiff=kIn-kOut=2*kIN-kTotal
    intrak=intramodularConnectivity(adjacency, colors=net$colors)
    assign(paste("intrak", genome, sep=""),intrak)
    assign(paste("net", genome, sep=""),net)
    # Boxplot connectivity
    pdf(paste(genome,"_modules.pdf",sep=""))
    par(mfrow=c(2,2))
    boxplot(intrak[,2]~net$colors, main = "Global connectivity",xlab="module",ylab="k")
    boxplot(intrak[,2]~net$colors, main = "Intramodular connectivity",xlab="module",ylab="k")
    boxplot(intrak[,3]~net$colors, main = "Out connectivity",xlab="module",ylab="k")
    boxplot(intrak[,4]~net$colors, main = "diff connectivity",xlab="module",ylab="k")
    # Study how connectivity is related to mean gene expression or variance of gene expression 
    mean1=function(x) mean(x,na.rm=T)
    var1=function(x) var(x,na.rm=T)
    meanExpr=apply( subDat,2,mean1)
    varExpr=apply( subDat,2,var1)
    plot(intrak$kTotal,meanExpr,  col=as.character(colorsa), main="Mean(Expression) vs kTotal",xlab="Connectivity")
    plot(intrak$kTotal, varExpr,  col=as.character(colorsa), main="Var(Expression) vs KTotal" ,xlab="Connectivity")
    plot(intrak$kWithin,meanExpr, col=as.character(colorsa), main="Mean(Expression) vs kWithin",xlab="Connectivity")
    plot(intrak$kWithin, varExpr, col=as.character(colorsa), main="Var(Expression) vs KWithin" ,xlab="Connectivity")
 
    # cluster coefficient measures the cliquishness of a gene, computational too much
    # cluster.coef= clusterCoef0(adjacency)
    # par(mfrow=c(1,1),mar=c(2,2,2,1))
    # plot(Connectivity, cluster.coef,col=as.character(colorh1),xlab="Connectivity",ylab="Cluster Coefficient")
    # clustering coefficient in a weighted network is roughly constant for highly connected genes inside of a given module. Across modules the clustering coefficient varies a lot.
    
    
    # Visualize the modules.
    dissTOM=as.matrix(1-TOM)
    # test this with a subset of the TOM, otherwise it will take forever
    select<-sample(dim(dissTOM)[1],2000)
    selectTOM = dissTOM[select, select];
    selectTree = flashClust(as.dist(selectTOM), method = "average")
    selectColors = labels2colors(net$colors[net$goodGenes])[select]
    # sizeGrWindow(9,9)
    # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
    # the color palette; setting the diagonal to NA also improves the clarity of the plot
    plotDiss = selectTOM^7
    diag(plotDiss) = NA
    TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, 2000 randomly selected genes")
    
    # use classical multi-dimensional scaling plots 
    # for visualizing the network. Here we chose 3 scaling dimensions
    # This also takes about 10 minutes...
    cmd1=cmdscale(as.dist(selectTOM),4)
    par(mfrow=c(2,3))
    plot(cmd1[,c(1,2)], col= as.character(selectColors) )
    plot(cmd1[,c(1,3)], col= as.character(selectColors) )
    plot(cmd1[,c(1,4)], col= as.character(selectColors) )
    plot(cmd1[,c(2,3)], col= as.character(selectColors) )
    plot(cmd1[,c(2,4)], col= as.character(selectColors) )
    plot(cmd1[,c(3,4)], col= as.character(selectColors) )
    # 3D plot, oh yeah
    par(mfrow=c(1,1))
    scatterplot3d(cmd1[,1:3], color=selectColors,  main="MDS plot, 2000 randomly selected genes",xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", zlab="Scaling Dimension 3",cex.axis=1.5,angle=320)
    
    # Eigengenes are the 1st principal component of modules in a given single dataset, which provide a summary profile for each module.
    # Displaying module heatmap and the eigengene
    # sizeGrWindow(8,7);
    MEs<-net$MEs
    plots <- list()  # new empty list
    dpa=as.factor(rep(seq(10,40,by=10),each=3))
    for(me in 0:(Nmodules-1)) {
       which.module=paste("ME",me,sep="")
       module.color=labels2colors(me)
       #heatmap
       par(mfrow=c(2,1), mar=c(0.3, 5.5, 4, 2))
       plotMat(t(scale(subDat[,net$colors==me ]) ),
               nrgcols=30,rlabels=T,rcols=module.color,
               main=paste(which.module, module.color, sep=": "), cex.main=2)
       #barplot
       par(mar=c(5, 4.2, 0, 0.7))
       barplot(MEs[,which.module], col=module.color, main="", cex.main=2,
       ylab="eigengene expression",xlab="seed development (dpa)", names.arg=rep(c(10,20,30,40), each=3) )
       #line, anova
       df<-data.frame(ME=MEs[,which.module], dpa, module = which.module )
       fit<-aov(ME~dpa,df)
       dfc<-summarySE(df, measurevar="ME", groupvars=c("dpa", "module"))
       plots[[me+1]]<- ggplot(dfc, aes(x=dpa, y=ME, group=module)) +
       geom_errorbar(aes(ymin=ME-se, ymax=ME+se), width=.1) +
       geom_line(colour=module.color) + geom_point( ) +
       ggtitle(paste(which.module," ",module.color,", p=", round(anova(fit)$"Pr(>F)"[1], 4), sep="") )+
       theme(plot.title=element_text( size=11))
    }
    for(page in 1:ceiling(Nmodules/9))
    {
        if(Nmodules>(9*page))
        {  multiplot(plotlist = plots[(9*page-8):(9*page)],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
        else
        {  multiplot(plotlist = plots[(9*page-8):Nmodules],  layout=matrix(1:9, nrow=3, byrow=TRUE) )  }
    }
       
    # Relate eigengenes to external traits or sample conditions
    dpa<-as.numeric(as.character(dpa))    
    moduleTraitCor = cor(MEs, dpa, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    # graphical representation for corelation 
    # sizeGrWindow(5,6)
    # Will display correlations and their p-values
    # pdf("A2_network20_modules.pdf")
    textMatrix = paste(signif(moduleTraitCor, 2), " (", signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    par(mfrow=c(1,1), mar = c(6, 8.5, 3, 3))
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = "dpa", yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = TRUE, colors = blueWhiteRed(50), 
               textMatrix = as.matrix(textMatrix), setStdMargins = FALSE,
               cex.text = 0.7,zlim = c(-1,1), main = paste("Module-development relationships"))
    # Another way to visualize the correlations with eigegene dendrogram and adjacency heatmap, BUT the heatmap color bar is NOT RIGHT, should be (-1, 1)
    MET = orderMEs(cbind(MEs, dpa))
    # sizeGrWindow(5,7.5);
    par(cex = 0.9)
    plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
    
    # end the big plot
    dev.off()
}

save(netA2, netD5, netTM1, netYuc, netAD3,
     intrakA2, intrakD5, intrakTM1, intrakYuc, intrakAD3,
     file = "concensus-04-networkTopology.RData")

