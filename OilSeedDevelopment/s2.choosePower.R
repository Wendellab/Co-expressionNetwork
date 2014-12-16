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


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=40, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for(set in 1:nSets){
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers, verbose = 2)[[2]])      }
collectGarbage()

# Plot the results:
colors=brewer.pal(5,"Set1")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
    for (col in 1:length(plotCols))
    {
        ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
        ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf("ChooseSoftThresholdPower.pdf")
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
    if (set==1)
    {
        plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
        main = colNames[col]);
        addGrid();
    }
    if (col==1)
    {
        text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
        labels=powers,cex=cex1,col=colors[set]);
    } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
    labels=powers,cex=cex1,col=colors[set]);
    if (col==1)
    {
        legend("bottomright", legend = shortLabels, col = colors, pch = 20) ;
    } else
    legend("topright", legend = shortLabels, col = colors, pch = 20) ;
}
dev.off()

save(powerTables, file = "concensus-02-choosePower.RData")
