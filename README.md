Co-expressionNetwork
====================

WGCNA R package
===============
[WGCNA](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/) is an R package for weighted correlation network analysis. Many examples and well documented tutorials are provided on the package website.

My goal is to generate an explicit and exact recipe how I use this tool to conduct co-expression network analysis, by providing:
* Raw data input files
* Analysis script that performs each of the analyses
* All output files/figures they generated

### General workflow
1. Basic data processing and cleaning
2. Choosing the soft-thresholding power: analysis of network topology
3. ...
4. ...

### About data analysis
I found the WGCNA package [FAQ page](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html) extremely helpful, and many key questions have been discussed.

## Normalization for RNA-Seq data
RNA-seq data need to be properly normalized, and below approaches are recommended:
* Use normalized counts (or RPKM/FPKM data) and log-transform them using log2(x+1).
* Use variance-stabilizing transformation from Deseq or EdgeR.
* If data come from different batches, check for batch effects and consider using ComBat for batch effect removal but other methods should also work.
* Check quantile scatterplots to make sure there are no systematic shifts between samples; if sample quantiles show correlations (which they usually do), quantile normalization can be used to remove this effect. 

## Sample size
At least 15 samples are recommended by WGCNA to construct robust networks. However, in our case, often 4 developmental time points with triplicates were analyzed, which gives us 12 samples to work with. Therefore, we need to be aware of several problems with small sample size, such as poorer fit to scale-free topology, and the resulted network may be too noisy to be biologically meaningful.

## Filtering genes
Genes with across-the-board low expressions tend to reflect noise and correlations based on counts that are mostly zero aren't really meaningful. These genes should be removed; for example, genes with read count < 10 in more than 90% of the samples. The actual thresholds should be based on experimental design, sequencing depth and sample counts.

“We do not recommend filtering genes by differential expression. WGCNA is designed to be an unsupervised analysis method that clusters genes based on their expression profiles. Filtering genes by differential expression will lead to a set of correlated genes that will essentially form a single (or a few highly correlated) modules. It also completely invalidates the scale-free topology assumption, so choosing soft thresholding power by scale-free topology fit will fail.” - WGCNA FAQ

## Big data problem
The size of dataset that can be analyzed in ONE step is limited by memory, SO use block-wise construction wisely.

"If the reader has access to a large workstation with more than 4 GB of memory, the parameter maxBlockSize (default 5000 probes) can be increased. A 16GB workstation should handle up to 20000 probes; a 32GB workstation should handle perhaps 30000. A 4GB standard desktop or a laptop may handle up to 8000-10000 probes, depending on operating system and ihow much memory is in use by other running program."

Cytoscape
====================
To be added...
* 1
* 2
* 3
* 






