library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)

# Read the data into R
seqdata <- read.delim("data.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("sampleinfo.txt")
head(seqdata)
dim(seqdata)


countdata <- seqdata
# Look at the output
head(countdata)
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
head(countdata)
colnames(countdata)
countdata <- as.matrix(countdata)


barplot(countdata,names=colnames(countdata))
# Add a title to the plot
title("Barplot of library sizes")
logcounts <- countdata
logcounts <- as.matrix(logcounts)
class(logcounts) <- "numeric"
logcounts <- logcounts[,-1]
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts")
abline(h=median(logcounts),col="red")
title("Boxplots of miRNAs across GBM and Normal controls")

#ggfortify#
library(ggplot2)
library(rlang)
library(devtools)
library(digest)
if (!require("devtools")) install.packages("devtools")
install_github('sinhrks/ggfortify')
library(ggfortify)

###loading data
pca1<- prcomp(clinicaltoGBMNor[,-c(1,2)])
as.factor(clinicaltoGBMNor$type)
autoplot(pca1, data = clinicaltoGBMNor, 
         colour = 'type', size = 5) + scale_color_brewer(palette='Set1')+ theme_bw()

# How many cell types and in what order are they stored?
levels(sampleinfo$type)
# Let's choose purple for basal and orange for luminal
col.type <- c(brewer.pal(2,"Set1"))[sampleinfo$type]
data.frame(sampleinfo$type,col.type)





##Hierarchical clustering with heatmaps
# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:100]
head(select_var)
# Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)
## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",
          main="Top 100 most variable genes across samples",ColSideColors=col.type,scale="row",
          distfun = function(x) dist(x,method = 'euclidean'))

coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$type), col = unique(col.type), lty = 1, lwd= 5, cex=.7)




# Save the heatmap
png(file="High_var_genes.heatmap.png")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 50 most variable genes across samples",ColSideColors=col.type,scale="row")
dev.off()



y <- DGEList(logcounts)
plotMDS(y)

# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,1))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$type)
plotMDS(y,col=col.type)


###############
##limma package##
#################
#Three matrix: expression matrix, group matrix, contrast matrix
#Three step: lmfit, eBayes, topTable
group <- paste(sampleinfo$type)
group <- factor(group)
##group matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) <- colnames(logcounts)
design

##contrast matrix
contrast.matrix <- makeContrasts(paste0(unique(group), collapse = "-"), levels = design)
contrast.matrix

##step1
fit <- lmFit(logcounts, design)
##step2
fit2<- contrasts.fit(fit, contrast.matrix)
fit2<- eBayes(fit2)
dim(fit2)
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
normalDEG = na.omit(tempOutput) 
#write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)
head(normalDEG)
write.csv(normalDEG, file = "normalmiR.csv")




# For the volcano plot we have to specify how many of the top genes to hightlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit2,coef=1,highlight=100)


par(mfrow=c(1,1))
with(normalDEG, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", xlim=c(-10,10)), cex = 0.1 )
with(subset(normalDEG, abs(logFC)>2 & abs(log10(adj.P.Val))>2), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
library(calibrate)
normalDEG$name <- rownames(normalDEG)
with(subset(normalDEG, abs(logFC)>2 & abs(log10(adj.P.Val))>2), textxy(logFC, -log10(adj.P.Val), labs=name, cex=.9))
abline(h = 2, v=  -2, col= 'red', lty =2)
abline(h = 2, v=  2, col= 'red', lty =2)


normal_select <- subset(normalDEG, abs(logFC)>2 & abs(log10(adj.P.Val))>2)
normalDEG_select <- normal_select$name

normalDEG_18miR <- logcounts[normalDEG_select,]
dim(normalDEG_18miR)
head(normalDEG_18miR)

heatmap.2(normalDEG_18miR,col=rev(morecols(50)),trace="none", 
          main="DE miRs GBM vs Normal",ColSideColors=col.type,scale="row", margins = c(5,9), dendrogram = "column")
coords <- locator(1) #click plot to get coordinates
legend(coords, legend = unique(sampleinfo$type), col = unique(col.type), lty = 1, lwd= 5, cex=.7)




