# install DESeq2
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
# load DESeq2
library("DESeq2")
#load count data, naming the columns for the sample
countData <- read.delim("gene_counts_vivo.tab", header=F, row.names=1)
colnames(countData) <- c("WT-24hr-2", "WT-24hr-3", "WT-72hr-1", "WT-72hr-2", "WT-72hr-3", "myd88-24hr-1", "myd88-24hr-2", "myd88-72hr-1", "myd88-72hr-2", "myd88-72hr-3")
head(countData)
#separate out the samples you want to work with, in this case, I'm including wt at 24hrs and all the myd88 samples
countData24 <- cbind(countData[1:2],countData[6:10])
colnames(countData24) <- c("WT-24hr-2", "WT-24hr-3", "myd88-24hr-1", "myd88-24hr-2", "myd88-72hr-1", "myd88-72hr-2", "myd88-72hr-3")
#make a data frame to give to DESeq to tell it about your experiment, row.names being the genes of interest and condition being the genotype/time
my.design <- data.frame(
  row.names = colnames( countData24 ),
  condition = c("wt", "wt", "myd", "myd", "myd", "myd", "myd" )
)
# create a DESeqDataSet and run the analysis
dds <- DESeqDataSetFromMatrix(countData = countData24, colData = my.design, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
# look at just the genes with low p-adjusted
resOrdered <- res[order(res$padj),]
head(resOrdered)
summary(res)

# three types of data transformation, and PCA plots of each
des <- DESeqTransform(dds)
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)

plotPCA(des, intgroup=c("condition"))
plotPCA(rld, intgroup=c("condition"))
plotPCA(vsd, intgroup=c("condition"))

#heat maps of the count data, and each of the 3 transformations for the top 50 genes in terms of average raw counts across all samples (top, commented out currently) 
#and in terms of most significantly differentially exressed(bottom) 
# select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
library("RColorBrewer")
library("gplots")
select <- order(res$padj)[1:100]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          trace="none", scale="none", dendrogram = "column", margin=c(10,6))
heatmap.2(assay(des)[select,], col = hmcol,
          trace="none", scale="none", dendrogram = "column", margin=c(10,6))
heatmap.2(assay(rld)[select,], col = hmcol,
          trace="none", scale="none", dendrogram = "column", margin=c(10,6))
heatmap.2(assay(vsd)[select,], col = hmcol,
          trace="none", scale="none", dendrogram = "column", margin=c(10,6))

#write a file with all genes that pass a certain padj filter cutoff
highconf <- subset(res[order(res$log2FoldChange),], padj < 0.05)
write.table(as.data.frame(highconf), file="24hr_results.txt", sep=',', quote=FALSE)

library("biomaRt")
IDs_uniprot <- read.delim("META_dbIDs_UniProt.txt", header=F, row.names=1)
colnames(IDs_uniprot) <- c("accession", "Db")
uniProt <- useMart("unimart", dataset="uniprot")
GO_IDs <- getBM(attributes=c("accession","name","go_id","db2go_p__dm_description"),filter="accession",values=IDs_uniprot$accession,mart=uniProt)
diff_IDs <- IDs_uniprot[rownames(highconf),]
diff_IDs <- diff_IDs[complete.cases(diff_IDs),]
GO_IDs_24 <- getBM(attributes=c("accession","name","go_id","go_name"),filter="accession",values=diff_IDs$accession,mart=uniProt)
GO_IDs_24_comp <- getBM(attributes=c("accession","name","go_id","go_name","db2go_p__dm_primary_id","db2go_p__dm_description","db2go_f__dm_primary_id","db2go_f__dm_description","db2go_c__dm_primary_id","db2go_c__dm_description"),filter="accession",values=diff_IDs$accession,mart=uniProt)

write.table(as.data.frame(GO_IDs), file="GO_IDs.txt", sep='\t', quote=FALSE, row.names=FALSE)
write.table(as.data.frame(GO_IDs_24_comp), file="GO_IDs_24.txt", sep='\t', quote=FALSE, row.names=FALSE)

source("http://bioconductor.org/biocLite.R")
biocLite("goseq")
library("goseq")

gene_len <- read.delim("gene_lengths.txt", header=F)
gl <- c(t(gene_len$V2))
names(gl) <- c(t(gene_len$V1))
DEG<-read.table("DEG.txt", skip=1 )
ALL<-as.data.frame(gene_len$V1)
DEG.vector <- c(t(DEG))
ALL.vector<-c(t(ALL))
gene.vector=as.integer(ALL.vector%in%DEG.vector)
names(gene.vector)=ALL.vector 
head(gene.vector)
tail(gene.vector)

pwf <- nullp(gene.vector, bias.data = gl)
plotPWF(pwf)

GO.wall <- goseq(pwf, gene2cat = )








