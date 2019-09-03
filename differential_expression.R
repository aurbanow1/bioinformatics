library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("regionReport")

path="PATH"
setwd(path)
directory<-path
group=rep(c("zm","ctrl"), each=3)
files<-list.files()

table <- data.frame(sampleName = files, fileName = files, condition = group)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table,
                                       directory = directory,
                                       design= ~condition)

dds <- DESeq(ddsHTSeq)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

res<-results( dds, contrast=c("condition","zm","ctrl"))
res<-res[order(res$padj),]
head(res)

de<-as.data.frame(res)
write.csv(de,file='de.csv')
write.csv(counts(dds, normalized=T),file='normalized_counts.csv')
write.csv(counts(dds, normalized=F),file='raw_counts.csv')

plotMA(dds,ylim=c(-10,10),main='DE pAdjValue < 0.005')
dev.copy(png,'DE MAplot.png')
dev.off()

vsd <- vst(dds, blind=FALSE)
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df <- as.data.frame(colData(dds)[,c("condition")])
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotate=df)
dev.copy(png,'sample-to-sample.png')
dev.off()

plotPCA(vsd, intgroup=c("condition"))
dev.copy(png,'pca.png')
dev.off()

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 15 )

pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )

dev.copy(png,'15-gene-distance.png')
dev.off()

count_de<-dim(de[de[5]<0.05,])[1]
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ),count_de )

pheatmap(assay(rld)[ topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )

dev.copy(png,'all-gene-distance.png')
dev.off()

EnhancedVolcano(res,
                lab = rownames(res),
                x = "log2FoldChange",
                y = "pvalue",
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)

dev.copy(png,'volcano_plot.png')
dev.off()

report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
                       intgroup = c('condition'), outdir='DE-Report',output = 'index', theme = theme_bw())
