library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")

path="FILE_LOCATION"
setwd(path)
directory <-path
group=rep(c("c","c_imc"), each=5)
sampleFiles<-list.files()

sampleCondition = group

sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)
 ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                      directory = directory,
                                        design= ~condition)

dds <- DESeq(ddsHTSeq)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

res<-results( dds, contrast=c("condition","c_imc","c"))
res<-res[order(res$padj),]
head(res)

write.csv(as.data.frame(res),file='sim_condition_treated_results_deseq2.csv')
write.csv(counts(dds, normalized=T),file='normalized_counts.csv')
write.csv(counts(dds, normalized=F),file='raw_counts.csv')

plotMA(dds,ylim=c(-10,10),main='DESeq2: DE pAdjValue < 0.005')
dev.copy(png,'deseq2_MAplot.png')
dev.off()

vsd <- vst(dds, blind=FALSE)
//vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
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

library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 15 )

pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
     trace="none", dendrogram="column", 
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
     annotate_col=df )

dev.copy(png,'15-gene-distance.png')
dev.off()

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ),1000 )

pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
     trace="none", dendrogram="column", 
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
      annotate_col=df  )

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

library(regionReport)

report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
                        intgroup = c('condition'), outdir='DESeq2-example',output = 'index', theme = theme_bw())


