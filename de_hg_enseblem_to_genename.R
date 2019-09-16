library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("regionReport")
library( "org.Hs.eg.db" )
library(EnsDb.Hsapiens.v75)
library(dplyr)
library(hgu95av2.db)
library(biomaRt)


convertIDs <- function( ids, fromKey, toKey, db ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey), multiVals ="first" ) )
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


path="/Volumes/MackDisk/UCLA\ -\ Lab\ /UCLA\ QCBio\ Collaboratory/ARISPE\ Lab/Endothelilal\ cells-Arispe\ Lab/de/dmso_vegf_vs_zm323_vegf_v2/"
setwd(path)
directory<-path
group=rep(c("dmso","zm323"), each=3)
files<-list.files()

table <- data.frame(sampleName = files, fileName = files, condition = group)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = table,
                                       directory = directory,
                                       design= ~condition)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- data.frame(rownames(ddsHTSeq))
names(genes) <- "ensembl_gene_id"
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                values = genes, 
                mart = mart)

merged_names <- merge(x = genes, 
                      y = symbol, 
                      by="ensembl_gene_id", all.x = TRUE)
merged_names[merged_names==""]<-NA
merged_names_unique=merged_names[!duplicated(merged_names[,c('ensembl_gene_id')]),]
merged_names_character=unlist(as.character(merged_names_unique$hgnc_symbo))

rownames(ddsHTSeq)<- coalesce(convertIDs( rownames(ddsHTSeq) , "ENSEMBL", "SYMBOL", org.Hs.eg.db ), 
                              convertIDs(rownames(ddsHTSeq) , "ENSEMBL", "SYMBOL", hgu95av2.db ),
                              merged_names_character,
                              rownames(ddsHTSeq))
  
dds <- DESeq(ddsHTSeq)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

res<-results( dds, contrast=c("condition","dmso","zm323"))
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

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 50 )
pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )
dev.copy(png,'50-gene-distance.png')
dev.off() 

topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 100 )
pheatmap( assay(rld)[ topVarGenes, ], scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )
dev.copy(png,'100-gene-distance.png')
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
