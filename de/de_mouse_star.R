library("DESeq2")
library("EnhancedVolcano")
library("pheatmap")
library("RColorBrewer")
library("genefilter")
library("regionReport")
library( "org.Mm.eg.db" )
library(dplyr)
library(biomaRt)

convertIDs <- function( ids, fromKey, toKey, db ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey), multiVals ="first" ) )
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


path="/Volumes/rna/mouse_aorta/count_without_outliar"
setwd(path)
directory<-path
group=rep(c("ctrl","ko"), each=10)
files<-list.files()
group=group[2:20]


#Read mappings from STAR alligner. 
ff <- list.files( path = ".", pattern = "*_Zfp*", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 2 ] ) )
ff <- gsub( "[.]/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1

(coldata <- data.frame(row.names=colnames(counts), group))
ddsHTSeq <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~group)

#Enseble to gene name mapping
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- data.frame(rownames(ddsHTSeq))
names(genes) <- "ensembl_gene_id"
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","mgi_symbol"),
                values = genes, 
                mart = mart)
merged_names <- merge(x = genes, 
                      y = symbol, 
                      by="ensembl_gene_id", all.x = TRUE)
merged_names[merged_names==""]<-NA
merged_names_unique=merged_names[!duplicated(merged_names[,c('ensembl_gene_id')]),]
merged_names_character=unlist(as.character(merged_names_unique$mgi_symbol))

rownames(ddsHTSeq)<- coalesce(convertIDs( rownames(ddsHTSeq) , "ENSEMBL", "SYMBOL", org.Mm.eg.db ), 
                              merged_names_character,
                              rownames(ddsHTSeq))

# DE
dds <- DESeq(ddsHTSeq)

keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

res<-results( dds, contrast=c("group","ko","ctrl"))
res<-res[order(res$padj),]
head(res)

de<-as.data.frame(res)
write.csv(de,file='de.csv')
write.csv(counts(dds, normalized=T),file='normalized_counts.csv')
write.csv(counts(dds, normalized=F),file='raw_counts.csv')

plotMA(dds,ylim=c(-10,10),main='DE pAdjValue < 0.1')
dev.copy(png,'DE MAplot.png')
dev.off()

vsd <- vst(dds, blind=FALSE)
#vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$group, vsd$type)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
df <- as.data.frame(colData(dds)[,c("group")])
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotate=df)
dev.copy(png,'sample-to-sample.png')
dev.off()

plotPCA(vsd, intgroup=c("group"))
dev.copy(png,'pca.png')
dev.off()

####  Since the clustering is only relevant for genes that actually carry signal, one usually carries it out only for a 
####  subset of most highly variable genes. Here, for demonstration, let us select the 35 genes with the highest variance 
####  across samples. We will work with the rlog transformed counts:

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

#### Top DE
topVarGenesDE <- rownames(head( res[order(res$padj),], 15 ))

count_de <- subset(assay(rld), rownames(rld) %in% topVarGenesDE)
pheatmap( count_de, scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )
dev.copy(png,'15-gene-distance_de.png')
dev.off()


topVarGenesDE <- rownames(head( res[order(res$padj),], 50 ))
count_de <- subset(assay(rld), rownames(rld) %in% topVarGenesDE)
pheatmap( count_de, scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )
dev.copy(png,'50-gene-distance_de.png')
dev.off()


topVarGenesDE <- rownames(head( res[order(res$padj),], 100 ))
count_de <- subset(assay(rld), rownames(rld) %in% topVarGenesDE)
pheatmap( count_de, scale="row", 
          trace="none", dendrogram="column", 
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
          annotate_col=df )
dev.copy(png,'100-gene-distance_de.png')
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
                       intgroup = c('group'), outdir='DE-Report',output = 'index', theme = theme_bw())

savehistory("de.R")
