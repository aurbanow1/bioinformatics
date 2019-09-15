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

# 1 way
rownames(ddsHTSeq)<- coalesce(convertIDs( rownames(ddsHTSeq) , "ENSEMBL", "SYMBOL", hgu95av2.db ), 
                              convertIDs(rownames(ddsHTSeq) , "ENSEMBL", "SYMBOL", hgu95av2.db ),
                              rownames(ddsHTSeq))
                              
# 2 way                              

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- data.frame(rownames(ddsHTSeq))
names(genes) <- "ensembl_gene_id"
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                values = genes, 
                mart = mart)

merged_names <- merge(x = symbol, 
                      y = genes, 
                      by.x="ensembl_gene_id",
                      by.y="ensembl_gene_id")

rownames(ddsHTSeq)<-merged_names[2]
