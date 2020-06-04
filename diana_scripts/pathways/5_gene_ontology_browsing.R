# use the libraries available to see the GO of our list of genes
library(data.table)
library(org.Hs.eg.db)
# ontologySimliarity comes with data objects encapsulating the GO (Gene Ontology) annotation of genes [1]:
# gene_GO_terms, a list of character vectors of term IDs of GO terms annotating each gene, named by gene,
# GO_IC, a numeric vector containing the information content of Gene Ontology terms based on frequencies of annotation in gene_GO_terms.
# https://cran.r-project.org/web/packages/ontologySimilarity/vignettes/ontologySimilarity-GO-example.html
library(ontologyIndex)
library(ontologySimilarity)



library(AnnotationDbi)
suppressPackageStartupMessages({library(hgu95av2.db)})
ls("package:hgu95av2.db")
hgu95av2.db
columns(hgu95av2.db) # If we want to know what kinds of data are retriveable via select
keytypes(hgu95av2.db) # If we are curious about what kinds of fields we could potentiall use as keys to query
# the database, we can use the keytypes method. In a perfect world, this method will
# return values very similar to what was returned by columns, but in reality, some kinds
# of values make poor keys and so this list is often shorter.
help('select') 

k <- head(keys(hgu95av2.db,keytype="PROBEID"))
mapIds(hgu95av2.db, keys=k, column=c("GENENAME"), keytype="PROBEID")


library(GO.db)
GOIDs <- c("GO:0042254","GO:0044183")
select(GO.db, keys=GOIDs, columns="DEFINITION", keytype="GOID")
help(GOpro)




# https://github.com/hbctraining/DGE_workshop_salmon/blob/master/lessons/genomic_annotation.md 

# AnnotationHub includes ENSEMBL, UCSC, ENCODE, Broad Institute, KEGG, NIH Pathway Interaction Database
library(AnnotationHub)
library(ensembldb)

unique(ah$dataprovider) # query Data Providers
unique(ah$rdataclass) # Explore the types of Data Objects available. The rdataclass allows you to see which kinds of R objects the hub will return to you. 

ah <- AnnotationHub() # Query AnnotationHub

human_ens <- query(ah, c("Homo sapiens", "EnsDb")) #  we would like to return the Ensembl EnsDb information for Human.
# The query retrieves all hits for the EnsDb objects, and you will see that they are listed by the release number.
# the most current release for GRCh38 is Ensembl98 and AnnotationHub offers that as an option to use. 
# if you were using an older genome build like hg19/GRCh37, you would need to load the EnsDb package 
# if available for that release or you might need to build your own with ensembldb.

human_ens <- human_ens[["AH53211"]] # Extract annotations of interest, for GRCh...?

# Extract gene-level information
genes=genes(human_ens, return.type = "data.frame")
genes[which(genes$gene_id == "ENSG00000223972"), ]$description

human_orgdb <- query(ah, c("Homo sapiens", "OrgDb"))
human_orgdb <- human_ens[["AH75742"]]
annotations_orgdb <- select(human_orgdb, res_tableOE_tb$gene, c("SYMBOL", "GENENAME", "ENTREZID"), "ENSEMBL")

# https://www.bioconductor.org/packages/release/bioc/vignettes/AnnotationHub/inst/doc/AnnotationHub.html
human <- query(ah, c("GeneOntology", "Homo sapiens")) 
human <- query(ah, c("ChainFile","UCSC", "Homo sapiens"))  # length(human)= 7507
human <- human[["AH14107"]] 

# https://seandavi.github.io/ITR/AnnotationHub.html
sub_ah = query(ah, c("OrgDb", "Homo sapiens"))
sub_ah
orgdb <- query(sub_ah, "OrgDb")[[1]]
orgdb
columns(orgdb)
keytypes(orgdb)


# https://www.bioconductor.org/packages/release/bioc/vignettes/GOexpress/inst/doc/GOexpress-UsersGuide.pdf
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GOexpress")

"""
# http://homer.ucsd.edu/homer/ngs/annotation.html
# annotate the genome with Homer

BED files tab separated contain
Column1: chromosome
Column2: starting position
Column3: ending position
Column4: Unique Peak ID
Column5: not used
Column6: Strand

 annotatePeaks.pl outputfile has gene description in col 18
 i.e. annotatePeaks.pl ERpeaks.txt hg18 -gtf gencode.gtf  > outputfile.txt
 i.e. annotatePeaks.pl ERpeaks.txt hg18   > outputfile.txt
 annotatePeaks.pl <peak/BED file> <genome>   > <output file>

"""

"""
https://www.disgenet.org/disgenet2r#installation-and-first-run
gene disease heatmaps
The gene2disease function retrieves the GDAs in DisGeNET for a given gene, or a for a list of genes. 
diseases associated to a single gene
or search by disease of interest
"""
library(devtools)
library(disgenet2r)
myListOfGenes <- c( "KCNE1", "KCNE2", "KCNH1", "KCNH2", "KCNG1")
# this work but how do we get the genes names

data2 <- gene2disease(gene     = myListOfGenes,score =c(0.2, 1),verbose  = TRUE)
plot( data2,class = "Network",prop = 10)
plot( data2,class  ="Heatmap",limit  = 100 )





