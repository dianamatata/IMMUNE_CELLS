# only first time
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("topGO")
BiocManager::install("ALL")

library(topGO)
library(ALL)

data(ALL) #  GO annotations is stored in the ALL object
data(geneList) # The geneList data is based on a differential expression analysis of the ALL(Acute Lymphoblastic Leukemia)
affyLib <- paste(annotation(ALL), "db", sep = ".") # The microarray used in the experiment is the hgu95av2 from Affymetrix, as we can see from the affyLib object.
library(package = affyLib, character.only = TRUE)


####### EXAMPLE

# Once we have an object of class topGOdata we can start with the enrichment analysis. 
# We will use two types of test statistics: Fisher’s exact test which is based on gene counts, 
# and a Kolmogorov-Smirnov like test which computes enrichment based on gene scores. 
# We can use both these tests since each gene has a score (representing how differentially expressed a gene is) 
# and by the means of topDiffGenes functions the genes are categorized into differentially expressed or not differentially expressed genes. 
# All these are stored into sampleGOdata object.
# 
# The function runTest is used to apply the specified test statistic and method to the data. It has three main arguments. 
# class topGOdata, method for dealing with the GO graph structure and the test statistic, respectively.
# First, we perform a classical enrichment analysis by testing the over-representation of GO terms within the group of differentially expressed genes. 
# For the method classic each GO category is tested independently.
# 
# annFUN.db this function is intended to be used as long as the chip used by the user has an annotation package available in Bioconductor.
# annFUN.org this function is using the mappings from the ”org.XX.XX” annotation packages. 
#     Currently, the function supports the following gene identifiers: Entrez, GenBank, Alias, Ensembl, Gene Symbol, GeneName and UniGene.
# annFUN.gene2GO this function is used when the annotations are provided as a gene-to-GOs mapping.
# annFUN.GO2gene this function is used when the annotations are provided as a GO-to-genes mapping.
# annFUN.file this function will read the annotationsof the type gene2GO or GO2genes from a text file.

sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db, affyLib = affyLib)

resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

#table
allRes <- GenTable(sampleGOdata, classicFisher = resultFisher,
   classicKS = resultKS, elimKS = resultKS.elim,
   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

#plot
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

####### ME

file="/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes/genesonly_70_cluster1.txt"
geneNames=list(read.table(file, header = FALSE, sep = "", dec = "."))

geneID2GO <- readMappings(file = system.file("examples/geneid2go.map", package = "topGO"))
str(head(geneID2GO))

geneNames <- names(geneID2GO)

## linking ensembl gene ID to GO term? https://www.biostars.org/p/102088/ 

library(biomaRt)
# select mart and data set
bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)

# Get ensembl gene ids and GO terms
EG2GO <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id'))

# examine result
head(EG2GO,15)
# Remove blank entries
EG2GO <- EG2GO[EG2GO$go_id != '',]

# convert from table format to list format
geneID2GO <- by(EG2GO$go_id,
                EG2GO$ensembl_gene_id,
                function(x) as.character(x))

# examine result
head(geneID2GO)

# terms can be accessed using gene ids in various ways
geneID2GO$ENSMUSG00000098488

####

# My genes are ensembl IDs and are not taken from a microarray, so feed "topGOdata" with a gene2GO list.
# http://thread.gmane.org/gmane.science.biology.informatics.conductor/14627)
# construct that list by mapping all ensembl IDs to GO IDs using the package "biomaRt".then

GOdata <- new("topGOdata", ontology = "MF", allGenes = selectedList,
  description = "Ensembl GO enrichment", annot = annFUN.gene2GO, gene2GO =gene2GO)

GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)


# allGenes: named vector of type numeric or factor. The names attribute contains the genes identifiers. The genes listed in this object define the gene universe.
sampleGOdata <- new("topGOdata",
                    description = "Simple session", ontology = "BP",
                    allGenes = geneList, geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db , affyLib = affyLib)




