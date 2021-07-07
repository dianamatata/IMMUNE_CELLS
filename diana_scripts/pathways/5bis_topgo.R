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
# fix the allgenes 

library(biomaRt) # biomaRt provides an interface to a growing collection of databases. The most prominent examples of BioMart databases 
# are maintain by Ensembl
library(org.Hs.eg.db) # Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.

file="/Users/dianaavalos/Programming/IMMUNE_CELLS/diana_scripts/pathways/trans_hubs_genes/genesonly_70_cluster3.txt"
geneNames=list(read.table(file, header = FALSE, sep = "", dec = "."))

#all_genes <- rep(c(0, 1, 1), each = 991) # length(unlist(geneNames))
all_genes <- rep(1,  length(unlist(geneNames)))
names(all_genes) <- unlist(geneNames, use.names=FALSE)

GOdata <- new("topGOdata", description = "Test", ontology = "BP", allGenes = all_genes, geneSel = function(p) p < 
                0.01,  annot = annFUN.org, mapping = "org.Hs.eg.db", 
              ID = "Ensembl")


resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   classicKS = resultKS, elimKS = resultKS.elim,
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)


par(cex = 0.6)
showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 5, useInfo = 'all')







