# RGCCA on BCN

load("filtered_data.RData")
source("helper_functions.R")

library("fgsea") # To test the distribution of the pathways
library("reactome.db") # To get the pathways
library("org.Hs.eg.db") # To get the names to entrez from the ensembl
library("data.table")

# Load previous data
load("sgcca.RData")
load("bootstrap.RData")

epithelium <- read.csv("epithelium.csv")
epithelium <- epithelium$Epithelium

pdf(paste0("Figures/", today, "_enrichments.pdf"))

# RNAseq ####
b <- STAB[["RNAseq"]]
d <- sgcca.centroid$a[["RNAseq"]][, 1]

# Remove duplicated if sgcca failed due to LAPACK subroutine
b <- b[!is.na(b[, 1]), ]

warning(sum(is.na(b[, 1])), " iterations failed.")

# Test if the gene is significant by comparing to how much times is different
# from 0 (because CCA tends to compensate itself)
count <- apply(b, 2, function(x) sum(x != 0, na.rm = TRUE))
freq <- count / nrow(b)

library(ggplot2)
ggplot(as.data.frame(cbind(freq, d))) +
    geom_point(aes(d, freq)) +
    xlab("Score") +
    ylab("Frequency (!= 0)") +
    ggtitle("RNAseq")

# Select those genes that are significant
significant <- names(d)[freq > 0.5]
significant <- gsub("(.*)\\..*", "\\1", significant)

loadings <- sgcca.centroid$a[["RNAseq"]]

# Convert the ids to the entrezIds
ensemblID <- rownames(loadings)
ensemblID <- gsub("(.*)\\..*", "\\1", ensemblID)
rownames(loadings) <- gsub("(.*)\\..*", "\\1", rownames(loadings))
entrezID <- mapIds(
  org.Hs.eg.db, keys = ensemblID, keytype = "ENSEMBL",
  column = "ENTREZID"
)
comp1 <- loadings[, 1]
names(comp1) <- entrezID


epitheliumE <- mapIds(
  org.Hs.eg.db, keys = as.character(epithelium),
  keytype = "SYMBOL", column = "ENTREZID"
)

epitheliumE <- unlist(epitheliumE, use.names = TRUE)
epitheliumE <- epitheliumE[!is.na(epitheliumE)]


# Extract the information of the pathways
genes2Pathways <- as.list(reactome.db::reactomeEXTID2PATHID)
pathways <- unlist(genes2Pathways, use.names = FALSE)
genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
paths2genes <- split(genes, pathways) # List of genes and the gene sets

# Subset to only human pathways
paths2genes <- paths2genes[grep("R-HSA-", names(paths2genes))]

## Compute the hypergeometric/enrichment analysis ####
library("ReactomePA")
enrich <- enrichPathway(
  gene = entrezID[significant], pvalueCutoff = 0.05,
  readable = TRUE, universe = unique(entrezID)
)
if (nrow(enrich) == 0) {
    warning("Enrichment didn't result in any pathway")
}
write.csv(as.data.frame(enrich), file = "RNAseq_enrichment.csv")

# Store the entrezid
entrezSig <- entrezID[significant]
entrezSig <- entrezSig[!is.na(entrezSig)]
paths2genes[["significant"]] <- entrezSig
paths2genes[["Epithelium"]] <- epitheliumE

## Compute the GSEA for the size effect ####
gseaSizeEffect <- fgsea(paths2genes[lengths(paths2genes) > 1],
                        comp1[comp1 != 0],
                        nperm = length(comp1))

# Get the name of the pathway
namesPaths <- select(
  reactome.db, keys = gseaSizeEffect$pathway,
  keytype = "PATHID", columns = "PATHNAME"
)
# Remove the homo sapiens part
namesPaths$PATHNAME <- gsub("Homo sapiens: (.*)", "\\1", namesPaths$PATHNAME)
# Add a column
gseaSizeEffect[, namesPaths := namesPaths$PATHNAME]
# Order the dataframe by size effect
data.table::setorder(gseaSizeEffect, -NES, padj, -size)
if (sum(gseaSizeEffect$padj < 0.05) == 0) {
    warning("GSEA didn't result in any pathway")
}
# Store the output
fwrite(gseaSizeEffect[pval < 0.05, ], file = "gsea_RNAseq_pathways.csv")

## 16S ####
b <- STAB[["16S"]]
d <- sgcca.centroid$a[["16S"]][, 1]

# Test if the gene is significant by comparing to how much times is different
# from 0 (because CCA tends to compensate itself)
count <- apply(b, 2, function(x) sum(x != 0, na.rm = TRUE))
freq <- count / nrow(b)

ggplot(as.data.frame(cbind(freq, d))) +
    geom_point(aes(d, freq)) +
    xlab("Score") +
    ylab("Frequency (!= 0)") +
    ggtitle("16S")

otus <- names(d)[freq > 0.5]

## Split the taxonomy
Taxon2Class <- as.list(as.data.frame(tax))

grouping <- split(Taxon2Class$Genus, Taxon2Class$Genus)
grouping <- sapply(grouping, names)

term2gene <- data.frame(
  "Gene" = tax[, "Genus"],
  "Term" = rownames(tax)
)
term2name <- data.frame(
  "Name" = tax[, "Genus"],
  "Term" = rownames(tax)
)

library("clusterProfiler")
enrich <- as.data.frame(enricher(
  gene = otus, universe = rownames(tax),
  minGSSize = 1, TERM2GENE = term2gene,
  TERM2NAME = term2name
))
if (nrow(enrich) == 0) {
    warning("Enrichment didn't result in any pathway")
}
write.csv(enrich, file = "Otus_genus_enrichment.csv")

# GSEA
comp1 <- sgcca.centroid$a[["16S"]][, 1]

gseaSizeEffect <- fgsea(grouping[lengths(grouping) > 1], comp1[comp1 != 0],
                        nperm = 20000)
data.table::setorder(gseaSizeEffect, -NES, padj, -size)
if (sum(gseaSizeEffect$padj < 0.05) == 0) {
    warning("GSEA didn't result in any pathway")
}
fwrite(gseaSizeEffect[pval < 0.05], file = "gsea_otus_genus.csv")
