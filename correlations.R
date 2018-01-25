# Load the helper file
source("helper_functions.R")
load("filtered_data.RData")

# Summarize to genus
library("metagenomeSeq")

tax <- tax[rownames(otus), ]
# Create the objects to summarize data
MR_i <- newMRexperiment(
  otus,
  # phenoData = AnnotatedDataFrame(meta),
  featureData = AnnotatedDataFrame(as.data.frame(tax))
)
genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")

# Subset if all the rows are 0 and if sd is 0
genus_i <- genus_i[apply(genus_i, 1, sd) != 0, ]
expr <- expr[apply(expr, 1, sd) != 0, ]

abundance <- 0.005 # 0.5%

## All samples ####
# Filter by abundance at 0.5%
a <- prop.table(genus_i, 2)
b <- rowSums(a >= abundance)

genus_i <- genus_i[b != 0, ]
# Correlate
p <- cor(t(genus_i), t(expr))
saveRDS(p, file = "correlations.RDS")

load("sgcca.RData")

# Find outliers/important genes
comp1 <- sgcca.centroid$a$RNAseq[, 1]
outliers <- comp1 != 0
p <- p[, names(comp1)]

library("org.Hs.eg.db")
names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))

symbol <- mapIds(
  org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL",
  column = "SYMBOL"
)

outliers[is.na(symbol)] <- FALSE

pdf(paste0("Figures/", today, "_plots.pdf"))
library("gplots")
a <- matrix(, ncol = ncol(p[, outliers]), nrow = nrow(p))
a[abs(p[, outliers]) >= 0.15] <- "*" # Significant threshold of 0.05
heatmap.2(
  p[, outliers], main = "Correlation heatmap all: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)

# To see if the weight have a relation with the direct correlation
# colCol = ifelse(comp1[outliers] > 0, "green", "black")

# To have a good feeling plot some random genes. Most of the time there are
# the double of 0 than in the selected heatmap
# sam <- sample(colnames(cors), 600)
# gplots::heatmap.2(cors[, sam], main = "Correlation heatmap: genes-genus",
#                   xlab = "Genes", ylab = "Genus", scale = "none",
#                   tracecol = "black", col = bluered(64), trace = "none",
#                   labCol = symbol[sam], margins = c(6, 9))




## All IBD ####
disease_i <- genus_i[, meta$IBD == "CD"]
disease_r <- expr[, meta$IBD == "CD"]

disease_i <- disease_i[apply(disease_i, 1, sd) != 0, ]
disease_r <- disease_r[apply(disease_r, 1, sd) != 0, ]

# Filter by abundance at 0.5%
a <- prop.table(disease_i, 2)
b <- rowSums(a > abundance)

disease_i <- disease_i[b != 0, ]

# Correlate
p <- cor(t(disease_i), t(disease_r))
saveRDS(p, file = "correlations_IBD.RDS")

disease <- p[, colnames(p) %in% names(outliers)]
a <- matrix(, ncol = ncol(disease[, outliers]), nrow = nrow(p))
a[abs(disease[, outliers]) >= 0.22] <- "*" # Significant threshold of 0.05
heatmap.2(
  disease[, outliers], main = "Correlation heatmap IBD: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)

## All Controls ####
disease_i <- genus_i[, meta$IBD == "CONTROL"]
disease_r <- expr[, meta$IBD == "CONTROL"]

disease_i <- disease_i[apply(disease_i, 1, sd) != 0, ]
disease_r <- disease_r[apply(disease_r, 1, sd) != 0, ]

# Filter by abundance at 0.5%
a <- prop.table(disease_i, 2)
b <- rowSums(a > abundance)

disease_i <- disease_i[b != 0, ]

# Correlate
p <- cor(t(disease_i), t(disease_r))
saveRDS(p, file = "correlations_C.RDS")

p <- p[, colnames(p) %in% names(outliers)]
outliers <- outliers[names(outliers) %in% colnames(p)]
a <- matrix(, ncol = ncol(p[, outliers]), nrow = nrow(p))
a[abs(p[, outliers]) >= 0.28] <- "*" # Significant threshold of 0.05
heatmap.2(
  p[, outliers], main = "Correlation heatmap controls: genes-genus",
  xlab = "Genes", ylab = "Genus", scale = "none",
  tracecol = "black", col = bluered(64), trace = "none",
  labCol = symbol[outliers], margins = c(6, 9), cellnote = a,
  notecol = "black", notecex = 0.5
)
