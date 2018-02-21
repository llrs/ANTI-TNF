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
abundance_genus <- prop.table(genus_i, 2)
abundance_genus_f <- rowSums(abundance_genus >= abundance)

genus_i_f <- abundance_genus[abundance_genus_f != 0, ]
# Correlate
corr_genus_f <- cor(t(genus_i_f), t(expr))
saveRDS(corr_genus_f, file = "correlations.RDS")

load("sgcca.RData")

# Find outliers/important genes
comp1 <- sgcca.centroid$a$RNAseq[, 1]
outliers <- comp1 != 0
corr_genus_f <- corr_genus_f[, names(comp1)]

library("org.Hs.eg.db")
names(comp1) <- gsub("(.*)\\..*", "\\1", names(comp1))

symbol <- mapIds(
  org.Hs.eg.db, keys = names(comp1), keytype = "ENSEMBL",
  column = "SYMBOL"
)

outliers[is.na(symbol)] <- FALSE

pdf(paste0("Figures/", today, "_correlations.pdf"))
library("gplots")
a <- matrix(, ncol = ncol(corr_genus_f[, outliers]), nrow = nrow(corr_genus_f))
a[abs(corr_genus_f[, outliers]) >= 0.15] <- "*" # Significant threshold of 0.05
heatmap.2(
    corr_genus_f[, outliers], main = "Correlation heatmap all: genes-genus",
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
CD_i <- genus_i[, meta$IBD == "CD"]
CD_r <- expr[, meta$IBD == "CD"]

CD_i <- CD_i[apply(CD_i, 1, sd) != 0, ]
CD_r <- CD_r[apply(CD_r, 1, sd) != 0, ]

# Filter by abundance at 0.5%
genus_CD <- prop.table(CD_i, 2)
genus_CD_f <- rowSums(genus_CD > abundance)

genus_CD <- genus_CD[genus_CD_f != 0, ]

# Correlate
corr_genus_CD <- cor(t(genus_CD), t(CD_r))
saveRDS(corr_genus_CD, file = "correlations_IBD.RDS")

int <- intersect(names(outliers), colnames(corr_genus_CD))
outliers_CD <- outliers[int]
corr_genus_CD <- corr_genus_CD[, int]

a <- matrix(, ncol = ncol(corr_genus_CD[, outliers_CD]), nrow = nrow(corr_genus_CD))
a[abs(corr_genus_CD[, outliers_CD]) >= 0.22] <- "*" # Significant threshold of 0.05
heatmap.2(
    corr_genus_CD[, outliers_CD], main = "Correlation heatmap IBD: genes-genus",
    xlab = "Genes", ylab = "Genus", scale = "none",
    tracecol = "black", col = bluered(64), trace = "none",
    labCol = symbol[outliers_CD], margins = c(6, 9), cellnote = a,
    notecol = "black", notecex = 0.5
)

## All Controls ####
C_i <- genus_i[, meta$IBD == "CONTROL"]
C_r <- expr[, meta$IBD == "CONTROL"]

C_i <- C_i[apply(C_i, 1, sd) != 0, ]
C_r <- C_r[apply(C_r, 1, sd) != 0, ]

# Filter by abundance at 0.5%
C_i <- prop.table(C_i, 2)
C_i_f <- rowSums(C_i > abundance)

C_i <- C_i[C_i_f != 0, ]

# Correlate
corr_controls <- cor(t(C_i), t(C_r))
saveRDS(corr_controls, file = "correlations_C.RDS")

int <- intersect(names(outliers), colnames(corr_controls))
outliers_C <- outliers[int]
corr_genus_C <- corr_controls[, int]

a <- matrix(, ncol = ncol(corr_genus_C[, outliers_C]), nrow = nrow(corr_genus_C))
a[abs(corr_genus_C[, outliers_C]) >= 0.28] <- "*" # Significant threshold of 0.05
heatmap.2(
    corr_genus_C[, outliers_C], main = "Correlation heatmap controls: genes-genus",
    xlab = "Genes", ylab = "Genus", scale = "none",
    tracecol = "black", col = bluered(64), trace = "none",
    labCol = symbol[outliers_C], margins = c(6, 9), cellnote = a,
    notecol = "black", notecex = 0.5
)
