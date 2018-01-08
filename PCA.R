# Load and create PCAs

# Load the mapping file
a <- read.csv("RNAseq/mapping_file.tab", sep = "\t", check.names = FALSE)
expr <- read.csv("RNAseq/matrix.tsv", sep = "\t", check.names = FALSE) # The RNAseq
db <- read.csv("RNAseq/db_biopsies_bcn_seq16S_noTRIM.txt", sep = "\t")

# Filter samples with metadata
sample_names <- gsub("^([0-9]+)", "00\\1", colnames(expr))
bcn_samples <- strsplit(sample_names[1:205], "-")
sapply(bcn_samples, function(x){
    paste0(substr(x[1], nchar(x[1])-3), "-",
           gsub("^w", "w0", x[2])) # TODO Continue here
})



# Remove low expressed genes
expr <- expr[rowSums(expr != 0) >= (0.25* ncol(expr)), ]
expr <- expr[rowMeans(expr) > quantile(rowMeans(expr), prob = 0.1), ]

# Filter by variance
SD <- apply(expr, 1, sd)
CV <- sqrt(exp(SD^2) - 1)
expr <- expr[CV > quantile(CV, probs = 0.1), ]

pca_i <- prcomp(t(expr), scale. = TRUE)
pca_i_x <- as.data.frame(pca_i$x)
pca_i_var <- round(summary(pca_i)$importance[2, ]*100, digits = 2)

pdf(paste0("Figures/", today, "PCA.pdf"))
plot(pca_i_x$PC1, pca_i_x$PC2, xlab = paste("PCA1", pca_i_var[1], "%"),
     ylab = paste("PCA1", pca_i_var[2], "%"), main = "PCA of RNAseq")





dev.off()
