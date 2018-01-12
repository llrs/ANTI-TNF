# Load and create PCAs

# Load the mapping file
a <- read.csv("RNAseq/mapping_file.tab", sep = "\t", check.names = FALSE)
expr <- read.csv("RNAseq/matrix.tsv", sep = "\t", check.names = FALSE) # The RNAseq
db <- read.csv("RNAseq/db_biopsies_bcn_seq16S_noTRIM.txt", sep = "\t",
               stringsAsFactors = FALSE)
otus <- read.csv("16S/OTU-tab-noTRIM.csv", check.names = FALSE, sep = ";",
                 stringsAsFactors = FALSE, row.names = 1)

tax <- otus[, ncol(otus)]
otus <- otus[, -ncol(otus)]

source("helper_functions.R")

# Extract the taxonomy and format it properly
tax <- taxonomy(tax, rownames(otus))

library("stringr")

# Filter samples with metadata
sample_names <- colnames(expr)
bcn_samples <- strsplit(sample_names[1:205], "-")
sample_names[1:205] <- sapply(bcn_samples, function(x) {
    patient <- str_pad(x[1], 3, side = "left", pad = "0")
    week <- gsub("^[w0]+([1-9]*)?$", "\\1", x[2])
    week <- paste0("w", str_pad(week, 3, side = "left", pad = "0"))
    paste(patient, week, sep = "-")

})

sample_names[206:length(sample_names)] <- gsub("(-T)?-T?TR-", "-T-DM-",
                                               sample_names[206:length(sample_names)])

# Missing samples are assumed to haven't been sequenced
# C2 wasn't sequenced
common <- intersect(db$Sample_Code, sample_names)
meta <- db[db$Sample_Code %in% common, ]
colnames(expr) <- sample_names

# We remove the samples for which we don't have the microorganism
expr <- expr[, common]
otus <- otus[, common]

# Remove low expressed genes
expr <- expr[rowSums(expr != 0) >= (0.25* ncol(expr)), ]
expr <- expr[rowMeans(expr) > quantile(rowMeans(expr), prob = 0.1), ]

# Filter by variance
SD <- apply(expr, 1, sd)
CV <- sqrt(exp(SD^2) - 1)
expr <- expr[CV > quantile(CV, probs = 0.1), ]

# PCA RNAseq
pca_i <- prcomp(t(expr), scale. = TRUE)
pca_i_x <- as.data.frame(pca_i$x)
pca_i_var <- round(summary(pca_i)$importance[2, ]*100, digits = 2)

# Filter otus
otus <- otus[rowSums(otus) != 0, ]

# PCA microbiota
pca_o <- prcomp(t(otus), scale. = TRUE)
pca_o_x <- as.data.frame(pca_o$x)
pca_o_var <- round(summary(pca_o)$importance[2, ]*100, digits = 2)

meta <- meta[match(meta$Sample_Code, colnames(expr)), ]
meta$region <- ifelse(grepl("COLON", meta$Aftected_area), "COLON",
                      meta$Aftected_area)

table(meta$Time[meta$ANTITNF_responder == 0])

# Save for other analysis
save(meta, otus, tax, expr, file = "filtered_data.RData")

ylabi <- paste("PCA2", pca_i_var[2], "%")
xlabi <- paste("PCA1", pca_i_var[1], "%")

ylabo <- paste("PCA2", pca_o_var[2], "%")
xlabo <- paste("PCA1", pca_o_var[1], "%")

pdf(paste0("Figures/", today, "_PCA.pdf"))

plot(pca_i_x$PC1, pca_i_x$PC2, xlab = xlabi, ylab = ylabi, main = "PCA of RNAseq",
     col = as.factor(meta$Aftected_area))
legend("bottomright", legend = levels(as.factor(meta$Aftected_area)),
       fill = as.numeric(as.factor(levels(as.factor(meta$Aftected_area)))))

samples <- cbind(pca_i_x, meta)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Aftected_area)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = region)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Patient_ID)) +
    guides(col = FALSE) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Time)) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = ANTITNF_responder)) +
    guides(col = guide_legend("Anti TNF responder?")) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Involved_Healthy)) +
    guides(col = guide_legend("Activity")) +
    ggtitle("PCA of RNAseq") +
    xlab(xlabi) +
    ylab(ylabi)

plot(pca_o_x$PC1, pca_o_x$PC2, xlab = xlabo, ylab = ylabo, main = "PCA of 16S",
     col = as.factor(meta$Aftected_area))
legend("bottomright", legend = levels(as.factor(meta$Aftected_area)),
       fill = as.numeric(as.factor(levels(as.factor(meta$Aftected_area)))))

samples <- cbind(pca_o_x, meta)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Aftected_area)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_text(aes(PC1, PC2, col = Aftected_area, label = Sample_Code)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

rm <- grep("^C11-.*-ILI", meta$Sample_Code)

# PCA microbiota 2
pca_o <- prcomp(t(otus[, -rm]), scale. = TRUE)
pca_o_x <- as.data.frame(pca_o$x)
pca_o_var <- round(summary(pca_o)$importance[2, ]*100, digits = 2)
samples <- cbind(pca_o_x, meta[-rm, ])
ylabo <- paste("PCA2", pca_o_var[2], "%")
xlabo <- paste("PCA1", pca_o_var[1], "%")

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Aftected_area)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = region)) +
    guides(col = guide_legend(title = "Afected Area")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Patient_ID)) +
    guides(col = FALSE) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Time)) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = ANTITNF_responder)) +
    guides(col = guide_legend("Anti TNF responder?")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

ggplot(samples) +
    geom_point(aes(PC1, PC2, col = Involved_Healthy)) +
    guides(col = guide_legend("Activity")) +
    ggtitle("PCA of microbiota") +
    xlab(xlabo) +
    ylab(ylabo)

dev.off()
