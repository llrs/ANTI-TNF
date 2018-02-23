# Load and create PCAs

# Load the mapping file
a <- read.csv("RNAseq/mapping_file.tab", sep = "\t", check.names = FALSE)
# The RNAseq
expr <- read.csv("RNAseq/matrix.tsv", sep = "\t", check.names = FALSE)
db <- data.table::fread(
  "RNAseq/db_biopsies_bcn_seq16S_noTRIM.txt", sep = "\t",
  stringsAsFactors = FALSE
)
otus <- read.delim("16S/OTUs-Table-BCN.tab", stringsAsFactors = FALSE,
                   row.names = 1,
                   check.names = FALSE)

tax <- otus[, ncol(otus)]
otus <- otus[, -ncol(otus)]

# PDF with the quality of the sequencing
pdf(paste0("Figures/", today, "_quality.pdf"))
# number of OTUs by sample
counts <- colSums(otus)
counts_ord <- counts[order(counts)]
barplot(counts_ord, col = ifelse(grepl("w", names(counts_ord)), "red", "black"),
        main = "Total otus", xlab = "Samples", ylab = "counts", border = NA)
abline(h = 100000, col = "green")
abline(h = mean(counts), col = "orange")
abline(h = 50000)

# Show the tendency between number of OTUs and the total amount of OTUs
o <- apply(otus, 2, function(x) sum(x != 0))
p <- data.frame(OTUs = o, Abundance = counts)
ggplot(p, aes(Abundance, OTUs)) +
    geom_point() +
    geom_smooth() +
    ggtitle("Biopsies") +
    xlim(c(0, max(p$Abundance)))

# Expression quality
counts <- colSums(expr)
counts_ord <- counts[order(counts)]
barplot(counts_ord, col = ifelse(grepl("w", names(counts_ord)), "red", "black"),
        main = "Total counts (RNAseq)", xlab = "Samples", ylab = "counts", border = NA)
abline(h = 100000, col = "green")
abline(h = mean(counts), col = "orange")
abline(h = 50000)

# RElationship between sequences and different number of genes
o <- apply(expr, 2, function(x) sum(x != 0))
p <- data.frame(Genes = o, Abundance = counts)
ggplot(p, aes(Abundance, Genes)) +
    geom_point() +
    geom_smooth() +
    ggtitle("Biopsies") +
    xlim(c(0, max(p$Abundance)))

dev.off()

# Load another metadata file
meta <- read.csv(
  "Metadata_BCN.csv", check.names = FALSE,
  na.strings = c("", "N/A")
)


colnames(meta) <- c(
  "Sample_Code", "Segmento", "Actividad", "Birth_date",
  "CDEIS parcial",
  "IBD", "Clinical_activity", "Endoscopic_Activity",
  "CDAI CD",
  "CDEIS CD", "Gender",
  "Date of diagnosis", "Date visit", "Ethnicity"
)

# Merge and format them to the same standard
db <- db[, -c(5, 6, 13, 14)]
meta <- meta[, -6]
a <- merge(db, meta, all.x = TRUE, by.x = "Sample_Code", by.y = "Sample_Code")
a <- as.data.frame(droplevels(a))
dates <- grep("date", colnames(a), value = TRUE, ignore.case = TRUE)
a[, dates] <- lapply(a[, dates], as.Date, format = "%m/%d/%Y")
a[, "TimeDiag"] <- as.numeric((a$`Date visit` - a$`Date of diagnosis`) / 365.25)
a[, "Age"] <- as.numeric((a$`Date visit` - a$Birth_date) / 365.25)
a[, "AgeDiag"] <- as.numeric((a$`Date of diagnosis` - a$Birth_date) / 365.25)

meta <- a

source("helper_functions.R")

# Extract the taxonomy and format it properly
tax <- taxonomy(tax, rownames(otus))

library("stringr")

# Filter samples with metadata and correct names to match between them
sample_names <- colnames(expr)
bcn_samples <- strsplit(sample_names[1:205], "-")
sample_names[1:205] <- sapply(bcn_samples, function(x) {
  patient <- str_pad(x[1], 3, side = "left", pad = "0")
  week <- gsub("^[w0]+([1-9]*)?$", "\\1", x[2])
  week <- paste0("w", str_pad(week, 3, side = "left", pad = "0"))
  paste(patient, week, sep = "-")
})

sample_names[206:length(sample_names)] <- gsub(
  "(-T)?-T?TR-", "-T-DM-",
  sample_names[206:length(sample_names)]
)

# Missing samples are assumed to haven't been sequenced
# C2 wasn't sequenced
colnames(expr) <- sample_names

# We remove the samples for which we don't have the microorganism
common <- intersect(colnames(otus), colnames(expr))
expr <- expr[, common]
otus <- otus[, common]
meta <- meta[meta$Sample_Code %in% common, ]
meta <- meta[match(common, meta$Sample_Code), ]
meta <- droplevels(meta)

# Normalize expression
expr_edge <- edgeR::DGEList(expr)
expr_edge <- edgeR::calcNormFactors(expr_edge, method = "TMM")

# Filter expression
expr_edge <- norm_RNAseq(expr_edge)

expr_norm <- edgeR::cpm(expr_edge, normalized.lib.sizes=TRUE, log = TRUE)
expr.centered <- expr_norm - apply(expr_norm, 1, median)

# Normalize OTUS
library("metagenomeSeq")
MR_i <- newMRexperiment(
    otus,
    featureData = AnnotatedDataFrame(as.data.frame(tax[rownames(tax), ]))
)
MR_i <- cumNorm(MR_i, metagenomeSeq::cumNormStat(MR_i))
otus <- MRcounts(MR_i, norm = TRUE, log = TRUE)

# Subset if all the rows are 0 and if sd is 0
otus <- otus[apply(otus, 1, sd) != 0, ]
otus <- otus[rowSums(otus) != 0, ]

# PCA RNAseq
pca_i <- prcomp(t(expr.centered), scale. = TRUE)
pca_i_x <- as.data.frame(pca_i$x)
pca_i_var <- round(summary(pca_i)$importance[2, ] * 100, digits = 2)

# PCA microbiota
pca_o <- prcomp(t(otus), scale. = TRUE)
pca_o_x <- as.data.frame(pca_o$x)
pca_o_var <- round(summary(pca_o)$importance[2, ] * 100, digits = 2)

meta$region <- ifelse(grepl("COLON", meta$Aftected_area), "COLON",
  meta$Aftected_area
)

# Save for other analysis
save(meta, otus, tax, expr, file = "filtered_data.RData")

ylabi <- paste("PCA2", pca_i_var[2], "%")
xlabi <- paste("PCA1", pca_i_var[1], "%")

ylabo <- paste("PCA2", pca_o_var[2], "%")
xlabo <- paste("PCA1", pca_o_var[1], "%")

pdf(paste0("Figures/", today, "_PCA.pdf"))

plot(
  pca_i_x$PC1, pca_i_x$PC2, xlab = xlabi, ylab = ylabi, main = "PCA of RNAseq",
  col = as.factor(meta$Aftected_area)
)
legend(
  "bottomright", legend = levels(as.factor(meta$Aftected_area)),
  fill = as.numeric(as.factor(levels(as.factor(meta$Aftected_area))))
)

samples <- cbind(pca_i_x, meta)

ggplot(samples) +
  geom_point(aes(PC1, PC2, col = Aftected_area)) +
  guides(col = guide_legend(title = "Afected Area")) +
  ggtitle("PCA of RNAseq") +
  xlab(xlabi) +
  ylab(ylabi)

ggplot(samples) +
  geom_point(aes(PC1, PC2, col = Segmento)) +
  guides(col = guide_legend(title = "Afected Area")) +
  ggtitle("PCA of RNAseq") +
  xlab(xlabi) +
  ylab(ylabi)

ggplot(samples) +
  geom_point(aes(PC1, PC2, col = Patient_ID)) +
  # guides(col = FALSE) +
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

plot(
  pca_o_x$PC1, pca_o_x$PC2, xlab = xlabo, ylab = ylabo, main = "PCA of 16S",
  col = as.factor(meta$Aftected_area)
)
legend(
  "bottomright", legend = levels(as.factor(meta$Aftected_area)),
  fill = as.numeric(as.factor(levels(as.factor(meta$Aftected_area))))
)

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

dev.off()
