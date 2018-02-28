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

# Load another metadata file
meta <- read.csv(
    "bd_BCN_tnf_biopsies_26022018.csv", check.names = FALSE,
    na.strings = c("", "N/A", "n.a."), row.names = NULL)
meta_bd <- read.csv("Metadata_BCN_28022018.csv")
source("helper_functions.R")

# Extract the taxonomy and format it properly
tax <- taxonomy(tax, rownames(otus))

library("stringr")

# Filter samples with metadata and correct names to match between them
sample_names <- colnames(expr)
bcn_samples <- strsplit(sample_names[1:205], "-")
rename <- function(x) {
    if (length(x) == 2){
        patient <- str_pad(x[1], 3, side = "left", pad = "0")
        week <- gsub("^[w0]+([1-9]*)?$", "\\1", x[2])
        week <- paste0("w", str_pad(week, 3, side = "left", pad = "0"))
        paste(patient, week, sep = "-")
    } else {
        id <- paste(x, sep = "-", collapse = "-")
        gsub("(-T)?-T?TR-", "-T-DM-", id)
    }
}
sample_names[1:205] <- sapply(bcn_samples, rename)

sample_names[206:length(sample_names)] <- gsub(
    "(-T)?-T?TR-", "-T-DM-",
    sample_names[206:length(sample_names)]
)

# Missing samples are assumed to haven't been sequenced
# C2 wasn't sequenced
colnames(expr) <- sample_names

# Columns to explain the datasets
expr_nam <- colnames(expr)[1:219]
# Remove sample from the other dataset
expr_nam <- grep("[0-9]{2}-T", expr_nam, invert = TRUE, value = TRUE)
nam <- unique(c(expr_nam, colnames(otus)))
expr_l <- as.numeric(nam %in% expr_nam)
otus_l <- as.numeric(nam %in% colnames(otus))
m <- data.frame("RNAseq" = expr_l, "Microbiome" = otus_l, row.names = nam)

# Columns to explain the time point of the sample
meta$Sample_id <- sapply(strsplit(as.character(meta$Sample_id), "-"), rename)
# One sample was sequenced in two batches
meta <- meta[!duplicated(meta$Sample_id), ]
rownames(meta) <- meta$Sample_id

meta2 <- meta[nam, ]
rownames(meta2) <- nam

incorporate <- function(vector){
    v <- unique(vector)
    v <- v[!is.na(v)]
    m2 <- sapply(v, function(x){as.numeric(vector %in% v)})
    colnames(m2) <- v
    m2
}

m_time <- incorporate(meta2$week)
time <- colnames(m_time)
m <- cbind(m, m_time)

# Columns to explain the region of the sample
m_region <- incorporate(meta2$sample_location)
region <- colnames(m_region)
m <- cbind(m, m_region)

# Columns to explain the disease of the sample
m_IBD <- incorporate(meta2$IBD)
IBD <- colnames(m_IBD)
m <- cbind(m, m_IBD)

# Some missing metadata
s <- rowSums(m)
missing <- s[s == 1 & m$RNAseq == 1]
if (length(missing) != 0) {
    stop("Missing data")
}

library("UpSetR")
pdf("Figures/samples_dataset.pdf")
st <- c("RNAseq", "Microbiome")
upset(m, order.by = "freq", sets = st)
upset(m, order.by = "freq", sets = IBD)
upset(m, order.by = "freq", sets = region)
upset(m, order.by = "freq", sets = time)
upset(m, order.by = "freq", sets = c(st, region))
upset(m, order.by = "freq", sets = c(st, IBD))
upset(m, order.by = "freq", sets = c(st, time))
upset(m, order.by = "freq", sets = c(st, region, time))
upset(m, order.by = "freq", sets = c(st, IBD, region))
upset(m, order.by = "freq", sets = c(st, time, IBD))
upset(m, order.by = "freq", sets = c(st, region, time, IBD))
dev.off()


# Create the same kind of plot but for the patients
# Summarize the information by patient
