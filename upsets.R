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
same_location <- seq_len(ncol(expr)) # Where the sequenced in the same place
bcn_samples <- strsplit(sample_names[same_location], "-")
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
sample_names[same_location] <- sapply(bcn_samples, rename)

# Missing samples are assumed to haven't been sequenced
# C2 wasn't sequenced
colnames(expr) <- sample_names

# Columns to explain the datasets
expr_nam <- sample_names[same_location]
# Remove sample from TRIM
otus_nam <- grep("(-TM?[0-9]+-)|(-S0-)", colnames(otus), invert = TRUE,
                 value = TRUE)

# Remove sample from the other dataset
nam <- unique(c(expr_nam, otus_nam))
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
    vector <- as.character(vector)
    v <- unique(vector)
    v <- v[!is.na(v)]
    m2 <- sapply(v, function(x){as.numeric(vector %in% x)})
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
if (any(rowSums(m) == 1)) {
    unmatched <- rownames(m)[rowSums(m) == 1]
    meta_bd <- meta_bd[meta_bd$muestra %in% unmatched, ]
    # Remove a duplicated (uncomplete) entry on the database
    meta_bd <- meta_bd[-4, ]
    rownames(meta_bd) <- meta_bd$muestra
    meta_bd2 <- meta_bd[unmatched, ]
    split <- strsplit(as.character(meta_bd2$muestra), "-")
    week <- sapply(split, function(x){
        if (length(x) == 2) {
            as.character(as.numeric(str_sub(x[2], 2, 4)))
        } else {
            NA
        }
    })
    meta_bd2$week <- week
    m_time <- incorporate(meta_bd2$week)
    rownames(m_time) <- unmatched
    m[unmatched, colnames(m_time)] <- m_time

    # Columns to explain the region of the sample
    region <- ifelse(meta_bd2$Segmento == "íleon", "ileum", "colon")
    names(region) <- rownames(meta_bd2)
    m_region <- incorporate(region)
    region <- unique(region)
    region <- region[!is.na(region)]
    rownames(m_region) <- unmatched
    m[unmatched, colnames(m_region)] <- m_region

    # Columns to explain the disease of the sample
    m_IBD <- incorporate(meta_bd2$IBD)
    rownames(m_IBD) <- unmatched
    m_IBD <- cbind(m_IBD, "ctrl" = ifelse(rowSums(m_IBD) == 0, 1, 0))
    m[unmatched, colnames(m_IBD)] <- m_IBD
}

# Remove the one that has less things being a control
m <- m[, -13]

# Assign the values to the controls
m$ileum[grepl("^C.*ILI$", rownames(m))] <- 1
m$colon[grepl("^C.*(CA|SIG)$", rownames(m))] <- 1
m$ctrl[grepl("^C", rownames(m))] <- 1

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
upset(m, order.by = "freq", sets = unique(c(st, time, IBD)))
upset(m, order.by = "freq", sets = unique(c(st, region, time, IBD)))
dev.off()


# Create the same kind of plot but for the patients
# Summarize the information by patient
m$Patient_ID <- as.character(sapply(strsplit(rownames(m), "-"), function(x){x[1]}))
library("data.table")
setDT(m)
cols <- colnames(m)[-ncol(m)]
m[, (cols) := lapply(.SD, "sum"), by = Patient_ID, .SDcols = cols]
m <- unique(m)
m3 <- apply(m[, , by = cols], 1:2, function(x){ifelse(x != 0, 1, 0)})
m3 <- as.data.frame(m3)
pdf("Figures/patients_dataset.pdf")
st <- c("RNAseq", "Microbiome")
barplot(table(m$RNAseq[m$RNAseq != 0]),
        main = "Samples per patient from RNASeq")
barplot(table(m$Microbiome[m$Microbiome != 0]),
        main = "Samples per patient from Microbiome")
upset(m3, order.by = "freq", sets = st)
upset(m3, order.by = "freq", sets = IBD)
upset(m3, order.by = "freq", sets = region)
upset(m3, order.by = "freq", sets = time)
upset(m3, order.by = "freq", sets = c(st, region))
upset(m3, order.by = "freq", sets = c(st, IBD))
upset(m3, order.by = "freq", sets = c(st, time))
upset(m3, order.by = "freq", sets = c(st, region, time))
upset(m3, order.by = "freq", sets = c(st, IBD, region))
upset(m3, order.by = "freq", sets = unique(c(st, time, IBD)))
upset(m3, order.by = "freq", sets = unique(c(st, region, time, IBD)))
dev.off()
