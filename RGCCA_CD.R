library("fgsea")

# Load the helper file
source("helper_functions.R")
source("helper_RGCCA.R")
load("filtered_data.RData", verbose = TRUE)
library("ggforce")

# Prepare meta for RGCCA
# Filter only the IBD patients
meta <- meta[meta$IBD == "CD", ]
meta <- droplevels(meta)

otus <- otus[, colnames(otus) %in% meta$Sample_Code]
expr <- expr[, colnames(expr) %in% meta$Sample_Code]

stopifnot(colnames(otus) == meta$Sample_Code)
stopifnot(colnames(expr) == meta$Sample_Code)

# Subset if all the rows are 0 and if sd is 0
otus <- otus[apply(otus, 1, sd) != 0, ]
# Filter expression
expr <- norm_RNAseq(expr)

# Select the features of metadata Time and Age_sample isn't the same?? perhaps removing them

variables <- c(
    "Patient_ID", "IBD", "Aftected_area",
    "ANTITNF_responder", "TimeDiag", "Age",
    "Clinical_activity", "Activity", "Gender",
    "Endoscopic_Activity"
)
meta_R <- meta[, variables]

meta_R[, variables] <- sapply(meta_R[, variables], function(x) {
    if (is.character(x) | is.factor(x)) {
        as.numeric(as.factor(x))
    } else if (is.numeric(x)) {
        x
    }
})
meta_R[is.na(meta_R)] <- 0


# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus), "metadata" = meta_R)
A <- sapply(A, function(x){
    x[, apply(x, 2, sd) != 0]
}, simplify = FALSE)

saveRDS(A, file = "BCN_CD.RDS")

# The design
C <- matrix(
    0, ncol = length(A), nrow = length(A),
    dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "metadata", 1)
C <- subSymm(C, "RNAseq", "metadata", 1)


# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.23390092, 0, 1) # We guess a 0.1 for the RNAseq expression
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
    1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))

ncomp <- rep(2, length(A))

sgcca.centroid <- sgcca(
    A, C, c1 = shrinkage,
    ncomp = ncomp,
    scheme = "centroid",
    scale = TRUE,
    verbose = FALSE
)
names(sgcca.centroid$Y) <- names(A)
names(sgcca.centroid$a) <- names(A)
names(sgcca.centroid$astar) <- names(A)
names(sgcca.centroid$AVE$AVE_X) <- names(A)
sgcca.centroid$AVE$AVE_X <- simplify2array(sgcca.centroid$AVE$AVE_X)

sgcca.factorial <- sgcca(
    A, C, c1 = shrinkage,
    ncomp = ncomp,
    scheme = "factorial",
    scale = TRUE,
    verbose = FALSE
)
names(sgcca.factorial$Y) <- names(A)
names(sgcca.factorial$a) <- names(A)
names(sgcca.factorial$astar) <- names(A)
names(sgcca.factorial$AVE$AVE_X) <- names(A)
sgcca.factorial$AVE$AVE_X <- simplify2array(sgcca.factorial$AVE$AVE_X)

sgcca.horst <- sgcca(
    A, C, c1 = shrinkage,
    ncomp = ncomp,
    scheme = "horst",
    scale = TRUE,
    verbose = FALSE
)
names(sgcca.horst$Y) <- names(A)
names(sgcca.horst$a) <- names(A)
names(sgcca.horst$astar) <- names(A)
names(sgcca.horst$AVE$AVE_X) <- names(A)
sgcca.horst$AVE$AVE_X <- simplify2array(sgcca.factorial$AVE$AVE_X)

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)
save(sgcca.centroid, file = "CD.RData")

samples <- data.frame(
    "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
    "microbiota" = sgcca.centroid$Y[["16S"]][, 1]
)


## Grouping of the variables ####
RNAseq1 <- samples$RNAseq
RNAseq2 <- sgcca.centroid$Y[["RNAseq"]][, 2]
microbiota2 <- sgcca.centroid$Y[["16S"]][, 2]
microbiota1 <- samples$microbiota

names(RNAseq1) <- rownames(samples)
names(microbiota1) <- rownames(samples)
names(RNAseq2) <- rownames(samples)
names(microbiota2) <- rownames(samples)
groups <- split(rownames(samples), as.factor(meta$ANTITNF_responder))
# First dimension seems to capture well the
fgsea(groups, RNAseq1, nperm = 1000)
fgsea(groups, microbiota1, nperm = 1000)
# Further dimensions
fgsea(groups, RNAseq2, nperm = 1000)
fgsea(groups, microbiota2, nperm = 1000)

km <- kmeans(samples[, c("RNAseq", "microbiota")], 2, nstart = 2)
plot(samples[, c("RNAseq", "microbiota")], col = km$cluster,
     main = "kmeans clustering\nTrying to predict responders")

## Plotting results ####
# Colors for the plots
samples <- cbind(samples, droplevels(meta))
samples$Patient_ID <- as.factor(samples$Patient_ID)
samples$Sample_Code <- as.character(samples$Sample_Code)

pdf(paste0("Figures/", today, "_RGCCA_plots_IBD.pdf"))

# Labels of the samples
samples$Time <- as.factor(samples$Time)
for (p in seq_along(levels(samples$Time))) {
    a <- ggplot(samples, aes(RNAseq, microbiota)) +
        geom_text(aes(label = Sample_Code)) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        ggtitle(paste0("Samples by time")) +
        xlab("RNAseq (component 1)") +
        ylab("16S (component 1)") +
        guides(col = guide_legend(title = "Patient")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = colors) +
        facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
    print(a)
}

for (p in seq_along(levels(samples$Patient_ID))) {
    a <- ggplot(samples, aes(RNAseq, microbiota)) +
        geom_text(aes(color = Patient_ID, label = Time)) +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0) +
        ggtitle(paste0("Samples by patient")) +
        xlab("RNAseq (component 1)") +
        ylab("16S (component 1)") +
        guides(col = guide_legend(title = "Patient")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = colors) +
        facet_wrap_paginate(~Patient_ID, ncol = 1, nrow = 1, page = p)
    print(a)
}
ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color = Patient_ID, label = Time)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle("All samples at all times ") +
    # xlab("RNAseq (component 1)") +
    # ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Patient")) +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(
        color = ANTITNF_responder,
        label = Sample_Code)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle("All samples at all times ") +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Responders")) +
    theme(plot.title = element_text(hjust = 0.5))


ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(
        color = Endoscopic_Activity,
        label =  Sample_Code)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle("All samples at all times ") +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Endoscopic Activity")) +
    theme(plot.title = element_text(hjust = 0.5))

ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color = Time, label = Sample_Code)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle("All samples at all times ") +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = guide_legend(title = "Time")) +
    theme(plot.title = element_text(hjust = 0.5))

variables <- data.frame(
    Origin = rep(names(A), sapply(A, ncol)),
    comp1 = unlist(sapply(
        sgcca.centroid$a,
        function(x) {
            x[, 1]
        }
    )),
    comp2 = unlist(sapply(
        sgcca.centroid$a,
        function(x) {
            x[, 2]
        }
    ))
)
variables$var <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(OTU_.*)$", "\\1", rownames(variables))
variables$var <- gsub("^RNAseq\\.(ENSG.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^.*\\.(ENSG.*)$", "\\1", rownames(variables))
rownames(variables) <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))
variables$var <- gsub("^metadata\\.(.*)$", "\\1", rownames(variables))

# Remove the variables that in both components are 0
keepComp1RNAseq <- mean(abs(variables$comp1)[variables$Origin == "RNAseq"])
keepComp1_16S <- mean(abs(variables$comp1)[variables$Origin == "16S"])
keepComp1_metadata <- mean(abs(variables$comp1)[variables$Origin == "metadata"])

keepComp2RNAseq <- mean(abs(variables$comp2)[variables$Origin == "RNAseq"])
keepComp2_16S <- mean(abs(variables$comp2)[variables$Origin == "16S"])
keepComp2_metadata <- mean(abs(variables$comp2)[variables$Origin == "metadata"])

keepComp1 <- c(
    abs(variables$comp1[variables$Origin == "RNAseq"]) > keepComp1RNAseq,
    abs(variables$comp1[variables$Origin == "16S"]) > keepComp1_16S,
    abs(variables$comp1[variables$Origin == "metadata"]) > keepComp1_metadata
)
keepComp2 <- c(
    abs(variables$comp2[variables$Origin == "RNAseq"]) > keepComp2RNAseq,
    abs(variables$comp2[variables$Origin == "16S"]) > keepComp2_16S,
    abs(variables$comp2[variables$Origin == "metadata"]) > keepComp2_metadata
)

subVariables <- variables[keepComp1 & keepComp2, ]

ggplot(subVariables, aes(comp1, comp2), color = Origin) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.1, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.2, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.3, npoints = 100)) +
    geom_path(aes(x, y), data = circleFun(c(0, 0), 0.4, npoints = 100)) +
    geom_text(aes(color = Origin, label = var)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    coord_cartesian() +
    ggtitle(
        "Variables important for the first two components",
        subtitle = "Integrating stools and mucosa samples"
    )

rnaseq_i <- subVariables$var[subVariables$Origin == "RNAseq"]
if (length(rnaseq_i) >= 2) {
    pr <- prcomp(t(expr[rnaseq_i, ]), scale. = TRUE)
    prS <- summary(pr)
    ggplot(as.data.frame(pr$x),
           aes(PC1, PC2,
               color = as.factor(meta$ANTITNF_responder))) +
        geom_point() +
        xlab(paste("PC1", prS$importance[2, "PC1"] * 100, "%")) +
        ylab(paste("PC2", prS$importance[2, "PC2"] * 100, "%")) +
        ggtitle("RNAseq PCA from the important variables") +
        guides(col = guide_legend(title = "Responders"))
}

micro_i <- subVariables$var[subVariables$Origin == "16S"]
if (length(micro_i) >= 2) {
    pr <- prcomp(t(otus[micro_i, ]), scale. = TRUE)
    prS <- summary(pr)
    ggplot(as.data.frame(pr$x),
           aes(PC1, PC2,
               color = as.factor(meta$ANTITNF_responder))) +
        geom_point() +
        xlab(paste("PC1", prS$importance[2, "PC1"] * 100, "%")) +
        ylab(paste("PC2", prS$importance[2, "PC2"] * 100, "%")) +
        ggtitle("16S PCA from the important variables") +
        guides(col = guide_legend(title = "Responders"))
}
# Plot for the same component the variables of each block
comp1 <- sapply(sgcca.centroid$a, function(x) {
    x[, 1]
})
variables_weight(comp1)

# Second component
comp2 <- sapply(sgcca.centroid$a, function(x) {
    x[, 2]
})
variables_weight(comp2)

# Bootstrap of sgcca
STAB <- boot_sgcca(A, C, shrinkage, 1000)

save(STAB, file = "bootstrap_CD.RData")

# Evaluate the boostrap effect and plot
boot_evaluate(STAB)

dev.off()
