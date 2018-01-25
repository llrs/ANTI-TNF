# RGCCA on BCN
load("filtered_data.RData")
source("helper_functions.R")
source("helper_RGCCA.R")
library("RGCCA")
library("ggforce")

# Prepare meta for RGCCA
variables <- c(
  "Patient_ID", "IBD", "Aftected_area", "Treatment",
  "ANTITNF_responder", "TimeDiag", "Age",
  "Clinical_activity", "Activity", "Gender",
  "Endoscopic_Activity"
)
meta_R <- meta[, variables]

meta_R[, variables] <- sapply(meta_R[, variables], function(x) {
  as.numeric(as.factor(x))
})
meta_R[is.na(meta_R)] <- 0

# Prepare input for the sgcca function
A <- list(RNAseq = t(expr), "16S" = t(otus), "metadata" = meta_R)

# The design
C <- matrix(
  0, ncol = length(A), nrow = length(A),
  dimnames = list(names(A), names(A))
)
C <- subSymm(C, "16S", "metadata", 1)
C <- subSymm(C, "RNAseq", "metadata", 1)

# We cannnot comput eht tau.estimate for A[[1]]
# (shrinkage <- sapply(A, tau.estimate))
shrinkage <- c(0.1034826, 0, 1) # We calculated the shrinkage for the RNAseq
# expression in the server
shrinkage[2] <- tau.estimate(A[[2]])
(min_shrinkage <- sapply(A, function(x) {
  1 / sqrt(ncol(x))
}))
# # Don't let the shrinkage go below the thershold allowed
(shrinkage <- ifelse(shrinkage < min_shrinkage, min_shrinkage, shrinkage))
# shrinkage <- rep(1, length(A))

ncomp <- c(2, 2, 2)

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

# list(sgcca.centroid = sgcca.centroid, sgcca.horst = sgcca.horst,
# sgcca.factorial = sgcca.factorial)
save(sgcca.centroid, file = "sgcca.RData")

samples <- data.frame(
  "RNAseq" = sgcca.centroid$Y[["RNAseq"]][, 1],
  "microbiota" = sgcca.centroid$Y[["16S"]][, 1]
)

# Colors for the plots

samples <- cbind(samples, meta)
samples$Patient_ID <- as.factor(samples$Patient_ID)
samples$Sample_Code <- as.character(samples$Sample_Code)

pdf(paste0("Figures/", today, "_RGCCA_plots.pdf"))


samples$Time <- factor(samples$Time, levels(as.factor(samples$Time))[c(4, 1, 2, 3)])
for (p in seq_along(levels(samples$Time))) {
  a <- ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color = Patient_ID, label = Patient_ID)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by time")) +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = FALSE) +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap_paginate(~Time, ncol = 1, nrow = 1, page = p)
  print(a)
}

for (p in seq_along(levels(samples$ID))) {
  a <- ggplot(samples, aes(RNAseq, microbiota)) +
    geom_text(aes(color = Patient_ID, label = Patient_ID)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    ggtitle(paste0("Samples by patient")) +
    xlab("RNAseq (component 1)") +
    ylab("16S (component 1)") +
    guides(col = FALSE) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = colors) +
    facet_wrap_paginate(~ID, ncol = 1, nrow = 1, page = p)
  print(a)
}

ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(
    color = Patient_ID,
    label = Sample_Code
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  # xlab("RNAseq (component 1)") +
  # ylab("16S (component 1)") +
  guides(col = FALSE) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(
    color = Time,
    label = Patient_ID
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Time")) +
  theme(plot.title = element_text(hjust = 0.5))


ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(
    color = ANTITNF_responder,
    label = Patient_ID
  )) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Responders")) +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(samples, aes(RNAseq, microbiota)) +
  geom_text(aes(color = Aftected_area, label = Patient_ID)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  ggtitle("All samples at all times ") +
  xlab("RNAseq (component 1)") +
  ylab("16S (component 1)") +
  guides(col = guide_legend(title = "Afected area")) +
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
if (any(subVariables$Origin == "RNAseq")) {
  pr <- prcomp(t(expr[subVariables$var[subVariables$Origin == "RNAseq"], ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta$ANTITNF_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"] * 100, "%")) +
    ylab(paste("PC2", prS$importance[2, "PC2"] * 100, "%")) +
    guides(col = guide_legend("Anti TNF responder")) +
    ggtitle("RNAseq PCA from the important variables")
}
if (any(subVariables$Origin == "16S")) {
  pr <- prcomp(t(otus[subVariables$var[subVariables$Origin == "16S"], ]), scale. = TRUE)
  prS <- summary(pr)
  ggplot(as.data.frame(pr$x), aes(PC1, PC2, color = as.factor(meta$ANTITNF_responder))) +
    geom_point() +
    xlab(paste("PC1", prS$importance[2, "PC1"] * 100, "%")) +
    ylab(paste("PC2", prS$importance[2, "PC2"] * 100, "%")) +
    guides(col = guide_legend("Anti TNF responder")) +
    ggtitle("16S PCA from the important variables")
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

# stop("Control Flow")
# To calculate the conficence interval on selecting the variable
# this interval should reduce as we fit a better model/relationship
STAB <- boot_sgcca(A, C, shrinkage, 1000)
save(STAB, file = "bootstrap.RData")

boot_evaluate(STAB)

dev.off()
