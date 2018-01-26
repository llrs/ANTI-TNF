# Summarize by taxa
library("metagenomeSeq")

load("filtered_data.RData")
source("helper_functions.R")
source("helper_prevalence.R")

# Match the name of meta and rownames
otus <- otus[, match(meta$Sample_Code, colnames(otus))]
tax <- tax[rownames(otus), ]

MR_i <- newMRexperiment(
  otus,
  # phenoData = AnnotatedDataFrame(meta),
  featureData = AnnotatedDataFrame(as.data.frame(tax))
)


genus_i <- aggTax(MR_i, lvl = "Genus", out = "matrix")
species_i <- aggTax(MR_i, lvl = "Species", out = "matrix")
family_i <- aggTax(MR_i, lvl = "Family", out = "matrix")
order_i <- aggTax(MR_i, lvl = "Order", out = "matrix")
class_i <- aggTax(MR_i, lvl = "Class", out = "matrix")
phylum_i <- aggTax(MR_i, lvl = "Phylum", out = "matrix")

rownames(meta) <- meta$Sample_Code
otus_i <- otus

## Time ####
otus <- comb_prevalence(otus_i, meta, c("Time"))
write.csv(otus, "prevalence_time_otus.csv")
genus <- comb_prevalence(genus_i, meta, c("Time"))
write.csv(genus, "prevalence_time_genus.csv")
species <- comb_prevalence(species_i, meta, c("Time"))
write.csv(species, "prevalence_time_species.csv")
family <- comb_prevalence(family_i, meta, c("Time"))
write.csv(family, "prevalence_time_family.csv")
order <- comb_prevalence(order_i, meta, c("Time"))
write.csv(order, "prevalence_time_order.csv")
class <- comb_prevalence(class_i, meta, c("Time"))
write.csv(class, "prevalence_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta, c("Time"))
write.csv(phylum, "prevalence_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_i, meta, c("ANTITNF_responder"))
write.csv(otus, "prevalence_otus.csv")
genus <- comb_prevalence(genus_i, meta, c("ANTITNF_responder"))
write.csv(genus, "prevalence_genus_i.csv")
species <- comb_prevalence(species_i, meta, c("ANTITNF_responder"))
write.csv(species, "prevalence_species.csv")
family <- comb_prevalence(family_i, meta, c("ANTITNF_responder"))
write.csv(family, "prevalence_family.csv")
order <- comb_prevalence(order_i, meta, c("ANTITNF_responder"))
write.csv(order, "prevalence_order.csv")
class <- comb_prevalence(class_i, meta, c("ANTITNF_responder"))
write.csv(class, "prevalence_class.csv")
phylum <- comb_prevalence(phylum_i, meta, c("ANTITNF_responder"))
write.csv(phylum, "prevalence_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_i, meta, c("ANTITNF_responder", "Time"))
write.csv(otus, "prevalence_response_time_otus.csv")
genus <- comb_prevalence(genus_i, meta, c("ANTITNF_responder", "Time"))
write.csv(genus, "prevalence_response_time_genus_i.csv")
species <- comb_prevalence(species_i, meta, c("ANTITNF_responder", "Time"))
write.csv(species, "prevalence_response_time_species.csv")
family <- comb_prevalence(family_i, meta, c("ANTITNF_responder", "Time"))
write.csv(family, "prevalence_response_time_family.csv")
order <- comb_prevalence(order_i, meta, c("ANTITNF_responder", "Time"))
write.csv(order, "prevalence_response_time_order.csv")
class <- comb_prevalence(class_i, meta, c("ANTITNF_responder", "Time"))
write.csv(class, "prevalence_response_time_class.csv")
phylum <- comb_prevalence(phylum_i, meta, c("ANTITNF_responder", "Time"))
write.csv(phylum, "prevalence_response_time_phylum.csv")

# COLON ######

colon <- meta$Aftected_area != "ILEUM"

## Time ####
otus <- comb_prevalence(otus_i[, colon], meta[colon, ], c("Time"))
write.csv(otus, "prevalence_colon_time_otus.csv")
genus <- comb_prevalence(genus_i[, colon], meta[colon, ], c("Time"))
write.csv(genus, "prevalence_colon_time_genus.csv")
species <- comb_prevalence(species_i[, colon], meta[colon, ], c("Time"))
write.csv(species, "prevalence_colon_time_species.csv")
family <- comb_prevalence(family_i[, colon], meta[colon, ], c("Time"))
write.csv(family, "prevalence_colon_time_family.csv")
order <- comb_prevalence(order_i[, colon], meta[colon, ], c("Time"))
write.csv(order, "prevalence_colon_time_order.csv")
class <- comb_prevalence(class_i[, colon], meta[colon, ], c("Time"))
write.csv(class, "prevalence_colon_time_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], meta[colon, ], c("Time"))
write.csv(phylum, "prevalence_colon_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(otus, "prevalence_colon_otus.csv")
genus <- comb_prevalence(genus_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(genus, "prevalence_colon_genus_i.csv")
species <- comb_prevalence(species_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(species, "prevalence_colon_species.csv")
family <- comb_prevalence(family_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(family, "prevalence_colon_family.csv")
order <- comb_prevalence(order_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(order, "prevalence_colon_order.csv")
class <- comb_prevalence(class_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(class, "prevalence_colon_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], meta[colon, ], c("ANTITNF_responder"))
write.csv(phylum, "prevalence_colon_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(otus, "prevalence_colon_response_time_otus.csv")
genus <- comb_prevalence(genus_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(genus, "prevalence_colon_response_time_genus_i.csv")
species <- comb_prevalence(species_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(species, "prevalence_colon_response_time_species.csv")
family <- comb_prevalence(family_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(family, "prevalence_colon_response_time_family.csv")
order <- comb_prevalence(order_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(order, "prevalence_colon_response_time_order.csv")
class <- comb_prevalence(class_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(class, "prevalence_colon_response_time_class.csv")
phylum <- comb_prevalence(phylum_i[, colon], meta[colon, ], c("ANTITNF_responder", "Time"))
write.csv(phylum, "prevalence_colon_response_time_phylum.csv")


# ILEUM ####

## Time ####
otus <- comb_prevalence(otus_i[, !colon], meta[!colon, ], c("Time"))
write.csv(otus, "prevalence_ileum_time_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], meta[!colon, ], c("Time"))
write.csv(genus, "prevalence_ileum_time_genus.csv")
species <- comb_prevalence(species_i[, !colon], meta[!colon, ], c("Time"))
write.csv(species, "prevalence_ileum_time_species.csv")
family <- comb_prevalence(family_i[, !colon], meta[!colon, ], c("Time"))
write.csv(family, "prevalence_ileum_time_family.csv")
order <- comb_prevalence(order_i[, !colon], meta[!colon, ], c("Time"))
write.csv(order, "prevalence_ileum_time_order.csv")
class <- comb_prevalence(class_i[, !colon], meta[!colon, ], c("Time"))
write.csv(class, "prevalence_ileum_time_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], meta[!colon, ], c("Time"))
write.csv(phylum, "prevalence_ileum_time_phylum.csv")

# Responders
otus <- comb_prevalence(otus_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(otus, "prevalence_ileum_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(genus, "prevalence_ileum_genus_i.csv")
species <- comb_prevalence(species_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(species, "prevalence_ileum_species.csv")
family <- comb_prevalence(family_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(family, "prevalence_ileum_family.csv")
order <- comb_prevalence(order_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(order, "prevalence_ileum_order.csv")
class <- comb_prevalence(class_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(class, "prevalence_ileum_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], meta[!colon, ], c("ANTITNF_responder"))
write.csv(phylum, "prevalence_ileum_phylum.csv")

# Responders Time
otus <- comb_prevalence(otus_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(otus, "prevalence_ileum_response_time_otus.csv")
genus <- comb_prevalence(genus_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(genus, "prevalence_ileum_response_time_genus_i.csv")
species <- comb_prevalence(species_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(species, "prevalence_ileum_response_time_species.csv")
family <- comb_prevalence(family_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(family, "prevalence_ileum_response_time_family.csv")
order <- comb_prevalence(order_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(order, "prevalence_ileum_response_time_order.csv")
class <- comb_prevalence(class_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(class, "prevalence_ileum_response_time_class.csv")
phylum <- comb_prevalence(phylum_i[, !colon], meta[!colon, ], c("ANTITNF_responder", "Time"))
write.csv(phylum, "prevalence_ileum_response_time_phylum.csv")
