#!/usr/bin/make

R_OPTS=--no-site-file --no-environ --no-restore --no-save

pre_files=16S/OTU-tab-noTRIM.csv \
          16S/OTUs-Table-BCN.tab \
          RNAseq/mapping_file.tab \
          RNAseq/mapping_file_noTRIM.tab \
          RNAseq/matrix.tsv \
          RNAseq/controls_RNAseq.RDS \
          RNAseq/db_biopsies_bcn_seq16S_noTRIM.txt

.PHONY: PCA

filtered_data.RData PCA: PCA.R $(pre_files) helper_functions.R
	@echo "Creating the PCAs"
	R CMD BATCH $(R_OPTS) $(<F)
