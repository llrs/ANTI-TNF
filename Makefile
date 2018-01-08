#!/usr/bin/make

R_OPTS=--vanilla

pre_files=16S/OTU-tab-noTRIM \
            16S/OTUs-Table-BCN.tab \
            RNAseq/mapping_file.tab \
            RNAseq/mapping_file_noTrim.tab \
            RNAseq/matrix.tsv \
            RNAseq/controls_RNAseq.RDS \
            RNAseq/db_biopsies_bcn_seq16_noTrim.txt

$(pre_files): PCA.R
	@echo "Creating the PCAs"
	R CMD BATCH $(R_OPTS) $(<F)
