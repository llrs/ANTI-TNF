#!/usr/bin/make

R_OPTS=--no-site-file --no-environ --no-restore --no-save

pre_files=16S/OTU-tab-noTRIM.csv \
          16S/OTUs-Table-BCN.tab \
          RNAseq/mapping_file.tab \
          RNAseq/mapping_file_noTRIM.tab \
          RNAseq/matrix.tsv \
          RNAseq/controls_RNAseq.RDS \
          RNAseq/db_biopsies_bcn_seq16S_noTRIM.txt \
          Metadata_BCN.csv

.PHONY: clean all correlations PCA RGCCA biological CD

all: PCA RGCCA biological correlations

filtered_data.RData PCA: PCA.R $(pre_files) helper_functions.R
	@echo "Creating the PCAs"
	R CMD BATCH $(R_OPTS) $(<F)

bootstrap.RData sgcca.RData: RGCCA.R helper_functions.R filtered_data.RData helper_RGCCA.R
	@echo "Integrating the datasets"
	R CMD BATCH $(R_OPTS) $(<F)

bootstrap_CD.RData CD.RData BCN_CD.RDS: RGCCA_CD.R helper_functions.R filtered_data.RData helper_RGCCA.R
	@echo "Integrating the datasets"
	R CMD BATCH $(R_OPTS) $(<F)

RNAseq_enrichment.csv Otus_genus_enrichment.csv: biological_relationships.R helper_functions.R sgcca.RData bootstrap.RData
	@echo "Interpreting biologically the loadings"
	R CMD BATCH $(R_OPTS) $(<F)

RNAseq_enrichment_CD.csv Otus_genus_enrichment_CD.csv: br_CD.R helper_functions.R bootstrap_CD.RData CD.RData BCN_CD.RDS
	@echo "Interpreting biologically the loadings"
	R CMD BATCH $(R_OPTS) $(<F)

bootstrap_superblock.RData superblock.RData superblock: superblock.R helper_functions.R filtered_data.RData
	@echo "Integrating the datasets in one big block"
	R CMD BATCH $(R_OPTS) $(<F)

correlations: correlations.R helper_functions.R filtered_data.RData sgcca.RData
	@echo "Heatmaps of correlations"
	R CMD BATCH $(R_OPTS) $(<F)

RGCCA: bootstrap.RData bootstrap_CD.RData

biological: RNAseq_enrichment.csv RNAseq_enrichment_CD.csv

CD: bootstrap_CD.RData RNAseq_enrichment_CD.csv

clean:
	@echo "Deleting all the outputs of R"
	rm *.Rout
	@echo "Deleting all the prevalence tables"
	rm prevalence_*.csv
