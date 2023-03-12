## Required libraries
library(logger)
suppressPackageStartupMessages(library(biomaRt))
library(foreach)
library(doSNOW)
suppressPackageStartupMessages(library(tidyverse))

## Load additional helper functions
source("helperFunctions.R")

## Logger options
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)
logger::log_threshold(logger::INFO)
log_file <- paste0("Extract_TPM.log")
logger::log_appender(logger::appender_tee(log_file, append = T))

## Files
metadata_path <- "../Metadata_results/metadata_samples.tsv"
metadata_output <-"../Metadata_results/metadata_TPM_kEff.tsv"
metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()

## Define the algorithm variables
main_path = "RBPs/"
download_only = F
download_cores = 16
overwrite_results = F

## Generate the variables
target_RBPs <- metadata %>%
  filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  filter(!is.na(gene_quantification_id)) %>%
  pull(target_gene) %>%
  unique()
ensembl_target_RBPs <- translateGenes(target_RBPs)

metadata_filtered <- metadata %>%
  filter(target_gene %in% target_RBPs, !is.na(gene_quantification_id)) %>%
  dplyr::select(target_gene, cell_line, experiment_type, experiment_id, gene_quantification_id) %>%
  mutate(path = paste0(main_path, target_gene, "/", experiment_type, "/")) %>%
  left_join(ensembl_target_RBPs, by = c("target_gene" = "hgnc_symbol")) %>%
  relocate(ensembl_gene_id, .before = cell_line)

## Create the directories
createDirectories(target_RBPs, metadata_filtered)

## Download the files and add a column with their path
metadata_quantifications <- downloadGeneQuantifications(metadata_filtered,
                                                        download_cores, 
                                                        overwrite_results)

## Extract the TPMs
metadata_TPM <- extractTPM(metadata_quantifications)

## Generate and save the knockdown efficiencies
metadata_kEff <- generateKnockdownEfficiency(metadata_TPM,
                                             output_file = metadata_output)
