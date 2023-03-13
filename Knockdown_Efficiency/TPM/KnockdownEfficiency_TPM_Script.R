# 1. Load libraries and variables ----
## Required libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(biomaRt))
shhh(library(here))
shhh(library(logger))
shhh(library(foreach))
shhh(library(doSNOW))
shhh(library(tidyverse))
options(dplyr.summarise.inform = FALSE)

## Load additional helper functions ----
source(here::here("Helper_Functions/hf_KnockdownEfficiency_TPM.R"))

## Logger options ----
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)
logger::log_threshold(logger::INFO)
log_file <- paste0(here::here("Knockdown_Efficiency/TPM/Extract_TPM.log"))
logger::log_appender(logger::appender_tee(log_file, append = T))

## Input Files ----
main_metadata_path <- here::here("Metadata_Results/")
metadata_path <- paste0(main_metadata_path, "metadata_samples.tsv")
metadata_TPM_output <- paste0(main_metadata_path, "metadata_TPM_kEff.tsv")

## Define the algorithm variables ----
main_path = here::here("Knockdown_Efficiency/TPM/RBPs/")
download_only = F
download_cores = 16
overwrite_results = F

# 2. Pipeline ----
## Generate the variables ----
metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()

target_RBPs <- metadata %>%
  dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  dplyr::filter(!is.na(gene_quantification_id)) %>%
  dplyr::pull(target_gene) %>%
  unique()
ensembl_target_RBPs <- translateGenes(target_RBPs)

metadata_filtered <- metadata %>%
  dplyr::filter(target_gene %in% target_RBPs, !is.na(gene_quantification_id)) %>%
  dplyr::select(target_gene, cell_line, experiment_type, experiment_id, gene_quantification_id) %>%
  dplyr::mutate(path = paste0(main_path, target_gene, "/", experiment_type, "/")) %>%
  dplyr::left_join(ensembl_target_RBPs, by = c("target_gene" = "hgnc_symbol"), multiple = "all") %>%
  dplyr::relocate(ensembl_gene_id, .before = cell_line)

## Create the directories ----
createDirectories(target_RBPs, metadata_filtered)

## Download the files and add a column with their path ----
metadata_quantifications <- downloadGeneQuantifications(metadata_filtered,
                                                        download_cores, 
                                                        overwrite_results)

## Extract the TPMs ----
metadata_TPM <- extractTPM(metadata_quantifications)

## Generate and save the knockdown efficiencies ----
metadata_kEff <- generateKnockdownEfficiency(metadata_TPM,
                                             output_file = metadata_TPM_output)
