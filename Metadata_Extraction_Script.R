# 1. Load libraries and variables ----
## Required libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(httr))
shhh(library(jsonlite))
shhh(library(foreach))
shhh(library(tibble))
shhh(library(logger))
shhh(library(tidyverse))

## Load additional helper functions ----
source("HelperFunctions/hf_MetadataDownloadExtraction.R")

## Logger options ----
logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout)
logger::log_threshold(logger::INFO)
log_file <- paste0(main_path, "ENCODE_API.log")
logger::log_appender(logger::appender_tee(log_file, append = T))

## Valid or required values ----
required_cell_lines <- c("HepG2", "K562")
valid_target_genes <- c() 
valid_genome_annotation <- "V29"
valid_file_format <- "bam"
valid_output_type <- "alignments"
valid_nucleic_acid_type <- "polyadenylated mRNA"

## Other parameters ----
overwrite_db <- F
download_method = "gene_silencing_series"

## Files ----
main_path <- "Metadata_results_experiments/"
if(!dir.exists(main_path)) dir.create(main_path)

### Output files
output_json <- paste0(main_path, "response.json")
output_search <- paste0(main_path, "all_experiments.tsv")
output_metadata <- paste0(main_path, "metadata_samples.tsv")
#output_metadata_combined <- paste0(main_path, "metadata_combined.tsv")

### Input files
input_target_gene_categories <- paste0("Additional_files/Target_gene_categories.tsv")
input_target_gene_NMD <- paste0("Additional_files/NMD.txt")

# 2. Pipeline ----

## Download the search data ----
URL = "https://www.encodeproject.org/search/?type=Experiment&assay_title=shRNA+RNA-seq&target.investigated_as=RNA+binding+protein&control_type!=*&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
#URL = "https://www.encodeproject.org/search/?status=released&type=Experiment&target.investigated_as=RNA+binding+protein&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=CRISPR+RNA-seq&limit=all&format=json"

## Download the json data from ENCODE
response_data <- getUrlResponse(URL, output_json)

## Summarize the search ----
summary_df <- generateSummary(response_data, 
                              valid_target_genes = valid_target_genes,
                              output_file = output_search)

## Extract the metadata ----
metadata_df <- generateMetadata(summary_df, 
                                download_method = download_method,
                                required_cell_lines = required_cell_lines,
                                valid_file_format = valid_file_format,
                                valid_genome_annotation = valid_genome_annotation,
                                valid_output_type = valid_output_type,
                                valid_nucleic_acid_type = valid_nucleic_acid_type,
                                output_file = output_metadata,
                                overwrite_db = overwrite_db)

## Add category information ----
metadata_df <- addTargetGeneCategory(metadata_df,
                                     input_target_gene_categories,
                                     input_target_gene_NMD,
                                     output_metadata)

# ## Combine with archived metadata
# metadata_archived_path <- path.expand("~/RytenLab-Research/02-ENCODE_API/Metadata_results_archived/metadata_samples.tsv")
# metadata_combined <- addArchivedMetadata(metadata_df,
#                                          metadata_archived_path,
#                                          input_target_gene_NMD,
#                                          output_file = output_metadata_combined)