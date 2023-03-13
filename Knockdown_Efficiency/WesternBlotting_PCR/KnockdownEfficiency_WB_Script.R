# 1. Load libraries and variables ----
## Required libraries ----
shhh <- suppressPackageStartupMessages
shhh(library(here))
shhh(library(logger))
shhh(library(foreach))
shhh(library(doSNOW))
shhh(library(tidyverse))
shhh(library(reticulate))

## Load additional helper functions ----
source(here::here("Helper_Functions/hf_KnockdownEfficiency_WB.R"))

## Load python interface ----
use_virtualenv("UCL")
pytesseract <- import("pytesseract")
PIL <- import("PIL")

## Input Files ----
main_metadata_path <- here::here("Metadata_Results/")
metadata_path <- paste0(main_metadata_path, "metadata_samples.tsv")
metadata_WB_output <- paste0(main_metadata_path, "metadata_WB_kEff.tsv")
metadata_PCR_output <- paste0(main_metadata_path, "metadata_PCR_kEff.tsv")

## Define the algorithm variables ----
main_path = here::here("Knockdown_Efficiency/WesternBlotting_PCR/RBPs/")
download_cores = 16
overwrite_results = F
resize_perc = 0.25
min_width = 600

# 2. Pipeline ----
## Generate the variables ----
metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()

target_RBPs <- metadata %>%
  dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  dplyr::filter(!is.na(document)) %>%
  dplyr::pull(target_gene) %>%
  unique()

metadata_filtered <- metadata %>%
  dplyr::filter(target_gene %in% target_RBPs, !is.na(document)) %>%
  dplyr::select(target_gene, cell_line, experiment_id, biosample, bio_rep, document, biosample_alias) %>%
  dplyr::mutate(path = paste0(main_path, target_gene, "/", experiment_id, "/", bio_rep, "/"))
 
## Create the directories ----
createDirectories(target_RBPs, metadata_filtered)

## Download the files and add a column with their path ----
metadata_documents <- downloadCharacterizationDocuments(metadata_filtered, 
                                                        download_cores, 
                                                        overwrite_results,
                                                        silent = F)

## Check the existence of the files ----
for(row_index in seq(nrow(metadata_documents))){
  row = metadata_documents[row_index, ]
  file_path <- row$file_path
  path <- row$path
  
  if(!file.exists(file_path) || file.info(file_path)$size < 10){
    logger::log_warn("Error for row ", row_index, "! File path ", path)
  }
}

## Extract the images of all files ----
metadata_images <- extractImages(metadata_documents, 
                                 overwrite_results = overwrite_results)

## Extract text from images ----
## This cannot be separated into an external function.
## Probably because of some incompatibility with the reticulate library to use
## python.
metadata_kEff <- foreach(row_index = seq(nrow(metadata_images))) %do%{
  row = metadata_images[row_index, ]
  path <- row$path
  image_path <- row$image_path
  
  ## If the image does not exists, return the unmodified row
  if(!file.exists(image_path)) return(row)
  
  image <- PIL$Image$open(image_path)
  image_cropped <- image$crop(list(0, image$height*0.9, 0.8*image$width, image$height))
  image_small <- resizeImage(image, resize_perc, min_width)
  
  image_small$save(paste0(path, "cropped_image.png"))
  text_df <- pytesseract$image_to_string(image_small) %>%
    str_replace_all("\\f", "") %>%
    str_split("\\n", simplify = T) %>%
    str_split(" ", simplify = T) %>%
    .[1:2, ] %>%
    as_tibble()
  
  if(ncol(text_df) == 5){
    kEff_df <- text_df %>% 
      `colnames<-`(c("method", "K562_1", "K562_2", "HepG2_1", "HepG2_2")) %>%
      dplyr::mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
      dplyr::mutate(K562 = rowMeans(dplyr::select(., K562_1, K562_2)),
                    HepG2 = rowMeans(dplyr::select(., HepG2_1, HepG2_2))) %>%
      dplyr::select(method, K562, HepG2)
    
    row$WB_HepG2 <- kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(HepG2)
    row$WB_K562 <- kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(K562)
    
    row$PCR_HepG2 <- kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(HepG2)
    row$PCR_K562 <- kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(K562)
  }else if(ncol(text_df) == 3){
    kEff_df <- text_df %>% 
      `colnames<-`(c("method", "cell_line_1", "cell_line_2")) %>%
      dplyr::mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
      dplyr::mutate(cell_line = rowMeans(dplyr::select(., cell_line_1, cell_line_2))) %>%
      dplyr::select(method, cell_line)
    
    row$WB_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
    row$WB_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method == "Western") %>% dplyr::pull(cell_line), NA)
    
    row$PCR_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
    row$PCR_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% dplyr::filter(method != "Western") %>% dplyr::pull(cell_line), NA)
  }else{
    logger::WARN("Error in row ", row_index, ". Columns are not valid")
    row$WB_HepG2 <- NA
    row$WB_K562 <- NA
    row$PCR_HepG2 <- NA
    row$PCR_K562 <- NA
  }
  
  return(row)
} %>% dplyr::bind_rows()

## Test for consistency between the cell lines ----
for(target_RBP in target_RBPs){
  metadata_RBP <- metadata_kEff %>%  filter(target_gene == target_RBP)
  WB_HepG2 <- metadata_RBP$WB_HepG2
  WB_K562 <- metadata_RBP$WB_K562
  PCR_HepG2 <- metadata_RBP$PCR_HepG2
  PCR_K562 <- metadata_RBP$PCR_K562
  
  if(length(unique(na.omit(WB_HepG2))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line HepG2, method WB.")
  if(length(unique(na.omit(WB_K562))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line K562, method WB.")
  
  if(length(unique(na.omit(PCR_HepG2))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line HepG2, method PCR.")
  if(length(unique(na.omit(PCR_K562))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line K562, method PCR.")
}

## Write the knockdown efficiency table to disk ----
writeEfficiencyTable(metadata_kEff, "WB", metadata_WB_output)
writeEfficiencyTable(metadata_kEff, "PCR", metadata_PCR_output)
