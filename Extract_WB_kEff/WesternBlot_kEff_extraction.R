## Required libraries
library(logger)
library(foreach)
library(doSNOW)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reticulate))


## Load additional helper functions
source("helperFunctions.R")

## Load python interface
use_virtualenv("UCL")
pytesseract <- import("pytesseract")
PIL <- import("PIL")

## Files
metadata_path <- path.expand("~/RytenLab-Research/02-ENCODE_API/Metadata_results/metadata_combined.tsv")
metadata_output <- path.expand("~/RytenLab-Research/02-ENCODE_API/Metadata_results/metadata_kEff.tsv")
metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()

## Define the algorithm variables
main_path = "RBPs/"
download_only = F
download_cores = 8
overwrite_results = F
resize_perc = 0.25

## Generate the variables
target_RBPs <- metadata %>%
  filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  filter(!is.na(document)) %>%
  pull(target_gene) %>%
  unique()

metadata_filtered <- metadata %>%
  filter(target_gene %in% target_RBPs, !is.na(document)) %>%
  select(target_gene, cell_line, experiment_id, biosample, bio_rep, document, biosample_alias) %>%
  mutate(path = paste0(main_path, target_gene, "/", experiment_id, "/", bio_rep, "/"))

## Create the directories
createDirectories(target_RBPs, metadata_filtered)

## Download the files and add a column with their path
metadata_documents <- downloadDocuments(metadata_filtered, 
                                        download_cores, 
                                        overwrite_results)

## Check the existence of the files
for(row_index in seq(nrow(metadata_documents))){
  row = metadata_documents[row_index, ]
  file_path <- row$file_path
  path <- row$path
  
  if(!file.exists(file_path) || file.info(file_path)$size < 10) logger::log_warn("Error for row ", row_index, "! File path ", path)
}

## Extract the images of all files
metadata_images <- extractImages(metadata_documents, 
                                 overwrite_results = overwrite_results)

## Extract text from images. This cannot be separated into an external function.
## Probably because of some incompatibility with the reticulate library to use
## python.
metadata_kEff <- foreach(row_index = seq(nrow(metadata_images)), .combine = dplyr::bind_rows) %do%{
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
      mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
      mutate(K562 = rowMeans(select(., K562_1, K562_2)),
             HepG2 = rowMeans(select(., HepG2_1, HepG2_2))) %>%
      select(method, K562, HepG2)
    
    row$WB_HepG2 <- kEff_df %>% filter(method == "Western") %>% pull(HepG2)
    row$WB_K562 <- kEff_df %>% filter(method == "Western") %>% pull(K562)
  }else if(ncol(text_df) == 3){
    kEff_df <- text_df %>% 
      `colnames<-`(c("method", "cell_line_1", "cell_line_2")) %>%
      mutate(across(-method, function(x) as.numeric(sub("%", "", x)))) %>%
      mutate(cell_line = rowMeans(select(., cell_line_1, cell_line_2))) %>%
      select(method, cell_line)
    
    row$WB_HepG2 <- ifelse(row$cell_line == "HepG2", kEff_df %>% filter(method == "Western") %>% pull(cell_line), NA)
    row$WB_K562 <- ifelse(row$cell_line == "K562", kEff_df %>% filter(method == "Western") %>% pull(cell_line), NA)
  }else{
    logger::WARN("Error in row ", row_index, ". Columns are not valid")
    row$WB_HepG2 <- NA
    row$WB_K562 <- NA
  }
  
  return(row)
}

## Test for consistency between the cell lines
for(target_RBP in target_RBPs){
  metadata_RBP <- metadata_kEff %>%  filter(target_gene == target_RBP)
  WB_HepG2 <- metadata_RBP$WB_HepG2
  WB_K562 <- metadata_RBP$WB_K562
  
  if(length(unique(na.omit(WB_HepG2))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line HepG2")
  if(length(unique(na.omit(WB_K562))) > 1) logger::ERROR("Error in RBP ", target_RBP, " cell line K562")
}

## Add to Metadata
metadata_kEff_output <- metadata_kEff %>%
  select(target_gene, WB_HepG2, WB_K562) %>%
  pivot_longer(c(WB_HepG2, WB_K562)) %>%
  mutate(name = str_replace_all(name, "WB_", "")) %>%
  group_by(target_gene, name) %>%
  summarize(kEff = ifelse(all(is.na(value)), NA, max(value, na.rm = T))) %>%
  mutate(kEff_avg = mean(kEff, na.rm = T)) %>%
  mutate_all(~ifelse(is.nan(.), NA, .)) %>%
  pivot_wider(c(target_gene, kEff_avg), values_from = kEff, names_from = name, names_prefix = "kEff_")

## Write to disk
write.table(metadata_kEff_output, metadata_output, sep = "\t", row.names = F, quote = FALSE)
