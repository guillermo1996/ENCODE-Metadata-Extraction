#' Creates the required subdirectories
#'
#' @param target_RBPs List of target genes for which to create the
#'   subdirectories.
#' @param metadata_filtered Data.frame containing the information required to
#'   download the gene quantifications per sample.
#'
#' @return
#' @export
createDirectories <- function(target_RBPs, 
                              metadata_filtered){
  for(i in seq(length(target_RBPs))){
    target_RBP <- target_RBPs[i]
    RBP_metadata <- metadata_filtered %>% dplyr::filter(target_gene == target_RBP)
    for (row in seq(nrow(RBP_metadata))) {
      document.path <- RBP_metadata[row, "path", T]
      dir.create(document.path, recursive = T, showWarnings = F)
    }
  }
}

#' Download the gene quantifications files
#'
#' From the samples' metadata previously generated, it download the gene
#' quantification files per sample. This file contains all the gene counts (and
#' TPM) found in the particular sample.
#'
#' @param metadata_filtered Data.frame containing the information required to
#'   download the gene quantifications per sample.
#' @param download_cores (Optional) Number of cores to use to download. Defaults to 1.
#' @param overwrite_results (Optional) Whether to redownload already found files. Defaults to FALSE.
#' @param silent (Optinal) Whether to print progress bar of the download process. Defaults to FALSE.
#'
#' @return Data.frame with the information of the downloaded files.
#' @export
downloadGeneQuantifications <- function(metadata_filtered,
                                        download_cores = 1,
                                        overwrite_results = F,
                                        silent = F){
  if(!silent){
    pb <- txtProgressBar(max=nrow(metadata_filtered), style=3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  }else{
    opts <- list()
  }
  
  cl <- makeCluster(download_cores)
  registerDoSNOW(cl)
  metadata_gene_quant <- foreach(row_index = seq(nrow(metadata_filtered)), .combine = "rbind", .export = "getDownloadLinkGeneQuantification", .options.snow=opts) %dopar%{
    row = metadata_filtered[row_index, ]
    
    ## Two possible download paths:
    download_link <- getDownloadLinkGeneQuantification(row)
    
    ## Where to save the file
    file_path <- paste0(row$path, row$gene_quantification_id, ".tsv")
    row$file_path <- file_path
    
    ## If overwrite results is set to TRUE or if the file does not exists or if
    ## the file is too small, try to download it from the different links. If
    ## none it success, remove the resulted file. We also modify the column
    ## "file_path" if the file was successfully created.
    if(overwrite_results || !file.exists(file_path) || file.info(file_path)$size < 10){
      tryCatch({
        download.file(download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
      }, error = function(e){
        file.remove(file_path)
      })
    }
    
    if(!file.exists(file_path) || file.info(file_path)$size < 10){
      row$file_path <- NA
    }
    
    return(row)
  }
  if(!silent) close(pb)
  stopCluster(cl)
  
  return(metadata_gene_quant)
}

#' Generates the download link
#'
#' @param row Row of the data.frame with the gene quantification ID.
#'
#' @return Download link of the gene quantification file.
#' @export
getDownloadLinkGeneQuantification <- function(row){
  gene_quantification_id <- row$gene_quantification_id
  
  return(paste0("https://www.encodeproject.org/files/", gene_quantification_id, "/@@download/", gene_quantification_id, ".tsv"))
}

#' Converts HGNC gene names to ENSEMBL ID
#'
#' @param gene_list List of genes in HGCN nomenclature
#'
#' @return Data.frame containing two columns: HGCN gene names and ENSEMBL gene
#'   names. In some situations, two or more ENSEMBL names are associated with
#'   the same HGCN name.
#' @export
translateGenes <- function(gene_list){
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
  
  
  gene_df <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filter = "hgnc_symbol",
                   mart = ensembl,
                   values = gene_list)
  
  return(gene_df)
}

#' Extract TPM of a particular
#'
#' Given a Data.frame with the samples' metadata, it reads the downloaded gene
#' quantification files and extract the TPM of the target gene.
#'
#' @param metadata_quantifications Data.frame with the information of the
#'   downloaded files.
#'
#' @return Data.frame with the TPM of the target gene for each sample.
#' @export
extractTPM <- function(metadata_quantifications){
  metadata_TPM <- foreach(row_index = seq(nrow(metadata_quantifications)), .combine = dplyr::bind_rows) %do%{
    suppressPackageStartupMessages(library(tidyverse))
    
    row = metadata_quantifications[row_index, ]
    
    target_gene <- row$target_gene
    enseml_target_gene <- row$ensembl_gene_id
    
    file_path <- row$file_path
    read_tpm <- read_delim(file_path, delim = "\t", show_col_types = F) %>%
      filter(str_extract(gene_id, "[^.]+") %in% enseml_target_gene) %>% pull(TPM)
    
    if(length(read_tpm) == 0){
      logger::log_warn("No ", target_gene, " (", enseml_target_gene, ") gene found in: ", file_path)
      row$TPM <- NA
    }else{
      row$TPM <- read_tpm
    }
    
    return(row)
  }
  
  metadata_TPM <- metadata_TPM %>%
    filter(!is.na(TPM))
  
  return(metadata_TPM)
}

#' Calculates the knockdown efficiency (kEff)
#'
#' From the metadata data.frame with the TPM information, it groups by target
#' gene and by experiment type (control/case) and calculates the knockdown
#' efficiency first by calculating the averages in both clusters, and then by
#' applying the following equation:
#'
#' $$ kEff=(1-\frac{TPM_case}{TPM_control})*100% $$
#'
#' The same procedure is also executed but grouped by cell line, to also report
#' the efficiency for different cell lines.
#'
#' @param metadata_TPM Data.frame with the TPM of the target gene for each
#'   sample.
#' @param output_file (Optional) Path to where the resulted data.frame will be
#'   stored
#'
#' @return Data.frame with the target gene and the estimated knockdown
#'   efficiency. The efficiency is also calculated by cell line.
#' @export
generateKnockdownEfficiency <- function(metadata_TPM,
                                        output_file = ""){
  kEff_global <- metadata_TPM %>% 
    group_by(target_gene, experiment_type) %>%
    summarize(TPM_avg = mean(TPM, na.rm = T)) %>%
    pivot_wider(target_gene, names_from = experiment_type, values_from = TPM_avg) %>%
    mutate(kEff = (1 - case/control)*100)
  
  kEff_cell_line <- metadata_TPM %>% 
    group_by(target_gene, cell_line, experiment_type) %>%
    summarize(TPM_avg = mean(TPM, na.rm = T)) %>%
    pivot_wider(c(target_gene, cell_line), names_from = experiment_type, values_from = TPM_avg) %>%
    mutate(kEff = (1 - case/control)*100) %>%
    pivot_wider(target_gene, names_from = cell_line, values_from = kEff)
  
  metadata_kEff <- kEff_global %>% 
    left_join(kEff_cell_line, by = "target_gene") %>%
    `colnames<-`(c("target_gene", "Avg_TPM_case", "Avg_TPM_control", "kEff", "kEff_HepG2", "kEff_K562"))
  
  if(output_file != ""){
    write.table(metadata_kEff, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_kEff)
}
