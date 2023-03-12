#' Creates the required subdirectories
#'
#' @param target_RBPs List of target genes for which to create the
#'   subdirectories.
#' @param metadata_filtered Data.frame containing the information required to
#'   download the biosample characterization document per sample.
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

#' Download the biosample preparation and characterization documents
#'
#' From the samples' metadata previously generated, it download the biosample
#' preparation and characterization documents per sample. This file contains the
#' reported efficiencies in both qRT-PCR and Western Blot.
#'
#' @param metadata_filtered Data.frame containing the information required to
#'   download the biosample preparation and characterization documents per sample.
#' @param download_cores (Optional) Number of cores to use to download. Defaults
#'   to 1.
#' @param overwrite_results (Optional) Whether to download already found
#'   files. Defaults to FALSE.
#' @param silent (Optional) Whether to print progress bar of the download
#'   process. Defaults to FALSE.
#'
#' @return Data.frame with the information of the downloaded files.
#' @export
downloadCharacterizationDocuments <- function(metadata_filtered,
                                              download_cores = 1,
                                              overwrite_results = FALSE,
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
  metadata_documents <- foreach(row_index = seq(nrow(metadata_filtered)), .combine = "rbind", .options.snow=opts) %dopar%{
    row = metadata_filtered[row_index, ]
    
    ## Two possible download paths:
    download_link <- paste0("https://www.encodeproject.org", row$document, "@@download/attachment/", row$biosample_alias, ".pdf")
    alt_download_link <- paste0("https://www.encodeproject.org", row$document, "@@download/attachment/", row$biosample_alias, "_update.pdf")
    
    ## Where to save the file
    file_path <- paste0(row$path, row$biosample_alias, ".pdf")
    row$file_path <- file_path
    
    ## If overwrite results is set to TRUE or if the file does not exists or if
    ## the file is too small, try to download it from the different links. If
    ## none it success, remove the resulted file. We also modify the column
    ## "file_path" if the file was successfully created.
    if(overwrite_results || !file.exists(file_path) || file.info(file_path)$size < 10){
      tryCatch({
        download.file(download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
      }, error = function(e){
        tryCatch({
          download.file(alt_download_link, file_path, method = "wget", extra = "--quiet --no-check-certificate")
        }, error = function(k){
          file.remove(file_path)
        })
      })
    }
    
    if(!file.exists(file_path) || file.info(file_path)$size < 10){
      row$file_path <- NA
    }
    
    return(row)
  }
  if(!silent) close(pb)
  stopCluster(cl)
  
  return(metadata_documents)
}

#' Extract the images from documents
#'
#' Given the data.frame with the sample's information and their characterization
#' document path, this function extracts the images from the documents and keeps
#' the last found image (which always contains the reported knockdown
#' efficiency).
#'
#' @param metadata_documents Data.frame with the information of the downloaded
#'   files.
#' @param overwrite_results (Optional) Whether to download already found files.
#'   Defaults to FALSE.
#'
#' @return Data.frame with the information of the extracted images.
#' @export
extractImages <- function(metadata_documents,
                          overwrite_results = FALSE){
  metadata_images <- foreach(row_index = seq(nrow(metadata_documents)), .combine = dplyr::bind_rows) %do%{
    row = metadata_documents[row_index, ]
    
    biosample = row$biosample
    file_path <- row$file_path
    image_path <- paste0(row$path, biosample, "_Western_Blot_Analysis.png")
    
    ## If the file does not exists return the unmodified row
    if(!file.exists(file_path)) return(row)
    
    ## If overwrite results is set to TRUE or the image does not exists, execute
    ## a python command to extract the figures. All the figures must have a .png
    ## or .jpeg extension. We remove all images but the last one, which contains
    ## the Western Blotting results.
    if(overwrite_results || !file.exists(image_path)){
      system2(command = path.expand("~/.virtualenvs/UCL/bin/python"),
              args = c("-m fitz extract -images", file_path, "-output", row$path))
      
      short_images <- list.files(row$path, pattern = ".png|.jpeg")
      long_images <- list.files(row$path, pattern = ".png|.jpeg", full.names = T)
      
      last_image <- sapply(short_images, function(x) str_extract(x, "\\d+") %>% as.numeric()) %>% max(na.rm = T)
      
      for(i in seq(length(short_images))){
        image_name <- short_images[i]
        file <- long_images[i]
        
        if(!grepl(last_image, image_name)){
          file.remove(file)
        }else{
          file.rename(file, image_path)
        }
      }
    }
    
    row$image_path <- image_path 
    row
  }
  
  return(metadata_images)
}

#' Resize an image
#'
#' @param image PIL image object.
#' @param resize_perc Final size of the output image in percentage. Defaults to
#'   0.25.
#' @param min_width Minimum width to apply the resize. If the image is too
#'   small, the text recognition algorithm will fail. Defaults to 600px.
#'
#' @return PIL image object.
#' @export
resizeImage <- function(image, 
                        resize_perc = 0.25,
                        min_width = 600){
  if(image$width > 600){
    image_small <- image_cropped$resize(list((image_cropped$width*resize_perc) %>% as.integer(), 
                                             (image_cropped$height*resize_perc) %>% as.integer()), 
                                        PIL$Image$Resampling$LANCZOS)
  }else{
    image_small <- image
  }
  
  return(image_small)
}

#' Summarize and write to disk the Knockdown Efficiency
#'
#' @param metadata_kEff Data.frame containing the extracted efficiencies for
#'   both methods.
#' @param method Whether to extract the Western blotting or qRT-PCR results.
#'   Valid values are: "WB" and "PCR". Defaults to "WB".
#' @param output_file (Optional) Path to where the resulted data.frame will be
#'   stored
#'
#' @return Data.frame with the summarized knockdown efficiencies reported by
#'   ENCODE with either qRT-PCR or Western blotting.
#' @export
writeEfficiencyTable <- function(metadata_kEff,
                                 method = "WB",
                                 output_file = ""){
  if(method == "WB"){
    metadata_kEff_output <- metadata_kEff %>%
      dplyr::select(target_gene, WB_HepG2, WB_K562) %>%
      tidyr::pivot_longer(c(WB_HepG2, WB_K562)) %>%
      dplyr::mutate(name = str_replace_all(name, "WB_", ""))
  }else if(method == "PCR"){
    metadata_kEff_output <- metadata_kEff %>%
      dplyr::select(target_gene, PCR_HepG2, PCR_K562) %>%
      tidyr::pivot_longer(c(PCR_HepG2, PCR_K562)) %>%
      dplyr::mutate(name = str_replace_all(name, "PCR_", ""))
  }
  
  metadata_kEff_output <- metadata_kEff_output %>%
    dplyr::group_by(target_gene, name) %>%
    dplyr::summarize(kEff = ifelse(all(is.na(value)), NA, max(value, na.rm = T))) %>%
    dplyr::mutate(kEff_avg = mean(kEff, na.rm = T)) %>%
    dplyr::mutate_all(~ifelse(is.nan(.), NA, .)) %>%
    tidyr::pivot_wider(c(target_gene, kEff_avg), values_from = kEff, names_from = name, names_prefix = "kEff_")
  
  ## Write to disk
  if(output_file != ""){
    write.table(metadata_kEff_output, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_kEff_output)
}
