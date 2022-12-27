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

downloadDocuments <- function(metadata_filtered,
                              download_cores,
                              overwrite_results){
  pb <- txtProgressBar(max=nrow(metadata_filtered), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
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
  close(pb)
  stopCluster(cl)
  
  return(metadata_documents)
}

extractImages <- function(metadata_documents,
                          overwrite_results){
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

extractText <- function(metadata_images,
                        resize_perc = 0.25,
                        min_width = 600){
  metadata_kEff <- foreach(row_index = seq(nrow(metadata_images[1:8, ])), .combine = dplyr::bind_rows) %do%{
    print(row_index)
    row = metadata_images[row_index, ]
    
    path <- row$path
    image_path <- row$image_path
    
    ## If the image does not exists, return the unmodified row
    if(!file.exists(image_path)) return(row)
    
    print(image_path)
    image <- PIL$Image$open(image_path)
    image_cropped <- image$crop(list(0, image$height*0.9, 0.8*image$width, image$height))
    image_small <- resizeImage(image, resize_perc, min_width)
    
    image_small$save(paste0(path, "cropped_image.png"))
    
    print(pytesseract$image_to_string(image_small))
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
  
  print(metadata_kEff)
  return(metadata_kEff)
}

resizeImage <- function(image, 
                        resize_perc,
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
