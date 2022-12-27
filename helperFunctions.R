#' Extract API response
#'
#' Given an URL, it requests the information from the ENCODE API.
#'
#' @param url Link to the ENCODE page to extract the information.
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#'
#' @return A list object containing the json downloaded from the ENCODE portal.
#' @export
getUrlResponse <- function(url, output_file = NULL) {
  response <- GET(url)
  r <- content(response, as = "text", encoding = "UTF-8")
  
  if(!is.null(output_file)){
    write(prettify(r, indent = 4), output_file)
  }
  
  return(fromJSON(r, flatten = T))
}

#' Generates a summary of the ENCODE experiment search
#'
#' From the generated list object built in [getUrlResponse()], it reads every
#' experiment found and extract their experiment ID, Gene Silencing Series (gss)
#' ID and cell line. All this information is summarized in a dataframe, where
#' every row is an experiment.
#'
#' @param response_data List object build from [getUrlResponse()] and a search
#'   URL.
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#'
#' @return Data.frame containing a summary of every experiment found in the
#'   ENCODE search.
#' @export
generateSummary <- function(response_data,
                            valid_target_genes = c(),
                            output_file = NULL) {
  summary_df <- tibble()
  for(i in seq(nrow(response_data$`@graph`))) {
    experiment <- response_data$`@graph`[i, ]
    
    biosample_ontology <- experiment$biosample_ontology.term_name
    target_gene <- experiment$target.label
    sample_id <- experiment$accession
    
    ## Extract the Gene Silencing Series
    gene_silencing_series <- sapply(experiment$related_series, function(x) x$accession)
    if(length(gene_silencing_series) > 1){
      logger::log_warn("More than one gene silencing series found for target gene ", target_gene, " and biosample ", biosample_ontology)
      gene_silencing_series <- gene_silencing_series[which(sapply(experiment$related_series[[1]]$`@type`, function(x) "GeneSilencingSeries" %in% x))] 
    }
    
    tmp_df <- tibble(target_gene = target_gene,
                     experiment_id = sample_id,
                     cell_line = biosample_ontology,
                     gene_silencing_series = gene_silencing_series)
    summary_df <- rbind(summary_df, tmp_df)
  }
  
  summary_df <- summary_df %>% arrange(target_gene, cell_line)
  
  if(!is.null(output_file)){
    write.table(summary_df, output_file, sep = "\t", row.names = F)
  }
  
  if(!is.null(valid_target_genes)){
    summary_df <- summary_df %>% filter(target_gene %in% valid_target_genes)
  }
  
  return(summary_df)
}


#' Generates the metadata data.frame
#'
#' Given a summary data.frame with the experiments information and the filters
#' to apply in the extraction process, it generates a data.frame with each
#' sample as a row and the sample's metadata as the columns.
#'
#' @param summary_df Data.frame containing a summary of every experiment found
#'   in the ENCODE search.
#' @param download_method (Optional) Whether to use the Gene Silencing Series
#'   (gss) ID or experiment ID to access the experiment's information. If
#'   available, the gss method is recomended since it requires less calls to the
#'   API. Valid inputs: "gene_silencing_series", "experiments". Defaults to
#'   "gene_silencing_series".
#' @param required_cell_lines (Optional) Required cell lines to extract
#'   information about the experiments. It requires that a least one experiment
#'   for each cell line provided is present in the summary data.frame. Defaults
#'   to "HepG2" and "K562".
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param output_file (Optional) Path to where the json downloaded file will be
#'   stored.
#' @param overwrite_db (Optional) If set to TRUE, it will read first the
#'   metadata found in "output_file" and will try to not redownload the metadata
#'   of previously successfull runs.
#'
#' @return Data.frame where every row is an ENCODE sample and the columns
#'   contain their metadata.
#' @export
generateMetadata <- function(summary_df,
                             download_method = "gene_silencing_series",
                             required_cell_lines = c("HepG2", "K562"), 
                             valid_nucleic_acid_type = "polyadenylated mRNA",
                             valid_genome_annotation = "V29",
                             valid_file_format = "bam",
                             valid_output_type = "alignments",
                             output_file = "",
                             overwrite_db = T){
  ## Check the overwrite_db configuration. If set to FALSE, we will try to read the
  ## output file and extract the previous information from it.
  if(output_file == "") overwrite_db = T
  if(!overwrite_db & file.exists(output_file)){
    df_previous <- read.csv(output_file, sep = "\t") %>% as_tibble() %>% select(-any_of(c("Splicing_regulation", "Spliceosome", "Exon_junction_complex", "NMD")))
  }else{
    df_previous <- tibble() 
    overwrite_db = T
  }
  
  ## Get the valid target genes by cell lines
  valid_target_genes <- summary_df %>% 
    group_by(target_gene) %>%
    summarise(cell_lines = list(cell_line))%>%
    rowwise %>%
    filter(all(required_cell_lines %in% cell_lines)) %>%
    pull(target_gene)
  
  metadata_df <- tibble()
  ## Loop through every target gene
  for(iter_tg in seq(length(valid_target_genes))){
    iter_target_gene <- valid_target_genes[iter_tg]
    logger::log_info("Starting target gene ", iter_target_gene, ":")
    
    ## Local target gene information
    summary_tg_df <- summary_df %>% filter(target_gene == iter_target_gene)
    
    ## If we find a total of 8 entries in the previous search, we skip this
    ## target gene. 
    if(!overwrite_db){
      if (df_previous %>% filter(target_gene == iter_target_gene) %>% nrow() == 8){
        logger::log_info("\t Already found in output file.")
        metadata_df <- rbind(metadata_df, df_previous %>% filter(target_gene == iter_target_gene))
        next
      }
    }
    
    ## Download the target gene metadata
    if(download_method == "gene_silencing_series"){
      tg_df <- LoopGeneSilencingSeries(summary_tg_df,
                                       valid_nucleic_acid_type,
                                       valid_file_format,
                                       valid_output_type,
                                       valid_genome_annotation)
    }else if(download_method == "experiments"){
      tg_df_2 <- LoopExperiments(summary_tg_df,
                                 valid_nucleic_acid_type,
                                 valid_file_format,
                                 valid_output_type,
                                 valid_genome_annotation)
    }else{
      logger::ERROR("No valid download method provided. Only gene_silencing_series or experiments are allowed.")
    }
    
    ## Add final information and sort the columns
    tg_df <- tg_df %>%
      mutate(target_gene = iter_target_gene, .before = cell_line) %>%
      filter(nucleic_acid_type == valid_nucleic_acid_type)
    
    column_order <- c("target_gene", "experiment_type", "cell_line", 
                      "gene_silencing_series", "experiment_id", 
                      "experiment_doi", "sample_id", "rin", "read_depth", 
                      "bio_rep", "tech_rep", "sex", "age", "life_stage", 
                      "gene_quantification_id", "file_format", 
                      "output_type", "genome_annotation", 
                      "mapped_run_type", "lab", "assay", "cellosaurus")
    
    tg_df <- tg_df %>% 
      select(c(intersect(column_order, names(.)), setdiff(names(.), column_order)))
    
    ## If the required cell lines are not found, ignore the target gene
    if(!all.equal(tg_df$cell_line %>% unique %>% sort, required_cell_lines %>% sort)){
      next
    }
    
    ## Appends the data and store to disk
    metadata_df <- rbind(metadata_df, tg_df)
    if(output_file != ""){
      write.table(metadata_df, output_file, sep = "\t", row.names = F, quote = FALSE)
    }
  }
  
  return(metadata_df)
}

#' Extract metadata from the Gene Silencing Series (gss)
#'
#' @param summary_tg_df Local target gene information from the summary
#'   data.frame generated with [generateSummary()].
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the metadata for a particular gss.
#' @export
LoopGeneSilencingSeries <- function(summary_tg_df,
                                    valid_nucleic_acid_type,
                                    valid_file_format,
                                    valid_output_type,
                                    valid_genome_annotation){
  ## Loop through every gene silencing series found for the target gene
  tg_df <- foreach(iter_gss = seq(nrow(summary_tg_df)), .combine = dplyr::bind_rows) %do%{
    gss <- summary_tg_df[iter_gss, ]
    
    ## Gene silencing series information
    gss_id <- gss$gene_silencing_series
    gss_cell_line <- gss$cell_line
    logger::log_info("\t Starting Gene Silencing Series ", gss_id, " (cell line = ", gss_cell_line, "):")
    
    ## Get the API response for the gene silencing series
    response_gss <- getUrlResponse(paste0("https://www.encodeproject.org/gene-silencing-series/", gss_id, "?format=json"))
    
    ## Information about the case and control samples
    gss_experiments <- response_gss$related_datasets
    
    ## Loop through every experiment found in the gene silencing series
    gss_df <- foreach(iter_experiment = seq(nrow(gss_experiments)), .combine = dplyr::bind_rows) %do%{
      experiment <- gss_experiments[iter_experiment, ]
      
      experiment_id <- experiment$accession
      experiment_type <- getExperimentType(experiment, gss$target_gene, gss_cell_line) 
      
      logger::log_info("\t\t Starting experiment ", experiment_id, " (type = ", experiment_type, ").")
      
      ## Requirements metadata
      experiment_additional_info <- getAdditionalInformation(experiment)
      
      ## Main metadata
      experiment_sample_files <- getSampleFiles(experiment$files[[1]], valid_file_format, valid_output_type, valid_genome_annotation)
      if(experiment_additional_info$nucleic_acid_type != valid_nucleic_acid_type) return(cbind(experiment_sample_files, experiment_additional_info))
      
      ## Other metadata
      experiment_rin <- getRin(experiment$replicates[[1]])
      experiment_read_depth <- getReadDepth(experiment$analyses[[1]], valid_genome_annotation)
      experiment_donor_info <- getDonorInfo(experiment$replicates[[1]])
      experiment_documents <- getDocumentFiles(experiment$replicates[[1]], gss$target_gene, experiment_type)
      experiment_doi <- experiment$doi
      experiment_gene_quantifications <- getGeneQuantificationFiles(experiment$files[[1]])
      
      ## Combine all information
      experiment_combined <- experiment_sample_files %>% 
        left_join(experiment_rin, by = "bio_rep") %>%
        left_join(experiment_read_depth, by = "bio_rep") %>%
        left_join(experiment_donor_info, by = "bio_rep") %>% 
        left_join(experiment_gene_quantifications %>% select(gene_quantification_id, bio_rep), by = "bio_rep") %>%
        left_join(experiment_documents, by = "bio_rep") %>%
        mutate(experiment_type = experiment_type, experiment_id = experiment_id, experiment_doi = experiment_doi, .before = "bio_rep") %>%
        cbind(experiment_additional_info)
      
      return(experiment_combined)
    }
    
    ## Add the relevant information
    gss_df <- gss_df %>%
      mutate(cell_line = gss_cell_line, 
             gene_silencing_series = gss_id, .before = sample_id)
  }
  
  return(tg_df)
}


#' Gets the experiment type (e.g. control/case)
#'
#' @param related_dataset Data.frame with the information about the experiment
#'   as provided by the ENCODE API.
#' @param target_gene Target gene of study.
#' @param gss_cell_line Cell line of the gss.
#'
#' @return The experiment type (e.g. control/case)
#' @export
getExperimentType <- function(related_dataset,
                              target_gene,
                              gss_cell_line){
  if(!"control_type" %in% names(related_dataset)){
    logger::log_warn(" No control experiment found for target gene ", target_gene, " and cell line ", gss_cell_line, ".")
    experiment_type = "case"
  }else{
    experiment_type <- ifelse(is.na(related_dataset$control_type), "case", "control")
  }
  
  return(experiment_type)
}

#' Gets the experiment additional information.
#'
#' @param related_dataset Data.frame with the information about the experiment
#'   as provided by the ENCODE API.
#'
#' @return A data.frame with the metadata about the laboratory, the assay, the
#'   cellosaurus, the nucleic acid type, the extraction method, the
#'   fragmentation method, the size selection method and the strand specificity.
#' @export
getAdditionalInformation <- function(related_dataset){
  sample_lab <- related_dataset$lab.title
  if(is.null(sample_lab)) 
    sample_lab <- related_dataset$lab$title
  sample_assay <- related_dataset$assay_term_name
  sample_cellosaurus <- related_dataset$biosample_ontology.dbxrefs %>% unlist %>% unique
  if(is.null(sample_cellosaurus)) 
    sample_cellosaurus <- related_dataset$biosample_ontology$dbxrefs %>% unlist %>% unique
  
  sample_nucleic_acid_type <- related_dataset$replicates %>% data.frame() %>% as_tibble() %>% pull(library.nucleic_acid_term_name) %>% unique
  sample_extraction_method <- related_dataset$replicates %>% data.frame() %>% as_tibble() %>% pull(library.extraction_method) %>% unique
  sample_fragmentation_method <- related_dataset$replicates %>% data.frame() %>% as_tibble() %>% pull(library.fragmentation_methods) %>% unlist %>% unique
  sample_size_selection_method <- related_dataset$replicates %>% data.frame() %>% as_tibble() %>% pull(library.library_size_selection_method) %>% unique
  sample_strand_specificity <- related_dataset$replicates %>% data.frame() %>% as_tibble() %>% pull(library.strand_specificity) %>% unique
  
  tibble(lab = sample_lab,
         assay = sample_assay,
         cellosaurus = sample_cellosaurus,
         nucleic_acid_type = sample_nucleic_acid_type,
         extraction_method = sample_extraction_method,
         fragmentation_method = sample_fragmentation_method,
         size_selection_method = sample_size_selection_method,
         strand_specificity = sample_strand_specificity) %>%
    return()
}

#' Gets the experiment's sample files
#'
#' @param files Data.frame containing all the files found for the experiment.
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return A data.frame with the metadata about the sample files, biological
#'   replicate, output type, etc.
#' @export
getSampleFiles <- function(files,
                           valid_file_format = "bam",
                           valid_output_type = "alignments",
                           valid_genome_annotation = "V29") {
  sample_files_info <- files %>% 
    filter(file_format == valid_file_format) %>%
    filter(output_type == valid_output_type) %>%
    filter(genome_annotation == valid_genome_annotation) %>%
    select(accession, biological_replicates, file_format, output_type, genome_annotation, technical_replicates, mapped_run_type) %>%
    tidyr::unnest(c(biological_replicates, technical_replicates)) %>%
    dplyr::rename("bio_rep" = "biological_replicates",
                  "tech_rep" = "technical_replicates",
                  "sample_id" = "accession")
  
  return(sample_files_info)
}

#' Gets the experiment's RIN information
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#'
#' @return Data.frame with the RIN information for each experiment's isogenic
#'   replicates.
#' @export
getRin <- function(replicates) {
  if(!"library.rna_integrity_number" %in% names(replicates)){
    logger::log_warn("\t\t\t No RIN found. Defaulted to NA.")
    rin_info <- replicates %>% select(biological_replicate_number) %>% mutate(rin = NA)
  }else{
    rin_info <- replicates %>% select(biological_replicate_number, library.rna_integrity_number)
  }
  names(rin_info) <- c("bio_rep", "rin")
  
  return(rin_info)
}

#' Gets the experiment's sample read depth
#'
#' @param analyses Data.frame containing information about the different
#'   analyses executed within the experiment.
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the read depth for each isogenic replicate.
#' @export
getReadDepth <- function(analyses,
                         valid_genome_annotation = "V29") {
  read_depth_info <- analyses %>% 
    filter(genome_annotation == valid_genome_annotation) %>%
    pull(`quality_metrics_report.Read depth`) %>%
    magrittr::extract2(1) %>% 
    select(biological_replicates, metric) %>%
    tidyr::unnest(biological_replicates)
  names(read_depth_info) <- c("bio_rep", "read_depth")
  
  return(read_depth_info)
}

#' Gets the sample donor information
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#'
#' @return Data.frame containing information about the sample donor.
#' @export
getDonorInfo <- function(replicates) {
  donor_info <- replicates %>% 
    select(biological_replicate_number, library.biosample.sex, library.biosample.age, library.biosample.life_stage)
  names(donor_info) <- c("bio_rep", "sex", "age", "life_stage")
  
  return(donor_info)
}

#' Get documents related to the experiment
#'
#' @param replicates Data.frame containing information about the isogenic
#'   replicates within the experiment.
#' @param target_gene Target gene of study.
#' @param experiment_type Experiment type (e.g. control/case)
#'
#' @return Data.frame containing the information about the biosample and characterization documents.
#' @export
getDocumentFiles <- function(replicates, 
                             target_gene, 
                             experiment_type){
  if(experiment_type == "case"){
    documents_info <- foreach(iter_replicate = seq(nrow(replicates)), .combine = "rbind") %do%{
      replicate <- replicates[iter_replicate, ]
      
      bio_rep <- replicate$biological_replicate_number
      accession <- replicate$library.accession
      documents <- replicate$library.biosample.documents %>% unlist()
      if(is.null(documents)) documents <- NA
      
      aliases <- replicate$library.biosample.aliases %>% unlist()
      if(length(aliases) > 1){
        aliases <- aliases[which(!grepl(":BG[a-zA-Z]LV", aliases))]
      }
      aliases <- str_split(aliases, ":|,", simplify = T)[, 2]
      
      tibble(bio_rep = bio_rep, 
             biosample = accession, 
             document = documents, 
             biosample_alias = aliases)
    }
  }else{
    documents_info <- tibble(bio_rep = c(1, 2))
  }
  
  return(documents_info)
}

#' Gets the gene quantification files from the experiment
#'
#' @param files Data.frame containing all the files found for the experiment.
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "TSV".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "gene quantifications".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the information of the gene quantifications
#'   files found within the experiment.
#' @export
getGeneQuantificationFiles <- function(files,
                                       valid_file_format = "tsv",
                                       valid_output_type = "gene quantifications",
                                       valid_genome_annotation = "V29") {
  tsv_files_info <- files %>% 
    filter(file_format == valid_file_format) %>%
    filter(output_type == valid_output_type) %>%
    filter(genome_annotation == valid_genome_annotation) %>%
    select(accession, biological_replicates, file_format, output_type, genome_annotation, technical_replicates) %>%
    tidyr::unnest(c(biological_replicates, technical_replicates)) %>%
    dplyr::rename("bio_rep" = "biological_replicates",
                  "tech_rep" = "technical_replicates",
                  "gene_quantification_id" = "accession")
  
  return(tsv_files_info)
}


#' Extract metadata from the case/control experiments
#'
#' @param summary_tg_df Local target gene information from the summary
#'   data.frame generated with [generateSummary()].
#' @param valid_nucleic_acid_type (Optional) Required nucleic acid type of the
#'   experiment. Defaults to "polyadenylated mRNA".
#' @param valid_file_format (Optional) Required output file format of the
#'   sample. Defaults to "BAM".
#' @param valid_output_type (Optional) Required output type of the sample.
#'   Defaults to "alignments".
#' @param valid_genome_annotation (Optional) Required gene annotation version of
#'   the sample to extract its metadata. Defaults to "V29".
#'
#' @return Data.frame containing the metadata for case and control experiments.
#' @export
LoopExperiments <- function(summary_tg_df,
                            valid_nucleic_acid_type,
                            valid_file_format,
                            valid_output_type,
                            valid_genome_annotation){
  ## Loop through every gene silencing series found for the target gene
  tg_df <- foreach(iter_exp = seq(nrow(summary_tg_df)), .combine = dplyr::bind_rows) %do%{
    exp <- summary_tg_df[iter_exp, ]
    
    ## Experiment information
    exp_id <- exp$experiment_id
    gss_id <- exp$gene_silencing_series
    exp_cell_line <- exp$cell_line
    
    logger::log_info("\t Starting Experiment ", exp_id, " (cell line = ", exp_cell_line, "):")
    
    ## Get the API response for the experiment
    response_case <- getUrlResponse(paste0("https://www.encodeproject.org/experiment/", exp_id, "?format=json"))
    control_id = response_case$possible_controls$accession
    response_control = getUrlResponse(paste0("https://www.encodeproject.org/experiment/", control_id, "?format=json"))
    
    case_df <- extractMetadataExperiment(response_case, 
                                         "case", 
                                         valid_nucleic_acid_type, 
                                         valid_file_format,
                                         valid_output_type,
                                         valid_genome_annotation)
    control_df <- extractMetadataExperiment(response_control, 
                                         "control", 
                                         valid_nucleic_acid_type, 
                                         valid_file_format,
                                         valid_output_type,
                                         valid_genome_annotation)
    
    ## Combine cases and controls
    exp_df <- dplyr::bind_rows(case_df, control_df)
    
    ## Add the relevant information
    exp_df <- exp_df %>% mutate(cell_line = exp_cell_line, gene_silencing_series = gss_id, .before = sample_id)
  }
  
  return(tg_df)
}


extractMetadataExperiment <- function(experiment, 
                                      experiment_type,
                                      valid_nucleic_acid_type,
                                      valid_file_format,
                                      valid_output_type,
                                      valid_genome_annotation){
  experiment_id = experiment$accession
  
  logger::log_info("\t\t Starting experiment ", experiment_id, " (type = ", experiment_type, ").")
  
  ## Requirements metadata
  experiment_additional_info <- getAdditionalInformation(experiment)
  
  ## Main metadata
  experiment_sample_files <- getSampleFiles(experiment$files, valid_file_format, valid_output_type, valid_genome_annotation)
  if(experiment_additional_info$nucleic_acid_type != valid_nucleic_acid_type) return(cbind(experiment_sample_files, experiment_additional_info))
  
  ## Other metadata
  experiment_rin <- getRin(experiment$replicates)
  experiment_read_depth <- getReadDepth(experiment$analyses, valid_genome_annotation)
  experiment_donor_info <- getDonorInfo(experiment$replicates)
  experiment_documents <- getDocumentFiles(experiment$replicates, gss$target_gene, experiment_type)
  experiment_doi <- experiment$doi
  experiment_gene_quantifications <- getGeneQuantificationFiles(experiment$files)
  
  ## Combine all information
  experiment_combined <- experiment_sample_files %>% 
    left_join(experiment_rin, by = "bio_rep") %>%
    left_join(experiment_read_depth, by = "bio_rep") %>%
    left_join(experiment_donor_info, by = "bio_rep") %>% 
    left_join(experiment_gene_quantifications %>% select(gene_quantification_id, bio_rep), by = "bio_rep") %>%
    left_join(experiment_documents, by = "bio_rep") %>%
    mutate(experiment_type = experiment_type, experiment_id = experiment_id, experiment_doi = experiment_doi, .before = "bio_rep") %>%
    cbind(experiment_additional_info)
  
}

#' Adds the functional category to the target genes
#'
#' Given the metadata dataframe, it adds a column to register the functional
#' category of the different target genes. A column is added for each category
#' being considered: Splicing regulation, Spliceosome, Exon junction complex or
#' Nononsense-mediated decay.
#'
#' @param metadata_df Data.frame where every row is an ENCODE sample and the
#'   columns contain their metadata.
#' @param input_Category Path to the file containing the functional category
#'   information obtained from Nostrand et. al. publication
#'   (https://www.nature.com/articles/s41586-020-2077-3).
#' @param input_NMD Path to the file containing the NMD functional category
#'   information as a list of genes from this category.
#' @param output_file Path to the output file where the metadata data.frame will
#'   be stored.
#'
#' @return Data.frame where every row is an ENCODE sample and the columns
#'   contain their metadata, including their functional category.
#' @export
addTargetGeneCategory <- function(metadata_df,
                                  input_Category = "",
                                  input_NMD = "",
                                  output_file = ""){
  target_RBPs_metadata <- readr::read_delim(input_Category, show_col_types = F)
  
  metadata_df <- metadata_df %>% 
    left_join(target_RBPs_metadata %>% 
                select(name, `Splicing regulation`, Spliceosome, `Exon Junction Complex`) %>%
                rename("Splicing_regulation" = "Splicing regulation",
                       "Exon_junction_complex" = "Exon Junction Complex"),
              by = c("target_gene" = "name"))
  
  if(input_NMD != ""){
    NMD_list <- readr::read_delim(input_NMD, show_col_types = F, delim = "\n") %>% dplyr::rename("Name" = "Name\t") %>% mutate(across(Name, str_replace, "\t", "")) %>% pull(Name)
    metadata_df <- metadata_df %>%
      mutate(NMD = ifelse(target_gene %in% NMD_list, 1, 0))
  }
  
  if(output_file != ""){
    write.table(metadata_df, output_file, sep = "\t", row.names = F, quote = FALSE)
  }
  
  return(metadata_df)
}

# addArchivedMetadata <- function(metadata_df,
#                                 metadata_archived_path,
#                                 input_target_gene_NMD,
#                                 output_file = ""){
#   metadata_archived <- readr::read_delim(metadata_archived_path, show_col_types = F)
#   NMD_list <- readr::read_delim(input_target_gene_NMD, show_col_types = F, delim = "\n") %>% dplyr::rename("Name" = "Name\t") %>% mutate(across(Name, str_replace, "\t", "")) %>% pull(Name)
#   
#   metadata_archived <- metadata_archived %>%
#     mutate(NMD = ifelse(target_gene %in% NMD_list, 1, 0))
#   
#   removed_genes <- setdiff(metadata_archived$target_gene, metadata_df$target_gene)
#   metadata_combined <- dplyr::bind_rows(metadata_df %>%
#                                           mutate(age = as.numeric(age)),
#                                         metadata_archived %>%
#                                           filter(target_gene %in% removed_genes))
#   if(output_file != ""){
#     write.table(metadata_combined, output_file, sep = "\t", row.names = F, quote = FALSE)
#   }
#   return(metadata_combined)
# }