-   <a href="#introduction" id="toc-introduction">1 Introduction</a>
-   <a href="#methods" id="toc-methods">2 Methods</a>
    -   <a href="#prerequisites" id="toc-prerequisites">2.1 Prerequisites</a>
    -   <a href="#steps" id="toc-steps">2.2 Steps</a>
        -   <a href="#step-0-set-the-extraction-parameters"
            id="toc-step-0-set-the-extraction-parameters">Step 0: Set the extraction
            parameters</a>
        -   <a href="#step-1-call-the-encode-api"
            id="toc-step-1-call-the-encode-api">Step 1: Call the ENCODE API</a>
        -   <a href="#step-2-generate-the-search-summary"
            id="toc-step-2-generate-the-search-summary">Step 2: Generate the search
            summary</a>
        -   <a href="#step-3-retrieval-of-metadata"
            id="toc-step-3-retrieval-of-metadata">Step 3: Retrieval of metadata</a>
        -   <a href="#AddCat" id="toc-AddCat">Step 4: Add target gene category</a>
-   <a href="#results" id="toc-results">3 Results</a>
    -   <a href="#common-metadata" id="toc-common-metadata">3.1 Common
        metadata</a>
    -   <a href="#statistics" id="toc-statistics">3.2 Statistics</a>
-   <a href="#conclusions" id="toc-conclusions">4 Conclusions</a>
-   <a href="#appendix1" id="toc-appendix1">5 Appendix 1: CRISPR metadata
    download</a>
-   <a href="#session-info" id="toc-session-info">6 Session info</a>

<style type="text/css">
.dataTables_scrollHeadInner{
  width:100% !important;
}
.dataTables_scrollHeadInner table{
  width:100% !important;
}
.code-folding-btn {
  display: none;
}
</style>

# 1 Introduction

Print markdown_github-yaml_metadata_block

In this tutorial, we will use [the ENCODE REST
API](https://www.encodeproject.org/help/rest-api/) to automate the
metadata extraction of experiments related to shRNA knockdown against
specific target genes. The functions provided can also be used to extend
its functionality (see [appendix 1](#appendix1)).

# 2 Methods

## 2.1 Prerequisites

For this tutorial, we first need the following libraries:

``` r
## Required libraries
library(httr)
library(jsonlite)
library(foreach)
library(tibble)
library(logger)
library(tidyverse)
library(here)

## Helper functions
source(here::here("Helper_Functions/hf_MetadataDownloadExtraction.R"))

## Logger options
logger::log_threshold(logger::WARN, index = 1)
logger::log_threshold(logger::INFO, index = 2)

logger_layout <- logger::layout_glue_generator(format = '[{time}] [{level}] {msg}')
logger::log_layout(logger_layout, index = 1)
logger::log_layout(logger_layout, index = 2)
```

## 2.2 Steps

### Step 0: Set the extraction parameters

In this step, we need to set the parameters that will determine the
program’s behaviour. The different options are:

-   Valid or required values:
    -   **required_cell_lines**: set to `c("HepG2", "K562")`. It
        requires each target gene to contain samples for both HepG2 and
        K562 cell lines.
    -   **valid_target_genes**: the specific target genes to focus the
        search on. Leave as an empty vector `c()` to not filter.
    -   **valid_genome annotation**: set to `V29`.
    -   **valid_file_format**: set to `bam`.
    -   **valid_output_type**: set to `alignments`.
    -   **valid_nucleic_acid_type**: set to `polyadenylated mRNA`.
-   Output files:
    -   **output_search**: location to store the results of the ENCODE
        search. Set to `Metadata_Results/all_experiments.tsv`
    -   **output_metadata**: location to store the metadata extraction
        results. Set to `Metadata_Results/metadata_samples.tsv`
-   Input files:
    -   **input_target_gene_categories**: path to file extracted from
        *Van Nostrand et. al.* publication with information about the
        functional category of the target genes. More information
        [here](#AddCat).
    -   **input_target_gene_NMD**: path to file which contain the list
        of genes categorized as relevant for nononsense-mediated decay
        processes.
-   Other parameters:
    -   **overwrite_db**: whether to remove the previous metadata found.
        Set to `TRUE` if some API call failed during the process so that
        it only download what the program needs. Defaults to `FALSE`.
    -   **download_method**: it can only be set to
        `gene_silencing_series` or `experiments`. It controls how the
        script access the API to download the information. Using the
        gene silencing series, we can access every control and case
        sample for each cell line with just one call, while using the
        experiments themselves needs two calls (one for cases and one
        for controls). It is recommended to set to gene silencing
        series, unless the ENCODE portal has set the access as
        restricted for them.

``` r
## Valid or required values
required_cell_lines <- c("HepG2", "K562")
valid_target_genes <- c() 
valid_genome_annotation <- "V29"
valid_file_format <- "bam"
valid_output_type <- "alignments"
valid_nucleic_acid_type <- "polyadenylated mRNA"

## Other parameters:
overwrite_db <- F
download_method <- "gene_silencing_series" # Other valid option is "experiments"

## Files
main_path <- here::here("Metadata_Results/")
if(!dir.exists(main_path)) dir.create(main_path)

#### Output files
output_json <- paste0(main_path, "response.json")
output_search <- paste0(main_path, "all_experiments.tsv")
output_metadata <- paste0(main_path, "metadata_samples.tsv")

#### Input files
input_target_gene_categories <- here::here("Additional_Files/Target_gene_categories.tsv")
input_target_gene_NMD <- here::here("Additional_Files/NMD.txt")

## Additional logger info
log_file <- paste0(main_path, "Metadata_Extraction.log")
logger::log_appender(logger::appender_file(log_file, append = T), index = 2)
```

### Step 1: Call the ENCODE API

First, we need to generate an intermediary dataframe containing the
information of an ENCODE search. We head to the [ENCODE experiment
search portal](https://www.encodeproject.org/search/?type=Experiment),
and input the different search filters we are interested in. The
following filters are what the program has been tested with (additional
search executed in [appendix 1](#appendix1)):

-   **Assay**:
    -   **Assay title**: shRNA RNA-seq
    -   **Target category**: RNA binding protein
    -   **Hide control experiments**: yes
-   **Biosample**:
    -   **Organism**: *Homo sapiens*
-   **Quality**:
    -   **Status**: released
-   **Other filters**:
    -   **Data Type**: Experiment

Once we have set the filters, we copy the URL and add
`&limit=all&format=json` to return all search results in `json` format.

``` r
URL = "https://www.encodeproject.org/search/?type=Experiment&assay_title=shRNA+RNA-seq&target.investigated_as=RNA+binding+protein&control_type!=*&status=released&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&limit=all&format=json"
```

Once we have decided on the URL, we execute the `getUrlResponse()`
function to call the ENCODE API and convert the returned json object
into an `R` list.

``` r
response_data <- getUrlResponse(URL, output_json)
```

The returned object contains all the information about the search
results:

``` r
names(response_data)
```

    ##  [1] "@context"      "@graph"        "@id"           "@type"        
    ##  [5] "clear_filters" "columns"       "facet_groups"  "facets"       
    ##  [9] "filters"       "notification"  "sort"          "title"        
    ## [13] "total"

### Step 2: Generate the search summary

Once we have the json object from the ENCODE API, we can extract the
relevant information and summarize it in a dataframe. To do so, we
execute the function `generateSummary()`, where we specify an output
path to store the results in a .tsv file format.

``` r
summary_df <- generateSummary(response_data, 
                              valid_target_genes = valid_target_genes,
                              output_file = output_search)
```

The summary dataframe contains a target gene and gene silencing series
per row, and have the case experiment IDs and cell line for that given
target gene. We store the gene silencing series because it contains the
information of both the case and control samples, while the experiment
IDs only provide information for the case samples. Here is an example of
the 10 first elements of the dataframe:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Experiment ID
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cell line
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Gene silencing series
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
AATF
</td>
<td style="text-align:left;">
ENCSR424YSV
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR634FMS
</td>
</tr>
<tr>
<td style="text-align:left;">
AATF
</td>
<td style="text-align:left;">
ENCSR973QSV
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR146UCW
</td>
</tr>
<tr>
<td style="text-align:left;">
ABCF1
</td>
<td style="text-align:left;">
ENCSR610VTA
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR266ZIC
</td>
</tr>
<tr>
<td style="text-align:left;">
ABCF1
</td>
<td style="text-align:left;">
ENCSR721MXZ
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR242RKD
</td>
</tr>
<tr>
<td style="text-align:left;">
ABT1
</td>
<td style="text-align:left;">
ENCSR756VLW
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR456IPG
</td>
</tr>
<tr>
<td style="text-align:left;">
ABT1
</td>
<td style="text-align:left;">
ENCSR233UVM
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR954JGH
</td>
</tr>
<tr>
<td style="text-align:left;">
ACO1
</td>
<td style="text-align:left;">
ENCSR511SYK
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR339RBD
</td>
</tr>
<tr>
<td style="text-align:left;">
ADAR
</td>
<td style="text-align:left;">
ENCSR104OLN
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR553SDV
</td>
</tr>
<tr>
<td style="text-align:left;">
ADAR
</td>
<td style="text-align:left;">
ENCSR164TLB
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR415ETS
</td>
</tr>
<tr>
<td style="text-align:left;">
AGO1
</td>
<td style="text-align:left;">
ENCSR533HXS
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR213SDO
</td>
</tr>
</tbody>
</table>

### Step 3: Retrieval of metadata

In this last step, we loop through every row of the summary dataframe.
Since we are using the *Gene silencing series*, only one ENCODE API call
for cell line is require to extract all the information of both case and
control samples.

``` r
metadata_df <- generateMetadata(summary_df, 
                                download_method = download_method,
                                required_cell_lines = required_cell_lines,
                                valid_file_format = valid_file_format,
                                valid_genome_annotation = valid_genome_annotation,
                                valid_output_type = valid_output_type,
                                valid_nucleic_acid_type = valid_nucleic_acid_type,
                                output_file = output_metadata,
                                overwrite_db = overwrite_db)
```

The argument `overwrite_db` of the function determines whether the
previous file will be updated or overwritten. If set to `False`, the
function will only retrieve the information for the missing target
genes, or those which do not have exactly 8 entries. This is because
sometimes the API will return bad responses, and we may not want to
repeat the process for all target genes.

The other parameters regulate which files or samples we are interested
in.

The returned dataframe contains the following columns:

-   **target_gene**: the target gene of the shRNA knockdown.
-   **experiment_type**: whether the experiment is case or control for
    that particular cell line and target gene.
-   **cell_line**: the cell type to which the knockdown was executed.
-   **gene_silencing_series**: the ENCODE ID of the *Gene Silencing
    Series* (*GSS*). Information about the *GSS* can be found in the
    link
    `https://www.encodeproject.org/gene-silencing-series/[Gene_silencing_series]/`
    (i.e. [ENCSR406ZEB](https://www.encodeproject.org/gene-silencing-series/ENCSR406ZEB/))
-   **experiment_id**: the ENCODE ID of the experiment in particular.
    Information about the experiment can be found in the link
    `https://www.encodeproject.org/experiments/[Experiment_ID]/`
    (i.e. [ENCSR047QHX](https://www.encodeproject.org/experiments/ENCSR047QHX/))
-   **sample_id**: the ENCODE ID of the sample in particular.
    Information about the sample can be found in the link
    `https://www.encodeproject.org/files/[Sample_ID]/`
    (i.e. [ENCFF946HGK](https://www.encodeproject.org/files/ENCFF946HGK/))
-   **RIN**: the RNA integrity number of the sample.
-   **read_depth**: number of reads per sample.
-   **bio_rep**: the isogenic replicate of the sample.
-   **tech_rep**: the technical replicate of the sample.
-   **sex**: the sex of the sample donor.
-   **age**: the age of the sample donor.
-   **life_stage**: the life stage of the sample donor.
-   **gene_quantification_id**: the ENCODE ID of the associated gene
    quantification file. It is employed to download the gene expressions
    (in TPM) before and after the knockdown, to study the efficiency.
-   **file_format**: format of the selected sample.
-   **output_type**: output type of the selected sample.
-   **genome_annotation**: genome annotation of the selected sample.
-   **mapped_run_type**: mapped run type of the selected sample.
-   **lab**: laboratory where the selected sample was produced.
-   **assay**: assay of the selected sample.
-   **cellosaurus**: cellosaurus ID of the sample. Usually describes the
    cell line.
-   **biosample**: the ENCODE ID of the biosample from which a case
    sample was generated.
-   **document**: ID of the biosample preparation and characterization
    document from which the reported knockdown efficiency from ENCODE is
    extracted.

More information about the metadata provided by ENCODE can be found in
their [*Terms and
Definitions*](https://www.encodeproject.org/data-standards/terms/)
portal.

### Step 4: Add target gene category

Additionally, we can add a target gene category based on *Van Nostrand
et. al.* [A large-scale binding and functional map of human RNA-binding
proteins](https://www.nature.com/articles/s41586-020-2077-3)
publication, which adds a category for whether it is relevant for
*Splicing regulation*, the *Spliceosome* or an *Exon Junction Complex*
(not mutually exclusive). If provided with a list of NMD genes, it will
also classify the target genes in this category.

``` r
metadata_df <- addTargetGeneCategory(metadata_df,
                                     input_target_gene_categories,
                                     input_target_gene_NMD,
                                     output_metadata)
```

The final dataframe is stores in a .tsv file if provided in the
`output_file` argument.

# 3 Results

An example of the final results can be seen in the following table:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
target_gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
experiment_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
cell_line
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
gene_silencing_series
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
experiment_id
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
experiment_doi
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
sample_id
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
rin
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
read_depth
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
bio_rep
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
tech_rep
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
sex
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
age
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
life_stage
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
gene_quantification_id
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
file_format
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
output_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
genome_annotation
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
mapped_run_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
lab
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
assay
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
cellosaurus
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
biosample
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
document
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
biosample_alias
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
nucleic_acid_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
extraction_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
fragmentation_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
size_selection_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
strand_specificity
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
Splicing_regulation
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
Spliceosome
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
Exon_junction_complex
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
NMD
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
HNRNPC
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR308YXN
</td>
<td style="text-align:left;">
ENCSR305XWT
</td>
<td style="text-align:left;">
10.17989/ENCSR305XWT
</td>
<td style="text-align:left;">
ENCFF857QEU
</td>
<td style="text-align:right;">
9.9
</td>
<td style="text-align:right;">
22724294
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
child
</td>
<td style="text-align:left;">
ENCFF570CRU
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0027
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
FAM120A
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR712BXM
</td>
<td style="text-align:left;">
ENCSR661HEL
</td>
<td style="text-align:left;">
10.17989/ENCSR661HEL
</td>
<td style="text-align:left;">
ENCFF114NYX
</td>
<td style="text-align:right;">
9.3
</td>
<td style="text-align:right;">
20037853
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
2_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF146WBD
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
GPKOW
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR888VLT
</td>
<td style="text-align:left;">
ENCSR424QCW
</td>
<td style="text-align:left;">
10.17989/ENCSR424QCW
</td>
<td style="text-align:left;">
ENCFF525ALF
</td>
<td style="text-align:right;">
9.4
</td>
<td style="text-align:right;">
27735686
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF812VKX
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
ESF1
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR484XYM
</td>
<td style="text-align:left;">
ENCSR032YMP
</td>
<td style="text-align:left;">
10.17989/ENCSR032YMP
</td>
<td style="text-align:left;">
ENCFF162FQZ
</td>
<td style="text-align:right;">
9.9
</td>
<td style="text-align:right;">
32549372
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
2_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF396VZW
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
CPSF7
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR664RDU
</td>
<td style="text-align:left;">
ENCSR667PLJ
</td>
<td style="text-align:left;">
10.17989/ENCSR667PLJ
</td>
<td style="text-align:left;">
ENCFF554WDW
</td>
<td style="text-align:right;">
9.7
</td>
<td style="text-align:right;">
16099743
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF096IBO
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RRP9
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR420UVJ
</td>
<td style="text-align:left;">
ENCSR661HEL
</td>
<td style="text-align:left;">
10.17989/ENCSR661HEL
</td>
<td style="text-align:left;">
ENCFF114NYX
</td>
<td style="text-align:right;">
9.3
</td>
<td style="text-align:right;">
20037853
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
2_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF146WBD
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS3A
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR803YTK
</td>
<td style="text-align:left;">
ENCSR419JMU
</td>
<td style="text-align:left;">
10.17989/ENCSR419JMU
</td>
<td style="text-align:left;">
ENCFF448GQT
</td>
<td style="text-align:right;">
9.8
</td>
<td style="text-align:right;">
40534571
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
2_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF240MKR
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
PPIG
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
ENCSR919QPP
</td>
<td style="text-align:left;">
ENCSR529MBZ
</td>
<td style="text-align:left;">
10.17989/ENCSR529MBZ
</td>
<td style="text-align:left;">
ENCFF999LFP
</td>
<td style="text-align:right;">
10.0
</td>
<td style="text-align:right;">
37116348
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
adult
</td>
<td style="text-align:left;">
ENCFF742KEN
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
<td style="text-align:left;">
ENCLB479ITK
</td>
<td style="text-align:left;">
/documents/f2485a7c-5fa3-4253-bfba-16cea47b6cc9/
</td>
<td style="text-align:left;">
PPIG_BGKLV29-47
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
XRN1
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR786TDX
</td>
<td style="text-align:left;">
ENCSR225PRV
</td>
<td style="text-align:left;">
10.17989/ENCSR225PRV
</td>
<td style="text-align:left;">
ENCFF085MQH
</td>
<td style="text-align:right;">
9.6
</td>
<td style="text-align:right;">
21601625
</td>
<td style="text-align:right;">
2
</td>
<td style="text-align:left;">
2_1
</td>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
child
</td>
<td style="text-align:left;">
ENCFF881HQM
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0027
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
PRPF6
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
ENCSR196YEP
</td>
<td style="text-align:left;">
ENCSR491FOC
</td>
<td style="text-align:left;">
10.17989/ENCSR491FOC
</td>
<td style="text-align:left;">
ENCFF055FJV
</td>
<td style="text-align:right;">
9.4
</td>
<td style="text-align:right;">
27706727
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:left;">
1_1
</td>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
child
</td>
<td style="text-align:left;">
ENCFF803IXT
</td>
<td style="text-align:left;">
bam
</td>
<td style="text-align:left;">
alignments
</td>
<td style="text-align:left;">
V29
</td>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0027
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0
</td>
</tr>
</tbody>
</table>

## 3.1 Common metadata

As seen in the table, every row corresponds to a sample related to both
a target gene, a cell line and an experiment type. Given our
requirements, all rows have a file format of `bam`, an output type of
`alignments` and a genome annotation of `V29`. When we study common
aspects between the different samples, we observe a clear difference
between the two cell lines: **all samples extracted from a same cell
line comes from the same donor (sex and age), and from the same tissue
(explained in the cellosaurus).**

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cell line
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Sex
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Age
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cellosaurus
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
male
</td>
<td style="text-align:left;">
15
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0027
</td>
</tr>
<tr>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
female
</td>
<td style="text-align:left;">
53
</td>
<td style="text-align:left;">
Cellosaurus:CVCL_0004
</td>
</tr>
</tbody>
</table>

Each cellosaurus corresponds to a different tissue:

-   **[CVCL_0004](https://www.cellosaurus.org/CVCL_0004)**: derived from
    *Pleural effusion*.
-   **[CVCL_0027](https://www.cellosaurus.org/CVCL_0027)**: derived from
    *Liver*.

Other parameters that are constant across the different samples
(independently of the cell line) are: the mapped run type, the lab, the
assay, the nucleic acid type (required), the extraction method, the
fragmentation method, the size selection method and the strand
specificity.

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
mapped_run_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
lab
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
assay
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
nucleic_acid_type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
extraction_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
fragmentation_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
size_selection_method
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
strand_specificity
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
paired-ended
</td>
<td style="text-align:left;">
Brenton Graveley, UConn
</td>
<td style="text-align:left;">
shRNA knockdown followed by RNA-seq
</td>
<td style="text-align:left;">
polyadenylated mRNA
</td>
<td style="text-align:left;">
Maxwell 16 LEV simpleRNA Cells Kit (Promega cat#: AS1270)
</td>
<td style="text-align:left;">
chemical (Illumina TruSeq)
</td>
<td style="text-align:left;">
SPRI beads
</td>
<td style="text-align:left;">
reverse
</td>
</tr>
</tbody>
</table>

## 3.2 Statistics

The total number of target gene founds under our conditions is 160 (for
a total of 1284 samples).

The average RIN is 9.62 (ranging between 7 and 10), considerably high.
We do find an experiment without a RIN (experiment ENCSR438UOT, case
samples for target gene XPO1 and cell line HepG2), which we define as
`NA`.

If we added the target gene categories, we have a total of 38 marked as
*Splicing regulation*, 19 marked as *Spliceosome*, 1 marked as *Exon
junction complex* and 7 marked as *NMD*. If we combine all these
categories, we obtained a total of 51 target genes:

    ##  [1] "ADAR"     "AQR"      "BUD13"    "CELF1"    "DAZAP1"   "EFTUD2"  
    ##  [7] "EIF4G1"   "EWSR1"    "FUBP1"    "GEMIN5"   "GPKOW"    "HNRNPC"  
    ## [13] "HNRNPU"   "KHDRBS1"  "KHSRP"    "MAGOH"    "MATR3"    "NCBP2"   
    ## [19] "NONO"     "PABPC1"   "PCBP1"    "PCBP2"    "PPIG"     "PRPF4"   
    ## [25] "PRPF6"    "PSIP1"    "PTBP1"    "PUF60"    "QKI"      "RAVER1"  
    ## [31] "RBM15"    "RBM22"    "RBM39"    "RPS10"    "RPS19"    "RPS3A"   
    ## [37] "SART3"    "SF1"      "SF3A3"    "SF3B4"    "SMN1"     "SMNDC1"  
    ## [43] "SND1"     "SNRNP200" "SRSF1"    "SUGP2"    "TARDBP"   "TIAL1"   
    ## [49] "U2AF1"    "U2AF2"    "ZRANB2"

# 4 Conclusions

With the developed script, it is possible to automate the ENCODE RBP
metadata extraction using the [the ENCODE REST
API](https://www.encodeproject.org/help/rest-api/). It is also proven
that the samples share relevant aspects to allow their direct
comparison.

# 5 Appendix 1: CRISPR metadata download

Additionally, the software was tested to download the metadata from
CRISPR knockdown. There are a few considerations we would have to keep
in mind:

-   Some gene silencing series from the CRISPR projects are not public.
    As such, we have to set the `download_method` to `experiments`, so
    that the metadata is extracting from the experiment page, and not
    from the gene silencing series.
-   We need to provide a different URL for the search. More precisely,
    we need to modify the “Assay title” from `shRNA RNA-seq` to
    `CRISPR RNA-seq`. Other than that, everything else is the same.

With these modifications in mind, the following script will download all
metadata found for CRISPR projects where two cell lines are found, each
with 4 samples (2 case and 2 control).

``` r
## Valid or required values
required_cell_lines <- c("HepG2", "K562")
valid_target_genes <- c() 
valid_genome_annotation <- "V29"
valid_file_format <- "bam"
valid_output_type <- "alignments"
valid_nucleic_acid_type <- "polyadenylated mRNA"

## Other parameters:
overwrite_db <- F
download_method <- "experiments"

## Files
main_path <- here::here("Metadata_Results_CRISPR/")
if(!dir.exists(main_path)) dir.create(main_path)

#### Output files
output_json <- paste0(main_path, "response.json")
output_search <- paste0(main_path, "all_experiments.tsv")
output_metadata <- paste0(main_path, "metadata_samples.tsv")

#### Input files
input_target_gene_categories <- here::here("Additional_Files/Target_gene_categories.tsv")
input_target_gene_NMD <- here::here("Additional_Files/NMD.txt")

## URL to CRISPR experiments
URL = "https://www.encodeproject.org/search/?status=released&type=Experiment&target.investigated_as=RNA+binding+protein&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=CRISPR+RNA-seq&limit=all&format=json"

## Pipeline
response_data <- getUrlResponse(URL, output_json)
summary_df <- generateSummary(response_data, 
                              valid_target_genes = valid_target_genes,
                              output_file = output_search)
metadata_df <- generateMetadata(summary_df, 
                                download_method = download_method,
                                required_cell_lines = required_cell_lines,
                                valid_file_format = valid_file_format,
                                valid_genome_annotation = valid_genome_annotation,
                                valid_output_type = valid_output_type,
                                valid_nucleic_acid_type = valid_nucleic_acid_type,
                                output_file = output_metadata, 
                                overwrite_db = overwrite_db)
metadata_df <- addTargetGeneCategory(metadata_df,
                                     input_target_gene_categories,
                                     input_target_gene_NMD,
                                     output_metadata)
```

# 6 Session info

<details>
<summary>
Show/hide
</summary>

    ## ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23)
    ##  os       Ubuntu 20.04.4 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       Etc/UTC
    ##  date     2023-03-13
    ##  pandoc   2.18 @ /usr/lib/rstudio-server/bin/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
    ##  package     * version    date (UTC) lib source
    ##  bit           4.0.5      2022-11-15 [1] RSPM (R 4.2.0)
    ##  bit64         4.0.5      2020-08-30 [1] RSPM (R 4.2.0)
    ##  cli           3.6.0      2023-01-09 [1] RSPM (R 4.2.0)
    ##  codetools     0.2-18     2020-11-04 [2] CRAN (R 4.2.1)
    ##  colorspace    2.1-0      2023-01-23 [1] RSPM (R 4.2.0)
    ##  crayon        1.5.2      2022-09-29 [1] RSPM (R 4.2.0)
    ##  curl          5.0.0      2023-01-12 [1] RSPM (R 4.2.0)
    ##  digest        0.6.31     2022-12-11 [1] RSPM (R 4.2.0)
    ##  dplyr       * 1.1.0      2023-01-29 [1] RSPM (R 4.2.0)
    ##  ellipsis      0.3.2      2021-04-29 [1] CRAN (R 4.2.0)
    ##  evaluate      0.20       2023-01-17 [1] RSPM (R 4.2.0)
    ##  fansi         1.0.4      2023-01-22 [1] RSPM (R 4.2.0)
    ##  fastmap       1.1.1      2023-02-24 [1] RSPM (R 4.2.0)
    ##  forcats     * 1.0.0      2023-01-29 [1] RSPM (R 4.2.0)
    ##  foreach     * 1.5.2      2022-02-02 [1] RSPM (R 4.2.0)
    ##  generics      0.1.3      2022-07-05 [1] RSPM (R 4.2.0)
    ##  ggplot2     * 3.4.1      2023-02-10 [1] RSPM (R 4.2.0)
    ##  glue          1.6.2      2022-02-24 [1] CRAN (R 4.2.0)
    ##  gtable        0.3.1      2022-09-01 [1] RSPM (R 4.2.0)
    ##  here        * 1.0.1      2020-12-13 [1] RSPM (R 4.2.0)
    ##  highr         0.10       2022-12-22 [1] RSPM (R 4.2.0)
    ##  hms           1.1.2      2022-08-19 [1] RSPM (R 4.2.0)
    ##  htmltools     0.5.4      2022-12-07 [1] RSPM (R 4.2.0)
    ##  httr        * 1.4.5      2023-02-24 [1] RSPM (R 4.2.0)
    ##  iterators     1.0.14     2022-02-05 [1] RSPM (R 4.2.0)
    ##  jsonlite    * 1.8.4      2022-12-06 [1] RSPM (R 4.2.0)
    ##  kableExtra    1.3.4.9000 2023-01-30 [1] Github (haozhu233/kableExtra@292f607)
    ##  knitr         1.42       2023-01-25 [1] RSPM (R 4.2.0)
    ##  lifecycle     1.0.3      2022-10-07 [1] RSPM (R 4.2.0)
    ##  logger      * 0.2.2      2021-10-19 [1] RSPM (R 4.2.0)
    ##  lubridate   * 1.9.2      2023-02-10 [1] RSPM (R 4.2.0)
    ##  magrittr      2.0.3      2022-03-30 [1] CRAN (R 4.2.0)
    ##  munsell       0.5.0      2018-06-12 [1] RSPM (R 4.2.0)
    ##  pillar        1.8.1      2022-08-19 [1] RSPM (R 4.2.0)
    ##  pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 4.2.0)
    ##  purrr       * 1.0.1      2023-01-10 [1] RSPM (R 4.2.0)
    ##  R6            2.5.1      2021-08-19 [1] CRAN (R 4.2.0)
    ##  readr       * 2.1.4      2023-02-10 [1] RSPM (R 4.2.0)
    ##  rlang         1.0.6      2022-09-24 [1] RSPM (R 4.2.0)
    ##  rmarkdown     2.20       2023-01-19 [1] RSPM (R 4.2.0)
    ##  rprojroot     2.0.3      2022-04-02 [1] CRAN (R 4.2.0)
    ##  rstudioapi    0.14       2022-08-22 [1] RSPM (R 4.2.0)
    ##  rvest         1.0.3      2022-08-19 [1] RSPM (R 4.2.0)
    ##  scales        1.2.1      2022-08-20 [1] RSPM (R 4.2.0)
    ##  sessioninfo * 1.2.2      2021-12-06 [1] RSPM (R 4.2.0)
    ##  stringi       1.7.12     2023-01-11 [1] RSPM (R 4.2.0)
    ##  stringr     * 1.5.0      2022-12-02 [1] RSPM (R 4.2.0)
    ##  svglite       2.1.1      2023-01-10 [1] RSPM (R 4.2.0)
    ##  systemfonts   1.0.4      2022-02-11 [1] RSPM (R 4.2.0)
    ##  tibble      * 3.1.8      2022-07-22 [1] RSPM (R 4.2.0)
    ##  tidyr       * 1.3.0      2023-01-24 [1] RSPM (R 4.2.0)
    ##  tidyselect    1.2.0      2022-10-10 [1] RSPM (R 4.2.0)
    ##  tidyverse   * 2.0.0      2023-02-22 [1] RSPM (R 4.2.0)
    ##  timechange    0.2.0      2023-01-11 [1] RSPM (R 4.2.0)
    ##  tzdb          0.3.0      2022-03-28 [1] RSPM (R 4.2.0)
    ##  utf8          1.2.3      2023-01-31 [1] RSPM (R 4.2.0)
    ##  vctrs         0.5.2      2023-01-23 [1] RSPM (R 4.2.0)
    ##  viridisLite   0.4.1      2022-08-22 [1] RSPM (R 4.2.0)
    ##  vroom         1.6.1      2023-01-22 [1] RSPM (R 4.2.0)
    ##  webshot       0.5.4      2022-09-26 [1] RSPM (R 4.2.0)
    ##  withr         2.5.0      2022-03-03 [1] CRAN (R 4.2.0)
    ##  xfun          0.37       2023-01-31 [1] RSPM (R 4.2.0)
    ##  xml2          1.3.3      2021-11-30 [1] CRAN (R 4.2.0)
    ##  yaml          2.3.7      2023-01-23 [1] RSPM (R 4.2.0)
    ## 
    ##  [1] /usr/local/lib/R/site-library
    ##  [2] /usr/local/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

</details>
