-   <a href="#objective" id="toc-objective">1 Objective</a>
-   <a href="#methods" id="toc-methods">2 Methods</a>
-   <a href="#steps" id="toc-steps">3 Steps</a>
    -   <a href="#prerequisites" id="toc-prerequisites">3.1 Prerequisites</a>
    -   <a href="#step-1-load-previous-metadata-and-set-script-parameters"
        id="toc-step-1-load-previous-metadata-and-set-script-parameters">Step 1:
        Load previous metadata and set script parameters</a>
    -   <a href="#step-2-download-the-biosample-characterization-documents"
        id="toc-step-2-download-the-biosample-characterization-documents">Step
        2: Download the biosample characterization documents</a>
    -   <a href="#step-3-extract-the-images-from-the-documents"
        id="toc-step-3-extract-the-images-from-the-documents">Step 3: Extract
        the images from the documents</a>
    -   <a href="#step-4-extract-the-text-from-the-images"
        id="toc-step-4-extract-the-text-from-the-images">Step 4: Extract the
        text from the images</a>
    -   <a href="#step-5-summarize-and-store-to-disk"
        id="toc-step-5-summarize-and-store-to-disk">Step 5: Summarize and store
        to disk</a>
-   <a href="#results" id="toc-results">4 Results</a>
    -   <a href="#kEff" id="toc-kEff">4.1 Knockdown Efficiency</a>
    -   <a href="#cell-line-influence" id="toc-cell-line-influence">4.2 Cell
        line influence</a>

# 1 Objective

In this tutorial, we will focus on downloading and extracting the
reported knockdown efficiency in the ENCODE shRNA knockdown experiments.

# 2 Methods

First, we needed a list of target RBPs to study. As with the previous
studies, we will employ the list of 356 RBPs published by *Van Nostrand
et. al.* [A large-scale binding and functional map of human RNA-binding
proteins](https://www.nature.com/articles/s41586-020-2077-3) categorized
as splicing regulation, spliceosome or exon-junction complex. We also
added a list of 118 genes classified as involved in nononsense-mediated
decay processes. In total, 56 of those genes were found as ENCODE
projects with the same number of experiments and similar metadata.

For most experiments (being an experiment a collection of 2 isogenic
replicates for one cell line and experiment type), the knockdown
efficiency is reported in their *Graveley Lab shRNA knockdown followed
by RNA-seq Biosample Preparation and Characterization Document* (see
[example](https://www.encodeproject.org/documents/e5391dd5-cf87-451c-97e0-fe9967170285/@@download/attachment/U2AF2-LV08-64.pdf)).
They applied two alternative techniques to evaluate the efficiency:
qRT-PCR and Western blotting.

In both situations, they measure the depletion level in the shRNA
treated sample against a control non-target shRNA sample. The depletion
is normalized using the gene GAPDH as a control. The technique we are
most interested in is the Western blotting, a powerful and important
procedure to analyse the detection of proteins, particularly in low
abundance. As such, it represents the actual concentration of protein in
the sample and not the RNA abundance.

The main procedure to obtain the reported knockdown efficiency for a
particular experiment is the following:

1.  Identify if the report is available for the particular experiment.
    See for example the [ADAR K562
    experiment](https://www.encodeproject.org/experiments/ENCSR164TLB/).
    At the bottom of the page, we need to find any document with the
    description: *biosample protocol for shRNA followed by RNA-seq*. The
    procedure is executed in the ENCODE metadata extraction framework,
    using [the ENCODE REST
    API](https://www.encodeproject.org/help/rest-api/). If the document
    is found, we store the information in the metadata from the samples.

2.  Download the documents for every interest RBP/NMD experiment. Even
    if the document is found, there are many inconsistencies inside the
    ENCODE experiments that make it hard to identify the download link
    (e.g. [some](https://www.encodeproject.org/experiments/ENCSR220TBR/)
    documents has the suffix `_update`,
    [others](https://www.encodeproject.org/experiments/ENCSR605MFS/)
    have an alternative HGCN name for the gene). Because of that, not
    all possible reports are guaranteed to be found.

3.  Since the knockdown efficiency is written in an image inside the
    pdf, we employed the `extract -images`
    [script](https://pymupdf.readthedocs.io/en/latest/module.html#extracting-fonts-and-images)
    from the *fitz* module within the python package *PyMuPDF* version
    1.21.1. There are several figures inside the pdf, but fortunately
    the Western blotting results are always the last figure. In figure
    @ref(fig:wbimages) we can see an example of the figures that contain
    the knockdown efficiency inside the reports.

<img src="../Knockdown_Efficiency/WesternBlotting_PCR/images/ENCLB155EFP_Western_Blot_Analysis.png" alt="Different images in which the knockdown efficiency is reported." width="49%" height="20%" /><img src="../Knockdown_Efficiency/WesternBlotting_PCR/images/ENCLB867RJA_Western_Blot_Analysis.png" alt="Different images in which the knockdown efficiency is reported." width="49%" height="20%" />
<p class="caption">
Different images in which the knockdown efficiency is reported.
</p>

1.  Once the images are extracted, we need to employ an OCR (Optical
    Character Recognition) algorithm to identify the text within the
    figures. To do so, we employed the open-source
    [Tesseract-OCR](https://github.com/tesseract-ocr/tesseract#tesseract-ocr)
    engine currently developed by Google. Its main advantages are the
    simplicity of use, overall efficiency and easy implementation in a
    python pipeline using the package
    [pytesseract](https://pypi.org/project/pytesseract/) version 0.3.10.
    Using the R library
    [reticulate](https://rstudio.github.io/reticulate/) version 1.26, we
    called the python functions to extract the text from the image. Some
    considerations about the procedure:

    -   The images contain a lot of other information and text not
        related to the Western blot results. This hinders the
        identification of the knockdown efficiencies from the output of
        the OCR algorithm. Thus, we decided to crop the image to only
        show the percentages and the method to measure the efficiency.
        See figure @ref(fig:croppedimage).
    -   The images extracted from the ENCODE report are usually in high
        resolution. However, the tesseract-ocr training probably used
        lower resolution images. Therefore, the accuracy was quite low
        in our initial tests. To solve the issue we had to resize the
        images to a lower resolution. We found that the most accurate
        resolution was 25% of the original.
    -   Because of the large amount of total reports from ENCODE, the
        accuracy could not be manually estimated. As such, a 100%
        accuracy is not guaranteed.

<img src="../Knockdown_Efficiency/WesternBlotting_PCR/images/cropped_image.png" alt="Cropped image from which we extract the text." width="49%" height="20%" />
<p class="caption">
Cropped image from which we extract the text.
</p>

1.  The last step was to summarize the results. Since two knockdown
    efficiencies were provided for each experiment, we estimated that
    the knockdown efficiency for that particular cell line was the
    average of the two. In the data-analysis pipeline we made no
    distinction between the two cell lines, thus the final knockdown
    efficiency was estimated as the average of the two cell lines.

# 3 Steps

## 3.1 Prerequisites

For this tutorial, we need to install the
[Tesseract-OCR](https://github.com/tesseract-ocr/tesseract#tesseract-ocr)
engine. We also need the library
[reticulate](https://rstudio.github.io/reticulate/) version 1.26 to use
the python modules [Pillow](https://pillow.readthedocs.io/en/stable/)
and [pytesseract](https://pypi.org/project/pytesseract/).

Therefore, we need to load the necessary libraries and helper functions:

``` r
# Required libraries
library(logger)
library(foreach)
library(doSNOW)
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reticulate))

# Load additional helper functions
source(here::here("Helper_Functions/hf_KnockdownEfficiency_WB.R"))

# Load python interface
use_virtualenv("UCL")
pytesseract <- import("pytesseract")
PIL <- import("PIL")
```

## Step 1: Load previous metadata and set script parameters

First, we specify where the input metadata is located and where the
output files will be stored. We also specify some parameters for the
scripts, like the path to store the downloaded files, the number of
download cores to use, the resize percentage to apply to the figures and
the minimum width to apply the resize (more information on that later).

## Step 2: Download the biosample characterization documents

Next step is to generate the subfolder structures where the biosample
characterization files will be stored (by default
`main_path/target_gene/experiment_id/biosample_replicate/`). We use the
function `downloadCharacterizationDocuments()` to download the documents
into the specified paths.

``` r
# Create the directories
createDirectories(target_RBPs, metadata_filtered)

# Download the files and add a column with their path
metadata_documents <- downloadCharacterizationDocuments(metadata_filtered, 
                                                        download_cores, 
                                                        overwrite_results,
                                                        silent = T)
```

It is recommended to test for the existence of the documents, since it
is possible that some experiments have their documents file names
modified, or they might not be available to download.

``` r
# Check the existence of the files
for(row_index in seq(nrow(metadata_documents))){
  row = metadata_documents[row_index, ]
  file_path <- row$file_path
  path <- row$path
  
  if(!file.exists(file_path) || file.info(file_path)$size < 10){
    logger::log_warn("Error for row ", row_index, "! File path ", path)
  }
}
```

## Step 3: Extract the images from the documents

Once we have all characterization documents available, we use the python
package [PyMuPDF](https://pymupdf.readthedocs.io/en/latest/) version
1.21.1. The function `extractImages()` extract the images and keeps only
the last one, which contains the reported efficiency for both methods.

``` r
# Extract the images of all files
metadata_images <- extractImages(metadata_documents, 
                                 overwrite_results = overwrite_results)
```

## Step 4: Extract the text from the images

The next step is to process the images from the step above and to
extract the written text in them to obtain the reported knockdown
efficiencies. As mentioned, we use the python modules *pytesseract* and
*Pillow*. First, we crop the image and reduce its dimensions
(tesseract-OCR works better with lower resolution images). We apply the
basic tesseract-OCR configuration to the image and simply modify the
returned string to obtain the results in data.frame format.

Since there are two types of images, we need to consider both scenarios,
differenciated by the number of columns that tesseract-OCR reads from
the image. The efficiencies are extracted and combined per cell line
(i.e. two efficiencies are provided per cell line, but the average is
reported). If no valid efficiency is extracted from the image, an `NA`
is returned instead:

``` r
# Extract text from images. This cannot be separated into an external function.
# Probably because of some incompatibility with the reticulate library to use
# python.
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
```

An important step is to test for consistency between the returned
efficiencies. In many experiments, each sample’s document reports the
efficiency for the whole target gene. As such, the same image is
processed up to 4 times per target gene. It is important to test that
all iterations of the same image returns the same output:

``` r
# Test for consistency between the cell lines
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
```

## Step 5: Summarize and store to disk

From the reported efficiencies for each cell line, we calculate the
average of them as the reported knockdown efficiency for a particular
target gene. These reported efficiencies are stored for both Western
blot and qRT-PCR.

``` r
# Write the knockdown efficiency table to disk
metadata_WB_kEff <- writeEfficiencyTable(metadata_kEff, "WB", metadata_WB_output)
metadata_PCR_kEff <- writeEfficiencyTable(metadata_kEff, "PCR", metadata_PCR_output)
```

# 4 Results

## 4.1 Knockdown Efficiency

In the following table, we can observe the obtained knockdown
efficiencies for every studied target gene:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Category
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
kEff \[%\]
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PTBP1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
99.00
</td>
</tr>
<tr>
<td style="text-align:left;">
PSIP1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
95.50
</td>
</tr>
<tr>
<td style="text-align:left;">
KHSRP
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
94.00
</td>
</tr>
<tr>
<td style="text-align:left;">
EIF4G1
</td>
<td style="text-align:left;">
NMD
</td>
<td style="text-align:right;">
93.50
</td>
</tr>
<tr>
<td style="text-align:left;">
SART3
</td>
<td style="text-align:left;">
Spliceosome
</td>
<td style="text-align:right;">
92.75
</td>
</tr>
<tr>
<td style="text-align:left;">
SUGP2
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
91.50
</td>
</tr>
<tr>
<td style="text-align:left;">
EWSR1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
90.75
</td>
</tr>
<tr>
<td style="text-align:left;">
GEMIN5
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
87.75
</td>
</tr>
<tr>
<td style="text-align:left;">
PPIG
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
86.00
</td>
</tr>
<tr>
<td style="text-align:left;">
SRSF1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
85.50
</td>
</tr>
<tr>
<td style="text-align:left;">
HNRNPC
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
85.25
</td>
</tr>
<tr>
<td style="text-align:left;">
FUBP1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
84.00
</td>
</tr>
<tr>
<td style="text-align:left;">
BUD13
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
82.50
</td>
</tr>
<tr>
<td style="text-align:left;">
DAZAP1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
81.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RBM22
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
81.00
</td>
</tr>
<tr>
<td style="text-align:left;">
ZRANB2
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
81.00
</td>
</tr>
<tr>
<td style="text-align:left;">
SND1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
78.75
</td>
</tr>
<tr>
<td style="text-align:left;">
RAVER1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
77.50
</td>
</tr>
<tr>
<td style="text-align:left;">
U2AF2
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
76.75
</td>
</tr>
<tr>
<td style="text-align:left;">
ADAR
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
74.75
</td>
</tr>
</tbody>
</table>

From the 51 studied projects, 5 of them had no available information for
any of the two cell lines. The average knockdown efficiency is 72.9%,
ranging from 46.8 to 99% in some cases. In the following visualization,
we represent in the X-axis the knockdown efficiency for every studied
target RBP/NMD, split by their functional category:

<img src="Tutorial_KnockdownEfficiency_WB_files/figure-markdown_github/kEff-category-graph-1.png" width="85%" style="display: block; margin: auto;" />

## 4.2 Cell line influence

In some cases, the difference in `kEff` between the two cell lines is
not irrelevant. In the following table and graph, we show the `kEff` for
each cell line and the difference between the two:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Category
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
HepG2 kEff \[%\]
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
K562 kEff \[%\]
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
Difference abs(HepG2 - K562)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PUF60
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
98.5
</td>
<td style="text-align:right;">
50.0
</td>
<td style="text-align:right;">
48.5
</td>
</tr>
<tr>
<td style="text-align:left;">
CELF1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
91.0
</td>
<td style="text-align:right;">
51.5
</td>
<td style="text-align:right;">
39.5
</td>
</tr>
<tr>
<td style="text-align:left;">
EFTUD2
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
40.5
</td>
<td style="text-align:right;">
78.5
</td>
<td style="text-align:right;">
38.0
</td>
</tr>
<tr>
<td style="text-align:left;">
PRPF6
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
51.5
</td>
<td style="text-align:right;">
86.5
</td>
<td style="text-align:right;">
35.0
</td>
</tr>
<tr>
<td style="text-align:left;">
PCBP2
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
70.0
</td>
<td style="text-align:right;">
44.5
</td>
<td style="text-align:right;">
25.5
</td>
</tr>
<tr>
<td style="text-align:left;">
TIAL1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
86.0
</td>
<td style="text-align:right;">
61.0
</td>
<td style="text-align:right;">
25.0
</td>
</tr>
<tr>
<td style="text-align:left;">
HNRNPC
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
76.5
</td>
<td style="text-align:right;">
94.0
</td>
<td style="text-align:right;">
17.5
</td>
</tr>
<tr>
<td style="text-align:left;">
NONO
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
61.5
</td>
<td style="text-align:right;">
78.5
</td>
<td style="text-align:right;">
17.0
</td>
</tr>
<tr>
<td style="text-align:left;">
MAGOH
</td>
<td style="text-align:left;">
Exon_junction_complex
</td>
<td style="text-align:right;">
77.0
</td>
<td style="text-align:right;">
61.5
</td>
<td style="text-align:right;">
15.5
</td>
</tr>
<tr>
<td style="text-align:left;">
U2AF1
</td>
<td style="text-align:left;">
Splicing_regulation
</td>
<td style="text-align:right;">
73.5
</td>
<td style="text-align:right;">
58.5
</td>
<td style="text-align:right;">
15.0
</td>
</tr>
</tbody>
</table>

<img src="Tutorial_KnockdownEfficiency_WB_files/figure-markdown_github/kEff-cell_line-graph-1.png" width="85%" style="display: block; margin: auto;" />
