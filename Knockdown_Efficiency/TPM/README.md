-   <a href="#objective" id="toc-objective">1 Objective</a>
-   <a href="#methods" id="toc-methods">2 Methods</a>
-   <a href="#steps" id="toc-steps">3 Steps</a>
    -   <a href="#step-0-load-previous-metadata-and-set-script-parameters"
        id="toc-step-0-load-previous-metadata-and-set-script-parameters">Step 0:
        Load previous metadata and set script parameters</a>
    -   <a href="#step-1-download-the-gene_quantification-files"
        id="toc-step-1-download-the-gene_quantification-files">Step 1: Download
        the <code>gene_quantification</code> files</a>
    -   <a href="#step-2-tpm-extraction" id="toc-step-2-tpm-extraction">Step 2:
        TPM extraction</a>
    -   <a href="#step-3-estimating-the-knockdown-efficiency"
        id="toc-step-3-estimating-the-knockdown-efficiency">3.1 Step 3:
        Estimating the Knockdown Efficiency</a>
-   <a href="#results" id="toc-results">4 Results</a>
    -   <a href="#knockdown-efficiency" id="toc-knockdown-efficiency">4.1
        Knockdown Efficiency</a>
    -   <a href="#cell-line-influence" id="toc-cell-line-influence">4.2 Cell
        line influence</a>

# 1 Objective

In this tutorial, we aim to estimate the knockdown efficiency of the
shRNA knockdown targeted genes by comparing the TPM values between
control and cases. This method is not the recommended approach to
measure the knockdown effiency, please refer to the Western Blotting
method extracted from ENCODE.

# 2 Methods

First, we needed a list of target RBPs to study. As with the previous
studies, we will employ the list of 356 RBPs published by *Van Nostrand
et. al.* [A large-scale binding and functional map of human RNA-binding
proteins](https://www.nature.com/articles/s41586-020-2077-3) categorized
as splicing regulation, spliceosome or exon-junction complex. We also
added a list of 118 genes classified as involved in nononsense-mediated
decay processes. In total, 56 of those genes were found as ENCODE
projects with the same number of experiments and similar metadata.

Once we have set the list of ENCODE projects to focus our studies, we
need to download the gene expressions for all samples. We have to
consider:

-   Each sample studied in the ENCODE RBP analysis has an associated
    `gene_quantification.tsv` file to download.
-   The files are the byproduct of using the software *RSEM* over both
    `FASTQ` files. In all ENCODE projects, we find two sets of `FASTQ`
    files which contain the first and second mates of the paired-end
    reads. Since we have 2 isogenic replicates in all experiments, a
    total of 4 `FASTQ` files are processed. To produce these files, one
    can follow the *RMSE*
    [tutorial](https://www.encodeproject.org/experiments/ENCSR896CFV/)
    and the ENCODE [pipeline
    information](https://www.encodeproject.org/documents/6354169f-86f6-4b59-8322-141005ea44eb/@@download/attachment/Long%20RNA-seq%20pipeline%20overview.pdf).
-   Therefore, a total of 4 files for cases and 4 files for controls are
    downloaded. Each with 2 for K562 cell line and 2 for hepG2 cell
    line.
-   Out of the 56 initial target genes, only 51 have their gene
    expressions available to download.

Then, we need to extract the TPM value of the target gene in that
particular experiment.

-   Since the `gene_quantification.tsv` files use the gene ENSEMBL
    notation, we use the `biomaRt` package to translate from the HGNC
    nomenclature to their ENSEMBL Id.
-   The provided tables in `gene_quantification` contain many useful
    columns: TPM, posterior mean estimates of the TPM (using Gibbs
    sampling), lower and upper CI bounds and TPM coefficient of quartile
    variation. We opted to use the TPM column as is.

Once we have the TPM values before and after the knockdown, we evaluate
the knockdown efficiency as defined in equation (1). For that, we first
need to calculate the average TPM values for each cluster
(case/control), independent of the cell line. A knockdown efficiency of
100% implies that the target gene was not found in the case sample,
while a 0% indicates that no change in TPM was observed.

$$
kEff\\;\\% =\left(1-\frac{TPM\_{avg,\\;case}}{TPM\_{avg,\\;control}}\right)\*100\\%\tag{1}
$$

# 3 Steps

## Step 0: Load previous metadata and set script parameters

First, we need to indicate the file from which to obtain the metadata of
the samples. In our case, it is the path to the file generated from the
metadata extraction tutorial. We also need to specify certain parameters
to the algorithm, like the directory to save the
`gene_quantification.tsv` files, or the number of cores to use.

``` r
## Input Files
main_metadata_path <- here::here("Metadata_Results/")
metadata_path <- paste0(main_metadata_path, "metadata_samples.tsv")
metadata_TPM_output <- paste0(main_metadata_path, "metadata_TPM_kEff.tsv")

## Define the algorithm variables
main_path = here::here("Knockdown_Efficiency/TPM/RBPs/")
download_only = F
download_cores = 16
overwrite_results = F

## Generate the variables ----
metadata <- readr::read_delim(metadata_path, show_col_types = F) %>% as_tibble()

target_RBPs <- metadata %>%
  dplyr::filter(if_any(c(Splicing_regulation, Spliceosome, Exon_junction_complex, NMD), ~ . != 0)) %>%
  dplyr::filter(!is.na(gene_quantification_id)) %>%
  dplyr::pull(target_gene) %>%
  unique()

## Convert HGCN nomenclature to ENSEMBL ID
ensembl_target_RBPs <- translateGenes(target_RBPs)

metadata_filtered <- metadata %>%
  dplyr::filter(target_gene %in% target_RBPs, !is.na(gene_quantification_id)) %>%
  dplyr::select(target_gene, cell_line, experiment_type, experiment_id, gene_quantification_id) %>%
  dplyr::mutate(path = paste0(main_path, target_gene, "/", experiment_type, "/")) %>%
  dplyr::left_join(ensembl_target_RBPs, by = c("target_gene" = "hgnc_symbol"), multiple = "all") %>%
  dplyr::relocate(ensembl_gene_id, .before = cell_line)
```

## Step 1: Download the `gene_quantification` files

Next, we use the functions `createDirectories()` and
`downloadGeneQuantifications()` to generate the folder structure where
the gene quantifications will be stored and to download them from the
ENCODE portal. The download process can be parallelized, since most of
the download times is spent stablishing the connection to the server and
not downloading the file itself (meaning that it is unlikely to saturate
the download bandwidth).

``` r
# Create the directories
createDirectories(target_RBPs, metadata_filtered)

# Download the files and add a column with their path
metadata_quantifications <- downloadGeneQuantifications(metadata_filtered,
                                                        download_cores, 
                                                        overwrite_results,
                                                        silent = T)
```

## Step 2: TPM extraction

Using the package `biomaRt`, we converted the HGCN nomenclature of the
genes into their ENSEMBL ID. Then, we use the function `extractTPM()` to
read the gene quantification file and extract the TPM for the target
gene of that particular samples.

``` r
# Extract the TPMs
metadata_TPM <- extractTPM(metadata_quantifications)
```

Results are shown in the following table:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
ENSEMBL ID
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cell line
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Experiment type
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Gene quantification ID
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
TPM
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
PTBP1
</td>
<td style="text-align:left;">
ENSG00000011304
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF597DBP
</td>
<td style="text-align:right;">
124.62
</td>
</tr>
<tr>
<td style="text-align:left;">
KHSRP
</td>
<td style="text-align:left;">
ENSG00000088247
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF774RTC
</td>
<td style="text-align:right;">
135.71
</td>
</tr>
<tr>
<td style="text-align:left;">
QKI
</td>
<td style="text-align:left;">
ENSG00000112531
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF597DBP
</td>
<td style="text-align:right;">
57.14
</td>
</tr>
<tr>
<td style="text-align:left;">
RBM15
</td>
<td style="text-align:left;">
ENSG00000162775
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:left;">
ENCFF773VSW
</td>
<td style="text-align:right;">
3.68
</td>
</tr>
<tr>
<td style="text-align:left;">
BUD13
</td>
<td style="text-align:left;">
ENSG00000137656
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF528UQK
</td>
<td style="text-align:right;">
13.58
</td>
</tr>
<tr>
<td style="text-align:left;">
SF1
</td>
<td style="text-align:left;">
ENSG00000168066
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:left;">
ENCFF452IUK
</td>
<td style="text-align:right;">
46.91
</td>
</tr>
<tr>
<td style="text-align:left;">
TARDBP
</td>
<td style="text-align:left;">
ENSG00000120948
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF587DSG
</td>
<td style="text-align:right;">
211.51
</td>
</tr>
<tr>
<td style="text-align:left;">
EFTUD2
</td>
<td style="text-align:left;">
ENSG00000108883
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:left;">
ENCFF446CYK
</td>
<td style="text-align:right;">
47.17
</td>
</tr>
<tr>
<td style="text-align:left;">
PCBP1
</td>
<td style="text-align:left;">
ENSG00000169564
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:left;">
ENCFF389GXH
</td>
<td style="text-align:right;">
37.27
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS19
</td>
<td style="text-align:left;">
ENSG00000105372
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
control
</td>
<td style="text-align:left;">
ENCFF146WBD
</td>
<td style="text-align:right;">
2625.84
</td>
</tr>
</tbody>
</table>

## 3.1 Step 3: Estimating the Knockdown Efficiency

Lastly, we use the function `generateKnockdownEfficiency()` to apply the
knockdown efficiency formula. Results are obtained for both all samples
at the same time and divided by cell lines. In most cases, the overall
knockdown efficiency is the average of the efficiencies for each cell
line (with some notable exceptions).

``` r
# Generate and save the knockdown efficiencies
metadata_kEff <- generateKnockdownEfficiency(metadata_TPM,
                                             output_file = metadata_TPM_output)
```

Final results are stored in disk:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
kEff \[%\]
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
kEff HepG2 \[%\]
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
kEff K562 \[%\]
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
RPS3A
</td>
<td style="text-align:right;">
86.2
</td>
<td style="text-align:right;">
84.2
</td>
<td style="text-align:right;">
89.5
</td>
</tr>
<tr>
<td style="text-align:left;">
SND1
</td>
<td style="text-align:right;">
86.7
</td>
<td style="text-align:right;">
88.1
</td>
<td style="text-align:right;">
85.9
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS10
</td>
<td style="text-align:right;">
75.3
</td>
<td style="text-align:right;">
82.9
</td>
<td style="text-align:right;">
70.2
</td>
</tr>
<tr>
<td style="text-align:left;">
SUGP2
</td>
<td style="text-align:right;">
50.3
</td>
<td style="text-align:right;">
42.7
</td>
<td style="text-align:right;">
54.2
</td>
</tr>
<tr>
<td style="text-align:left;">
PUF60
</td>
<td style="text-align:right;">
65.1
</td>
<td style="text-align:right;">
69.9
</td>
<td style="text-align:right;">
62.1
</td>
</tr>
<tr>
<td style="text-align:left;">
CELF1
</td>
<td style="text-align:right;">
30.6
</td>
<td style="text-align:right;">
42.0
</td>
<td style="text-align:right;">
22.9
</td>
</tr>
<tr>
<td style="text-align:left;">
RBM39
</td>
<td style="text-align:right;">
70.1
</td>
<td style="text-align:right;">
69.1
</td>
<td style="text-align:right;">
73.0
</td>
</tr>
<tr>
<td style="text-align:left;">
PABPC1
</td>
<td style="text-align:right;">
92.3
</td>
<td style="text-align:right;">
90.2
</td>
<td style="text-align:right;">
93.1
</td>
</tr>
<tr>
<td style="text-align:left;">
PSIP1
</td>
<td style="text-align:right;">
79.5
</td>
<td style="text-align:right;">
75.5
</td>
<td style="text-align:right;">
81.4
</td>
</tr>
<tr>
<td style="text-align:left;">
PRPF4
</td>
<td style="text-align:right;">
73.0
</td>
<td style="text-align:right;">
66.3
</td>
<td style="text-align:right;">
76.5
</td>
</tr>
</tbody>
</table>

# 4 Results

## 4.1 Knockdown Efficiency

The average knockdown efficiency is 65%, ranging from 31 to 92% in some
cases. In the following visualization, we represent in the X-axis the
estimated knockdown efficiency for every studied target RBP/NMD, split
by their functional category:

<img src="Tutorial_KnockdownEfficiency_TPM_files/figure-markdown_github/app1-kEff-category-graph-1.png" width="85%" style="display: block; margin: auto;" />

## 4.2 Cell line influence

If we study with further details the table provided with the results, we
can observe that, in some instances, there are RBPs in which the
knockdown efficiency is significantly different between the two
different cell lines.

Let’s study the gene `MATR3` for example, where we observe that the
*T**P**M*<sub>*a**v**g*, *c**o**n**t**r**o**l*</sub> of <span
style="background-color: #c2ffa1">K562</span> is around 332 TPM while
for <span style="background-color: #ffaca1">HepG2</span> is only 152
TPM.

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cell line
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cluster
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
TPM
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Gene quantification ID
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:right;">
51.33
</td>
<td style="text-align:left;">
ENCFF363DHS
</td>
</tr>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:right;">
49.57
</td>
<td style="text-align:left;">
ENCFF221ZUL
</td>
</tr>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:right;">
68.67
</td>
<td style="text-align:left;">
ENCFF140UQA
</td>
</tr>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:left;">
case
</td>
<td style="text-align:right;">
71.04
</td>
<td style="text-align:left;">
ENCFF326HCV
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
MATR3
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
HepG2
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
control
</td>
<td style="text-align:right;color: black !important;background-color: #ffaca1 !important;">
135.96
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
ENCFF548ZDN
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
MATR3
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
HepG2
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
control
</td>
<td style="text-align:right;color: black !important;background-color: #ffaca1 !important;">
168.13
</td>
<td style="text-align:left;color: black !important;background-color: #ffaca1 !important;">
ENCFF785WTG
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
MATR3
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
K562
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
control
</td>
<td style="text-align:right;color: black !important;background-color: #c2ffa1 !important;">
314.99
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
ENCFF587DSG
</td>
</tr>
<tr>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
MATR3
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
K562
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
control
</td>
<td style="text-align:right;color: black !important;background-color: #c2ffa1 !important;">
350.41
</td>
<td style="text-align:left;color: black !important;background-color: #c2ffa1 !important;">
ENCFF774RTC
</td>
</tr>
</tbody>
</table>

If we study the knockdown efficiency for each cell line individually, we
observe a *k**E**f**f* of 66.8% for HepG2 and 79% for K562:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Cell line
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
TPM Cases
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
TPM Controls
</th>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
kEff
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
HepG2
</td>
<td style="text-align:right;">
50.450
</td>
<td style="text-align:right;">
152.045
</td>
<td style="text-align:left;">
66.8%
</td>
</tr>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:left;">
K562
</td>
<td style="text-align:right;">
69.855
</td>
<td style="text-align:right;">
332.700
</td>
<td style="text-align:left;">
79%
</td>
</tr>
</tbody>
</table>

In some cases, the difference in *k**E**f**f* between the two cell lines
is significantly greater. In the following table and graph, we show the
*k**E**f**f* for each cell line and the difference between the two:

<table class=" lightable-classic lightable-hover" style="font-size: 14px; font-family: Cambria; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;font-weight: bold;font-size: 16px;">
Target gene
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
HepG2 kEff %
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
K562 kEff %
</th>
<th style="text-align:right;font-weight: bold;font-size: 16px;">
Difference abs(HepG2 - K562)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
SRSF1
</td>
<td style="text-align:right;">
12.3
</td>
<td style="text-align:right;">
82.1
</td>
<td style="text-align:right;">
69.7
</td>
</tr>
<tr>
<td style="text-align:left;">
U2AF2
</td>
<td style="text-align:right;">
1.7
</td>
<td style="text-align:right;">
68.7
</td>
<td style="text-align:right;">
67.1
</td>
</tr>
<tr>
<td style="text-align:left;">
SNRNP200
</td>
<td style="text-align:right;">
61.6
</td>
<td style="text-align:right;">
30.0
</td>
<td style="text-align:right;">
31.5
</td>
</tr>
<tr>
<td style="text-align:left;">
SF3A3
</td>
<td style="text-align:right;">
56.4
</td>
<td style="text-align:right;">
84.3
</td>
<td style="text-align:right;">
27.9
</td>
</tr>
<tr>
<td style="text-align:left;">
GPKOW
</td>
<td style="text-align:right;">
76.0
</td>
<td style="text-align:right;">
51.1
</td>
<td style="text-align:right;">
24.9
</td>
</tr>
<tr>
<td style="text-align:left;">
PRPF6
</td>
<td style="text-align:right;">
42.1
</td>
<td style="text-align:right;">
66.3
</td>
<td style="text-align:right;">
24.2
</td>
</tr>
<tr>
<td style="text-align:left;">
U2AF1
</td>
<td style="text-align:right;">
66.4
</td>
<td style="text-align:right;">
42.9
</td>
<td style="text-align:right;">
23.5
</td>
</tr>
<tr>
<td style="text-align:left;">
NCBP2
</td>
<td style="text-align:right;">
77.7
</td>
<td style="text-align:right;">
55.8
</td>
<td style="text-align:right;">
21.9
</td>
</tr>
<tr>
<td style="text-align:left;">
QKI
</td>
<td style="text-align:right;">
70.4
</td>
<td style="text-align:right;">
50.2
</td>
<td style="text-align:right;">
20.2
</td>
</tr>
<tr>
<td style="text-align:left;">
CELF1
</td>
<td style="text-align:right;">
42.0
</td>
<td style="text-align:right;">
22.9
</td>
<td style="text-align:right;">
19.1
</td>
</tr>
<tr>
<td style="text-align:left;">
ZRANB2
</td>
<td style="text-align:right;">
40.9
</td>
<td style="text-align:right;">
60.1
</td>
<td style="text-align:right;">
19.1
</td>
</tr>
<tr>
<td style="text-align:left;">
FUBP1
</td>
<td style="text-align:right;">
61.6
</td>
<td style="text-align:right;">
77.8
</td>
<td style="text-align:right;">
16.2
</td>
</tr>
<tr>
<td style="text-align:left;">
GEMIN5
</td>
<td style="text-align:right;">
66.6
</td>
<td style="text-align:right;">
82.8
</td>
<td style="text-align:right;">
16.1
</td>
</tr>
<tr>
<td style="text-align:left;">
TIAL1
</td>
<td style="text-align:right;">
69.7
</td>
<td style="text-align:right;">
54.8
</td>
<td style="text-align:right;">
14.9
</td>
</tr>
<tr>
<td style="text-align:left;">
PTBP1
</td>
<td style="text-align:right;">
55.7
</td>
<td style="text-align:right;">
70.4
</td>
<td style="text-align:right;">
14.7
</td>
</tr>
<tr>
<td style="text-align:left;">
EFTUD2
</td>
<td style="text-align:right;">
64.8
</td>
<td style="text-align:right;">
78.8
</td>
<td style="text-align:right;">
13.9
</td>
</tr>
<tr>
<td style="text-align:left;">
RPS10
</td>
<td style="text-align:right;">
82.9
</td>
<td style="text-align:right;">
70.2
</td>
<td style="text-align:right;">
12.7
</td>
</tr>
<tr>
<td style="text-align:left;">
MATR3
</td>
<td style="text-align:right;">
66.8
</td>
<td style="text-align:right;">
79.0
</td>
<td style="text-align:right;">
12.2
</td>
</tr>
<tr>
<td style="text-align:left;">
RBM22
</td>
<td style="text-align:right;">
45.7
</td>
<td style="text-align:right;">
57.5
</td>
<td style="text-align:right;">
11.8
</td>
</tr>
<tr>
<td style="text-align:left;">
SUGP2
</td>
<td style="text-align:right;">
42.7
</td>
<td style="text-align:right;">
54.2
</td>
<td style="text-align:right;">
11.5
</td>
</tr>
</tbody>
</table>

<img src="Tutorial_KnockdownEfficiency_TPM_files/figure-markdown_github/app1-kEff-cell_line-graph-1.png" width="85%" style="display: block; margin: auto;" />

We represent in the X and Y-axis the same information as before, but
divided by the cell line (black for HepG2 and blue for K562). The green
line represents the distance between the two measurements.

For our analysis, we grouped both cell lines independently without any
considerations. The previous graph shows that the knockdown effects are
not always consistent in both cell lines. These results could have a
relevant effect in the mis-splicing ratio of introns found in the more
affected genes (e.g. U2AF2, SRSF1).

<!-- For example, for target gene U2AF2, the variation in TPM from <span style="background-color: #c2ffa1">control</span> to <span style="background-color: #ffaca1">cases</span> in the HepG2 cell line produces a knockdown efficiency of 1.7%, while nearly a 70% efficiency for the K562 cell line. -->
<!-- ```{r app1-kEff-U2AF2, echo = F} -->
<!-- metadata_TPM %>% -->
<!--   dplyr::filter(target_gene == "U2AF2") %>% -->
<!--   dplyr::select(target_gene, cell_line, experiment_type, TPM, gene_quantification_id) %>% -->
<!--   dplyr::arrange(cell_line, experiment_type) %>% -->
<!--   `colnames<-`(c("Target gene", "Cell line", "Cluster", "TPM", "Gene quantification ID")) %>% -->
<!--   kableExtra::kbl(booktabs = T, linesep = "") %>% -->
<!--   kableExtra::kable_classic(full_width = F, "hover", "striped", html_font = "Cambria", font_size = 14) %>% -->
<!--   kableExtra::row_spec(1:2, color = "black", background = "#ffaca1") %>% -->
<!--   kableExtra::row_spec(3:4, color = "black", background = "#c2ffa1") -->
<!-- ``` -->
