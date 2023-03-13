# ENCODE - Metadata and knockdown efficiency extraction

This repository contains the code necessary to download and extract the ENCODE metadata for different shRNA knockdown experiments. It is designed to be easily expanded to any other experiments' categories (tested with CRISPR knockdown). 

Provided are a standalone script ([Metadata_Extraction_Script.R](Metadata_Extraction_Script.R)) and a detailed tutorial ([Tutorial_ENCODE_metadata.html](https://guillermo1996.github.io/ENCODE_Metadata_Extraction/RMarkdown/Tutorial_ENCODE_metadata.html)).

Additionally, two different tutorials are also available to estimate/extract the knockdown efficiency:

* [Reported knockdown efficiency by ENCODE (Western Blotting and qRT-PCR)](https://guillermo1996.github.io/ENCODE_Metadata_Extraction/RMarkdown/Tutorial_KnockdownEfficiency_WB.html): automates the search and download of the experiment characterization document from ENCODE to extract the reported knockdown efficiency by two methods: Western Blotting and qRT-PCR. The former is the most relevant method, since it reports the concentration of the protein in the sample and not the RNA abundance.

* [Estimated knockdown efficiency from TPM](https://guillermo1996.github.io/ENCODE_Metadata_Extraction/RMarkdown/Tutorial_KnockdownEfficiency_TPM.html): automates the search and download of the `gene_quantifications.tsv` file provided by ENCODE, which includes the gene expression levels in TPM for the different samples. We compare the average gene expressions in case and control samples to estimate the knockdown efficiency. This method reports results related to the qRT-PCR results, but the previously mentioned Western Blotting method is still preferred.
