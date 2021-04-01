


# Accessing public data

<br><br>


## Introduction
Thanks to research becoming more open in the last decade, there is now a ***huge amount*** of freely available data online.  The push towards repeatable analysis and open access to data has had enormous impact on medical science, and will continue to do so in the future. Clearly, there is too much data to properly cite in a single document, but hopefully the main sections have been covered here. 

This document is a ***beginners guide*** for how to access this data. 



Data has been divided into sections as best as possible. <br>In each section, the following is covered:
* Who uses the data? (what type of analysis / field)
* How to access (download or use)
* Format of the data 

<br>

## Contents
* [Summary of common data repositories](#Common-data-Repositories)
    * [NCBI](#NCBI)
    * [EMBL-EBI](#EMBL-EBI)
    * [DDBJ](#DDBJ)

* [Reads](#Reads)
    * [DNA](#DNA)
    * [RNA](#RNA)
    * [Epigenetics](#Epigenetics)

* [Genomes](#Genomes)
    * [Assemblies](#Assemblies)
    * [Annotations](#Annotations)
    * [Genes](#Genes)
    * [Taxonomy](#Taxonomy)

* [Transcriptomics](#Transcriptomics)
    * [Expression data](#Expression-data)
    * [Transcript isoforms](#Transcript-isoforms)
    * [Networks & pathways](#Networks-&-pathways)

* [Variation](#Variation)
    * [Sequence variation (SNVs/SNPs, Indels)](#Sequence-variation-(SNVs/SNPs,-Indels))
    * [Structural variation](#Structural-variation)

* [Proteomics](#Proteomics)
    * [Protein structures](#Protein-structures)
    * [Protein interactions](#Protein-interactions)
    * [Protein expression](#Protein-expression)
    * [Protein domains and functions](#Protein-domains-and-functions)

<br><br>

## Common data Repositories 
Jump to
* [Overview](#Overview)
* [Data sharing between organisations](#Shared-Data)
* [Unique archives per organisation](#Unique-Archives-per-Organisation)
* [Non open-access archives](#Non-open-access-archives)

<br>

### Overview
There are 3 major organisations which house publically available data. While other independent resources exist, the bulk of open-access data is stored within these 3 groups:

* The National Center for Biotechnology Information (NCBI)
* The European Bioinformatics Institute (EMBL-EBI)
* The DNA Data Bank of Japan (DDBJ)

<br>

### Data Sharing

<br>NCBI, EMBL-EBI, and DDBJ are all part of the International Nucleotide Sequence Database Collaboration (INSDC). The INSDC is an initiative which *encourages the sharing of nucleotide sequence data between organisations*. As a result, data submitted to the archives below will be shared across relevant NCBI, EMBL-EBI, and DDBJ databases.  A large amount of information - from raw read data, to alignments, assemblies, functional annotations, and sample information - are shared. 

The following is a table of  ***archives which share data***

<br>

| Data | NCBI | EMBL-EBI | DDBJ |
| ------ | ------ | ------ | ------ |
| NGS reads | [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/?term=all%5Bfilter%5D) | [ENA](https://www.ebi.ac.uk/ena/browser/home) | [DDBJ Sequence Read Archive (DRA)](https://ddbj.nig.ac.jp/DRASearch/) |
| Capillary electrophoresis reads | [Trace Archive](https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi?view=search) | [ENA](https://www.ebi.ac.uk/ena/browser/home) | [DDBJ Trace Archive (DTA)](https://www.ddbj.nig.ac.jp/dta/index-e.html) |
| Annotated sequences (genes, assemblies etc) | [GenBank / RefSeq](https://www.ncbi.nlm.nih.gov/nuccore/?term=all%5Bfilter%5D) | [ENA](https://www.ebi.ac.uk/ena/browser/home) | [DDBJ](http://ddbj.nig.ac.jp/arsa/?lang=en) |
| Samples | [BioSample](https://www.ncbi.nlm.nih.gov/biosample/?term=all%5Bfilter%5D) | [BioSamples](https://www.ebi.ac.uk/biosamples/) | [BioSample](https://ddbj.nig.ac.jp/BSSearch/) |
| Studies | [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/?term=all%5Bfilter%5D) | [ENA](https://www.ebi.ac.uk/ena/browser/home) | [BioProject](https://ddbj.nig.ac.jp/BPSearch/) |

<br>

Similar databases are available at each of these organisations, and they mirror (hold a copy of) each others relevant data. This is helpful as it reduces data redundancy (NCBI, EMBL-EBI, and DDBJ accessions for a single piece of data are linked), improves download/upload speeds due to geographic closeness, and encourages international sharing of data. 

<br>

### Projects, Samples, and data

Understanding the hierarchy between archives is one of the most tricky aspects when navigating public data. Anyone who has worked with databases will know that the relationships between data are often hard to express in a standard way.   The 3 main organisations (NCBI, EMBL-EBI and DDBJ) arrange information into **BioProjects**, **BioSamples**, and Data, which is a good solution given the challenge. 

**BioProjects** are containers which store links.  They are like folders which hold links to all the data and metadata associated with some project. The links can be directly to data, or can be to descriptions of the data (metadata). 

Side note: EMBL-EBI call these ‘BioStudies’ instead of **BioProjects** for some unknown reason. We will use the term BioProject from here.  

**BioSamples** are actually just descriptions of biological material. They do not relate to the data which was generated, but they can link to data which was derived from the particular biological sample / material.  For example, if you isolated a colony of bacteria for whole genome sequencing (WGS), a BioSample entry would be created to describe the bacterial isolate. The BioSample would then have a link to the WGS data, specifying “the WGS dataset was generated from this biological material!”.

<p align="center">
    <img src="media/project_organisation.png" width="70%">
</p>

<br>

### NCBI
American 

Read more here:<br>
https://www.ncbi.nlm.nih.gov/home/documentation/


<p align="center">
    <img src="media/nucleotide_databases.png" width="70%">
</p>

<p align="center">
    <img src="media/variants.png" width="70%">
</p>

<br>

### EMBL-EBI
ELIXIR data archives:<br>
https://elixir-europe.org/services/list?field_scientific_domain_tid=All&field_elixir_badge_tid=All&field_type_of_service_tid=1134&field_elixir_node_target_id=All&combine=

https://www.ebi.ac.uk/ena/browser/about/content




<br>

### DDBJ


<br>

### Conclusion

Many academic journals now insist on data being publically accessible. They may even specifically ask for your data to be housed in the above, and will request accession numbers to confirm this before publishing. #TODO confim this


<br><br>

## Reads

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| [European nucleotide archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide data | All | 🟢 | 🟢 | 🔴 |
| [Sequence read archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) | High-throughput sequence data - capillary, short & long NGS reads | All | 🟡 | 🟢 | 🔴 |
| [DDBJ Sequence Read Archive (DRA)](https://www.ddbj.nig.ac.jp/dra/index-e.html) | High-throughput sequence data - capillary, short & long NGS reads | All | 🔴 | 🟢 | 🔴 |
| [NCBI Trace Archive](https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi) | Capillary sequencing only (replaced by SRA)  | All | 🔴 | 🟢 | 🔴 |
| [DDBJ Trace Archive (DTA)](https://www.ddbj.nig.ac.jp/dta/index-e.html) | Capillary sequencing only (replaced by DRA)  | All | 🔴 | 🟢 | 🔴 |


<br>


### Who uses the data

Reads are one of the most commonly used forms of bioinformatics data. They often are represented as DNA sequence (using ATGC) regardless of whether they are DNA or RNA sequence. 

Commonly used for: 
* Differential gene expression (DE) analysis
* Transcript isoform detection
* RNA base modification
* Transcript Regulation (miRNAs etc ) 

### How to access

### Data format


<br><br>

## Genomics
Jump to
* [Assemblies](#Assemblies)
* [Annotations](#Annotations)
* [Genes](#Genes)
* [Taxonomy](#Taxonomy)

🔴 - Low

🟡 - Med

🟢 - High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
|  <br><br>**All nucleotide sequences** |
| [European nucleotide archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide sequences | All | 🟢 | 🟢 | 🔴 |
| [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) | All nucleotide sequences | All | 🟡 | 🟢 | 🔴 |
| [DNA Data Bank of Japan (DDBJ)](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide sequences | All | 🔴 | 🟢 | 🔴 |
|  <br><br>**Genome Assemblies** | 
| [European nucleotide archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide sequences | All | 🟢 | 🟢 | 🔴 |
| [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly/) | Genome Assemblies | All | 🟡 | 🟢 | 🔴 |
| [DNA Data Bank of Japan (DDBJ)](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide sequences | All | 🔴 | 🟢 | 🔴 |
| <br><br>**Taxonomy** | 
| [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1423) | The standard taxonomy system | All | 🟡 | 🟢 | 🔴 |
| <br><br>**Functional Elements** | 
| [Encyclopedia of DNA Elements (ENCODE)](https://www.encodeproject.org/) | Annotations for human functional DNA elements | Human + select model organisms | 🟢 | 🟢 | 🟢 |
| [Encyclopedia of genes and gene variants (GENCODE)](https://www.gencodegenes.org/) | Annotations for human (and mouse) genes | Human, Mouse | 🟢 | 🟢 | 🟢 |
| [GeneCards](https://www.genecards.org/) | Aggregator for all gene-centric data. Each gene listed once. | Human | 🟢 | 🟢 | 🟡 |
| [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) | Genes and links to data/metadata | All | 🟡 | 🟢 | 🔴 |
| [European nucleotide archive (ENA)]() | All nucleotide information | All | 🟢 | 🟢 | 🔴 |
| [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) | All nucleotide information | All | 🟡 | 🟢 | 🔴 |
| [DNA Data Bank of Japan (DDBJ)](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide information | All | 🔴 | 🟢 | 🔴 |




<br>

### Assemblies
* NCBI Assembly
    * GenBank vs RefSeq 
* ENA
<br>

### Annotations
* ENCODE
* GENCODE 
<br>

### Genes
* ENA
* SILVA


<br>

### Taxonomy

<br>





<br><br>

## Transcriptomics
Jump to
* [Expression data](#Expression-data)
* [Transcript isoforms](#Transcript-isoforms)
* [Networks & pathways](#Networks-&-pathways)

<br>

🔴 - Low

🟡 - Med

🟢 - High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Bulk RNA** | 
| [Genotype-Tissue Expression (GTEx)](https://gtexportal.org/home/) | Tissue-specific gene expression and regulation | Human | 🟢 | 🟢 | 🟢 |
| [All of expression (AOE)](https://aoe.dbcls.jp/en) | Aggregates publicly available gene expression data | All | 🟢 | 🟢 | 🟡 |
| [Expression Atlas](https://www.ebi.ac.uk/gxa/home) | Abundance and localisation of RNA | All | 🟢 | 🟡 | 🟢 |
| [Gene Expression Omnibus (GEO) datasets](https://www.ncbi.nlm.nih.gov/gds/?term=all%5Bfilter%5D) | Functional Genomics Data (from NGS, Arrays etc) | All | 🟡 | 🟢 | 🟡 |
| [Gene Expression Omnibus (GEO) profiles](https://www.ncbi.nlm.nih.gov/geoprofiles/?term=all%5Bfilter%5D) | Expression profiles for a specific condition | All | 🟡 | 🟢 | 🟡 |
| [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) | Functional Genomics Data (NGS, Arrays etc) | All | 🟡 | 🟢 | 🟡 |
| <br><br>**Single Cell** |
| [The Human Cell Atlas (HCA)](https://www.humancellatlas.org/learn-more/human-cell-atlas/) | Single cell studies | Human | 🟢 | 🟢 | 🟡 |
| [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) | Single cell studies | All | 🟢 | 🟡 | 🟢 |
| [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell) | Single cell studies | All | 🟡 | 🟢 | 🟡 |
| [Tabula Muris](https://tabula-muris.ds.czbiohub.org/) | Single-cell transcriptome data | Mouse | 🟡 | 🔴 | 🟢 |
| [Human cell landscape (HCL)](https://db.cngb.org/HCL/) | Cell types and localisations | Human | 🔴 | 🔴 | 🟢 |
| <br><br>**Networks** |
| [Connectivity Map (CMap)](https://clue.io/data) | Transcriptional responses to chemical, genetic, and disease perturbation | Human | 🔴 | 🟢 | 🟡 |
| <br><br>**Transcript Isoforms** |
| [Genotype-Tissue Expression (GTEx)](https://gtexportal.org/home/) | Tissue-specific gene expression and regulation | Human | 🟢 | 🟢 | 🟢 |
| <br><br>**Noncoding RNA** |
| [RNAcentral ](https://rnacentral.org/) | All RNA information | :-: | :-: | :-: | :-: |

<br>

### Expression data
* GEO
* GTEx
* Tabula Muris
* 

<br>

### Transcript isoforms

<br>


### Networks & pathways
* REACTOME
* KEGG






<br><br>

## Variation
Jump to
* [Sequence variation (SNVs/SNPs, Indels)](#Sequence-variation-(SNVs/SNPs,-Indels))
* [Structural variation](#Structural-variation)

### Sequence variation (SNVs/SNPs, Indels)
* ClinVar
* OMIM
* GeneCards

<br>

### Structural variation






<br><br>

## Proteomics
jump to
* [Protein structures](#Protein-structures)
* [Protein interactions](#Protein-interactions)
* [Protein expression](#Protein-expression)
* [Protein domains and functions](#Protein-domains-and-functions)



### Protein structures
* PDB 
* PDBe 
<br>

### Protein interactions
* STRING
* BRENDA

<br>

### Protein expression
* PRIDE (mass spec data)
* Human Protein Atlas
<br>

### Protein domains and functions
* InterPro

<br><br>

