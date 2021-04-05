


# Accessing public data

<br><br>


## Introduction
Thanks to research becoming more open in the last decade, there is now a ***huge amount*** of freely available data online.  The push towards repeatable analysis and open access to data has had enormous impact on medical science, and will continue to do so in the future. Clearly, there is too much data to properly cite in a single document, but hopefully the main sections have been covered here. 

This document is a ***beginners guide*** for how to access this data. 

<br>

## Contents
* [Summary of common data repositories](#Common-data-Repositories)
    * [NCBI](#NCBI)
    * [EMBL-EBI](#EMBL-EBI)
    * [DDBJ](#DDBJ)

* [Sequence Reads](#Reads)
    * [DNA](#DNA)
    * [RNA](#RNA)
    * [Epigenetics](#Epigenetics)

* [Genomics](#Genomics)
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
    <img src="media/project_organisation.png" width="90%">
</p>

<br>

### NCBI
American 

You can also build complex queries using fields and logial operators. For example to search for horse dopamine receptor D2:

("Equus caballus"[Organism] OR horse[All Fields]) AND (dopamine 
receptor D2[Protein Name] OR (dopamine[All Fields] AND receptor[All Fields] AND D2[All Fields]))

Here is a list of all the fields which you can search:
https://www.ncbi.nlm.nih.gov/books/NBK49540/

Read more here:<br>
https://www.ncbi.nlm.nih.gov/home/documentation/


<br>

### EMBL-EBI
ELIXIR data archives:<br>
https://elixir-europe.org/services/list?field_scientific_domain_tid=All&field_elixir_badge_tid=All&field_type_of_service_tid=1134&field_elixir_node_target_id=All&combine=

https://www.ebi.ac.uk/ena/browser/about/content

<br>

### DDBJ

Is bad

<br>

### Conclusion

Many academic journals now insist on data being publically accessible. They may even specifically ask for your data to be housed in the above, and will request accession numbers to confirm this before publishing. #TODO confim this


<br><br>

## Reads

Jump to:
<br>[Next Gen Sequencing](#Next-Gen-Sequencing)
<br>[Capillary Electrophoresis](#Capillary-Electrophoresis)

<br><br>

🔴 Low 🟡 Med 🟢 High

<br>


| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| :-: | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Next Gen Sequencing** |
| [ENA](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide data | All | 🟢 | 🟢 | 🔴 |
| [SRA](https://www.ncbi.nlm.nih.gov/sra) | High-throughput sequence data - capillary, short & long NGS reads | All | 🟡 | 🟢 | 🔴 |
| [DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html) | High-throughput sequence data - capillary, short & long NGS reads | All | 🔴 | 🟢 | 🔴 |
| <br><br>**Capillary Electrophoresis** |
| [NCBI Trace Archive](https://trace.ncbi.nlm.nih.gov/Traces/trace.cgi) | Capillary sequencing only (replaced by SRA)  | All | 🔴 | 🟢 | 🔴 |
| [DDBJ Trace Archive (DTA)](https://www.ddbj.nig.ac.jp/dta/index-e.html) | Capillary sequencing only (replaced by DRA)  | All | 🔴 | 🟢 | 🔴 |


<br><br>

### Next Gen Sequencing

<br>

**Who uses these archives**

Sequence reads are one of the most commonly used forms of bioinformatics data. They often are represented as DNA sequence (using ATGC) in FASTQ format regardless of whether they are DNA or RNA sequence. 

Reads form the basis of modern-day bioinformatics, and are used in nearly every field. Common uses of reads include: 
* Differential gene expression (DE) analysis
* Variant calling & GWAS
* Genome assembly
* Microbiome characterisation
* Epigenetic analysis
* Transcript isoform detection

<br>

**How to access**

ENA

The [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/search) is perhaps the easiest platform to find and download read data sets. This said, it sometimes yields less results than NCBI SRA. 

Search using the bar in the top right hand corner of the page. Anything relevant should appear in the results. There are filters on the left side of the page which you can use to restrict the results to only sequence reads (click 'Run' under 'Reads' for raw reads). 

After the read sets have appeared, click on the accession (SRR10390707 in this case)

<br>

<p align="center">
    <img src="media/ENA1.PNG" width="90%">
</p>

<p align="center">
    <img src="media/ENA2.PNG" width="90%">
</p>

This will take you to a page of metadata and download links for the sequencing run. We can see the organism, accession, instrument, and other information. 

At the bottom there is a section called 'Read Files'. Highlighted in pink we see a .fastq.gz file which can be clicked on to download the read set. For sequencing runs with multiple read sets, they will all appear in this section. 

For extra information about the read set, expand the 'Show Column Selection' dropdown.

<br>

NCBI SRA

The process for NCBI SRA is similar to ENA, but is a little more involved. It often returns more results than ENA, and there are ***more filtering options*** available in the left-hand panel, and in the 'Results by taxon' section on the right.

Clicking one of the results takes you to a SRA experiment page, which contains metadata and a link to the actual sequencing run. At the bottom is the 'Runs' section, which contains the run accession we want (in this case SRR14102829).  

<br>

<p align="center">
    <img src="media/SRA1.PNG" width="90%">
</p>

<p align="center">
    <img src="media/SRA2.PNG" width="90%">
</p>

<p align="center">
    <img src="media/SRA3.PNG" width="90%">
</p>

<br>

After clicking on the run accession (SRR14102829), you are transferred to the SRA archive which contains the run information and download links. 

Go to the 'Data access' tab to find the data. 

Once here, the read set can be downloaded via the links in the 'Original format' section on the bottom. 

<br><br>

**Data format**

Read data is almost always in FASTQ format. This format contains 4 lines per sequence:
* Read identifier
* Sequence
* Spacer (the '+' character)
* Quality string

<br>
<p align="left">
    <img src="media/read1.PNG" width="90%">
</p>


The actual sequence of the read is line 2, and the base-accuracy quality score associated with each base is line 4. Base qualities are seen as characters, but actually represent numbers (convert from ASCII character to decimal then deduct 33 usually to get the Phred quality score)

<br><br>

### Capillary Electrophoresis 

<br>

**Who uses these archives**

Capillary Electrophoresis (Sanger) sequencing is no longer commonly used, but is still viable data. 

**How to access**

Capillary Electrophoresis data has been migrated to NCBI SRA, ENA, and DRA, so can be searched as above. It included or excluded by advanced searches and filters. 

**Data format**

FASTQ format as above. 


<br><br>

## Genomics
Jump to
<br>[Assemblies](#Assemblies)
<br>[Annotations](#Annotations)
<br>[Genes](#Genes)
<br>[Taxonomy](#Taxonomy)

<br><br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
|  <br><br>**All nucleotide sequences** |
| [ENA](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide sequences | All | 🟢 | 🟢 | 🔴 |
| [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) | All nucleotide sequences | All | 🟡 | 🟢 | 🔴 |
| [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide sequences | All | 🔴 | 🟢 | 🔴 |
|  <br><br>**Genome Assemblies** | 
| [ENA](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide sequences | All | 🟢 | 🟢 | 🔴 |
| [NCBI Assembly](https://www.ncbi.nlm.nih.gov/assembly/) | Genome Assemblies | All | 🟡 | 🟢 | 🔴 |
| [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide sequences | All | 🔴 | 🟢 | 🔴 |
| <br><br>**Taxonomy** | 
| [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1423) | The standard taxonomy system | All | 🟡 | 🟢 | 🔴 |
| <br><br>**Functional Elements** | 
| [ENCODE](https://www.encodeproject.org/) | Annotations for human functional DNA elements | Human + select model organisms | 🟢 | 🟢 | 🟢 |
| [GENCODE](https://www.gencodegenes.org/) | Annotations for human (and mouse) genes | Human, Mouse | 🟢 | 🟢 | 🟢 |
| [GeneCards](https://www.genecards.org/) | Aggregator for all gene-centric data. Each gene listed once. | Human | 🟢 | 🟢 | 🟡 |
| [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/) | Genes and links to data/metadata | All | 🟡 | 🟢 | 🔴 |
| [ENA](https://www.ebi.ac.uk/ena/browser/home) | All nucleotide information | All | 🟢 | 🟢 | 🔴 |
| [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) | All nucleotide information | All | 🟡 | 🟢 | 🔴 |
| [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) | All nucleotide information | All | 🔴 | 🟢 | 🔴 |

<br><br>

### All nucleotide sequences

<br>

**Who uses these archives**

NCBI Nucleotide, ENA and DDBJ are locations where you can search for nucleotide data of any kind. These act as central repositories for data search and storage. Many other databases link to these repositories. 

NCBI nucleotide is a search which pulls data from sources which actually store the data, such as GenBank, RefSeq, and the TPA (see below).  

ENA and DDBJ are more straightforward, and appear as a single entity to the user.  

The data in these repositories appear in more specific archives.  For example, [NCBI gene]() contains gene sequences deposited in GenBank and RefSeq, and contains additional information for each listed gene. 



<br>

<p align="center">
    <img src="media/nucleotide_databases.png" width="70%">
</p>

<br>

People doing a wide search may find these archives helpful:
* If your biological question is narrow, you may just want to know what sequence data exists (if any)
* Searching these large central archives will inform you on the entirety of the data you may find in any location


<br>

**How to access**

As these repositories are for *all nucleotide data*, they aren't always great at retreiving the data you want. The only archive I would suggest for beginners is ENA, as the filters are very straightforward and helpful. NCBI can be quite powerful, but requires some experience. 

[NCBI Nucleotide advanced search](https://www.ncbi.nlm.nih.gov/nuccore/advanced) and [ENA advanced search](https://www.ebi.ac.uk/ena/browser/advanced-search) allow more specific queries and are the best option for finding data. 

If advanced searches are hard to understand at first, simple filtering is a good alternative. Filtering will also teach you the fields which can be used to narrow your results. 

In the example below, I am using NCBI Nucleotide to search for the AIRE gene. 2,228 results are given of varying types, including gene sequences, predicted gene sequences, genome assemblies, and read sets.  

<p align="center">
    <img src="media/nucleotide1.PNG" width="80%">
</p>

On the right we can select a specific taxonomic group to narrow our search, and below we can view the search details which NCBI nucleotide used. You can enter an advanced search here for more specific results (see this list of valid fields)

Filters can be seen on the left side of the screen, and can be useful if you know how to use them. By setting the 'Sequence length' filter to 1700 - 2000 bp, we can narrow our search to only sequences which are orthologues of the AIRE gene. This narrows our results to 146 items. 

<br>

<p align="center">
    <img src="media/nucleotide2.PNG" width="80%">
</p>

<br>

These sequences can then be downloaded by clicking 'Send to' in the top part of the results page. A number of options are available, and in this case I have selected FASTA to download the actual nucleotide sequences in FASTA format. 

The same search can be performed with the following query on [NCBI Nucleotide advanced search](https://www.ncbi.nlm.nih.gov/nuccore/advanced):

<p align="center">
    <img src="media/nucleotide3.PNG" width="80%">
</p>

<br>

**Data format**

* Depends on the data - keep in mind, these are collections of any nucleotide sequence, so the actual format may differ. Nucleotide sequences are generally stored in FASTA format. 

* FASTA files 

<br><br>

### Genome Assemblies

<br>

**Who uses the data**

**How to access**

**Data format**

<br><br>

### Taxonomy

<br>

**Who uses the data**

**How to access**

**Data format**

<br><br>

### Functional Elements

<br>

**Who uses the data**

**How to access**

**Data format**


<br><br>

## Transcriptomics
Jump to
* [Expression data](#Expression-data)
* [Transcript isoforms](#Transcript-isoforms)
* [Networks & pathways](#Networks-&-pathways)

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Bulk RNA** | 
| [GTEx](https://gtexportal.org/home/) | Tissue-specific gene expression and regulation | Human | 🟢 | 🟢 | 🟢 |
| [AOE](https://aoe.dbcls.jp/en) | Aggregates publicly available gene expression data | All | 🟢 | 🟢 | 🟡 |
| [Expression Atlas](https://www.ebi.ac.uk/gxa/home) | Abundance and localisation of RNA | All | 🟢 | 🟡 | 🟢 |
| [GEO datasets](https://www.ncbi.nlm.nih.gov/gds/?term=all%5Bfilter%5D) | Functional Genomics Data (from NGS, Arrays etc) | All | 🟡 | 🟢 | 🟡 |
| [GEO profiles](https://www.ncbi.nlm.nih.gov/geoprofiles/?term=all%5Bfilter%5D) | Expression profiles for a specific condition | All | 🟡 | 🟢 | 🟡 |
| [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/) | Functional Genomics Data (NGS, Arrays etc) | All | 🟡 | 🟢 | 🟡 |
| <br><br>**Single Cell** |
| [The Human Cell Atlas](https://www.humancellatlas.org/learn-more/human-cell-atlas/) | Single cell studies | Human | 🟢 | 🟢 | 🟡 |
| [Single Cell Expression Atlas](https://www.ebi.ac.uk/gxa/sc/home) | Single cell studies | All | 🟢 | 🟡 | 🟢 |
| [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell) | Single cell studies | All | 🟡 | 🟢 | 🟡 |
| [Tabula Muris](https://tabula-muris.ds.czbiohub.org/) | Single-cell transcriptome data | Mouse | 🟡 | 🔴 | 🟢 |
| [Human cell landscape](https://db.cngb.org/HCL/) | Cell types and localisations | Human | 🔴 | 🔴 | 🟢 |
| <br><br>**Networks** |
| [Connectivity Map (CMap)](https://clue.io/data) | Transcriptional responses to chemical, genetic, and disease perturbation | Human | 🔴 | 🟢 | 🟡 |
| <br><br>**Transcript Isoforms** |
| [GTEx](https://gtexportal.org/home/) | Tissue-specific gene expression and regulation | Human | 🟢 | 🟢 | 🟢 |
| <br><br>**Noncoding RNA** |
| [RNAcentral ](https://rnacentral.org/) | All RNA information | All | 🟢 | 🟢 | 🟡 |

<br>


<br><br>

## Variation
Jump to
* [Sequence variation (SNVs/SNPs, Indels)](#Sequence-variation-(SNVs/SNPs,-Indels))
* [Structural variation](#Structural-variation)

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Sequence Variation (SNPs, Indels etc)** | 
| [EVA](https://www.ebi.ac.uk/eva/) | All variant data | All | 🟡 | 🟢 | 🔴 |
| [NCBI dbSNP](https://www.ncbi.nlm.nih.gov/snp/) | All sequence variant data | Human | 🔴 | 🟢 | 🔴 |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Variant-phenotype relationship (health) | Human | 🔴 | 🟢 | 🟡 |
| [OMIM](https://www.omim.org/) | Gene-phenotype relationship | Human | 🔴 | 🟡 | 🟢 |
| <br><br>**Structural Variation (SVs)** | 
| [EVA](https://www.ebi.ac.uk/eva/) | All variant data | All | 🟡 | 🟢 | 🔴 |
| [NCBI dbVar](https://www.ncbi.nlm.nih.gov/dbvar/) | All structural variant data | Human | 🔴 | 🟢 | 🔴 |
| [DGV](http://dgv.tcag.ca/dgv/app/home) | Structural variation in healthy control samples (archived) | Human | 🔴 | 🟡 | 🟡 |

<br>

### Sequence variation (SNP/SNV, Indel etc)

<p align="center">
    <img src="media/variants.png" width="90%">
</p>

### Structural variation


<br><br>

## Proteomics
jump to
* [Protein structures](#Protein-structures)
* [Protein interactions](#Protein-interactions)
* [Protein expression](#Protein-expression)
* [Protein domains and functions](#Protein-domains-and-functions)

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Protein Sequences** | 
| [UniProt](https://www.uniprot.org/help/about) | Protein sequences and annotations | All | 🟢 | 🟢 | 🟢 |
| [Enzyme portal](https://www.ebi.ac.uk/enzymeportal/) | Concise summary of enzymes | All | 🟢 | 🟡 | 🟢 |
| [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/) | Protein sequences and annotations | All | 🔴 | 🟢 | 🔴 |
| <br><br>**Protein Domains & Families** | 
| [InterPro](https://www.ebi.ac.uk/interpro/about/interpro/) | Protein domains & families | All | 🟢 | 🟢 | 🟡 |
| [Pfam](http://pfam.xfam.org/) | Protein families | All | 🔴 | 🟢 | 🟡 |
| <br><br>**Protein Structures** | 
| [PDB](https://www.rcsb.org/) | Protein structures & associated data | All | 🟢 | 🟢 | 🟢 |
| [PDBe](https://www.ebi.ac.uk/pdbe/node/1) | Protein structures & associated data | All | 🟡 | 🟢 | 🟢 |
| [PDBJ](https://pdbj.org/) | Protein structures & associated data | All | 🔴 | 🟢 | 🟢 |
| <br><br>**Protein Expression** | 
| [The Human Protein Atlas](https://www.proteinatlas.org/) | Antibody-based imaging, mass spectrometry, transcriptomics data  | Human | 🟢 | 🟢 | 🟢 |
| [PRIDE](https://www.ebi.ac.uk/pride/) | Mass spectrometry data | All | 🟡 | 🟢 | 🟢 |
| <br><br>**EM, XRay, NMR Data** | 
| [EMDB](https://wwwdev.ebi.ac.uk/emdb/) | 3D EM density maps | All | 🟡 | 🟢 | 🟡 |
| [EMDataResource](https://www.emdataresource.org/index.html) | 3D EM density maps, models & metadata | All | 🔴 | 🟢 | 🟡 |
| [EMPIRE](https://www.ebi.ac.uk/pdbe/emdb/empiar/) | Raw electron microscopy images | All | 🟡 | 🟡 | 🟡 |
| [BMRB](https://bmrb.io/) | NMR data | All | 🔴 | 🟢 | 🟡 |

<br>

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

## Pathways & Reactions

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| [Reactome](https://reactome.org/) | Biological pathways | All | 🟢 | 🟢 | 🟢 |
| [Rhea](https://www.rhea-db.org/) | Reactions of biological interest | All | 🟡 | 🟢 | 🟢 |
| [KEGG](https://www.genome.jp/kegg/kegg2.html) | Biological systems | All | 🔴 | 🟢 | 🟢 |

<br>

<br><br>

## Metagenomics / Microbiome

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| [SILVA](https://www.arb-silva.de/) | ribosomal RNA sequences | All | 🟡 | 🟢 | 🟢 |
| [Ribosomal database project (RDP)](https://rdp.cme.msu.edu/index.jsp) | ribosomal RNA sequences | Bacteria, Archaea, Fungi | 🟡 | 🟢 | 🟡 |
| [MGnify](https://www.ebi.ac.uk/metagenomics/) | Microbiome experiments & data | All | 🟢 | 🟢 | 🔴 |
| [BacDrive](https://bacdive.dsmz.de/) | Bacterial information (Geographical, biochemical) | Bacteria | 🟢 | 🟢 | 🟢 |

<br>

<br><br>

## Metabolomics

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| [ChEMBL](https://www.ebi.ac.uk/chembl/)  | Bioactive molecules | All | 🟢 | 🟢 | 🟢 |
| [MetaboLights](https://www.ebi.ac.uk/metabolights)  | Studies of Metabolites | All | 🟡 | 🟢 | 🟡 |

<br>

<br><br>

## Domain Specific

<br>

🔴 Low 🟡 Med 🟢 High

<br>

| Name | Data stored | Organisms | Ease<br>of<br>Access | Amount<br>of<br>data | Data<br>curation<br>/quality | 
| --- | :-: | :-: | :-: | :-: | :-: |
| <br><br>**Cancer** |   
| [COSMIC](https://cancer.sanger.ac.uk/cosmic) | Somatic mutations in human cancer | Human | 🟢 | 🟢 | 🟢 |  
| <br><br>**Biological Images** | 
| [BioImage archive](https://www.ebi.ac.uk/biostudies/BioImages/studies) | All biological image data  | All | 🟢 | 🟢 | 🟢 |
| [Image Data Resource (IDR)](https://idr.openmicroscopy.org/about/) | Image datasets from published studies | All | 🟢 | 🟡 | 🟢 |   
| [Cell Image Library ](http://www.cellimagelibrary.org/home) | Images, videos, and animations of cells | All | 🟢 | 🟢 | 🟡 |
| <br><br>**Neuroscience** |   
| [Allen Brain Map](https://portal.brain-map.org/) | Data and analysis related to the brain | Human, Mouse | 🟢 | 🟢 | 🟢 |  
| <br><br>**Immunology** |   
| [ImmGen](https://www.immgen.org/) | Microarray gene expression & regulation | Mouse | 🟢 | 🟢 | 🟢 |  
| <br><br>**Biodiversity** |   
| [GBIF](https://www.gbif.org/) | Biodiversity data | All | 🟡 | 🟢 | 🟢 | 
| <br><br>**Disease Biomarkers** |   
| [BIONDA](http://bionda.mpc.ruhr-uni-bochum.de/start.php) | Biomarker candidates published in PubMed articles | Human | 🔴 | 🟢 | 🟡 |   
| <br><br>**Fruit flies** |   
| [FlyBase](https://flybase.org/) | All data types | Fruit flies | 🔴 | 🟢 | 🔴 |  
| <br><br>**Epigenomics** |   
| [MethBase](http://smithlabresearch.org/software/methbase/) | Reference methylomes (bisulfide-seq) | Selected model organisms | 🔴 | 🟡 | 🟡 | 

<br>

<br><br>


## Graveyard



Data has been divided into sections as best as possible. <br>In each section, the following is covered:
* Who uses the data? (what type of analysis / field)
* How to access (download or use)
* Format of the data 