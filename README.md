# RNAseq_Kallisto_workflow
Kallisto-based gene level RNA-seq analysis with Nexflow.
This repository contains a workflow description and Nextflow pipeline designed to obtain raw transcript and gene counts starting from Illumina sequencing reads in FASTQ format. 
The workflow includes following steps:
1. Downloading transcript sequences and creating Kallisto index;
2. Infer strandness of the experiment using RSEQC;
3. Apply Kallisto-based pipeline to process the reads and obtain raw transcript counts;
4. Summarize transcript counts at the gene level to obtain raw counts matrix that can be used with DESeq2, edgeR or similar.

## Environment management
The analysis was done in the separate conda environment designed for bulk RNA-seq. Nextflow workflow management software was intalled system-wide. All of the software within conda environment was installed using mamba package manager. A   

## Section 1. 
Download Rat transcriptome from Ensemble Rnor6.0 and unzip
```
wget ftp://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
gunzip Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
```
Create kallisto index with transcripts fasta file
```
kallisto index -i RNOR6_transcripts.idx Rattus_norvegicus.Rnor_6.0.cdna.all.fa
```



