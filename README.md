# RNAseq_Kallisto_workflow
Kallisto-based gene level RNA-seq analysis with Nexflow.
This repository contains a workflow description and Nextflow pipeline designed to obtain raw transcript and gene counts starting from Illumina sequencing reads in FASTQ format. 
The workflow includes following steps:
1. Downloading transcript sequences and creating Kallisto index;
2. Infer strandness of the experiment using RSEQC;
3. Apply Kallisto-based pipeline to process the reads and obtain raw transcript counts;
4. Summarize transcript counts at the gene level to obtain raw counts matrix that can be used with DESeq2, edgeR or similar.

## Environment management
The analysis was done in the separate conda environment designed for bulk RNA-seq. Nextflow workflow management software was intalled system-wide. All of the software within conda environment was installed using mamba package manager. I had to create a separate environement for RSEQC due to conflicts that could not be resolved trough conda.  

## Section 1. Create kallisto index.
Download Rat transcriptome from Ensemble Rnor6.0 and unzip
```
wget ftp://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
gunzip Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
```
Create kallisto index with transcripts fasta file
```
kallisto index -i RNOR6_transcripts.idx Rattus_norvegicus.Rnor_6.0.cdna.all.fa
```
## Section 2. Infer strandness of the experiment.
Download Rat genomes (Ensembl, Rnor6.0) form Illumine iGenome website
```
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Rattus_norvegicus/Ensembl/Rnor_6.0/Rattus_norvegicus_Ensembl_Rnor_6.0.tar.gz
```
Create HISAT2 index and map a small sub-sample of reads to the genome
```
hisat2-build -p <N_THREADS> <GENOME_FASTA> <INDEX>
hisat2 -p <N_THREADS> -x <INDEX> -1 <READ1> -2 <READ2> | samtools sort -@ <N_THREADS> -o <SORTED_BAM>
```
Convert GTF annotation file to BED using bedops
```
gtf2bed < <ANNOTATION_GTF> > <ANNOTATION_BED>
infer_experiment.py -r <ANNOTATION_BED> -i <SORTED_BAM>
```

The output of infer_experiment.py was:
This is PairEnd Data
Fraction of reads failed to determine: 0.0214
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0176
Fraction of reads explained by "1+-,1-+,2++,2--": 0.9610

These output is interpreted as follows:
Fraction of reads explained by "1++,1--,2+-,2-+" --> Fraction of forward stranded reads
Fraction of reads explained by "1+-,1-+,2++,2--" --> Fraction of reverse stranded reads
The strandness of the library is inferred based on the majority fraction. For example, in our case, 0.961 of the reads are reverse stranded indicating a reverse stranded library. In unstanded library the fraction should be similar to 0.5 in both cases. 
For reverse stranded libraries we have to use --rf-stranded option for Kallisto.

## Section 3. Apply Kallisto-based Nextflow pipeline to counts of transcripts.
















