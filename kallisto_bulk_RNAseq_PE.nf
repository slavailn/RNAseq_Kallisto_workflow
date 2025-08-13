/*
 * Step‑wise integration: Raw FastQC -->  fastp trimming --> FastQC (trimmed) --> Kallisto
 */

nextflow.enable.dsl = 2

//-------------------------------------
// Parameters
//-------------------------------------
params.rawdir     = 'raw_fastq'          // input PE FASTQs
params.threads    = 24                   // CPUs per task
params.qc_dir     = 'qc_raw'             // FastQC (raw)
params.trim_dir   = 'trimmed_reads'      // fastp outputs
params.qc_trimdir   = 'qc_trimmed'         // FastQC (trimmed)
params.kallisto_index = '/mnt/bioinf_data/GENOMES/Rattus_norvegicus/Ensembl/Rnor_6.0/Transcripts/RNOR6_transcripts.idx'   // Kallisto index
params.kallisto_dir   = 'kallisto_quant'    // Kallisto outputs

//-------------------------------------
// Channel: paired‑end tuples (sample_id, [R1,R2])
//-------------------------------------
reads_pattern = "${params.rawdir}/*_R{1,2}.fastq.gz"
Channel.fromFilePairs(reads_pattern, sampleRegex:/\.([A-Za-z0-9]+)_R[12]/)
       .set { read_pairs }

//-------------------------------------
// Workflow
//-------------------------------------
workflow {
    fastqc_raw(read_pairs)
    trimmed_pairs = trim_reads(read_pairs)
    fastqc_trimmed(trimmed_pairs)
    kallisto_quant(trimmed_pairs)
}

//-------------------------------------
// Processes
//-------------------------------------


//-------------------------------------
// Run fastQC on raw reads
//-------------------------------------
process fastqc_raw {
    tag { sample_id }
    cpus params.threads
    publishDir params.qc_dir, mode:'copy'
    conda 'bioconda::fastqc=0.12.1'

    input:
        tuple val(sample_id), path(reads)

    output:
        path '*_fastqc.html'
        path '*_fastqc.zip'

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    """
}

//-------------------------------
// Trim reads with fastp
//-------------------------------
process trim_reads {
    tag { id }
    cpus params.threads
    publishDir params.trim_dir, mode:'copy'
    conda 'bioconda::fastp=0.23.4'

    input:
        tuple val(id), path(reads)

    output:
        tuple val(id), path("${id}_R{1,2}_trimmed.fastq.gz")

    script:
    """
    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} \
        -o ${id}_R1_trimmed.fastq.gz -O ${id}_R2_trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --dont_eval_duplication \
        --trim_poly_g \
        --cut_front \
        --cut_front_mean_quality 30 \
        -w 16
    """
}


//-------------------------------
// Run fastQC on trimmed reads
//------------------------------
process fastqc_trimmed {
    tag { id }
    cpus params.threads
    publishDir params.qc_trimdir, mode:'copy'
    conda 'bioconda::fastqc=0.12.1'

    input:
        tuple val(id), path(reads)

    output:
        path '*_fastqc.html'
        path '*_fastqc.zip'

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${reads.join(' ')}
    """
}

//-------------------------------------
// Kallisto quantification
//-------------------------------------
process kallisto_quant {
    tag { id }
    cpus params.threads
    publishDir "${params.kallisto_dir}/${id}", mode:'copy', overwrite:false
    conda 'bioconda::kallisto=0.46.2'

    input:
        tuple val(id), path(reads)

    output:
        path 'abundance.tsv'
        path 'run_info.json'

    script:
    """
    kallisto quant \
        -i ${params.kallisto_index} \
        -o . \
        -t ${task.cpus} \
        --rf-stranded \
        -b 100 \
        --bias \
        ${reads.join(' ')}
    """
}
