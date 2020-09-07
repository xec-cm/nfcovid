#!/usr/bin/env nextflow
/*
========================================================================================
                                    xec-cm/nfcovid
========================================================================================
    nfcovid Analysis Pipeline.
    #### Homepage / Documentation
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nfcovid --input samplesheet.csv --genome 'MN908947.3' -profile docker

    Mandatory arguments
        --input [file]                    Comma-separated file containing information about the samples in the experiment (see docs/usage.md)
        --fasta [file]                    Path to fasta reference for viral genome. 
        --bedpe [file]                    Path to BED file containing amplicon positions. Mandatory when calling variants with --protocol amplicon
        --adapter [file]                  Parh of adapter file for trimmomatic triming. Mandatory when --fastp is false 
        -profile [str]                    Configuration profile to use. Can use multiple (comma separated) Available: conda and docker
    
    Generic 
        --single_end [bool]               Specifies that the input is single-end reads (Default: false)

    Read trimming
        --save_trimmed [bool]             Save the trimmed FastQ files in the results directory (Default: false)
        --fastp [bool]                    Use Fasp for adapter trimming (Defauls: false)

    Read alignment 
        --save_align_intermeds [bool]     Save the bowtie2 file intermediates in the results directory (Default: false)
        --filter_unmapped [bool]          Remove unmapped reads from alignments (Default: false) 
        --min_mapped_reads [int]          Minimum number of mapped reads below which samples are removed from further processing (Default: 1000)

    QC
        --skip_fastqc [bool]              Skip FastQC (Default: false)

    Other options:
        --outdir [file]                   The output directory where the results will be saved
        -name [str]                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
        --awsqueue [str]                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion [str]                 The AWS Region for your AWS Batch job to run on
        --awscli [str]                    Path to the AWS CLI tool
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                SET UP CONFIGURATION VARIABLES                       -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

if (params.input) { ch_input = file(params.input, checkIfExists: true) } else { exit 1, "Input samplesheet file not specified!" }

if (params.fasta) {
    file(params.fasta, checkIfExists: true)
    lastPath = params.fasta.lastIndexOf(File.separator)
    lastExt = params.fasta.lastIndexOf(".")
    fasta_base = params.fasta.substring(lastPath+1)
    index_base = params.fasta.substring(lastPath+1,lastExt)
} else {
    exit 1, "Viral genome fasta file not specified!"
}
ch_fasta = file(params.fasta)

if (!params.fastp) {
    file(params.adapter, checkIfExists: true)
    ch_adapter = file(params.adapter)
}


if (params.bedpe) { ch_bedpe = file(params.bedpe, checkIfExists: true) } else { exit 1, 'Path to BED file containing amplicon positionsnot specified!'}

// Print warning if viral genome fasta has more than one sequence
def count = 0
ch_fasta.withReader { reader ->
    while (line = reader.readLine()) {
        if (line.contains('>')) {
            count++
            if (count > 1) {
                log.info "[nf-core/viralrecon] WARNING: This pipeline does not support multi-fasta genome files. Please amend the '--fasta' parameter."
                break
            }
        }
    }
}

////////////////////////////////////////////////////
/* --                   AWS                    -- */
////////////////////////////////////////////////////

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                       HEADER LOG INFO                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name'] = custom_runName ?: workflow.runName
summary['Samplesheet'] = params.input
summary['Viral Fasta File'] = params.fasta
summary['Amplicon bed'] = params.bedpe
if (params.adapter) summary['Adapter file'] = params.adapter
summary['Data Type'] = params.single_end ? 'Single-End' : 'Paired-End'
summary['Trimming tool'] = params.fastp ? 'FASTP' : 'TRIMMOMATIC'          
if (params.save_trimmed)             summary['Save Trimmed'] = 'Yes'
summary['Min Mapped Reads']      = params.min_mapped_reads
if (params.filter_unmapped)      summary['Remove Unmapped Reads']  = 'Yes'
if (params.save_align_intermeds) summary['Save Align Intermeds'] =  'Yes'
if (params.skip_fastqc)              summary['Skip FastQC'] = 'Yes'
summary['Output dir']                = params.outdir
summary['Publish dir mode']          = params.publish_dir_mode
summary['Launch dir']                = workflow.launchDir
summary['Working dir']               = workflow.workDir
summary['Script dir']                = workflow.projectDir
summary['User']                      = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']            = params.awsregion
    summary['AWS Queue']             = params.awsqueue
    summary['AWS CLI']               = params.awscli
}
summary['Config Profile']            = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(22)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                     PARSE DESIGN FILE                               -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * PREPROCESSING: Reformat design file and check validitiy
 */
process CHECK_DESIGN {
    tag "$design"
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    path design from ch_input

    output:
    path '*.csv' into ch_design_reads_csv

    script:  // This script is bundled with the pipeline, in nf-core/atacseq/bin/
    """
    check_design.py $design design_reads.csv
    """
}

/*
 * Create channels for input fastq files
 */
if (params.single_end) {
    ch_design_reads_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1) ] ] }
        .into { ch_raw_reads_fastqc;
                ch_raw_reads_trim }
} else {
    ch_design_reads_csv
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1), file(row.fastq_2) ] ] }
        .into { ch_raw_reads_fastqc;
                ch_raw_reads_trim }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        FASTQ QC                                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 1: FASTQC
 */
process FASTQC {
    tag "$name"
    label 'process_medium'
    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename -> 
                    filename.endsWith('.zip') ? "zips/$filename" : filename
                }

    when:
    !params.skip_fastqc

    input:
    tuple val(name), path(reads) from ch_raw_reads_fastqc

    output:
    path '*.{zip,html}' into ch_fastqc_reports_mqc

    script:
    // Added soft-links to original fastqs for consistent naming in MultiQC
    if (params.single_end) {
        """
        [ ! -f  ${name}.fastq.gz ] && ln -s $reads ${name}.fastq.gz
        fastqc -q -t $task.cpus ${name}.fastq.gz
        """
    } else {
        """
        [ ! -f  ${name}_1.fastq.gz ] && ln -s ${reads[0]} ${name}_1.fastq.gz
        [ ! -f  ${name}_2.fastq.gz ] && ln -s ${reads[1]} ${name}_2.fastq.gz
        fastqc -q -t $task.cpus ${name}_1.fastq.gz
        fastqc -q -t $task.cpus ${name}_2.fastq.gz
        """
    }
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ADAPTER TRIMMING                             -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 2: ADAPTER TRIMMING 
 */
if (params.fastp) { 
    process FASTP {
        tag "$sample"
        label 'process_medium'
        publishDir "${params.outdir}/fastp", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        if (filename.endsWith(".json")) filename
                        else if (filename.endsWith(".fastp.html")) filename
                        else if (filename.endsWith("_fastqc.html")) "fastqc/$filename"
                        else if (filename.endsWith(".zip")) "fastqc/zips/$filename"
                        else if (filename.endsWith(".log")) "log/$filename"
                        else params.save_trimmed ? filename : null
                    }

        input:
        tuple val(sample), path(reads) from ch_raw_reads_trim

        output:
        tuple val(sample), path("*.trim.fastq.gz") into ch_trimmed_bowtie2
        path "*.json" into ch_fastp_mqc
        path "*_fastqc.{zip,html}" into ch_fastp_fastqc_mqc
        path "*.{log,fastp.html}"
        path "*.fail.fastq.gz"

        script:
        // Added soft-links to original fastqs for consistent naming in MultiQC
        if (params.single_end) {
            """
            [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
            IN_READS='--in1 ${sample}.fastq.gz'
            OUT_READS='--out1 ${sample}.trim.fastq.gz --failed_out ${sample}.fail.fastq.gz'
            
            fastp \\
                \$IN_READS \\
                \$OUT_READS \\
                --cut_front \\
                --cut_tail \\
                --length_required 50 \\
                --trim_poly_x \\
                --json ${sample}.fastp.json \\
                --html ${sample}.fastp.html \\
                2> ${sample}.fastp.log
            fastqc --quiet *.trim.fastq.gz
            """
        } else {
            """ 
            [ ! -f  ${sample}_1.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.fastq.gz
            [ ! -f  ${sample}_2.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.fastq.gz
            IN_READS='--in1 ${sample}_1.fastq.gz --in2 ${sample}_2.fastq.gz'
            OUT_READS='--out1 ${sample}_1.trim.fastq.gz --out2 ${sample}_2.trim.fastq.gz --unpaired1 ${sample}_1.fail.fastq.gz --unpaired2 ${sample}_2.fail.fastq.gz'

            fastp \\
                \$IN_READS \\
                \$OUT_READS \\
                --detect_adapter_for_pe \\
                --cut_front \\
                --cut_tail \\
                --length_required 50 \\
                --trim_poly_x \\
                --json ${sample}.fastp.json \\
                --html ${sample}.fastp.html \\
                2> ${sample}.fastp.log
            fastqc --quiet *.trim.fastq.gz
            """
        }
    }
}

if (!params.fastp) { 
    process TRIMMOMATIC {
        tag "$sample"
        label "process_medium"
        publishDir "${params.outdir}/trimmomatic", mode: params.publish_dir_mode,
            saveAs: { filename ->
                        if (filename.endsWith(".log")) "log/$filename"
                        else params.save_trimmed ? filename : null
                    }

        input:
        tuple val(sample), path(reads) from ch_raw_reads_trim
        path adapterFile from ch_adapter

        output:
        tuple val(sample), path("*_trim.fastq.gz") into ch_trimmed_bowtie2
        if (!params.single_end) { path "*_fail.fastq.gz" }
        path "*.log"

        script:
        // Added soft-links to original fastqs for consistent naming in MultiQC
        if (params.single_end) {
            """ 
            [ ! -f  ${sample}.fastq.gz ] && ln -s $reads ${sample}.fastq.gz
            IN_READS='${sample}.fastq.gz'
            OUT_READS='${sample}_trim.fastq.gz'
            trimmomatic SE -threads $task.cpus -phred33 \${IN_READS} \${OUT_READS}  \\
                ILLUMINACLIP:${adapterFile}:2:30:10 LEADING:30 TRAILING:30 MINLEN:75 SLIDINGWINDOW:30:20 2> ${sample}.trimmomatic.log
            """
        } else {
            """
            [ ! -f  ${sample}_1.fastq.gz ] && ln -s ${reads[0]} ${sample}_1.fastq.gz
            [ ! -f  ${sample}_2.fastq.gz ] && ln -s ${reads[1]} ${sample}_2.fastq.gz
            IN_READS='${sample}_1.fastq.gz ${sample}_2.fastq.gz'
            OUT_READS='${sample}_1_trim.fastq.gz ${sample}_1_fail.fastq.gz ${sample}_2_trim.fastq.gz ${sample}_2_fail.fastq.gz'
            
            trimmomatic PE -threads $task.cpus -phred33 \${IN_READS} \${OUT_READS}  \\
                ILLUMINACLIP:${adapterFile}:2:30:10 LEADING:30 TRAILING:30 MINLEN:75 SLIDINGWINDOW:30:20 2> ${sample}.trimmomatic.log
            """
        }
    } 
}  

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                        ALIGN                                        -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 3.1: BUILD BOWTIE2 INDEX FROM FASTA
 */
process BOWTIE2_INDEX {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}/genome", mode: params.publish_dir_mode

    input:
    path fasta from ch_fasta

    output:
    path "Bowtie2Index" into ch_index

    script:
    """
    bowtie2-build \\
        --seed 1 \\
        $fasta \\
        $index_base
    mkdir Bowtie2Index && mv ${index_base}* Bowtie2Index
    """
}

/*
 * STEP 3.2: MAP READS WITH BOWTIE2
 */
process BOWTIE2 {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".log")) "log/$filename"
                    else params.save_align_intermeds ? filename : null
                }

    input:
    tuple val(sample), path(reads) from ch_trimmed_bowtie2
    path index from ch_index

    output:
    tuple val(sample), path("*.bam") into ch_bowtie2_bam
    path "*.log" into ch_bowtie2_mqc

    script:
    input_reads = params.single_end ? "-U $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    filter = params.filter_unmapped ? "-F4" : ""
    """
    bowtie2 \\
        --threads $task.cpus \\
        --local \\
        --very-sensitive-local \\
        -x ${index}/${index_base} \\
        $input_reads \\
        2> ${sample}.bowtie2.log \\
        | samtools view -@ $task.cpus -b -h -O BAM -o ${sample}.bam $filter -
    """
    }

/*
 * STEP 3.3: SORT BAM FILE
 */
process SORT_BAM {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/variants/bam", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".flagstat")) "samtools_stats/$filename"
                    else if (filename.endsWith(".idxstats")) "samtools_stats/$filename"
                    else if (filename.endsWith(".stats")) "samtools_stats/$filename"
                }

    input:
    tuple val(sample), path(bam) from ch_bowtie2_bam

    output:
    tuple val(sample), path("*.sorted.{bam,bam.bai}"), path("*.flagstat") into ch_sort_bam
    path "*.{flagstat,idxstats,stats}" into ch_sort_bam_flagstat_mqc

    script:
    """
    samtools sort -@ $task.cpus -o ${sample}.sorted.bam -T $sample $bam
    samtools index ${sample}.sorted.bam
    samtools flagstat ${sample}.sorted.bam > ${sample}.sorted.bam.flagstat
    samtools idxstats ${sample}.sorted.bam > ${sample}.sorted.bam.idxstats
    samtools stats ${sample}.sorted.bam > ${sample}.sorted.bam.stats
    """
}

// Get total number of mapped reads from flagstat file
def get_mapped_from_flagstat(flagstat) {
    def mapped = 0
    flagstat.eachLine { line ->
        if (line.contains(' mapped (')) {
            mapped = line.tokenize().first().toInteger()
        }
    }
    return mapped
}

// Function that checks the number of mapped reads from flagstat output
// and returns true if > params.min_mapped_reads and otherwise false
pass_mapped_reads = [:]
fail_mapped_reads = [:]
def check_mapped(sample,flagstat,min_mapped_reads=500) {
    mapped = get_mapped_from_flagstat(flagstat)
    c_reset = "\033[0m";
    c_green = "\033[0;32m";
    c_red = "\033[0;31m";
    if (mapped < min_mapped_reads.toInteger()) {
        log.info ">${c_red}>>>> $sample FAILED MAPPED READ THRESHOLD: ${mapped} < ${params.min_mapped_reads}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! <<<<${c_reset}<"
        fail_mapped_reads[sample] = mapped
        return false
    } else {
        //log.info "-${c_green}           Passed mapped read threshold > bowtie2 ($sample)   >> ${mapped} <<${c_reset}"
        pass_mapped_reads[sample] = mapped
        return true
    }
}

// Remove samples that failed mapped read threshold
ch_sort_bam
    .filter { sample, bam, flagstat -> check_mapped(sample,flagstat,params.min_mapped_reads) }
    .map { it[0..2] }
    .set { ch_sort_bam }

/*
 * STEP 3.4: TRIM SEQUENCES
 */
process TRIM {
    tag "$sample"
    label "process_medium"
    publishDir "${params.outdir}/trim", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".log")) "log/$filename"
                    else filename
                }
    
    input:
    tuple val(sample), path(bam), path(flagstat) from ch_sort_bam
    path bedpe from ch_bedpe

    output:
    tuple val(sample), path("*.bam") into ch_trim_bam
    path "*.log" into ch_trim_log_mqc
    
    script:
    input_bam = "${bam[0]}"
    """
    bamclipper.sh -b $input_bam -p $bedpe -n $task.cpus 2> ${sample}.bamclipper.log
    """
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/* --                                                                     -- */
/* --                   VARIANT CALLING AND CONSENSUS                     -- */
/* --                                                                     -- */
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*
 * STEP 4.1: GENERATE CONSENSUS FASTQ
 */
process CONSENSUS_FASTQ {
    tag "$sample"
    label "process_medium"
    publishDir "${params.outdir}/consensus", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".fastq")) "fastqs/$filename"
                    else filename
                }

    input: 
    tuple val(sample), path(bam) from ch_trim_bam
    path index from ch_index

    output:
    tuple val(sample), path("*.fastq") into ch_consensus_fastq

    script:
    """ 
    samtools mpileup -q 20 -Q 30 -uf ${index}/${index_base}.fasta $bam \\
        | bcftools call -c \\
        | vcfutils.pl vcf2fq -Q 20 -d 20 > ${sample}_consensus.fastq
    """  
}

/*
 * STEP 4.2: CONVERT FASTQ TO FASTA
 */
process FASTQ_TO_FASTA {
    tag "$sample"
    label "process_medium"
    publishDir "${params.outdir}/consensus", mode: params.publish_dir_mode,
        saveAs: { filename ->
                    if (filename.endsWith(".fasta")) "fastas/$filename"
                    else filename
                }

    input:
    tuple val(sample), path(cons_fastq) from ch_consensus_fastq

    output:
    tuple val(sample), path("*.fasta") into ch_consensus_fasta

    script:
    """
    sed -n '1~4s/^@/>/p;2~4p' $cons_fastq > ${sample}_consensus.fasta
    """
}