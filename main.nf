#!/usr/bin/env nextflow

nextflow.enable.dsl=2

project_dir = projectDir

def helpMessage() {
    log.info"""
    Description:
    An easy to use pipeline to separate endosymbiont genomes from their host's

    Pipeline summary:
    1. Trimming using Trimmomatic
    2. Quality Control using FastQC
    3. Mapping of endosymbiont reads on reference genome using bowtie2
    4. Filtering the mpped reads using samtools
    5. Assembly of endosymbiont genome using ABySS
    6. Assembly of host genome using ABySS
    7. Assembly quality assesment using BUSCO

    Usage:
        nextflow run main.nf --reads '*_R{1,2}\\.fastq.gz' --endosymbiont_reference '*_endosymRef\\.fna' --host_reference '*_hostRef\\.fna'
        
    Mandatory arguments:
        --reads             path to one or more sets of paired-ended reads (valid
                            file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq')
        --endosymbiont_reference
                            path to one or more reference genomes for the endosymbiont
                            assembly (valid file types: '.fna')
        --host_reference
                            path to one or more reference genomes for the host assembly
                            (valid file types: '.fna')

    Input/output options:
        --output            path to a directory which the results are written to
                            (default: $params.output)

    Resource allocation:
        --memory            memory limit for the assembly step in GB (default:
                            $params.memory)
        --threads           maximum number of threads to be used by the pipeline
                            (default: '$params.threads')

    Flow control:
        --endosymbiont_only
                            skip processing of reads not belonging to the endosymbiont (default: $params.endosymbiont_only)
        --skip_coverage     skip coverage estimate step (default: $params.skip_coverage)
        --skip_trimming     skip trimming step (default: $params.skip_trimming)
        --skip_qc           skip reads quality assessment (default: $params.skip_qc)
        --skip_endosymbiont_assembly
                            skip endosymbiont assembly step (default: $params.skip_endosymbiont_assembly)
        --skip_host_assembly
                            skip host assembly step (default: $params.skip_host_assembly)
                           
        --skip_assembly_quality
                            skip assembly quality assessment (default: $params.skip_assembly_quality)

    Miscellaneous:
        --help              display this help message and exit
        --version           display the pipeline's version number and exit
    """.stripIndent()
}

def versionNumber() {
    log.info"symbiontDivider ~ version $workflow.manifest.version"
}

// Display the version number on request
if ( params.version ) exit 0, versionNumber()

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

// Input validation
if ( params.reads == null) {
    exit 1, "Missing mandatory argument '--reads'\n" +
            "Launch this workflow with '--help' for more info"
}

if ( params.endosymbiont_reference == null) {
    exit 1, "Missing mandatory argument '--endosymbiont_reference'\n" +
            "Launch this workflow with '--help' for more info"
}

if ( params.host_reference == null) {
    exit 1, "Missing mandatory argument '--host_reference'\n" +
            "Launch this workflow with '--help' for more info"
}

rawReads = Channel
    .fromFilePairs( params.reads, size: 2, type: 'file' )
    .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
    .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.reads}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }

endosymbiont_reference = Channel
    .fromPath( params.endosymbiont_reference, type: 'file')

host_reference = Channel
    .fromPath( params.host_reference, type: 'file')


process raw_qc {

    tag "$name"

    input:
    tuple val(name), file(reads)

    when:
    ! skip_qc

    script:
    """

    fastqc --threads ${params.threads/2} --quiet ${reads[0]} ${reads[1]}

    """

}

process trimming {

    tag "$name"

    input:
    tuple val(name), file(reads)

    output:
    tuple val(name), file('*.fq'), emit: trimmed_reads

    script:
    """

    trim_galore --paired --cores ${params.threads/2} ${reads[0]} ${reads[1]}

    """

}

process trimmed_qc {

    tag "$name"

    input:
    tuple val(name), file(reads)

    when:
    ! skip_qc

    script:
    """

    fastqc --threads ${params.threads/2} --quiet ${reads[0]} ${reads[1]}

    """

}

process endosymbiont_mapping {

    tag "$name"

    input:
    tuple val(name), file(reads)
    file endosymbiont_reference

    output:
    tuple val(name), file('*_endosym.sam'), emit: endosym_mapped
    path 'log.txt', emit: alignment_stats

    script:
    """
    bowtie2-build ${endosymbiont_reference} endosymbiont
    bowtie2 -x endosymbiont -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive
    cat .command.log > log.txt
    """
}

process host_mapping {

    tag "$name"

    input:
    tuple val(name), file(reads)
    file host_reference

    output:
    tuple val(name), file('*_host.sam'), emit: host_mapped

    when:
    ! params.endosymbiont_only

    script:
    """
    bowtie2-build ${host_reference} host
    bowtie2 -x host -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_host.sam --very-sensitive
    """

}


process endosymbiont_read_filtering {
    
    tag "$name"

    input:
    tuple val(name), file(endosym_mapped)

    output:
    tuple val(name), file('*_endosym.sam'), emit: endosym_filtered


    script:
    """

    samtools view -@ ${params.threads/2} -h -F 4 $endosym_mapped > mapped_${endosym_mapped}

    """
}

process host_read_filtering {
    
    tag "$name"

    input:
    tuple val(name), file(host_mapped)

    output:
    tuple val(name), file('*_host.sam'), emit: host_filtered

    when:
    ! params.endosymbiont_only


    script:
    """

    samtools view -@ ${params.threads/2} -h -F 4 $host_mapped > mapped_${host_mapped}

    """
}

process endosymbiont_assembly {
    publishDir "${params.output}/endosymbiont"

    tag "$name"

    input:
    tuple val(name), file(filtered)

    output:
    tuple val(name), file('*scaffolds.fa'), emit: endosym_assembled

    when:
    ! params.skip_endosymbiont_assembly

    script:
    """

    abyss-pe np=${params.threads/2} name=marta2_endosym k=96 in='$filtered' B=${params.memory/2}G H=3 kc=3 v=-v

    """
}


process host_assembly {
    publishDir "${params.output}/host"

    tag "$name"

    input:
    tuple val(name), file(filtered)

    output:
    tuple val(name), file('*scaffolds.fa'), emit: host_assembled

    when:
    ! params.skip_host_assembly && ! params.endosymbiont_only

    script:
    """
    abyss-pe np=${params.threads/2} name=${name}_host k=96 in='$filtered' B=${params.memory/2}G H=3 kc=3 v=-v
    """

}


process endosymbiont_assembly_quality {
    publishDir "${params.output}/endosymbiont"

    tag "$name"

    input:
    tuple val(name), file(endosym)

    output:
    file '*'

    when:
    ! params.skip_assembly_quality && ! params.skip_endosymbiont_assembly

    script:
    """
    busco -c ${params.threads/2} -i $endosym -m genome -o $name --auto-lineage-prok

    """

}

process host_assembly_quality {
    publishDir "${params.output}/host"


    tag "$name"

    input:
    tuple val(name), file(host)

    output:
    file '*'

    when:
    ! params.skip_assembly_quality && ! params.skip_host_assembly && ! params.endosymbiont_only

    script:
    """
    busco -c ${params.threads/2} -i $host -m genome -o $name --auto-lineage-euk

    """

}

process coverage_estimate {

    input:
    path stats
    tuple val(name), file(reads)
    file endosymbiont_reference

    output:
    file 'coverage.txt'

    script:
    """
    grep -v ">" $endosymbiont_reference | tr -d "\n" | wc -c > host_count.txt
    cat $stats > alignment_rate.txt
    cat ${reads[0]} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > base_count.txt
    python3 $project_dir/bin/coverage_estimate.py


    """
}
/*
process compression {

    input:
    file sample from assembled

    output:
    file '*.gz' into compressed

    script:
    """
    """
}

process visualise_quality {

    input:
    file sample from quality

    script:
    """
    """
}






*/
workflow {

    raw_qc(rawReads)

    trimming(rawReads)

    trimmed_qc(trimming.out)

    endosymbiont_mapping(trimming.out, endosymbiont_reference)

    coverage_estimate(endosymbiont_mapping.out.alignment_stats, trimming.out, endosymbiont_reference)

    host_mapping(trimming.out, host_reference)

    endosymbiont_read_filtering(endosymbiont_mapping.out.endosym_mapped)

    host_read_filtering(host_mapping.out.host_mapped)

    endosymbiont_assembly(endosymbiont_read_filtering.out.endosym_filtered)
    
    host_assembly(host_read_filtering.out.host_filtered)
    
    endosymbiont_assembly_quality(endosymbiont_assembly.out.endosym_assembled)
    
    host_assembly_quality(host_assembly.out.host_assembled)
}

workflow.onComplete {
    // Display complete message
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
    // Display error message
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}

