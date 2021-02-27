#!/usr/bin/env nextflow

nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    Empty
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

if ( params.symbiont_reference == null) {
    exit 1, "Missing mandatory argument '--symbiont_reference'\n" +
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

endosymbionReference = Channel
    .fromPath( params.symbiont_reference, type: 'file')

hostReference = Channel
    .fromPath( params.host_reference, type: 'file')


process trimming_and_qc {

    tag "$name"

    input:
    tuple val(name), file(reads)

    output:
    tuple val(name), file('*.fq')

    script:
    """

    trim_galore --paired --fastqc --cores 8 ${reads[0]} ${reads[1]}

    """
}

process symbiont_mapping {

    tag "$name"

    input:
    tuple val(name), file(reads)
    file symbiont_reference

    output:
    tuple val(name), file('*_endosym.sam'), emit: endosym_mapped

    script:
    """
    bowtie2-build ${symbiont_reference} symbiont
    bowtie2-inspect -n symbiont
    bowtie2 -x symbiont -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive
    """
}

process host_mapping {

    tag "$name"

    input:
    tuple val(name), file(reads)
    file host_reference

    output:
    tuple val(name), file('*_host.sam'), emit: host_mapped

    script:
    """
    bowtie2-build ${host_reference} host
    bowtie2-inspect -n host
    bowtie2 -x host -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_host.sam --very-sensitive
    """

}


process read_filtering {
    
    tag "$name"

    input:
    tuple val(name), file(endosym_mapped)
    tuple val(name_host), file(host_mapped)

    output:
    tuple val(name), file('*_endosym.sam'), emit: endosym_filtered
    tuple val(name), file('*_host.sam'), emit: host_filtered


    script:
    """

    samtools view -h -F 4 $endosym_mapped > mapped_${endosym_mapped}
    samtools view -h -F 4 $host_mapped > mapped_${host_mapped}

    """
}


process endosymbiont_assembly {

    tag "$name"

    input:
    tuple val(name), file(filtered)

    output:
    tuple val(name), file('*scaffolds.fa'), emit: endosym_assembled

    script:
    """

    abyss-pe np=4 name=marta2_endosym k=96 in='$filtered' B=5G H=3 kc=3 v=-v

    """
}


process host_assembly {

    tag "$name"

    input:
    tuple val(name), file(filtered)

    output:
    tuple val(name), file('*scaffolds.fa'), emit: host_assembled

    script:
    """

    abyss-pe np=4 name=marta2_host k=96 in='$filtered' B=5G H=3 kc=3 v=-v

    """

}


process endosymbiont_assembly_quality {

    tag "$name"

    input:
    tuple val(name), file(endosym)

    output:
    file '*'

    script:
    """
    busco -i $endosym -m genome -o $name --auto-lineage-prok

    """

}

process host_assembly_quality {

    tag "$name"

    input:
    tuple val(name), file(host)

    output:
    file '*'

    script:
    """
    busco -i $host -m genome -o $name --auto-lineage

    """

}
/*
process coverage_estimate {

    input:
    file genome

    script:
    """
    grep -v ">" genome | perl -pe "s/\n//g" | wc -c | bin/coverage_estimate.py

    """
}
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
    trimming_and_qc(rawReads)
    symbiont_mapping(trimming_and_qc.out, endosymbionReference)
    host_mapping(trimming_and_qc.out, hostReference)
    read_filtering(symbiont_mapping.out.endosym_mapped, host_mapping.out.host_mapped)
    endosymbiont_assembly(read_filtering.out.endosym_filtered)
    host_assembly(read_filtering.out.host_filtered)
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
