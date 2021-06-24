#!/usr/bin/env nextflow

nextflow.enable.dsl=2

project_dir = projectDir

def helpMessage() {
    log.info"""
    Description:
        An easy to use pipeline to separate endosymbiont genomes from their host's
    Pipeline summary:
        1. Trimming using TrimGalore!
        2. Read quality control using FastQC
        3. De Novo assembly using megahit
        4. Filtering endosymbiont genome using blastn
        5. Filtering host mitogenome using blastn
        6. Read mapping for coverage estimate using bowtie2
        7. Coverage estimate
        8. Assembly quality assessment using QUAST
    Usage:
        nextflow run main.nf --reads '*_R{1,2}\\.fastq.gz' --endosymbiont_reference '*_endosymRef\\.fna' --host_reference '*_hostRef\\.fna'
        
    Mandatory arguments:
        --reads             path to one or more sets of paired-ended reads (valid
                            file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq')
        --endosymbiont_reference
                            path to one or more reference genomes for the endosymbiont
                            assembly (valid file type extensions: '.fa', '.fna', '.fasta', '.faa')
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


process raw_qc {

    publishDir "${params.output}/$name/quality_control", mode: 'copy'

    tag "$name"

    input:
    tuple val(name), file(reads)

    output:
    file '*.html'

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
    tuple val(name), file('*.fq*'), emit: trimmed_reads

    when:
    ! params.skip_trimming

    script:
    """
    trim_galore --paired --cores ${params.threads/2} ${reads[0]} ${reads[1]}
    """

}

process trimmed_qc {

    publishDir "${params.output}/$name/quality_control", mode: 'copy'

    tag "$name"

    input:
    tuple val(name), file(reads)

    output:
    file '*.html'

    when:
    ! skip_qc && ! params.skip_trimming

    script:
    """
    fastqc --threads ${params.threads/2} --quiet ${reads[0]} ${reads[1]}
    """

}

process first_assembly {

    tag "$name"

    input:
    tuple val(name), file(reads)

    output:
    tuple val(name), file('megahit_out/final.contigs.fa')

    script:
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -t ${params.threads} -m ${params.memory}
    """
}

process endosymbiont_mapping {

    publishDir "${params.output}/$name/endosymbiont_assembly", mode: 'copy'
    
    tag "$name"

    input:
    tuple val(name), file(contigs)
    file endosymbiont_reference

    output:
    tuple val(name), file('endosymbiont_genome.fa'), emit: endosym_mapped

    script:
    """
    makeblastdb -in $endosymbiont_reference -title endosymbiont -parse_seqids -dbtype nucl -hash_index -out db
    blastn -query $contigs -db db -outfmt "10 qseqid" > seqid.txt
    grep -F -f seqid.txt $contigs -A 1 > blasted_contigs.fa
    grep -v "-" blasted_contigs.fa > endosymbiont_genome.fa
    """
}

process host_read_filtering {

    publishDir "${params.output}/$name/host_assembly", mode: 'copy'
    
    tag "$name"

    input:
    tuple val(name), file(host_assembled)

    output:
    tuple val(name), file('mitogenome.fa'), emit: host_filtered

    when:
    ! params.endosymbiont_only 


    script:
    """
    makeblastdb -in $project_dir/seqs/cox1.fa -title cox1 -parse_seqids -dbtype nucl -hash_index -out db
    for i in {11..25..1}
      do
        blastn -query $host_assembled -db db -outfmt "10 qseqid" -word_size \$i > seqid.txt
        grep -F -f seqid.txt $host_assembled -A 1 > mitogenome.fa
        if [[ \$(wc -l mitogenome.fa) = "2 mitogenome.fa" ]];
        then
          break
        fi
      done
    """
}

process endosymbiont_assembly_quality {

    publishDir "${params.output}/$name/endosymbiont_assembly", mode: 'copy'

    tag "$name"

    input:
    tuple val(name), file(endosym)
    file endosymbiont_reference

    output:
    file '*'

    when:
    ! params.skip_assembly_quality && ! params.skip_endosymbiont_assembly

    script:
    """
    quast.py $endosym -r $endosymbiont_reference
    """

}

process host_assembly_quality {

    publishDir "${params.output}/$name/host_assembly", mode: 'copy'


    tag "$name"

    input:
    tuple val(name), file(host)

    output:
    file '*'

    when:
    ! params.skip_assembly_quality && ! params.endosymbiont_only

    script:
    """
    quast.py $host
    """

}

process mapping_for_coverage_estimate{

    tag "$name"

    input:
    tuple val(name), file(reads)
    tuple val(name), file(assembled_endosymbiont)

    output:
    path 'log.txt', emit: alignment_stats

    script:
    """
    bowtie2-build ${assembled_endosymbiont} endosymbiont
    bowtie2 -x endosymbiont -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive
    cat .command.log > log.txt
    """

}

process coverage_estimate {

    publishDir "${params.output}/$name"

    tag "$name"

    input:
    path stats
    tuple val(name), file(reads)
    file assembled_endosymbiont

    output:
    file 'coverage.txt'

    script:
    """
    grep -v ">" $assembled_endosymbiont | tr -d "\n" | wc -c > host_count.txt
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
    if (params.skip_trimming) {
        first_assembly(rawReads) }
    else {
        first_assembly(trimming.out) }
    endosymbiont_mapping(first_assembly.out, endosymbiont_reference)
    endosymbiont_assembly_quality(endosymbiont_mapping.out, endosymbiont_reference)
    if (params.skip_trimming) { 
        mapping_for_coverage_estimate(rawReads, endosymbiont_mapping.out)
        coverage_estimate(mapping_for_coverage_estimate.out, rawReads, endosymbiont_mapping.out)
    }
    else {
        mapping_for_coverage_estimate(trimming.out, endosymbiont_mapping.out)
        coverage_estimate(mapping_for_coverage_estimate.out, trimming.out, endosymbiont_mapping.out)
    }
    host_read_filtering(first_assembly.out)
    host_assembly_quality(host_read_filtering.out.host_filtered)

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
