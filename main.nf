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

// Input validation reads
if ( params.reads == null) {
    exit 1, "Missing mandatory argument '--reads'\n" +
            "Launch this workflow with '--help' for more info"
}
// Input validation endosymbiont reference
if ( params.endosymbiont_reference == null) {
    exit 1, "Missing mandatory argument '--endosymbiont_reference'\n" +
            "Launch this workflow with '--help' for more info"
}

// Creation of read pair channel with file extension filter and check if empty
rawReads = Channel
    .fromFilePairs( params.reads, size: 2, type: 'file' )
    .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
    .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.reads}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }

// Creation of endosymbiont reference channel
endosymbiont_reference = Channel
    .fromPath( params.endosymbiont_reference, type: 'file')


process RAWQC {

    /* 
        Process Description:
        Quality control of raw reads with FastQC
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/quality_control", mode: 'copy'

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw read files and the files themself
    tuple val(name), file(reads)

    output:
    // FastQC outputs the results as .html files -> those are output
    file '*.html'

    when:
    // This process only executes when QC is not skipped
    ! skip_qc

    script:
    // FastQC command with half of input threads and in quiet mode (see FastQC documentation for details)
    """
    fastqc --threads ${params.threads/2} --quiet ${reads[0]} ${reads[1]}
    """

}

process TRIMMING {

    /* 
        Process Description:
        Trimming of read files using TrimGalore!
    */
    
    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw read files and the files themself
    tuple val(name), file(reads)

    output:
    // TrimGalore! outputs the trimmed reads in .fg files that are output in a tuple comined with the name of the files
    tuple val(name), file('*.fq*'), emit: trimmed_reads

    when:
    // This process is only executed when trimming is not skipped
    ! params.skip_trimming

    script:
    // TrimGalore! command with half of input threads and paired read options (see TrimGalore! documentation for details)
    """
    trim_galore --paired --cores ${params.threads/2} ${reads[0]} ${reads[1]}
    """

}

process TRIMMEDQC {

    /* 
        Process Description:
        Quality control of trimmed reads with FastQC
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/quality_control", mode: 'copy'

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the trimmed read files and the files themself
    tuple val(name), file(reads)

    output:
    // FastQC outputs the results as .html files -> those are output
    file '*.html'

    when:
    // This process is only executed if QC and trimming are not skipped
    ! skip_qc && ! params.skip_trimming

    script:
    // FastQC command with half of input threads and in quiet mode (see FastQC documentation for details)
    """
    fastqc --threads ${params.threads/2} --quiet ${reads[0]} ${reads[1]}
    """

}

process DENOVOASSEMBLY {

    /* 
        Process Description:
        De novo assembly of reads using megahit assembler
    */

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name), file(reads)

    output:
    // megahit outputs the assembled genome as a .fa file called "final.contigs.fa" -> output of the process
    tuple val(name), file('megahit_out/final.contigs.fa')

    script:
    // megahit command using input threads and memory (for details see megahit documentation)
    """
    megahit -1 ${reads[0]} -2 ${reads[1]} -t ${params.threads} -m ${params.memory}
    """
}

process ENDOSYMBIONTREADFILTERING {

    /* 
        Process Description:
        Extraction of contigs belonging to the endosymbiont using blastn and grep
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/endosymbiont_assembly", mode: 'copy'
    
    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), file(contigs)
    // A fasta file containing the endosymbiont reference genome
    file endosymbiont_reference

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the endosymbiont genome
    tuple val(name), file('endosymbiont_genome.fa'), emit: endosym_mapped

    script:
    /*
    Script description:
    1. A blast database from the reference genome is created
    2. A blastn search with the de novo assembled contigs is performed and found contig ids are saved in a .txt
    3. Based on the contig ids, contigs are grepped from the de novo assembled contigs
    4. Dashes incorporated by grep are removed
    */
    """
    makeblastdb -in $endosymbiont_reference -title endosymbiont -parse_seqids -dbtype nucl -hash_index -out db
    blastn -query $contigs -db db -outfmt "10 qseqid" > seqid.txt
    grep -F -f seqid.txt $contigs -A 1 > blasted_contigs.fa
    grep -v "-" blasted_contigs.fa > endosymbiont_genome.fa
    """
}

process HOSTMITOGENOMEFILTERING {

    /* 
        Process Description:
        Extraction of contig/s belonging to the host mitogenome using blastn and grep
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/host_assembly", mode: 'copy'
    
    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), file(host_assembled)

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the host mitogenome
    tuple val(name), file('mitogenome.fa'), emit: host_filtered

    when:
    // This process is only executed if the endosymbont only mode is not selected
    ! params.endosymbiont_only 


    script:
    /*
    Script description:
    1. Create empty files for intermediate storage
    2. Create a blast database from a sequence for the cox1 gene (mitogenome exclusive gene)
    3. Iterate from 11 to 25
        3.1. Concatenate the previous found reads to prev_seqid.txt
        3.2. Blastn search with word size determined by iteration number, using de novo assembled reads as query
        3.3. Determined seqids are made unique and cat nto unique_seqid.txt
        3.4. If the unique_seqid.txt is empty -> grep contigs based on previously found seqids -> remove dashes implemented by grep -> break iteration
        3.5. If the unique_seqid.txt has 1 entry -> grep corrisponding contig -> break iteration
    */
    """
    touch mitogenome.fa
    touch prev_seqid.txt
    touch unique_seqid.txt
    makeblastdb -in $project_dir/seqs/cox1.fa -title cox1 -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {11..25..1}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query $host_assembled -db db -outfmt "10 qseqid" -word_size \$i > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        echo "made seqids unique"
        if [[ \$(wc -l unique_seqid.txt) = "0 unique_seqid.txt" ]];
        then
          grep -F -f prev_seqid.txt $host_assembled -A 1 > mitogenome.fa
          echo "multiple possible mitogenomes found"
          grep -v "--" mitogenome.fa > mitogenome_dashes_removed.fa
          echo "removed dashes"
          rm mitogenome.fa
          echo "removed old mitogenome"
          cat mitogenome_dashes_removed.fa > mitogenome.fa
          break
        fi
        if [[ \$(wc -l unique_seqid.txt) = "1 unique_seqid.txt" ]];
        then
          grep -F -f unique_seqid.txt $host_assembled -A 1 > mitogenome.fa
          echo "mitogenome found"
          break
        fi
      done
    echo "process successful"
    """
}

process ENDOSYMBIONTGENOMEQUALITY {

    /* 
        Process Description:
        Quality assessment of found endosymbiont genome using quast
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/endosymbiont_assembly", mode: 'copy'

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the endosymbiont genome
    tuple val(name), file(endosym)
    // A fasta file containing the endosymbiont reference genome
    file endosymbiont_reference

    output:
    // All files created by quast are output
    file '*'

    when:
    // This process is only executed if the last quality assessment is not skipped
    ! params.skip_assembly_quality

    script:
    // quast command (see quast documentation for details)
    """
    quast.py $endosym -r $endosymbiont_reference
    """

}

process HOSTMITOGENOMEQUALITY {

    /* 
        Process Description:
        Quality assessment of found host mitogenome using quast
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name/host_assembly", mode: 'copy'

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the host mitogenome
    tuple val(name), file(host)

    output:
    // All files created by quast are output
    file '*'

    when:
    // This process is only executed if the last quality assessment is not skipped and the endosymbiont only mode is not activated
    ! params.skip_assembly_quality && ! params.endosymbiont_only

    script:
    // quast command (see quast documentation for details)
    """
    quast.py $host
    """

}

process READMAPPINGFORCOVERAGE{

    /* 
        Process Description:
        Mapping of raw/trimmed reads onto the found endosymbiont genome using bowtie2
    */

    // Name of read files
    tag "$name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name), file(reads)
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the endosymbiont genome
    tuple val(name), file(assembled_endosymbiont)

    output:
    // The log file created by bowtie2 is output
    path 'log.txt', emit: alignment_stats

    script:
    /*
    Script description:
    1. Bowtie2 database is built from endosymbiont genome
    2. bowtie mapping with half the input threads (se bowtie2 documentation for details)
    3. command log is concatenated into log.txt
    */
    """
    bowtie2-build ${assembled_endosymbiont} endosymbiont
    bowtie2 -x endosymbiont -p ${params.threads/2} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive
    cat .command.log > log.txt
    """

}

process COVERAGEESTIMATE {

    /* 
        Process Description:
        Estimation of sequencing coverage using alignment stats from bowtie2, length of the assembled endosymbiont and length of the input reads
    */

    // Copies output file in output folder
    publishDir "${params.output}/$name"

    // Name of read files
    tag "$name"

    input:
    // log.txt file from previous process
    path stats
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name), file(reads)
    // A .fa file containing the assembled endosymbiont genome
    file assembled_endosymbiont

    output:
    // The process creates a file containign the coverage which is output
    file 'coverage.txt'

    script:
    /*
    Script description:
    1. Determining the number of bases in the endosymbiont genome
    2. Concatenating the alignment stats into a file
    3. Determining the number of bases in the read files
    4. Coverage estimate using a selfmade Python script
    */
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

    RAWQC(rawReads)
    TRIMMING(rawReads)
    TRIMMEDQC(TRIMMING.out)
    if (params.skip_trimming) {
        DENOVOASSEMBLY(rawReads) }
    else {
        DENOVOASSEMBLY(TRIMMING.out) }
    ENDOSYMBIONTREADFILTERING(DENOVOASSEMBLY.out, endosymbiont_reference)
    ENDOSYMBIONTGENOMEQUALITY(ENDOSYMBIONTREADFILTERING.out, endosymbiont_reference)
    if (params.skip_trimming) { 
        READMAPPINGFORCOVERAGE(rawReads, ENDOSYMBIONTREADFILTERING.out)
        COVERAGEESTIMATE(READMAPPINGFORCOVERAGE.out, rawReads, ENDOSYMBIONTREADFILTERING.out)
    }
    else {
        READMAPPINGFORCOVERAGE(TRIMMING.out, ENDOSYMBIONTREADFILTERING.out)
        COVERAGEESTIMATE(READMAPPINGFORCOVERAGE.out, TRIMMING.out, ENDOSYMBIONTREADFILTERING.out)
    }
    HOSTMITOGENOMEFILTERING(DENOVOASSEMBLY.out)
    HOSTMITOGENOMEQUALITY(HOSTMITOGENOMEFILTERING.out.host_filtered)

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
