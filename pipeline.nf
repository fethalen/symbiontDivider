#!/usr/bin/env nextflow

data = Channel.fromPath('/home/cmau/pipeline/sample_data/*.fastq')

process trimming_and_qc {

    input:
    file sample from data

    output:
    file '*.fq' into trimmed
    file '*.html' into quality

    script:
    """

    trim_galore $sample --fastqc --cores 4

    """
}

trimmed.view{it.name}
quality.view{it.name}

/*

process endosymbiont_mapping {

    input:
    file sample from trimmed

    output:
    file '*' into mapped

    script:
    """
    bowtie2-build GCA_003671365.1_ASM367136v1_cds_from_genomic.fna wolbachia
    bowtie2-inspect -n wolbachia
    bowtie2 -x wolbachia -1 NG-23851_Marta2_lib377150_6660_2_1.fastq -2 NG-23851_Marta2_lib377150_6660_2_2.fastq -S sample.sam --very-sensitive > bowtie_output.txt
    """

}

process coverage_size_est {

    input:
    file sample from mapped

    output:
    file '*' into quality

    script:
    """
    """

}

process assembly {

    input:
    file sample from mapped

    output:
    file '*' into assembled

    script:
    """
    """

}

process assembly_quality {

    input:
    file sample from assembled

    output:
    file '*' into quality

    script:
    """
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