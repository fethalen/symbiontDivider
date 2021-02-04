#!/usr/bin/env nextflow

params.reads = '/home/cmau/pipeline/sample_data/*_{1,2}.fastq'
params.genome = '/home/cmau/pipeline/sample_data/GCA_003671365.1_ASM367136v1_cds_from_genomic.fna'
params.host_genome= '/home/cmau/pipeline/sample_data/GCA_001015115.1_ASM101511v1_genomic.fna'  //better host genome desperately needed

Channel.fromFilePairs( params.reads )
       .ifEmpty{error "Cannot find files"}
       .set { data }


process trimming_and_qc {

    tag "$name"

    input:
    tuple val(name), file(reads) from data

    output:
    tuple val(name), file('*.fq') into trimmed
    file '*.html' into quality

    script:
    """

    trim_galore --paired --fastqc --cores 8 ${reads[0]} ${reads[1]}

    """
}


process mapping {

    tag "$name"

    input:
    tuple val(name), file(reads) from trimmed
    path endosym_genome from params.genome
    path host_genome from params.host_genome

    output:
    tuple val(name), file('*_endosym.sam') into endosym_mapped
    tuple val(name), file('*_host.sam') into host_mapped

    script:
    """
    bowtie2-build ${endosym_genome} wolbachia
    bowtie2-inspect -n wolbachia
    bowtie2 -x wolbachia -p 8 -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive


    bowtie2-build ${host_genome} host
    bowtie2-inspect -n host
    bowtie2 -x host -p 8 -1 ${reads[0]} -2 ${reads[1]} -S ${name}_host.sam --very-sensitive
    """

}


process read_filtering {
    
    tag "$name"

    input:
    tuple val(name), file(endosym_mapped) from endosym_mapped
    tuple val(name_host), file(host_mapped) from host_mapped

    output:
    tuple val(name), file('*_endosym.sam') into endosym_filtered
    tuple val(name), file('*_host.sam') into host_filtered


    script:
    """

    samtools view -h -F 4 $endosym_mapped > mapped_${endosym_mapped}
    samtools view -h -F 4 $host_mapped > mapped_${host_mapped}

    """
}


process endosymbiont_assembly {

    tag "$name"

    input:
    tuple val(name), file(filtered) from endosym_filtered

    output:
    file '*scaffolds.fa'

    script:
    """

    abyss-pe np=4 name=marta2_endosym k=96 in='$filtered' B=5G H=3 kc=3 v=-v

    """
}


process host_assembly {

    tag "$name"

    input:
    tuple val(name), file(filtered) from host_filtered

    output:
    file '*scaffolds.fa'

    script:
    """

    abyss-pe np=4 name=marta2_host k=96 in='$filtered' B=5G H=3 kc=3 v=-v

    """

}
/*
process coverage_size_est {

    input:
    file mapping_log from coverage_assess

    output:
    file '*.txt' into test

    script:
    """
    #!/usr/bin/env python3

    log = open('${mapping_log}', 'r')
    lines = log.readlines()
    f = open('test.txt', 'w')
    print('test')

    """

}
/*

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
