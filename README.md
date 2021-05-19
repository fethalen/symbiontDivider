# symbiontDivider
An easy to use pipeline to separate endosymbiont genomes from their host's

## Version

0.0.1

## Status

1. Trimming using TrimGalore!
2. Read quality control using FastQC
3. Read mapping using bowtie2
4. Coverage estimate
5. Read filtering using samtools
6. Assembly of endosymbiont genome using ABySS
7. Assembly of host mitogenome using megahit
8. Assembly quality assessment using QUAST

## Installation

```bash
git clone https://github.com/clemensma/symbiontDivider
cd symbiontDivider/Docker
./build.sh
cd ..
```

## Usage

```bash
nextflow run main.nf --reads '*_R{1,2}\\.fastq.gz' --endosymbiont_reference '*_endosymRef\\.fna' -with-docker
```

## ToDo

- (Assembling host genome de novo option)
- Visualisation of quality
- Compression for long time storage
- (CPU Core management)
- (Docker release)
