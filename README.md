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
7. De novo assembly using megahit
8. Finding of mitogenome using NCBI blastn
9. Assembly quality assessment using QUAST

## Installation

```bash
git clone https://github.com/clemensma/symbiontDivider
cd symbiontDivider/Docker
./build.sh
cd ..
```

## Usage

```bash
nextflow run main.nf --reads '*_R{1,2}\\.fastq' --endosymbiont_reference '*_endosymRef\\.fna' -with-docker
```
I highly recommend to unpack your files before starting the programm!

## ToDo

- Visualisation of quality
- Compression for long time storage
- Better strategy for endosymbiont assembly
- (CPU Core management)
