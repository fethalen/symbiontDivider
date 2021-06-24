# symbiontDivider
An easy to use pipeline to separate endosymbiont genomes from their host's

## Version

0.0.1

## Status

1. Trimming using TrimGalore!
2. Read quality control using FastQC
3. De Novo assembly using megahit
4. Filtering endosymbiont genome using blastn
5. Filtering host mitogenome using blastn
6. Read mapping for coverage estimate using bowtie2
7. Coverage estimate
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
I highly recommend unpacking your reads before using! 

## ToDo

- Visualisation of quality
- Compression for long time storage
- (CPU Core management)
