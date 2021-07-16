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

## Requirements

- Docker
- Nextflow

## Installation

```bash
git clone https://github.com/clemensma/symbiontDivider
cd symbiontDivider/Docker
chmod +x build.sh
./build.sh
cd ..
```

## Usage

For read files named something like 'name_of_sample_R1.fastq' and 'name_of_sample_R2.fastq' use the following.
```bash
nextflow run main.nf --reads '*_R{1,2}\.fastq' --endosymbiont_reference '*_endosymRef\.fna' -with-docker
```
If your read files are named differently adjust the pattern specified accordingly.

Also I highly recommend to unpack your files before starting the programm!


## ToDo

- Visualisation of quality
- Compression for long time storage
- Better strategy for endosymbiont assembly
- (CPU Core management)
