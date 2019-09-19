Novel-X: Novel sequence insertion detection using Linked-Reads
======

Novel-X detects and genotypes novel sequence insertions in 10X sequencing dataset using non-trivial read alignment signatures and barocde information.

# Table of contents
1. [Installation](#installation)
2. [Commands Options](#commands-options)
3. [Output Formats](#output-formats)
4. [Example Commands](#example-commands)
5. [Publications](#publications)
6. [Contact & Support](#contact)

## Installation

To start working with Novel-X please clone this repository recursively:

```
git clone --recursive git@github.com:1dayac/Novel-X.git
```

If you clone repository non-recursively Novel-X will not work. To fix this run:

```
git submodule update --init --recursive
```

Novel-X is a pipeline based on a popular Snakemake workflow management system and consists of multiple steps and requires a lot of external sofware.

Python dependencies are listed in requrements.txt file. The can be downloded and installed with following command:

```
pip install -r requirements.txt
```

Following software also should be installed (version numbers used for testing are shown in bracets):

* Longranger (version 2.15) - [Download Page](https://support.10xgenomics.com/genome-exome/software/downloads/latest)
* Velvet (commit 9adf09f) - [GitHub Page](https://github.com/dzerbino/velvet) - outdated but still useful assembler with minimal assumptions about the data
* BlastN (version 2.2.31) - [Download Page](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* Samtools (version 1.7) - [Project Page](http://www.htslib.org/)
* SPAdes (version 3.13) - [Project Page](http://cab.spbu.ru/software/spades/)
* Quast (version 4.4)- [Project Page](http://cab.spbu.ru/software/quast/)

Some of these programms can be installed with conda package. Highly recommended. 
Path to executables should be provided in path_to_executables_config.json file.
Inside bxtools folder run following commands (estimated execution time is around 2 minutes):

```
./configure
make
make install
```


We tested our tool using CentOS Linux 7 OS, but we suppose that it should work at any modern Unix-like system.
 
Then you are ready to go.


## Commands Options

Novel-X can be run with novel-x.py script with two modes:

* run - run pipeline from the scratch
* restart - if previous pipeline was not finished for some reason you can try to catch up with novel-x.py restart command.

A typical command to start Novel-X is 
```
python novel-x.py run --bam my_bam.bam --genome my_genome.fasta --outdir my_dir
```

Optional arguments are:
* --lr20 - needed if you run pipeline on a bam file obtained by LongRanger2.0 pipeline
* --nt - optinal filtering of non-human sequences from the orphan contigs

You can invoke help message by typing:

```
python novel-x.py run --help
```
or

```
python novel-x.py restart --help
```
## Output Formats

Novel-X write results into vcf-file. If your bam-file was named HM2KYBBXX_NA18509.bam, the resulting vcf-file will be named HM2KYBBXX_NA18509.vcf and will be stored inside the outdir folder.

## Example Commands

Run from the start:

```
python ~/Novel-X/novel-x.py run --bam /athena/ihlab/scratch/dmm2017/70_samples_data/HLF3WBBXX_NA12006_longranger.bam -t 8 -m 200 --nt /athena/ihlab/scratch/dmm2017/blast_database/  --genome /athena/ihlab/scratch/dmm2017/hg38/hg38.fa --outdir /athena/ihlab/scratch/dmm2017/70_samples/novelx_NA12006
```

Restart from the last stage:
```
python ~/Novel-X/novel-x.py restart --outdir novelx_NA12006
```

There is a problem on filter_target_contig stage at the moment. It can exit with non-zero exit code. We recommend to comment out the next line before using restart option:

```
parallel --jobs {THREADS} filter_target_contigs ::: {input.contigs}/*
```

## Demo command

We placed a toy dataset in demo folder to test that software is installed correctly. You can run command:

```
python ~/Novel-X/novel-x.py run --bam ~/Novel-X/demo/demo.bam -t 1 -m 1  --genome ~/Novel-X/demo/demo.fasta --outdir out
```

This command takes about 15 minutes to finish on our hardware. It produces a vcf-file with a single vcf record.
```
chr1_25500000_25535000  29503   .       T       TGTATTGTGTGTATGAGGGTTGTGTGCTGTGTGTTGTGTATATATTGTATGTGTTATGTGTATGTATGTCGTATGAGTGTATTCTGTATATGTGTTTTGTGTGGTCTATTATGTATGTGGCATGTGTTGTGTATGTGTGTTGTGTGTGATGTGTTGTATGTGTGTTGTGCATATATGTTGTTTCTGTGTATGTATGTTATGTGTATGTGTATGTTGTGTTGTATGTATGGGTTGTGCCTATGTGCTGTGTTGTGTGCTGCATGCATGTTTGTGTGGTGTGTGTATTTAGGTTGTGTGCTATTTATGTGTCTATATTGTATGTGTTGTATGTGTGTTGTATGTATGTGTAGTGTATGTGTGTTGTGTGTGATGTGTATATGTGGTGTGTGTATGTCTGTTATGTGTATGTATGAGTGTATGTGTGTTGTGTGTGTTGTGTATATGTGTTGTGTGTGTTGTGTATGTGTGTTGTGAGTTGTGTATATGTGGTGTGAGTTGTGTTGTGTCATGTATGTGTGCATTGTGTATAGGTGTTGCATGTGTGTTGTGTTGTGTGTATGTGTTGTGTTGTGTATATGTGGTATGTGAATGTGTATGTTGTATGTTGTGTTGTATGTATATGTGTTATGTATATGTGATGTGTGTGTTGTGTATATGCTGGGTGTGTGTGTACATGTGTGTATGTGTGTTGTATGTATGTGTGTATGCATGTGTGTTGCGTATATGTGGTATGTGTGCATGTGTGTTGTCATGTGTATGTGTGTTGTGTATATGTGTGTGTTGTGTATATGTGTTGTGTGTATGTGTATCATGTTGTGTGTATGTGTTATGTTGTGTATATGTGGTGTGTGAATGTGTGTTGTGTGTATGTGTATGTTGTCTGTTTTGTGTGTGTATACGTGGTGTGTGTGTGTTGTGTTGTGTATATGTGTTGTGTGTGTTGCGTGTATGTGTTGTGTGTT      .       PASS    DP=100  NODE_1_length_4180_cov_43.887399        2776    276     347     1306
```
Output may slightly differ based on your software versions.


## Publications

"Novel sequence insertion detection using Linked-Reads" preprint is available at [https://www.biorxiv.org/content/10.1101/551028v1](https://www.biorxiv.org/content/10.1101/551028v1).

## Contact & Support

Feel free to drop any inquiry to [meleshko.dmitrii@gmail.com]() 
