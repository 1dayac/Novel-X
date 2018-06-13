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

Novel-X is a pipeline based on a popular Snakemake workflow management system and consist of multiple steps and requires a lot of external sofware.

Python dependencies are listed in requrements.txt file. The can be downloded and installed with following command:

```
pip install -r reqiurements.txt
```

Following software also should be installed:

* Longranger (version >= 2.15) - [Download Page](https://support.10xgenomics.com/genome-exome/software/downloads/latest)
* Velvet - [GitHub Page](https://github.com/dzerbino/velvet) - outdated but still useful assembler with minimal assumptions about the data
* BlastN - [Download Page](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* Samtools - [Project Page](http://www.htslib.org/)
* SPAdes - [Project Page](http://cab.spbu.ru/software/spades/)
* Quast - [Project Page](http://cab.spbu.ru/software/quast/)

Some of these programms can be installed with conda package. Highly recommended. 
Path to executables should be provided in path_to_executables_config.json file. Then you are ready to go.


## Commands Options

Novel-X can be run with novel-x.py script with two modes:

* run - run pipeline from the scratch
* restart - if previous pipeline was not finished for some reason you can try to catch up with novel-x.py restart command.

A typical command to start Novel-X is 
```
python novel-x.py run --bam my_bam.bam --genome my_genome.fasta --out-dir my_dir
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

## Example Commands

## Publications

"Novel sequence insertion detection using Linked-Reads" is submitted to ECCB-2018 conference. Text is available by request.

## Contact & Support

Feel free to drop any inquiry to [dmitrii.meleshko@gmail.com]() 
