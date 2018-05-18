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

You can invoke message by typing:

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
