# HITAG data analysis pipeline

This repository contains the code and configuration files needed to analyze our NGS results to quantify the junctions produced by HITAG. 



[TOC]

#### Description

This pipeline intakes a configuration file in the form of a CSV that species the name of targets and the gRNA that targets that gene. 

01_find_ensembles.py - Queries ensembl server to identify the ENSG associated with the gene name. Some target names were found manually due to having different names than ENSG uses.



02_find_targets_from_ensembles.py -  Queries ensemble server for the genomic region associated with targets found in previous step, checks for PAM within their sequence, and saves the location +-500bp from the cut site into a fasta.



03_load_and_bin_reads.py - Reads fastq files, removes any reads that consist entirely of N's, concatenates remaining reads into a single fasta.



04_make_db.cmd - Runs blast commands to make databases from the genomic targets found in 02 and the linkers.



05_run_blast.cmd - Blasts the reads against the genomic targets database and the linkers, reporting top hit for each read in blastn tabular format.



06_analyze_blast.py - Processes alignments produced by blast, reports junction characteristics.

#### Configuration Requirements

Pipeline requires the following files:

- all_guides.csv, specifies target and gRNA sequence. Example provided.
- linkers.fasta, species linker sequence(s). Example provided.
- Fastq files of NGS reads. Example not provided, but available for download at https://www.ncbi.nlm.nih.gov/sra under BioProject accession number: PRJNA895413.

#### Installation Requirements

- Python Environment (3.8.10)
  - Pandas
  - Biopython
- NCBI-blast (2.13.0)