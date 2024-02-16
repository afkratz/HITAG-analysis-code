

## Chavez lab HITAG analysis code

This repo contains code and a few sequences related to the paper "High-throughput tagging of endogenous loci for rapid characterization of protein function"

##### Overview:

We provide a bash script, paper_analysis.sh, which provides examples of calling the other scripts in the repo to use them as was done in the paper.

The main script is hitag.py, which is the entry point to call bowtie2, which we use to align our reads to the genome and our libraries of linkers.

Additionally, we provide a script, process_hitag_results.py, which identifies off-targets and further categorizes reads into a few categories.



##### Requirements:

See requirements.txt
