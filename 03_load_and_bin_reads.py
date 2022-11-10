# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:44:41 2022

@author: Alexander Kratz

This file intakes 4 fastq's, filters out reads that consisted entirely of N's, and then concatenates the remaining reads into one larger file. 
Reads are available for download at https://www.ncbi.nlm.nih.gov/sra under BioProject accession number: PRJNA895413.

"""
from Bio import SeqIO
import pandas as pd

def revcom(ori):
    trans=str.maketrans('ATGC', 'TACG')
    return ori.upper().translate(trans)[::-1]

results=[]


print("EDIT THIS FILE TO LOAD READS, EXAMPLE PROVIDED BELOW")
"""
print("reading file 1...")
for record in SeqIO.parse('../reads/A252_R1.fastq','fastq'):
    fw_seq=revcom(str(record.seq))
    if fw_seq.count("N")==len(fw_seq): continue
    results.append(fw_seq)
print("reading file 2...")
for record in SeqIO.parse('../reads/A258_R1.fastq','fastq'):
    fw_seq=revcom(str(record.seq))
    if fw_seq.count("N")==len(fw_seq): continue
    results.append(fw_seq)
print("reading file 3...")
for record in SeqIO.parse('../reads/A272_R1.fastq','fastq'):
    fw_seq=revcom(str(record.seq))
    if fw_seq.count("N")==len(fw_seq): continue
    results.append(fw_seq)
print("reading file 4...")
for record in SeqIO.parse('../reads/A279_R1.fastq','fastq'):
    fw_seq=revcom(str(record.seq))
    if fw_seq.count("N")==len(fw_seq): continue
    results.append(fw_seq)

with open("../reads_combined.fasta","w") as oh:
    for i,sequence in enumerate (results):
        oh.write(">"+str(i)+"\n")
        oh.write(sequence+"\n")
"""
