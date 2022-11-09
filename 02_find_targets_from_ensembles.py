# -*- coding: utf-8 -*-
"""
Created on Wed May 18 13:05:50 2022

@author: Alexander Kratz
"""

import requests, sys
import pandas as pd

def get_ensg(ensg):
    server = "https://rest.ensembl.org"
    ext = "/sequence/id/"+ensg+"?"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    return(r.text)

def revcom(ori):
    trans=str.maketrans('ATGC', 'TACG')
    return ori.upper().translate(trans)[::-1]




guides_df=pd.read_csv("../all_guides_with_ens.csv")
guides_df["direction"]=""
target_record_dictionary={}
for i in range(0,len(guides_df)):
    ens=guides_df.at[i,"ens"]
    grna=guides_df.at[i,"Guide"]
    record=get_ensg(ens)
    if grna in str(record):
        location=record.find(grna)
        if(record[location+21:location+23]=="GG"):#Check for PAM downstream of guide
            target_record_dictionary[guides_df.at[i,"Target"]]=(record,location,'f')
            guides_df.at[i,"direction"]="f"
            
    elif revcom(grna) in str(record):
        location=record.find(revcom(grna))
        if(record[location-3:location-1]=="CC"):#Check for PAM on opposite strand
            target_record_dictionary[guides_df.at[i,"Target"]]=(record,location,'r')
            guides_df.at[i,"direction"]="r"
    else:
        print("Cut site not found for target:",guides_df.at[i,"Target"]) #Error, didn't find the cutsite
        quit()

    print(i,guides_df.at[i,"Target"],location,guides_df.at[i,"direction"])

target_sequence_area={}
for target in target_record_dictionary.keys():
    entry=target_record_dictionary[target]
    sequence_record=entry[0]
    location = entry[1]
    target_sequence_area[target]=sequence_record[location-500:location+500] # Grab the sequence from -500 bp to +500 bp around the cut site


from Bio import SeqIO
from Bio import Seq
records=[]
for target in target_sequence_area.keys():
    records.append(SeqIO.SeqRecord(Seq.Seq(target_sequence_area[target]),id=target,description=""))
    
with open ("../genomic_targets.fasta","w") as fh:
    SeqIO.write(records,fh,"fasta")
    