# -*- coding: utf-8 -*-
"""
Created on Wed May 18 12:04:53 2022

@author: Alexander Kratz

This script loads a csv from the folder that contains /scripts that has the following structure:
    
    Target	Guide	                Frame
    SEC24C	GTCATATGCACAAGGAGATT	0
    PXDNL	CGCTTCTCTGGGGAATCACT	0
    SETD2	TCTAATTCAGTGTCCTCTTT	0
    ...
    
It looks up the target name from the ENSEMBL server and extracts the ENSg tag for the genomic locus associated with that target.
For 12 of our targets, the name that we used internally does not match the name that the ENSEMBL database uses. 
These were manually filled in below.
The ENSG names are then added in to make a .csv that has the following structure:
    Target	Guide	                Frame	ens
    SEC24C	GTCATATGCACAAGGAGATT	0	    ENSG00000176986
    PXDNL	CGCTTCTCTGGGGAATCACT	0	    ENSG00000147485
    SETD2	TCTAATTCAGTGTCCTCTTT	0	    ENSG00000181555
    ...

"""

import requests, sys
import pandas as pd

input_df=pd.read_csv("../all_guides.csv")#load dataframe of targets/guides


def get_ens(gene_name):
    server = "https://rest.ensembl.org"
    ext = "/lookup/symbol/homo_sapiens/"+gene_name+"?expand=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    decoded = r.json()
    return(decoded['id'])


input_df["ens"]=""#add a blank column for ens tags
for i in range(0,len(input_df)):
    try:input_df.at[i,"ens"]=get_ens(input_df.at[i,"Target"])#try to load the ENS for a target
    except:
        input_df.at[i,"ens"]="UNKNOWN"#if can't load the target, set the target to "UNKNOWN"
    
    print(i,input_df.at[i,"Target"],input_df.at[i,"ens"])

Manually_Found={
        "H1F0":     "ENSG00000189060",  # Manually found ensemble Gene ID's for loci which have multiple names
        "FAM195B":  "ENSG00000225663",  #
        "C19orf43": "ENSG00000123144",  #
        "KIAA0355": "ENSG00000166398",  #
        "HN1L":     "ENSG00000007545",  #
        "H2AFV":    "ENSG00000105968",  # 
        "ATP5A1":   "ENSG00000152234",  #
        "C22orf28": "ENSG00000100220",  #
        "YARS":     "ENSG00000134684",  #
        "C15orf52": "ENSG00000188549",  #
        "SFRS3":    "ENSG00000112081",  #   
        "ORF1":     "ENSG00000213672",  #
        }

for i in range(0,len(input_df)):
    if input_df.at[i,"Target"] in Manually_Found.keys():
        input_df.at[i,"ens"]=Manually_Found[input_df.at[i,"Target"]]

assert "UNKNOWN" not in list(input_df["ens"])#check if there are any targets without an ENSg tag

input_df.to_csv("../all_guides_with_ens.csv",index=False)
