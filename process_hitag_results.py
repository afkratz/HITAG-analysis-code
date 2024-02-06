# -*- coding: utf-8 -*-
"""
Chavez Lab HITAG post-main-script-processing script
    Identifies off-target reads using the MANE transcript list
    Applies a series of summarizing functions to the results
author: Alexander F Kratz
"""

import pandas as pd
from progress.bar import Bar
import os
import gzip

def process_offtargets(file_name):
    """
    This code performs a few functions to analyze the off-targets.

    First, it groups the off-targets based on what base and chromosome their analysis occurs on.
    This is useful because some of the later code is relatively time-intensive so it's good to 
    minimize the number of times needed to perform that.

    Second, it searches the MANE.GTCh38 file for the locations of the off-targets. This file
    contains a canonical set of transcripted sequences. If an alignment occurs in a region 
    that is identifiable as exactly one transcript, it notes that off-target as being in
    that transcript. Due to the difficulty of telling if an arbitrary location in the genome
    is in an arbitrary frame, we do not attempt to identify sense or frame-ness.

    Importantly, this step is agnostic to what the "true targets" are. As a result, a deletion
    greater than 1kb but in a targeted gene is called an off-target at this stage, 
    even though it will end up getting converted to a "large deletion" in the last step.
    """
    
    df = pd.read_csv('{}_outputs.csv'.format(file_name))
    off_targets = df[df['categorization']=='offtarget']
    
    combined_off_targets = pd.DataFrame()
    for i in off_targets.index:
        chromosome,base = off_targets.at[i,'genomic_chromosome'],off_targets.at[i,'genomic_alignment_start']
        key = str(chromosome)+":"+str(base)
        if key not in combined_off_targets.index:
            combined_off_targets.at[key,'chromosome']=chromosome
            combined_off_targets.at[key,'base']=base
            combined_off_targets.at[key,'count']=0
        combined_off_targets.at[key,'count']+=1
    
    combined_off_targets.to_csv('{}_combined_offtargets.csv'.format(file_name))
            
    bar = Bar("Processing offtargets...",max=len(combined_off_targets),suffix='%(index)i / %(max)i - %(eta)ds')
    for i in combined_off_targets.index:
        bar.next()

        #For each off-target, we take its genomic location and compare it to each transcript
        chromosome,base = combined_off_targets.at[i,'chromosome'],combined_off_targets.at[i,'base']
        found = 0

        for transcript_name in transcript_df.index:
            name = transcript_df.at[transcript_name,'GRCh38_chr']
            if name in chromosome_record_names:
                rec_chrom = chromosome_record_names[name]
            else:
                continue
            chr_start = transcript_df.at[transcript_name,'chr_start']
            chr_end = transcript_df.at[transcript_name,'chr_end']
            
            if chromosome == rec_chrom and base > min(chr_start,chr_end) and base < max(chr_start,chr_end):
                found+=1
                target = transcript_name#save this for later, we will only use it if found==1

        #Done comparing to each transcript, lets see if we have exactly one match
        #More than one match is rare, but does occasionally occur.
        
        combined_off_targets.at[i,'number of transcripts matched']=found
        if found == 1:
            combined_off_targets.at[i,'off-target symbol'] = transcript_df.at[target,'symbol']
            combined_off_targets.at[i,'off-target Ensembl_Gene'] = transcript_df.at[target,'Ensembl_Gene']
            combined_off_targets.at[i,'off-target chr_strand'] = transcript_df.at[target,'chr_strand']
            chr_start = transcript_df.at[transcript_name,'chr_start']
            chr_end = transcript_df.at[transcript_name,'chr_end']
    bar.finish()

    combined_off_targets.to_csv('{}_offtargets_called.csv'.format(file_name))


def combine_called_offtargets(name,targets):
    """
    This code combines the off-targets which are at different bases within the same
    gene into a single count. Additionally, it checks a list of targeted genes (based on ensg)
    to identify if the previously-called off-targets are in fact large deletions. This is 
    encoded in the "in_screen" column of the output.
    """

    df = pd.read_csv('{}_offtargets_called.csv'.format(name),index_col='Unnamed: 0')
    genes = set(df['off-target symbol'])
    odf = pd.DataFrame()
    for gene in genes:
        if gene != gene : 
            continue
        subdf = df[df['off-target symbol']==gene]
        for col in ['chromosome', 'off-target symbol', 'off-target Ensembl_Gene','off-target chr_strand']:
            odf.at[gene,col]=subdf.iloc[0][col]
        ensg = odf.at[gene,'off-target Ensembl_Gene'][:-3]
        odf.at[gene,'in_screen'] = ensg in targets
        odf.at[gene,'count']=subdf['count'].sum()
    odf.to_csv('{}_offtargets_by_gene.csv'.format(name),index=False)

def summarize_by_target(name:str):
    """
    This code summarizes the results for a given condition by each target,
    breaking down the read count of each category of outputs for each target.
    It then integrates the off-target calling which was performed previously.
    """
    df=pd.read_csv('{}_outputs.csv'.format(name))
    odf = pd.DataFrame()
    targets = list(set(df['genomic_target']))
    categories = list(set(df['categorization']))
    for t in targets:
        for c in categories:
            odf.at[t,c]=0
    
    bar = Bar("Summarizing results...",max=len(df),suffix='%(index)i / %(max)i - %(eta)ds')
    for i in range(len(df)):
        bar.next()
        odf.at[
            df.at[i,'genomic_target'],
            df.at[i,'categorization']]+=1
    bar.finish()
       
    odf['thought_was_off_target']=0
    off_target_by_gene = pd.read_csv('{}_offtargets_by_gene.csv'.format(name))
    in_screen = off_target_by_gene[off_target_by_gene['in_screen']==True]
    for i in in_screen.index:
        gene = in_screen.at[i,'off-target symbol']
        odf.at[gene,'thought_was_off_target']=in_screen.at[i,'count']
    odf.at['true_off_targets_unaccounted_for',"offtarget"]=odf.at['Unknown','offtarget']-sum(in_screen['count'])
    odf.drop('Unknown',inplace=True)
    odf.to_csv('{}_target_summary.csv'.format(name))

def category_summarize(name):
    """
    This code summarizes the output of each run based on what category the runs fall into
    These include:
        perfect - a read with no insertions or deletions
        backwards_genomic - arread which aligned in the opposite frame from expected
        deletion_within_alignment - a read which did not align continuously to the genome. Hypothesized to be a spurious read
        large_deletion - a read which aligned to a targeted transcript, but not within 60 bp of the programmed cut-site. 
                         Includes two categories of reads, those which were within 1,000 bp but >60bp, and those which
                         Were more than 1,000 but within a targeted transcript, which were identified in process_offtargets
        inframe_indel - a read which aligned to a guide in frame and in the sense direction, but did have insertions or deletions
        outframe_indel - a read which aligned to a guide in the sense direction, but had out of frame insertions or deletions
        off-target transcripted - a read which aligned to exactly one non-targeted transcript. Sense/frame not identified
        off-target not transcripted - a read which did not align to exactly one transcript
            
    """
    
    df = pd.read_csv('{}_target_summary.csv'.format(name))
    perfect = df['perfect'].sum()
    backwards = df['backwards_genomic'].sum()
    deletion_within_alignment = df['deletion_within_alignment'].sum()
    inframe_indel = df.in_frame_indel.sum()
    outframe_indel = df.out_frame_indel.sum()
    large_del = df['thought_was_off_target'].sum()
    if 'large_deletion' in df.columns:
        large_del+=df['large_deletion'].sum()
    
    
    initial_unaligned = pd.read_csv('{}_offtargets_called.csv'.format(name))['count'].sum()

    transcripted_off_targets = pd.read_csv('{}_offtargets_by_gene.csv'.format(name))

    #This variable is the "initially thought to be off target but identified to a target" category
    fixed_unaligned = transcripted_off_targets[transcripted_off_targets['in_screen']]['count'].sum()

    #This variable is the "
    off_target_transcripted = transcripted_off_targets[~transcripted_off_targets['in_screen']]['count'].sum()

    off_target_not_transcripted = initial_unaligned - fixed_unaligned - off_target_transcripted
    
    odf = pd.DataFrame()
    odf.at['perfect','count']=perfect
    odf.at['backwards','count']=backwards
    odf.at['deletion_within_alignment','count']=deletion_within_alignment
    odf.at['inframe_indel','count']=inframe_indel
    odf.at['outframe_indel','count']=outframe_indel
    odf.at['large_del','count']=large_del
    odf.at['off-target transcripted','count'] = off_target_transcripted
    odf.at['off-target not transcripted','count'] = off_target_not_transcripted
    odf.to_csv('{}_category_summary_overall.csv'.format(name))

def main():

    process_offtargets('results/HCT_SG')
    process_offtargets('results/HEK_SG')
    process_offtargets('results/HEK_TF_F1_S1')
    process_offtargets('results/HEK_TF_F1_S2')
    process_offtargets('results/HEK_TF_F2_S1')
    process_offtargets('results/HEK_TF_F2_S2')


    sg = pd.read_csv('sg_grnas_validated.csv')
    stress_granuale_ens= list(sg['ensg'])
    
    combine_called_offtargets('results/HCT_SG',stress_granuale_ens)
    combine_called_offtargets('results/HEK_SG',stress_granuale_ens)

    tf1 = pd.read_csv('tf_grnas_f1_validated.csv')
    tf1_ens= list(tf1['ensg'])
    combine_called_offtargets('results/HEK_TF_F1_S1',tf1_ens)
    combine_called_offtargets('results/HEK_TF_F1_S2',tf1_ens)

    tf2 = pd.read_csv('tf_grnas_f2_validated.csv')
    tf2_ens= list(tf2['ensg'])
    combine_called_offtargets('results/HEK_TF_F2_S1',tf2_ens)
    combine_called_offtargets('results/HEK_TF_F2_S2',tf2_ens)

    summarize_by_target('results/HCT_SG')
    summarize_by_target('results/HEK_SG')
    summarize_by_target('results/HEK_TF_F1_S1')
    summarize_by_target('results/HEK_TF_F1_S2')
    summarize_by_target('results/HEK_TF_F2_S1')
    summarize_by_target('results/HEK_TF_F2_S2')

    category_summarize('results/HCT_SG')
    category_summarize('results/HEK_SG')
    category_summarize('results/HEK_TF_F1_S1')
    category_summarize('results/HEK_TF_F1_S2')
    category_summarize('results/HEK_TF_F2_S1')
    category_summarize('results/HEK_TF_F2_S2')
    




#we do this as a global variable to avoid re-doing it for each data-set
if not os.path.exists('MANE.GRCh38.v1.3.summary.txt.gz'):
    os.system('wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.3.summary.txt.gz')
with gzip.open('MANE.GRCh38.v1.3.summary.txt.gz') as fh:
    transcript_df = pd.read_csv(fh,sep='\t')

chromosome_record_names = {
    'NC_000001.11':1,
    'NC_000002.12':2,
    'NC_000003.12':3,
    'NC_000004.12':4,
    'NC_000005.10':5,
    'NC_000006.12':6,
    'NC_000007.14':7,
    'NC_000008.11':8,
    'NC_000009.12':9,
    'NC_000010.11':10,
    'NC_000011.10':11,
    'NC_000012.12':12,
    'NC_000013.11':13,
    'NC_000014.9' :14,
    'NC_000015.10':15,
    'NC_000016.10':16,
    'NC_000017.11':17,
    'NC_000018.10':18,
    'NC_000019.10':19,
    'NC_000020.11':20,
    'NC_000021.9' :21,
    'NC_000022.11':22,
    'NC_000023.11':23,#also known as X
    'NC_000024.10':24,#also known as Y
    }


if __name__ == "__main__":
    main()
