# -*- coding: utf-8 -*-
"""
Chavez Lab HITAG analysis code
author: Alexander F Kratz

"""
import numpy as np
from typing import List,Dict,Tuple
import pandas as pd
import gzip
from Bio import SeqIO
import os
import subprocess
from progress.bar import Bar

class GenomicLocation:
    def __init__(self,chromosome:int,base:int):
        if chromosome=='X':chromosome=23
        if chromosome not in range(1,25):
            print("Chromosome not found!")
        self.chromosome=chromosome
        self.base=int(base)

class GuideRNA:
    def __init__(self,chromosome:int,cutsite:int,sequence:str,strand:int,transcription_direction:int,target:str):
        if chromosome=='X':chromosome=23
        if chromosome not in range(1,25):
            print("Chromosome not found!")
        self.chromosome=chromosome
        self.sequence = sequence
        self.strand = strand
        self.target = target
        self.base = cutsite
        self.cutsite = cutsite

        self.transcription_direction = transcription_direction

def genomic_distance(a:GenomicLocation|GuideRNA,b:GenomicLocation|GuideRNA)->float:
    if a.chromosome!=b.chromosome:
        return float('inf')
    else:
        return abs(a.base-b.base)

def find_closest(a:GenomicLocation,b:List[GenomicLocation])->GenomicLocation:
    closest = np.inf
    best = None
    for other in b:
        distance = genomic_distance(a,other)
        if distance<closest:
            closest = distance
            best = other
    return best

def revcom(ori):
    trans=str.maketrans('ATGC', 'TACG')
    return ori.upper().translate(trans)[::-1]

def reverse_cigar(cigar:str)->str:
    return "".join(split_cigar(cigar)[::-1])

def split_cigar(cigar:str)->List:
    cigar_tokens = ('M','S','D','I')
    split = []
    growin = ""
    for i in range(0,len(cigar)):
        char = cigar[i]
        if char in cigar_tokens:
            growin = growin+char
            split.append(growin)
            growin=""
        else:
            growin = growin+char
    return split

def cigar_to_string(cigar:str)->str:
    return "".join(map(lambda x:x[-1]*int(x[:-1]),split_cigar(cigar)))

def download_genome():
    url = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
    command = ["wget", "-O",os.path.join("sequences","GRCh38_genomic.fna.gz"),url]
    subprocess.run(command)

chromosome_record_names:Dict[str,int] = {
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

def make_indexes():
    if not os.path.exists(os.path.join("indexes","genome_index.1.bt2")):
        if not os.path.exists(os.path.join('sequences','GRCh38_chromosomes.fna')):
            chromosomes = list()
            with gzip.open(os.path.join("sequences","GRCh38_genomic.fna.gz"),'rt') as handle:
                for record in SeqIO.parse(handle,'fasta'):
                    if record.name in chromosome_record_names:
                        chromosomes.append(record)
            SeqIO.write(chromosomes,os.path.join('sequences','GRCh38_chromosomes.fna'),'fasta')
        command = ['bowtie2-build', '--threads','8', '-f', 'sequences/GRCh38_chromosomes.fna','indexes/genome_index']
        print(command)
        subprocess.run(command)


    if not os.path.exists(os.path.join("indexes","sg_linker_index.1.bt2")):
        command = ['bowtie2-build', '--threads','8', '-f', 'sequences/sg_linkers.fasta','indexes/sg_linker_index']
        print(command)
        subprocess.run(command)
    if not os.path.exists(os.path.join("indexes","tf_linker_index.1.bt2")):
        command = ['bowtie2-build', '--threads','8', '-f', 'sequences/tf_linkers.fasta','indexes/tf_linker_index']
        print(command)
        subprocess.run(command)
    
    
def run_bowtie(input_file_name:str,linker_file_name:str,output_file_name:str):
    command = ['bowtie2','-x',os.path.join('indexes',linker_file_name),
               '-U',input_file_name,
               '-S','{}_aligned_linker.sam'.format(output_file_name),
               '-p','16',
               '--sam-nosq',
               '--local',
               '--rdg', '20,3',
               '--rfg', '20,3',
               ]
    print(" ".join(command))
    subprocess.run(command)
    command = ['bowtie2','-x',os.path.join('indexes','genome_index'),
               '-U',input_file_name,
               '-S','{}_aligned_genome.sam'.format(output_file_name),
               '-p','16',
               '--sam-nosq',
               '--local',
               '--rdg', '20,3',
               '--rfg', '20,3',
               ]
    print(" ".join(command))
    subprocess.run(command)  



def load_genome()->Dict[int,str]:
    chromosomes = dict()
    for record in SeqIO.parse(os.path.join('sequences','GRCh38_chromosomes.fna'),'fasta'):
        if record.name in chromosome_record_names:
            chromosome_number = chromosome_record_names[record.name]
            chromosomes[chromosome_number]=str(record.seq)
    return chromosomes


def load_guides(file_name:str)->List[GuideRNA]:
    df = pd.read_csv(file_name,index_col='Unnamed: 0')
    grnas = list()
    for i in df.index:
        sequence = df.at[i,'grna']
        chromosome = df.at[i,'chromosome']
        cutsite = df.at[i,'cutsite']
        strand = df.at[i,'grna_strand']
        transcription_direction = df.at[i,'ensg_strand']
        target = df.at[i,'target']
        grna =  GuideRNA(
            sequence = sequence,
            chromosome=chromosome,
            cutsite=cutsite,
            strand=strand,
            transcription_direction=transcription_direction,
            target = target,
        )
        grnas.append(grna)
    return grnas

class SAMResult():
    def __init__(self,sam_line:str):
        parsed = sam_line.split('\t')
        self.name:str = parsed[0]
        self.alignment_flags:int = int(parsed[1])
        self.target_name:str = parsed[2]
        self.target_base:int = int(parsed[3])
        self.alignment_quality = parsed[4]
        self.cigar = parsed[5]
        self.seq:str = parsed[9]

def load_sam(file_name:str)->Dict[str,SAMResult]:
    results = dict()
    with open(file_name,'r') as fh:
        fh.readline()
        fh.readline()
        rejected  = 0
        while True:
            line = fh.readline().strip()
            if line =="":return results
            parsed = SAMResult(line)
            key = parsed.name
            #if len(results)==500_000:return results
            if parsed.alignment_flags != 4:
                results[key]=parsed
            else:
                rejected+=1
            
def load_results(file_name:str)->List[Tuple[SAMResult,SAMResult]]:
    print('loading linker results')
    linker_results = load_sam("{}_aligned_linker.sam".format(file_name))
    print('loading genome results')
    genome_results = load_sam("{}_aligned_genome.sam".format(file_name))
    combined = list()
    keys = set(linker_results.keys()).union(set(genome_results.keys()))
    for k in keys:
        if k not in genome_results:continue
        if k not in linker_results:continue
        if linker_results[k].alignment_flags !=4:
            if genome_results[k].alignment_flags!=4:
                combined.append(
                    (linker_results.pop(k),genome_results.pop(k))
                )
    return combined



def process_read(res:Tuple[SAMResult,SAMResult],targets:List[GuideRNA])->Dict:
    result = dict()

    genomic_name = res[1].target_name
    chromosome = chromosome_record_names[genomic_name]
    genomic_base = int(res[1].target_base)-1
    location = GenomicLocation(chromosome,genomic_base)

    
    result['seq']=res[1].seq
    result['genomic_chromosome'] = location.chromosome
    result['genomic_alignment_start'] = res[1].target_base
    result['genomic_alignment_quality'] = res[1].alignment_quality
    result['genomic_direction']=res[1].alignment_flags
    result['genomic_cigar']=res[1].cigar

    result['linker_name']=res[0].target_name
    result['linker_alignment_start']=res[0].target_base
    result['linker_alignment_quality'] = res[0].alignment_quality
    result['linker_direction']=res[0].alignment_flags
    result['linker_cigar']=res[0].cigar

    closest_target = find_closest(location,targets)
    if closest_target == None:
        result['genomic_target'] = "Unknown"
        result['categorization'] = 'offtarget'#specifically due to aligning with Y chromosome
        return result
    
    else:
        distance = genomic_distance(location,closest_target)
        if distance<1000:
            result['cut_site'] = closest_target.cutsite
            result['genomic_target']=closest_target.target
            distance = genomic_distance(location,closest_target)
        else:
            result['genomic_target']="Unknown"
            result['categorization'] = 'offtarget'
            return result

    if result['linker_direction'] != 16:
        result['categorization'] = 'backwards_linker'
        return result
    
    if result['genomic_direction']==16:
        if closest_target.transcription_direction != 1:
            result['categorization'] = 'backwards_genomic'
            return result
        expanded_genomic_cigar = cigar_to_string(result['genomic_cigar'])
        expanded_linker_cigar = cigar_to_string(result['linker_cigar'])
        if 'D' in expanded_genomic_cigar or 'D' in expanded_linker_cigar:
            result['categorization'] = "deletion_within_alignment"
            return result
        assert len(expanded_genomic_cigar) == len(expanded_linker_cigar)

        query_to_target_last_base = expanded_genomic_cigar.rfind('M')
        query_to_linker_first_base = expanded_linker_cigar.find('M')
        target_to_query_last_base = result['genomic_alignment_start'] + expanded_genomic_cigar.count('M') - 1
        linker_to_query_first_base = result['linker_alignment_start'] - 1
    
    elif result['genomic_direction']==0:
        if closest_target.transcription_direction != -1:
            result['categorization'] = 'backwards_genomic'
            return result
        reversed_genomic_cigar = reverse_cigar(result['genomic_cigar'])
        expanded_genomic_cigar = cigar_to_string(reversed_genomic_cigar)
        expanded_linker_cigar = cigar_to_string(result['linker_cigar'])
        if 'D' in expanded_genomic_cigar or 'D' in expanded_linker_cigar:
            result['categorization'] = "deletion_within_alignment"
            return result
        assert len(expanded_genomic_cigar) == len(expanded_linker_cigar)
        query_to_target_last_base = expanded_genomic_cigar.rfind('M')
        query_to_linker_first_base = expanded_linker_cigar.find('M')
        target_to_query_last_base = result['genomic_alignment_start'] - 1
        linker_to_query_first_base = result['linker_alignment_start'] - 1
    else:
        raise(TypeError)
    

    number_of_overlapping_bases = max((
        query_to_target_last_base-query_to_linker_first_base+1,
        0))    
    query_to_target_last_base -= number_of_overlapping_bases
    if result['genomic_direction']==0:
        target_to_query_last_base +=number_of_overlapping_bases
        result['genomic_target_deleted_bases'] = -(result['cut_site'] - target_to_query_last_base)
    elif result['genomic_direction']==16:
        target_to_query_last_base -=number_of_overlapping_bases
        result['genomic_target_deleted_bases'] = result['cut_site'] - target_to_query_last_base
    
    
    
    result['linker_deleted_bases'] =  linker_to_query_first_base

    result['inserted_bases'] = query_to_linker_first_base - query_to_target_last_base - 1
    result['frame'] = (result['linker_deleted_bases'] + result['genomic_target_deleted_bases'] - result['inserted_bases'])%3

    """
    #debug section!
    result['linker_ali_seq']=res[0].seq
    result['genome_ali_seq']=res[1].seq
    result['expanded_genomic_cigar'] = expanded_genomic_cigar
    result['expanded_linker_cigar'] = expanded_linker_cigar
    result['nob'] = number_of_overlapping_bases
    result['qttlb'] = query_to_target_last_base
    result['qtlfb'] = query_to_linker_first_base
    result['ttqlb'] = target_to_query_last_base
    result['ltqfb'] = linker_to_query_first_base
    """



    if result['genomic_target_deleted_bases']==0 and result['linker_deleted_bases']==0 and result['inserted_bases'] == 0:
        result['categorization'] = 'perfect'
        return result
    
    if result['genomic_target_deleted_bases'] > 60:
        result['categorization'] = 'large_deletion'
        return result

    if result['frame']==0:
        result['categorization'] = 'in_frame_indel'
        return result

    if result['frame'] in (1,2):
        result['categorization'] = 'out_frame_indel'
        return result
    


def results_to_df(results,targets)->pd.DataFrame:
    processed_reads= dict()

    bar = Bar("Processing alignments...",max=len(results),suffix='%(index)i / %(max)i - %(eta)ds')
    for i in range(len(results)):
        if(i%100==0):bar.next(100)
        processed_read = process_read(results[i],targets)
        processed_reads[i] = processed_read
    bar.finish()
    df = pd.DataFrame.from_dict(processed_reads, orient='index')           
    return df


def run_analysis(input_file_name:str,linker_file_name:str,output:str,target_file_name:str):
    if not os.path.exists(os.path.join('sequences','GRCh38_genomic.fna.gz')):
        download_genome()
    make_indexes()
    #run_bowtie(input_file_name=input_file_name,linker_file_name=linker_file_name,output_file_name=output)  
    targets = load_guides(target_file_name)
    results = load_results(output)
    df=results_to_df(results,targets)
    df.to_csv('{}_outputs.csv'.format(output))

