# -*- coding: utf-8 -*-
"""
Chavez Lab HITAG analysis code
    src file of functions imported by hitag.py
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
            raise(ValueError("Chromosome {} invalid".format(chromosome)))
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

    if not os.path.exists(os.path.join('sequences','GRCh38_chromosomes.fna')):
            chromosomes = list()
            with gzip.open(os.path.join("sequences","GRCh38_genomic.fna.gz"),'rt') as handle:
                for record in SeqIO.parse(handle,'fasta'):
                    if record.name in chromosome_record_names:
                        chromosomes.append(record)
            SeqIO.write(chromosomes,os.path.join('sequences','GRCh38_chromosomes.fna'),'fasta')

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

def make_genome_index():
    if not os.path.exists('indexes'):os.mkdir('indexes')
    if not os.path.exists(os.path.join("indexes","genome_index.1.bt2")):
        command = ['bowtie2-build', '--threads','8', '-f', 'sequences/GRCh38_chromosomes.fna','indexes/genome_index']
        print(command)
        subprocess.run(command)

def make_other_index(index_input_file):
    if not os.path.exists('indexes'):os.mkdir('indexes')
    _,index_input_file_name = os.path.split(index_input_file)
    name,ext = os.path.splitext(index_input_file_name)
    if not os.path.exists(os.path.join("indexes","{}.1.bt2".format(name))):
        command = ['bowtie2-build', '--threads','8', '-f', index_input_file,'indexes/{}'.format(name)]
        print(command)
        subprocess.run(command)
    
    
def run_bowtie(input_file_name:str,linker_file_name:str,output_file_name:str):
    _,linker_file_name = os.path.split(linker_file_name)
    linker_name,ext = os.path.splitext(linker_file_name)
    command = ['bowtie2',
               '-x',os.path.join('indexes',linker_name), #Target library of sequences
               '-U',input_file_name, #Input fastq
               '-S','{}_aligned_linker.sam'.format(output_file_name),#Output file name
               '-p','16',     #Run on 16 threads
               '--sam-nosq',  #Skips a bunch of intro data in the output
               '--local',     #Generate a local alignment, ie an alignment within the sequence, not trying to align the whole sequence.
               '--rdg','20,3',#Assign a gap-open penalty of 20 and gap-extend of 3. This is to strongly penalize gapped alignments which otherwise cause many spurious alignments
               '--rfg','20,3',#Same as above but for the target instead of the query
               ]
    print(" ".join(command))
    subprocess.run(command)
    command = ['bowtie2',
               '-x',os.path.join('indexes','genome_index'), #target library of sequences
               '-U',input_file_name, #Input fastq
               '-S','{}_aligned_genome.sam'.format(output_file_name), #Output file name
               '-p','16',     #Run on 16 threads
               '--sam-nosq',  #Skips a bunch of intro data in the output
               '--local',     #Generate a local alignment, ie an alignment within the sequence, not trying to align the whole sequence.
               '--rdg','20,3',#Assign a gap-open penalty of 20 and gap-extend of 3. This is to strongly penalize gapped alignments which otherwise cause many spurious alignments
               '--rfg','20,3',#Same as above but for the target instead of the query
               ]
    print(" ".join(command))
    subprocess.run(command)  

def cleanup(output_file_name:str):
    if os.path.exists('{}_aligned_genome.sam'.format(output_file_name)):
        os.remove('{}_aligned_genome.sam'.format(output_file_name))
    if os.path.exists('{}_aligned_linker.sam'.format(output_file_name)):
        os.remove('{}_aligned_linker.sam'.format(output_file_name))

def load_genome()->Dict[int,str]:
    chromosomes = dict()
    for record in SeqIO.parse(os.path.join('sequences','GRCh38_chromosomes.fna'),'fasta'):
        if record.name in chromosome_record_names:
            chromosome_number = chromosome_record_names[record.name]
            chromosomes[chromosome_number]=str(record.seq).upper()
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

            #4 corresponds to no alignment. We discard those reads
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
                    #We store the filtered results as a list of tuples, 
                    #where the 0th item is the genome alignment and the 1st item is the linker alignment
                    (genome_results.pop(k),linker_results.pop(k))
                )
    return combined


def validate_guides(guides: List[GuideRNA])->bool:
    genome = load_genome()
    for guide in guides:
        sequence = guide.sequence
        chromosome = guide.chromosome
        cutsite = guide.cutsite 
        strand = guide.strand
        if strand == 1:
            genome_sequence  = genome[chromosome][cutsite-17:cutsite+3]
            if genome_sequence!=sequence:return False
            if genome[chromosome][cutsite+4:cutsite+6]!='GG':return False
        else:
            genome_sequence = genome[chromosome][cutsite-3:cutsite+17]
            if genome_sequence!=revcom(sequence):return False
            if genome[chromosome][cutsite-6:cutsite-4]!='CC':return False
    return True


def process_read(res:Tuple[SAMResult,SAMResult],targets:List[GuideRNA],debug=False)->Dict:
    result = dict()
    genomic_res = res[0]
    linker_res = res[1]

    genomic_name = genomic_res.target_name
    chromosome = chromosome_record_names[genomic_name]
    genomic_base = int(genomic_res.target_base)-1
    location = GenomicLocation(chromosome,genomic_base)

    
    result['seq']=genomic_res.seq#Could also get this from the linker alignment
    result['genomic_chromosome'] = location.chromosome
    result['genomic_alignment_start'] = genomic_res.target_base
    result['genomic_alignment_quality'] = genomic_res.alignment_quality
    result['genomic_direction']=genomic_res.alignment_flags
    result['genomic_cigar']=genomic_res.cigar

    result['linker_name']=linker_res.target_name
    result['linker_alignment_start']=linker_res.target_base
    result['linker_alignment_quality'] = linker_res.alignment_quality
    result['linker_direction']=linker_res.alignment_flags
    result['linker_cigar']=linker_res.cigar

    closest_target = find_closest(location,targets)
    if closest_target == None:
        #Closest target can fail to find any closest target if there are no targets on the chromosome the read aligns to
        result['genomic_target'] = "Unknown"
        result['categorization'] = 'offtarget'
        return result
    

    else:
        
        distance = genomic_distance(location,closest_target)
        if distance<1000:
            #If the cutsite is within 1kb of a cutsite, we assume that is the target this junction is related to
            result['cut_site'] = closest_target.cutsite
            result['genomic_target']=closest_target.target
            distance = genomic_distance(location,closest_target)

        else:
            #If we are greater than 1kb from any cutsite, we assume this is off-target.
            #We will check these later to see if they are in a transcript
            result['genomic_target']="Unknown"
            result['categorization'] = 'offtarget'
            return result

    #If the linker is, for some reason, aligning on the opposite strand as we expect, 
    # then something very strange is happening on the read. It should be discarded as useless
    if result['linker_direction'] != 16:
        result['categorization'] = 'backwards_linker'
        return result
    

    """
    There are two directions a gene can be transcribed in, and two directions that our read can align to the genome
    The correct orientation is where they're reversed from eachother, because base one of our read 
    
    FORWARD GENE, CORRRECT READ:
        base 1 of chromosome -------------------------- base N of chromosome
                        ---target gene transcription--->
                                                    <---our read---
                                                    ||||||  ||||
                                                    genome  ||||
                                                            linker
    FORWARD GENE, BACKWARDS READ:
        base 1 of chromosome -------------------------- base N of chromosome
                            ---target gene transcription--->
                                                    ---our read--->
                                                    ||||    ||||||
                                                    linker  genome

    REVERSE GENE, CORRECT READ:
        base 1 of chromosome -------------------------- base N of chromosome
                            <---target gene transcription---
                 ---our read--->
                 ||||    ||||||
                 linker  genome

    REVERSE GENE, BACKWARDS READ:
        base 1 of chromosome -------------------------- base N of chromosome
                        <---target gene transcription---
                        <---our read---
                        ||||||  ||||
                        genome  ||||
                                linker

    The "backwards reads" can occur because plasmids can enter the genome as multi-mers, 
    where (presumably) upstream of this there is a correct junction, but the cell then repairs the junction
    with additional copies of the opened up plasmid.
    Regardless, these reads are useless and cannot be evaluated as junctions. We therefore note them backwards genomic and continue.

    If the read is in the correct orientation, we calculate a number of values:

                                *
    Genome   -------------------a-------------------------------
                |||||||||||||||||
    Read    :-------------------bc-----------------------
                                 ||||||||||||||||||||||||
    Linker  :                    d--------------------------------------------------
    
    In this diagram, the points *,a,b,c & d correspond to:
        * : intended cut site of the grna
        a : last base of the genome (in the direction of transcription) that aligns to the read
        b : last base of the read that aligns to the genome
        c : first base of the read that aligns to the linker
        d : first base of the read that aligns to the linker

        if a == calculated cut site, b == c-1, and d==1, then this is a perfect junction    
    
    We do need to do one adjustment to account for a particular case:
    Genome  :-------------------------------------------------     \
                |||||||||||||||||||||||                            |In this condition, some bases of the read align to both the linker and genome. This is due to homology between the target downstream of the cut-site and the linker
    Read    :-------------------------------------------------     |We fix this by decreasing the alignment between the genome and the read by the number of bases that double-align
                                    ||||||||||||||||||||           |In effect, we assign those bases to the linker
    Linker  :-------------------------------------------------     /
    

    What various indels look like:

    Insertion: bases in read do not align to either the linker or the genome
                                *
    Genome   -------------------a-------------------------------
                |||||||||||||||||
    Read    :-------------------b---c-----------------------
                                    ||||||||||||||||||||||||
    Linker  :                       d--------------------------------------------------

    
    Deletion from genome: The read aligns to the genome but stops before the expected cut site
                                *
    Genome   ---------------a-------------------------------
                |||||||||||||
    Read    :---------------bc-----------------------
                             ||||||||||||||||||||||||
    Linker  :                d--------------------------------------------------

    
    Deletion from linker: The read aligns to the linker, but not at the first base of the linker.
                                *
    Genome   -------------------a-------------------------------
                |||||||||||||||||
    Read    :-------------------bc-----------------------
                                 ||||||||||||||||||||||||
    Linker  :               -----d----------------------------------------------



    """         
    


    if result['genomic_direction']==16:
        if closest_target.transcription_direction == -1:
            #If the direction of alignment is backwards as well as the direction of transcription, the read is backwards from what we expect
            result['categorization'] = 'backwards_genomic'
            return result
        
        expanded_genomic_cigar = cigar_to_string(result['genomic_cigar'])
        expanded_linker_cigar = cigar_to_string(result['linker_cigar'])

        if 'D' in expanded_genomic_cigar or 'D' in expanded_linker_cigar:
            #Rarely, reads will have alignments to the genome with internal deletions
            #This is despite our heavy penalization of insertions and deletions in BOWTIE.
            #It is difficult to explain how this could occur, and we believe they are likely spurious alignments
            #Regardless, we note them here and do not further analyze them. 
            result['categorization'] = "deletion_within_alignment"
            return result
        
        #Sanity check
        assert len(expanded_genomic_cigar) == len(expanded_linker_cigar)

        query_to_target_last_base = expanded_genomic_cigar.rfind('M')
        query_to_linker_first_base = expanded_linker_cigar.find('M')
        target_to_query_last_base = result['genomic_alignment_start'] + expanded_genomic_cigar.count('M') - 1
        linker_to_query_first_base = result['linker_alignment_start'] - 1
    
    elif result['genomic_direction']==0:
        if closest_target.transcription_direction == 1:
            #If the direction of alignment is forward as well as the direction of transcription, the read is backwards from what we expect
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

    #debug section!
    if debug:
        result['linker_aligned_seq']=linker_res.seq
        result['genome_aligned_seq']=genomic_res.seq
        result['expanded_genomic_cigar'] = expanded_genomic_cigar
        result['expanded_linker_cigar'] = expanded_linker_cigar
        result['number_of_overlapping_bases'] = number_of_overlapping_bases
        result['query_to_target_last_base'] = query_to_target_last_base
        result['query_to_linker_first_base'] = query_to_linker_first_base
        result['target_to_query_last_base'] = target_to_query_last_base
        result['linker_to_query_first_base'] = linker_to_query_first_base



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
    


def results_to_df(results,targets,debug=False)->pd.DataFrame:
    processed_reads= dict()
    bar = Bar("Processing alignments...",max=len(results),suffix='%(index)i / %(max)i - %(eta)ds')
    for i in range(len(results)):
        if(i%100==0):bar.next(100)
        processed_read = process_read(results[i],targets,debug=debug)
        processed_reads[i] = processed_read
    bar.finish()
    df = pd.DataFrame.from_dict(processed_reads, orient='index')           
    return df


def run_analysis(
        input_file_name:str,
        linker_file_name:str,
        output:str,
        guides:List[GuideRNA],
        debug:bool,
        ):
    if not os.path.exists('indexes'):
        os.mkdir('indexes')
    make_genome_index()
    make_other_index(linker_file_name)
    run_bowtie(input_file_name=input_file_name,linker_file_name=linker_file_name,output_file_name=output)  
    results = load_results(output)
    df=results_to_df(results,guides)
    df.to_csv(os.path.join('results','{}_outputs.csv'.format(output)))
