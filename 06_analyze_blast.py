import pandas as pd
import math
from Bio import SeqIO


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
            print()

def revcom(ori):#return the reverse complement of a sequence of DNA (ignores anything that is not a capital A,T,C, or G)
    trans=str.maketrans('ATGC', 'TACG')
    return ori.upper().translate(trans)[::-1]
    

def load_genomic_targets(): #loads the target loci, in the form of a dictionary where the key is the target name and the value is the cutsite +- 500 bp (sometimes truncated if the locus on emsemble does not extend that far)
    name_to_locus={}
    for record in SeqIO.parse("../genomic_targets.fasta","fasta"): 
        name_to_locus[record.name] = str(record.seq)
    return name_to_locus

def get_frames(guides_df):
    name_to_frame={}
    for i in range(0,len(guides_df)):
        name_to_frame[guides_df.at[i,"Target"]]=guides_df.at[i,"Frame"]
    return name_to_frame

def find_cs_fw_guide(seq,grna):  #return how many bases into a sequence is the cutsite of a forward guide in a sequence
    grna_location=seq.find(grna) #find what base the start of the n20 is
    return(grna_location+17)     #return 17 past where the start of the n20 is (forward direction cut site)
    """
    |--------n20---------|PAM
    |nnnnnnnnnnnnnnnnnnnn|NGG
    |----------------|---|NGG
                    CS
    """

def find_cs_rv_guide(seq,grna):          #return how many bases into a sequence is the cutsite of a reverse guide in a sequence
    grna_location=seq.find(revcom(grna)) #find what base the start of the n20 is
    return(grna_location+3)              #return 3 past where the start of the n20 is (reverse direction cut site)
    """
    PAM|--------n20---------|
    CCN|nnnnnnnnnnnnnnnnnnnn|
    CCN|---|----------------|
          CS
    """


def find_cutsites(targets_and_guides_df,target_sequences):
    target_to_cutsite={}
    for i in range(0,len(targets_and_guides_df)):
        grna=targets_and_guides_df.at[i,"Guide"]
        target_seq=target_sequences[targets_and_guides_df.at[i,"Target"]]
        if grna in target_seq:
            target_to_cutsite[targets_and_guides_df.at[i,"Target"]]=find_cs_fw_guide(target_seq,grna)
        elif revcom(grna) in target_seq:
            target_to_cutsite[targets_and_guides_df.at[i,"Target"]]=find_cs_rv_guide(target_seq,grna)
        else:
            assert False
    return target_to_cutsite


def load_blast_results(e_value = 1e-7, n_bases=20):
    target_results=[]
    with open("../genomic_target_alignments.txt","r") as target_fh:
        while(True):
            line = target_fh.readline().strip() #read a line from the results
            if(line==""):break                  #if at end of file, break out of loop
            
            parsed_line=line.split("\t")
            """
            [Query, match, percent match, length of match, Nmismatch, Ngaps, Qstart, Qend, Sstart, Send, Evalue, Bitscore]
            [0      1      2              3                4          5      6       7     8       9     10      11      ]
            
            The line is split into a list of 11 items which correspond to the fields that the BLAST call returned and stored in a file.
            Explanation of fields:
                0 - Query: What item was blasted against the database. Since we issue a number to each read, this is just a number that represents which read was aligned
                1 - Match: What database member this read aligned to. Will be a gene in this file, and a linker in the second
                2 - Percent match: What percent of the alignment is a perfect match. IE, %M=(length-mismatch)/length
                3 - Length of match: How many bases aligned
                4 - Nmismatch: How many mismatches were in the aligned region
                5 - Ngaps: How many gaps were in the aligned region
                6 - Qstart: What base of the query (IE, the read) the alignment begins at
                7 - Qend: What base of the query the alignment ends at
                8 - Sstart: What base of the subject (IE, the target or the linker) the alignment begins at
                9 - Send: What base of the subject the alignment ends at
                10 - Evalue: Odds of a hit of this quality by chance in this database
                11 - Bitscore: Absolute measure of the similarity between query and subject, regardless of database siz 
            """
            
            parsed_line=[int(parsed_line[0]),parsed_line[1],float(parsed_line[2])]+list(map(lambda x:int(x),parsed_line[3:10]))+[float(parsed_line[10]),float(parsed_line[11])]
            #this line converts the query, length, Nmismatch, Ngaps, Qstart, Qend, Sstart, Send values to integers and the percent match, Evalue and Bitscores to floats
            
            if(parsed_line[10]<=e_value and parsed_line[3]>=n_bases): #if this read aligned with at least 20 bases of the target
                target_results.append(parsed_line)                    #add it to the list of results to further process
    linker_results=[]
    with open("../linker_alignments.txt","r") as linker_fh:
        while(True):
            line = linker_fh.readline().strip() #read a line from the results
            if(line==""):break                  #if at end of file, break out of loop
            
            parsed_line=line.split("\t")
            parsed_line=[int(parsed_line[0]),parsed_line[1],float(parsed_line[2])]+list(map(lambda x:int(x),parsed_line[3:10]))+[float(parsed_line[10]),float(parsed_line[11])]

            #[Query, match, percent match, length of match, Nmismatch, Ngaps, Qstart, Qend, Sstart, Send, Evalue, Bitscore]
            #[0      1      2              3                4          5      6       7     8       9     10      11      ]
            
            if(parsed_line[10]<=e_value and parsed_line[3]>=n_bases):#if this read aligned with at least 20 bases of the linker
                linker_results.append(parsed_line)#add it to the list of results to further process
    return(target_results,linker_results)





def filter_target_and_linker_alignments(target_alignments,linker_alignments):    
    
    #first, we go through the alignments, removing any reads that do not align to both sets. IE, remove reads that align to linker but not target or align to target but not linker
    
    filtered_target_alignments=[]
    filtered_linker_alignments=[] #lists to contain the reads that pass the first filter, checking to see if a read aligns to both Target and Linker
    
    while(len(target_alignments)>0 and len (linker_alignments)>0):
        target_read_num=int(target_alignments[-1][0]) #get which read is at the end of the target alignment list (working from end of list for performance reasons)
        linker_read_num=int(linker_alignments[-1][0]) #get which read is at the end of the linker alignment list
       
        while(target_read_num<linker_read_num):#if the target read is a lower number than the linker read, then we need to pop linker reads to look for the pair to this read
            linker_alignments.pop()
            if(len(linker_alignments)==0):break
            linker_read_num=int(linker_alignments[-1][0])
        
        while(linker_read_num<target_read_num):#if the linker read is a lower number than the target read, then we need to pop target reads to look for the pair to this read
            target_alignments.pop()
            if(len(target_alignments)==0):break
            target_read_num=int(target_alignments[-1][0])
            
        if(len(target_alignments)==0 or len(linker_alignments)==0):#check if the lists are empty
            continue
        if(target_alignments[-1][0]==linker_alignments[-1][0]):#if the final for both have the same read count, then we have found a read that aligns to both.
        
            while(int(target_alignments[-1][0])==target_read_num):#For EVERY alignment that we get from this read, append it to the filtered reads
                    filtered_target_alignments.append(target_alignments.pop())
                    if(len(target_alignments)==0):break
                
            while(int(linker_alignments[-1][0])==linker_read_num):
                filtered_linker_alignments.append(linker_alignments.pop())
                if(len(linker_alignments)==0):break

    #next, we make lists of every read that aligns multiple times to one target or to one linker. These reads are ambiguous, and are very rare (~1/10,000 reads). 
    ambiguous_reads=0 #running count of ambiguous reads
    
    target_duplicates=set()
    for i in range(len(filtered_target_alignments)-1):
        if filtered_target_alignments[i][0]==filtered_target_alignments[i+1][0]:
           target_duplicates.add(filtered_target_alignments[i][0])
    
    linker_duplicates=set()
    for i in range(len(filtered_linker_alignments)-1):
        if filtered_linker_alignments[i][0]==filtered_linker_alignments[i+1][0]:
           linker_duplicates.add(filtered_linker_alignments[i][0])

    final_target_alignments=[]
    final_linker_alignments=[]#make list to contain final set of alignments to target and linker, exclusing ambigious alignments
    
    while(len(filtered_target_alignments)>0):
        index = filtered_target_alignments[-1][0]
        if(index not in target_duplicates and index not in linker_duplicates):
            final_target_alignments.append(filtered_target_alignments.pop())
        else:
            filtered_target_alignments.pop()
            ambiguous_reads+=1
    while(len(filtered_linker_alignments)>0):
        index = filtered_linker_alignments[-1][0]
        if(index not in target_duplicates and index not in linker_duplicates):
            final_linker_alignments.append(filtered_linker_alignments.pop())
        else:
            filtered_linker_alignments.pop()
            ambiguous_reads+=1
    return final_target_alignments,final_linker_alignments,ambiguous_reads


def analyze_junction(target_alignment,linker_alignment,name_to_cutsite,name_to_frame):
    #[Query, match, percent match, length of match, Nmismatch, Ngaps, Qstart, Qend, Sstart, Send, Evalue, Bitscore]
    #[0      1      2              3                4          5      6       7     8       9     10      11      ]
    
    result={}
    result["Read"]=target_alignment[0]
    result["Target"] = target_alignment[1]
    result["Linker"] = linker_alignment[1]
    
    #Step 1: Check if the target alignment is backwards. This can happen due to multiple plasmids inserting into the cut-site, with some happening to insert backwards into the locus. If this happens, this read is useless, although the locus probably IS tagged correctly, or it wouldn't have survived selection.
    if target_alignment[8]>target_alignment[9]:
        result["Direction"]='reverse'
        return result#since it's backwards, return from it. No point charachterizing it further
    else:
        result["Direction"]='forward'
    
    #Step 2: Determine which bases of the query aligns to the target and the linker
    query_to_target_last_base = target_alignment[7] #Qend of read alignment to target
    query_to_linker_first_base = linker_alignment[6] #Qstart of read alignment to linker
    
    #Step 3: Determine which bases of the target align to the query
    target_to_query_last_base = target_alignment[9] #Send of read alignment to target
    
    #Step 4: Determine which bases of the linker align to the query
    linker_to_query_first_base = linker_alignment[8] #Sstart of read alignment to linker
    
    
    #Step 5: Refine any "double alignments". By chance, some of the target cut sites have homology to the linker downstream of the cutsite. When we align the query to target and linker, some bases will therefore be aligned to both the target and the linker. Since this is confusing, we determined to assign any "double overlapping" alignments to the linker. This may overstate the disruption to the target, so is actually a conservative "pessimistic" approach
      
    """
    Target     :-------------------------------------------------     \
                   ||||||||||||||||||||                               |
    Read(query):-------------------------------------------------     |This condition does not need adjustment, there is no overlap due to homology between the target and the linker
                                       ||||||||||||||||||||           |
    Linker     :-------------------------------------------------     /
    
    
    
    Target     :-------------------------------------------------     \
                   |||||||||||||||||||||||||                          |
    Read(query):-------------------------------------------------     |This condition needs adjustment, there is a double-alignment overlap due to homology between the target and the linker
                                       ||||||||||||||||||||           |We can just decrease the alignment between the target and the query by the number of bases that double-align
    Linker     :-------------------------------------------------     /
    """
    n_overlap_bases = max([query_to_target_last_base-query_to_linker_first_base+1,0])
    if n_overlap_bases > 0:
        query_to_target_last_base -= n_overlap_bases #
        target_to_query_last_base -= n_overlap_bases #if there is an overlap between the two, remove the overlap from the query-target and target-query alignments
    
    #Step 6: determine number of bases lost from cut_site and linker
    target_bases_lost = name_to_cutsite[target_alignment[1]]-target_to_query_last_base
    result["Target deletions"] = target_bases_lost
    
    linker_bases_lost = linker_to_query_first_base-1
    result["Linker deletions"] = linker_bases_lost
    
    #Step 7: determine number of bases gained as insertions
    inserted_bases = query_to_linker_first_base-query_to_target_last_base-1
    result["Insertions"]=inserted_bases
    
    #Step 8, determine if the overall junction is in frame. If it is not in frame, return
    result["Frame"]=(linker_bases_lost+target_bases_lost-inserted_bases)%3
    if result["Frame"]!=0:
        return result
    
    #Step 9, check if we have a large deletion (over 60 base pairs) - these are reported as a category, since they are likely mechanistically distinct from small indels, and are more likely to seriously impede protein function
    if(abs(target_bases_lost) >60):
        result["Large deletions"]=1
        return result
    result["Large deletions"]=0
    
    #Step 10, determine the change in amino acids. Since different frames have different alignments to the codons, these are handled slightly differently depending on which base of the codon the cut-site is at. For a simpler aproximation, which is never off by more than 1, you could just say AA's changed = int(bases changed/3)
    guide_cut_frame=name_to_frame[target_alignment[1]]
    
    if(guide_cut_frame == 0):    
            """
            Target: ...|123|123|123|123| 
                                        Cut site 
            Linker:                             |123|123|123|123|...
            
            """
            target_amino_acids_lost = math.ceil(target_bases_lost/3)
            linker_amino_acids_lost = math.ceil(linker_bases_lost/3)
            inserted_amino_acids = math.ceil(inserted_bases/3)
    elif(guide_cut_frame == 1):  
            """
            Target: ...|123|123|123|123|12 
                                          Cut site 
            Linker:                               3|123|123|123|...
            
            """
            target_amino_acids_lost = math.ceil((target_bases_lost-2)/3) #The first two bases lost do NOT lose us any amino acids, since that codon was already disrupted by the cut
            linker_amino_acids_lost = math.ceil((linker_bases_lost-1)/3) #The first base does not lose us an amino acid, since that codon was part of the linker to get into frame
            inserted_amino_acids = math.ceil(inserted_bases/3)
    elif(guide_cut_frame == 2):  
            """
            Target: ...|123|123|123|123|1 
                                         Cut site 
            Linker:                              23|123|123|123|...
            
            """
            target_amino_acids_lost = math.ceil((target_bases_lost-1)/3) #The first base does not lose us an amino acid, since that codon was already disrupted by the cut
            linker_amino_acids_lost = math.ceil((linker_bases_lost-2)/3) #The first two bases lost do NOT lose us any amino acids, since that codon was part of the linker to get into frame
            inserted_amino_acids = math.ceil(inserted_bases/3)
    result["Target AA del"]=max((0,target_amino_acids_lost))
    result["Linker AA del"]=linker_amino_acids_lost
    result["AA ins"]=inserted_amino_acids
    return result


def process_filtered_alignments(filtered_target_alignments,filtered_linker_alignments):
    name_to_locus=load_genomic_targets()
    guides_df=pd.read_csv("../all_guides.csv")
    name_to_cutsite=find_cutsites(guides_df,name_to_locus)
    name_to_frame=get_frames(guides_df)
    
    results_by_target={}
    for target in name_to_locus.keys():
        results_by_target[target]={}
        results_by_target[target]["backwards"]=0
        results_by_target[target]["out of frame - 1"]=0
        results_by_target[target]["out of frame - 2"]=0
        results_by_target[target]["perfect"]=0
        results_by_target[target]["in frame indel"]=0
        results_by_target[target]["large deletions"]=0
        for change_type in ("target deletion(AA) ","linker deletion(AA) ","insertion(AA) "):
            for n in range(0,21):
                results_by_target[target][change_type+str(n)]=0
        results_by_target[target]["large insertion"]=0
    
    for i in range(0,len(filtered_target_alignments)):
        analyzed_junction = analyze_junction(filtered_target_alignments[i],filtered_linker_alignments[i],name_to_cutsite,name_to_frame)
        
        target = analyzed_junction["Target"]
        
        if analyzed_junction["Direction"] == "reverse":
            results_by_target[target]["backwards"]+=1
        else:
            if analyzed_junction["Frame"]==1:
                results_by_target[target]["out of frame - 1"]+=1
            elif analyzed_junction["Frame"]==2:
                results_by_target[target]["out of frame - 2"]+=1
            else:
                assert analyzed_junction["Frame"]==0 #sanity check
                
                if analyzed_junction["Large deletions"]==1:
                    results_by_target[target]["large deletions"]+=1
                else:
                    if analyzed_junction["Target deletions"]==analyzed_junction["Linker deletions"]==analyzed_junction["Insertions"]==0:
                        results_by_target[target]["perfect"]+=1
                    else:
                        results_by_target[target]["in frame indel"]+=1
                        
                    target_AA_del=analyzed_junction["Target AA del"]
                    if(target_AA_del)<=20:
                        results_by_target[analyzed_junction["Target"]]["target deletion(AA) "+str(target_AA_del)]+=1
                        
                    linker_AA_del=analyzed_junction["Linker AA del"]
                    if(linker_AA_del)<=20:
                        results_by_target[target]["linker deletion(AA) "+str(linker_AA_del)]+=1
                        
                    AA_ins=analyzed_junction["AA ins"]
                    if(AA_ins)<=20:
                        results_by_target[target]["insertion(AA) "+str(AA_ins)]+=1
                    else:
                        results_by_target[target]["large insertion"]+=1
    
    results_dataframe=pd.DataFrame.from_dict(results_by_target,orient="index")
    results_dataframe.to_csv("../results_by_target.csv")
    
    summary = {} #Sum traits over all targets
    traits = list(results_dataframe.columns)
    
    for trait in traits:
        summary[trait]=sum(results_dataframe[trait])
    summary_dataframe=pd.DataFrame.from_dict(summary,orient="index")
    summary_dataframe.rename({0:"count"},inplace=True,axis="columns")
    summary_dataframe.index.rename("Trait",inplace=True)
    summary_dataframe.to_csv("../summed_over_targets.csv")

print("Loading alignments...")
target_alignments,linker_alignments=load_blast_results()
print("Alignments loaded\nFiltering alignments")
filtered_target_alignments,filtered_linker_alignments,n_ambiguous_reads=filter_target_and_linker_alignments(target_alignments,linker_alignments)
print("Alignments filtered\nProcessing junctions")
process_filtered_alignments(filtered_target_alignments,filtered_linker_alignments)
