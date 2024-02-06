# -*- coding: utf-8 -*-
"""
Chavez Lab HITAG argument parse/main file to call
Args:
    -i,--input  : input file of a .fastq that contains the reads you want to search for linker and genomic alignments
    -o,--output : name of the run which will be used to record results
    -l,--linker : input file of a .fasta that contains the linker sequences you expect to find in the read
    -g,--guides : input file of a .csv with columns target,grna,frame,ensg,ensg_strand,grna_strand,chromosome,cutsite. Examples in grna/
    --validate  : flag, if true checks if the guide sequences are present in the genome with a NGG PAM in the expected locations, default false
    --cleanup   : flag to delete alignments after done, default is false
author: Alexander F Kratz
"""
import hitag_src as hs
import argparse
import os

def main():
    parser=argparse.ArgumentParser(
    )

    parser.add_argument("--validate",action='store_true')
    
    parser.add_argument("--cleanup",action='store_true')

    parser.add_argument("-g","--guides",
        type=str,required=True,help=('file path to guide rna file'))

    parser.add_argument('-l','--linker',
        type=str,required=True,help=('.fasta of the linkers to align to'))

    parser.add_argument('-i','--input',
        type=str,required=True,help=(".fastq of the reads"))

    parser.add_argument('-o','--output',
        type=str,required=True,help=(".csv of output")
    )

    #If true, then write results with a lot more info about the shape of alignments
    #Ultimately gets passed to hitag_source.process_read
    parser.add_argument("--debug",action='store_true')

    args = parser.parse_args()
    if not os.path.exists(args.guides):
        raise FileNotFoundError(args.guides)
    
    if not os.path.exists(args.linker):
        raise FileNotFoundError(args.linker)
    
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    guides = hs.load_guides(args.guides)
    if not os.path.exists('sequences'):os.mkdir('sequences')
    if not os.path.exists(os.path.join('sequences','GRCh38_genomic.fna.gz')):
        hs.download_genome()

    if args.validate:
        print('Validating {} guides...'.format(len(guides)))
        if hs.validate_guides(guides):
            print('Guides validated')
        else:
            print('Guides invalid')
            quit()
    
    hs.run_analysis(
        input_file_name=args.input,
        linker_file_name=args.linker,
        output=args.output,
        guides=guides,
        debug=args.debug,
    )
    if args.cleanup:
        hs.cleanup(output=args.output)
    
    

if __name__ == '__main__':
    main()