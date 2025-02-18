#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 11:14:37 2025

@author: hasan.yilmaz
"""

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pysam
import pandas as pd
import sys

def is_complete(read, refLengths, adapter5, adapter3, margin):

    #Check read start and read end. E.g. for adapter5=24, adapter=30 and margin=5
    #Read has to start before 19nt and end after tRNALength + 25th nucleotide with respect to reference positions
    if 'His' in str(read.reference_name):
        adapter5 = 44
    if read.reference_start <= adapter5-margin and \
        read.reference_end >= refLengths[read.reference_name]-adapter3+margin:
            return True
    return False

def filter_bam(readList, out_dir, is_complete, len5, len3, margin, qmap, complete=False, fragment=False, firststrand=False):
    
    print('tRNA Counting Initiated!')

    
    #Create empty final dataframe
    df = pd.DataFrame()
    
    #Iterate through samples
    for readName in list(readList.keys()):
        
        print(f'{readName}')
        
        read_path = readList[readName]
        
        bamfile = pysam.AlignmentFile(read_path, "rb", threads = 8)
        
        outputBam = pysam.AlignmentFile(f'{out_dir}/{readName}_filtered.bam', "wb", template=bamfile)
        
        #Create a dictionary for every tRNAs are keys and their lengths are values
        refLengths = dict(zip(bamfile.references, bamfile.lengths))
        
        # Create dictionary of tRNAs which every values are 0
        counts = {key: 0 for key in list(refLengths.keys())}
        counts.update({'antisense': 0})
        counts.update({'not_unique': 0})


        for read in bamfile:
            if complete:
                if not (is_complete(read, refLengths, adapter5=len5, adapter3=len3, margin=margin)):
                    continue
            if fragment:
                if (is_complete(read, refLengths, adapter5=len5, adapter3=len3, margin=margin)):
                    continue
            if read.is_secondary or read.is_supplementary:
                continue

            if read.mapq>=qmap:
                if not firststrand and read.is_reverse \
                   or firststrand and not read.is_reverse:
                    counts['antisense'] +=1
                else:
                    outputBam.write(read)
                    # Increase the count for tRNA of the read that passes filtering
                    counts[read.reference_name] +=1
            elif read.reference_name.startswith("oligo"):
                counts[read.reference_name] +=1
            else:
                counts['not_unique'] +=1
        
        bamfile.close()
        
        print(f'{readName} read and processed! next file..')
        
        tempdata = pd.DataFrame.from_dict(counts, orient='index')
        tempdata.columns = [readName]
        if df.empty:
            df = tempdata
        else:
            df = pd.concat([df, tempdata], axis=1)
    return(df)
    

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='''Filter .bams and count tRNAs''',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--inputList', help='tab separated file that includes sample name and directory information of file')
    parser.add_argument('-o', '--output', help='directory for outputs')
    parser.add_argument('-l5', '--length5', help="length of 5' adapter", type=int, default=24)
    parser.add_argument('-l3', '--length3', help="length of 3' adapter", type=int, default=30)
    parser.add_argument('-m', '--margin', help="margin for adapter length", type=int, default=5)
    parser.add_argument('-q', '--quality', help="mapq", type=int, default=10)
    parser.add_argument('-c', '--complete', help='count only full length tRNAs', action='store_true')
    parser.add_argument('-nc', '--notComplete', help='count only not full length tRNAs', action='store_true')
    parser.add_argument('-fs', '--firststrand', help="antisense strand sequencing",
                        action='store_true')

    
    arguments = parser.parse_args()  
    
    sampleList = dict(zip(pd.read_csv(arguments.inputList, sep='\t')['Sample'],
               pd.read_csv(arguments.inputList, sep='\t')['Path']))
    
    out_dir = arguments.output

    tRNA_counts = filter_bam(readList=sampleList, out_dir=out_dir, is_complete=is_complete,
                             len5=arguments.length5, len3=arguments.length3,
                             margin=arguments.margin, qmap=arguments.quality,
                             complete=arguments.complete, fragment=arguments.notComplete,
                             firststrand=arguments.firststrand)
    
    tRNA_counts.index.name = 'reference'
    tRNA_counts = tRNA_counts.reset_index()
    tRNA_counts.to_csv(f'{out_dir}/counts.tsv', sep='\t', header=True, index=False)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nMission Aborted!\nReturning to the Earth!     \n")