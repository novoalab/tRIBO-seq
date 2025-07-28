#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 11:39:07 2025

@author: hasan.yilmaz
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from collections import defaultdict, Counter
import pysam
import pandas as pd
import numpy as np
import gzip
import re
import sys

def is_complete(read, refLengths, adapter5, adapter3, margin):

    #Check read start and read end. E.g. for adapter5=24, adapter=30 and margin=0
    #Read has to start before 19nt and end after tRNALength + 25th nucleotide
    #with respect to reference positions
    if 'His' in str(read.reference_name):
        adapter5 = 44
    if read.reference_start <= adapter5 and \
        read.reference_end >= refLengths[read.reference_name]-adapter3+margin:
            return True
    return False

def three_prime_adapter_check(read, refLengths, margin):
    if read.reference_end >= refLengths[read.reference_name]-margin:
        return True
    return False


def parse_fasta(fasta_path):
        
    from Bio import SeqIO
            
    print("Fasta File is Reading...")
    # Create empty variables for iteration and empty dictionary for transcript sequences
    key = []
    seq = []
    fasta_dict = {}
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Take sequence ids
            key = record.id
            # Take sequence of that specific transcript
            seq = record.seq.upper()
            # Put the name and the sequence in dictionary
            temp_dict = {key : seq}
            fasta_dict.update(temp_dict)
    
    handle.close()
    
    return(fasta_dict)
    print("DONE!")
    
def readRef(fasta_aln_path):
    from Bio import SeqIO
    
    #Read fasta into a dictionary
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_aln_path, 'fasta'))
    
    #Save annotatio to ref variable
    ref = fasta['reference_annotation'].description
    ref = str.split(str.split(ref, ' ')[1], ',')
    ref[len(ref)-3] = 'C1'
    ref[len(ref)-2] = 'C2'
    

    references = dict()
    for item in list(fasta.keys()):
        if item in ['reference_annotation', 'secondary_structure']:
            continue
        sequence = list(fasta[item].seq.upper())
        
        #Create a dictionary that every sprinzl is key and
        # counterpart sequences are values
        refDict = dict(zip(ref, sequence))
        #Remove sprinzl numbers does not contain any base
        refDict = {k: v for k, v in refDict.items() if v != '-'}
        
        #Save sprinzl series for each tRNAs
        references[item] = list(refDict.keys())
    
    return(references)

def clean_string(basestring):
    """Remove +N and -N patterns, and collect insertions."""
    
    #Check how many reads have insertions not how long inserted sequences
    insertions = len(re.findall(r'\+\d+', basestring))
    
    # Remove the +N markers
    basestring = re.sub(r'\+\d+', '', basestring)

    # Remove the -N markers
    basestring = re.sub(r'\-\d+', '', basestring)
    
    return(basestring, insertions)
    

def filter_bam(readList, condition_list, out_dir, len5, len3,
               margin, qmap, complete=False, fragment=False, firststrand=False):
    
    print('tRNA Counting Initiated!')
    
    #Create empty final dataframe for counts
    df = pd.DataFrame()
    
    #Create empty final dataframe for coverages
    cov_df = pd.DataFrame()
    
    #Iterate through samples
    for readName in list(readList.keys()):
        
        print(f'{readName}')
        
        read_path = readList[readName]
        
        condition = condition_list[readName]
        
        bamfile = pysam.AlignmentFile(read_path, "rb", threads = 8)
        
        outputBam = pysam.AlignmentFile(f'{out_dir}/{readName}_filtered.bam', "wb", template=bamfile)
        
        #Create a dictionary for every tRNAs are keys and their lengths are values
        refLengths = dict(zip(bamfile.references, bamfile.lengths))
        
        # Create dictionary of tRNAs which every values are 0
        counts = {key: 0 for key in list(refLengths.keys())}
        counts.update({'antisense': 0})
        counts.update({'not_unique': 0})
        
        refs = list(bamfile.references)
        refs = [k for k in refs if k not in ['RDN5', 'RDN58', 'oligo3', 'oligo5', 'oligo5_oligo3']]
        proportions = {key: [] for key in refs}

        for read in bamfile:
            if not three_prime_adapter_check(read, refLengths, margin=margin):
                continue
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
                    
                    #From pysam; https://pysam.readthedocs.io/en/latest/api.html#api
                    #reference_length = aligned length of the read on the reference genome.
                    
                    if read.reference_name in refs:
                        alignment_length = read.reference_length
                        alignment_proportion = alignment_length / refLengths[read.reference_name]
                    
                        proportions[read.reference_name].extend([alignment_proportion])
                    
            elif read.reference_name.startswith("oligo"):
                counts[read.reference_name] +=1
            else:
                counts['not_unique'] +=1
        
        bamfile.close()
        outputBam.close()
        
        print(f'{readName} read and processed! next file..')
        
        tempdata = pd.DataFrame.from_dict(counts, orient='index')
        tempdata.columns = [readName]
        if df.empty:
            df = tempdata
        else:
            df = pd.concat([df, tempdata], axis=1)
            
        temp_cov = pd.DataFrame([(k, v) for k, values in proportions.items() for v in values],
                          columns=["reference", "proportion"])
        temp_cov['sample'] = readName
        temp_cov['condition'] = condition
        
        cov_df = pd.concat([cov_df, temp_cov], axis=0)    
        
    return(df, cov_df)


def get_counts(readList, out_dir, fasta_path, references, len5, len3, qmap):
    
    # Get a dictionary from fasta
    fasta = parse_fasta(fasta_path)
    
    print('Single Nucleotide Level Counting Initiated!')

    #Create final dataframe 
    data = pd.DataFrame()
    for sample in list(readList.keys()):
        
        print(f'{sample}')
        
        inputDir = f'{out_dir}/{sample}_filtered.bam'
        
        pysam.index(inputDir)
        
        bamfile = pysam.AlignmentFile(inputDir, "rb")
        
        refLengths = dict(zip(bamfile.references, bamfile.lengths))

        #Skip counts for the listed items
        refList = [item for item in list(bamfile.references)
                   if item not in(
                           'RDN5', 'RDN58', 'oligo3',
                           'oligo5', 'oligo5_oligo3', 'antisense')]
        
        tempdata = pd.DataFrame()
        
        for item in refList:
            
            reference = []
            position = []
            sprinzl = [''] * len5 + references[item] + [''] * len3
            A = []
            G = []
            C = []
            T = []
            insertion = []
            deletion = []
        
            pileupLine = ["--ff", "3840", "--no-output-ends",
                          "--no-output-del", "--no-output-ins",
                          "-q", str(qmap), "-Q", "0", inputDir,
                          "-Br", item, "-f", fasta_path]
            
            pile = pysam.mpileup(*pileupLine)
            
            seq = list(fasta[item].upper())
            
            # samtools mpileup is not reporting if there is no count in a position.
            # This extra data frame and left join step is to ensure to report all positions
            # This part continues with pd.merge after the for loop
            tRNA_temp = pd.DataFrame({'reference': [item] * refLengths[item],
                                     'position': list(range(1,refLengths[item]+1)),
                                     'sprinzl_number': sprinzl,
                                     'base': seq})
            
            if pile == '':
                continue
            
            
            for line in pile.rstrip().split('\n'):
                lineData = line.split('\t')
                ref, pos, base = lineData[:3]
                base = base.upper()
                pos = int(pos)
                basestring = lineData[4]
                
                # Clean the base string and extract insertions
                cleaned_bases, insertions = clean_string(basestring)
                
                # samtools mpileup reports "." if position matches with fasta (when you feed fasta)
                # To count directly the base those positions replaced with actual base
                baseCount = Counter(cleaned_bases.replace('.',base))
               
                reference.extend([ref])
                position.extend([pos])
                A.extend([baseCount['A']])
                T.extend([baseCount['T']])
                G.extend([baseCount['G']])
                C.extend([baseCount['C']])
                deletion.extend([baseCount['*']])
                insertion.extend([insertions])
                
            
            tRNA_temp2 = pd.DataFrame({'reference':reference,
                                     'position':position,
                                     f'{sample} A':A,
                                     f'{sample} T':T,
                                     f'{sample} G':G,
                                     f'{sample} C':C,
                                     f'{sample} insertion':insertion,
                                     f'{sample} deletion':deletion,})
            
            tRNA_temp = pd.merge(tRNA_temp, tRNA_temp2,
                                 on =['reference', 'position'], how ='left')
            tRNA_temp = tRNA_temp.fillna(0)
            
            if tempdata.empty:
                tempdata = tRNA_temp
            else:
                tempdata = pd.concat([tempdata, tRNA_temp])
            
        if data.empty:
            data = tempdata
        else:
            data = pd.merge(
                data,
                tempdata,
                on =['reference', 'position', 'sprinzl_number', 'base'],
                how ='outer'
                )
        
        print(f'{sample} read and processed! next file..')
    
    data = data.fillna(0)
    return(data)

def get_errors(data, readList):
    print('Errors are calculating..')
    tempData = data[['reference', 'position', 'sprinzl_number', 'base']].copy()
    
    for index, columns in data.iterrows():
        base = columns['base']
        for sample in list(readList.keys()):
            total = sum(columns[[f'{sample} A', f'{sample} T',
                                 f'{sample} G', f'{sample} C',
                                 f'{sample} deletion', f'{sample} insertion']])
            if total == 0:
                error = 0
            else:
                correct = sum(columns[[f'{sample} {base}', f'{sample} insertion']])
                error = round((total-correct)/total, 4)
            tempData.loc[index, f'{sample} sum_err'] = error
    return(tempData)

def createSummary(data, demux):
    
    data = data.set_index('reference')
    
    data = data.transpose()
    
    sumReads = data.sum(axis=1)
    
    uniqueReads = sumReads - data['not_unique']
    uniqueReadsPerc = round((uniqueReads * 100) / sumReads, 2)
    
    oligo3Perc = round((data['oligo3'] * 100) / sumReads, 2)
    oligo5Perc = round((data['oligo5'] * 100) / sumReads, 2)
    antisensePerc = round((data['antisense'] * 100) / sumReads,2)
    
    data['aligned'] = sumReads
    
    data['unique'] = uniqueReads
    data['% unique'] = uniqueReadsPerc
    
    data['% antisense'] = antisensePerc
    data['% oligo3'] = oligo3Perc
    data['% oligo5'] = oligo5Perc
    
    data = data.reset_index()
    data = data[['index', 'aligned', 'unique', '% unique', 'antisense',
                 '% antisense', 'oligo3', '% oligo3', 'oligo5', '% oligo5']]
    data = data.rename(columns= {'index': 'Sample'})
    
    data = pd.merge(data, demux, on='Sample', how='left')
    
    data['% aligned'] = round(data['aligned'] * 100 / data['Demux Reads'], 2)
    
    data = data[['Sample', 'Barcode', 'Demux Reads', 'aligned', '% aligned',
                 'unique', '% unique', 'antisense', '% antisense',
                 'oligo3', '% oligo3', 'oligo5', '% oligo5']]
    
    return(data)
    
def getDemuxStats(demuxPath, readList):
    
    experiments = readList['Experiment'].unique()
    
    final_df = pd.DataFrame({"Sample": [],
                             "Experiment": [],
                             "Barcode": [],
                             "Demux Reads": []})
    
    for experiment in experiments:
        temp_read_list = readList[readList['Experiment'] == experiment].copy()
        
        temp_path = f'{demuxPath}/{experiment}/demux/pod5.demux.tsv.gz'

        with gzip.open(temp_path, 'rb') as f:
            file_content = pd.read_csv(f, sep='\t')
        
        file_content = file_content.groupby('barcode').count()['read_id']
        file_content = file_content.reset_index()
        file_content.columns = ['Barcode', 'Demux Reads']

        file_content = pd.merge(temp_read_list,
                                file_content, on = 'Barcode', how = 'left')
        
        final_df = pd.concat([final_df, file_content], axis=0)
    
    final_df = final_df.reset_index(drop=True)
    
    return(final_df)
        

    

def main():    
    import argparse
    
    parser = argparse.ArgumentParser(description='''Count tRNAs in total and in single nucleotide level''',
                                     formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-i', '--input_list', help='tab separated file that includes sample name and directory information of file')
    parser.add_argument('-o', '--output', help='directory for outputs')
    parser.add_argument('-f', '--fasta', help='directory for fasta file')
    parser.add_argument('-a', '--alignments', help='directory for fasta_aln file')
    parser.add_argument('-l5', '--length5', help="length of 5' adapter", type=int, default=24)
    parser.add_argument('-l3', '--length3', help="length of 3' adapter", type=int, default=30)
    parser.add_argument('-m', '--margin', help="margin for adapter length", type=int, default=0)
    parser.add_argument('-q', '--quality', help="mapq", type=int, default=10)
    parser.add_argument('-c', '--complete', help='count only full length tRNAs', action='store_true')
    parser.add_argument('-nc', '--not_complete', help='count only not full length tRNAs', action='store_true')
    parser.add_argument('-fs', '--firststrand', help="antisense strand sequencing ie cDNA [sense strand sequencing ie DRS]",
                        action='store_true')
    parser.add_argument('-x', '--demuxStat', help="Path to demux stat.tsv.gz file")
    
    arguments = parser.parse_args()  
    
    # Transfer arguments to variables
    input_list = arguments.input_list
    out_dir    = arguments.output
    fasta      = arguments.fasta
    alignments = arguments.alignments
    length5    = arguments.length5
    length3    = arguments.length3
    margin     = arguments.margin
    quality    = arguments.quality
    complete   = arguments.complete
    fragment   = arguments.not_complete
    firststrand= arguments.firststrand
    demux_stat = arguments.demuxStat
    
    
    sampleList = dict(zip(pd.read_csv(input_list, sep='\t')['Sample'],
               pd.read_csv(input_list, sep='\t')['Path']))
    
    conditionList = dict(zip(pd.read_csv(input_list, sep='\t')['Sample'],
               pd.read_csv(input_list, sep='\t')['Condition']))
    
    barcodeList = pd.read_csv(
        input_list, sep='\t'
        )[['Sample', 'Experiment', 'Barcode']].copy()
    
    references = readRef(alignments)

    tRNA_counts, coverage_proportions = filter_bam(
        readList=sampleList, condition_list=conditionList, out_dir=out_dir,
        len5=length5, len3=length3, margin=margin, qmap=quality,
        complete=complete, fragment=fragment,
        firststrand=firststrand
        )
    
    tRNA_single_nucleotide = get_counts(
        readList=sampleList, out_dir=out_dir, fasta_path=fasta,
        references=references, len5 = length5,
        len3 = length3, qmap=quality
        )
    
    tRNA_sumErr = get_errors(data=tRNA_single_nucleotide, readList=sampleList)
    
    tRNA_sumErr.to_csv(f'{out_dir}/sumErr.tsv', sep='\t', header=True, index=False)
    
    tRNA_counts.index.name = 'reference'
    tRNA_counts = tRNA_counts.reset_index()
    tRNA_counts.to_csv(f'{out_dir}/counts.tsv', sep='\t', header=True, index=False)
    
    coverage_proportions.to_csv(f'{out_dir}/coverage_proportions.tsv',
                                sep='\t', header=True, index=False)
    
    tRNA_single_nucleotide.to_csv(f'{out_dir}/counts_single_nucleotide.tsv',
                                  sep='\t', header=True, index=False)
    
    demuxStats = getDemuxStats(demuxPath=demux_stat, readList=barcodeList)

    summary = createSummary(data=tRNA_counts, demux=demuxStats)
    
    summary.to_csv(f'{out_dir}/Stats.tsv', sep='\t', header=True, index=False)
    
    print('nanoCount landed on Mars!..')
    
    if __name__ == '__main__':
        try:
            main()
        except KeyboardInterrupt:
            sys.stderr.write("\nMission Aborted!\nReturning to the Earth!     \n")
