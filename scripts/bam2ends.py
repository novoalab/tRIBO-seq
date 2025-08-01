#!/usr/bin/env python
desc="""Report number of 5' read ends for all transcripts. 
"""
epilog="""Author: lpryszcz@crg.es
Torredembarra, 30/04/2025
"""

import os, sys
import pysam
import numpy as np
from datetime import datetime

def bam2ends(out, bams, mapq=10, skip_flag=3856, threads=1, verbose=False):
    """Report number of  5' read ends for all transcripts across BAM file(s)"""
    samples = [os.path.basename(bam) for bam in bams]
    # open alignments
    sams = [pysam.AlignmentFile(bam) for bam in bams]
    # get refs
    sam = sams[0]
    ref2len = {r: l for r, l in zip(sam.references, sam.lengths)}
    # filter algs
    read_callback = lambda a: not a.flag&skip_flag and a.mapq>=mapq
    # write header
    out.write("locus\t%s\n"%"\t".join(f"{s} no_ends\t{s} ends" for s in samples))
    # process refs
    for ref, rlen in ref2len.items():
        ends = np.zeros((len(bams), 2, rlen), dtype='int')
        for si, sam in enumerate(sams):
            # get cov
            cc = sam.count_coverage(ref, quality_threshold=0, read_callback=read_callback)
            cov = np.array(cc).sum(axis=0)
            #ends[si, 0] = cov
            # get number of reads with 5' ends for every position
            for a in sam.fetch(ref):
                if not read_callback(a): continue
                # cound ends
                ends[si, 1, a.pos] += 1
            # count number of reads that don't end there
            ends[si, 0] = cov - ends[si, 1]
        # vstact & report
        ends = np.vstack(ends)
        for i in range(rlen):
            out.write("%s:%s\t%s\n"%(ref, i+1, "\t".join(map(str, ends[:, i]))))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version="1.0a")
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"), 
                        help="output stream [stdout]")
    parser.add_argument("-i", "--input", nargs="+", default=[],
                        help="input BAM files(s)")
    parser.add_argument("-m", "--mapq", default=10, type=int,
                        help="min. mapping quality [%(default)s]")
    parser.add_argument("-F", "--skip_flag", default=3856, type=int,
                        help="SAM flag to skip [%(default)s]")
    parser.add_argument("-t", "--threads", default=1, type=int,
                        help="number of cores to use [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose: 
        sys.stderr.write("Options: %s\n"%str(o))

    bam2ends(o.output, o.input, o.mapq, o.skip_flag, o.threads, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
  
