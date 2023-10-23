#! /usr/bin/env python

import pysam
import sys

def main():
    bedfile = sys.argv[1]
    reference = sys.argv[2]
    
    fasta = pysam.FastaFile(reference)
    with open(bedfile, 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            contig = items[0]
            position = int(items[1])
            
            start = max(0, position - 100)
            end = min(position + 101, fasta.get_reference_length(contig))
            seq = fasta.fetch(contig, start, end)
            print(line.rstrip('\n'), seq, sep="\t")
    
    fasta.close()        

if __name__ == "__main__":
    main()
