#! /usr/bin/env python

import pysam
import sys

def calc(input_log: str) -> None:
    
    haplotyped = 0
    unassigned = 0
    s_unassigned = 0
    b_mixed = 0
    s_mixed = 0
    with open(input_log, 'r') as f:
        for line in f:
            if line[0] not in ['i', 'd', 'r']: continue
            
            items = line.rstrip('\n').split('\t')
            if items[3] == "identical": continue
            bp1_hap1 = int(items[1].split(',')[1])
            bp1_hap2 = int(items[1].split(',')[2])
            
            if bp1_hap1 > 0 and bp1_hap2 > 0:
                bp1_cat = "mixed"
            elif bp1_hap1 > 0 or bp1_hap2 > 0:
                bp1_cat = "haplotyped"
            else:
                bp1_cat = "unassigned"
            
            
            bp2_hap1 = int(items[2].split(',')[1])
            bp2_hap2 = int(items[2].split(',')[2])
            
            if bp2_hap1 > 0 and bp2_hap2 > 0:
                bp2_cat = "mixed"
            elif bp2_hap1 > 0 or bp2_hap2 > 0:
                bp2_cat = "haplotyped"
            else:
                bp2_cat = "unassigned"
            
            if bp1_cat == "haplotyped" and bp2_cat == "haplotyped":
                haplotyped += 1
            elif bp1_cat == "mixed" and bp2_cat == "mixed":
                b_mixed += 1
            elif bp1_cat == "mixed" or bp2_cat == "mixed":
                s_mixed += 1
            elif bp1_cat == "unassigned" and bp2_cat == "unassigned":
                unassigned += 1
            else:
                s_unassigned += 1 
    
    print("Both breakpoint assigned    :", haplotyped, sep="\t")
    print("Both breakpoint mixed       :", b_mixed, sep="\t")
    print("Single breakpoint mixed     :", s_mixed, sep="\t")
    print("Both breakpoint unassigned  :", unassigned, sep="\t")
    print("Single breakpoint unassigned:", s_unassigned, sep="\t")

def main():
    input = sys.argv[1]
    calc(input)

if __name__ == "__main__":
    main()
            
            
                
            