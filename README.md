# Nanomonsv_postprocess
[Nanomonsv](https://github.com/friend1ws/nanomonsv.git) postprocess filtering and haplotyping.

## How to install
```
git clone https://github.com/yos-sk/nanomonsv_postprocess.git
cd nanomonsv_postprocess
cargo build --release
```

## Preparation

1. Convert nanomonsv results to bed file.
```
awk '{if (NR != 1) print $1 "\t" $2 - 1 "\t" $2 "\t" $8 "\n" $4 "\t" $5 - 1 "\t" $5 "\t" $8 \
    nanomonsv.result.txt > nanomonsv_result.bed
```

2. Collect sequences around SV breakpoints \
[pysam](https://github.com/pysam-developers/pysam.git) library is needed.
```
samtools faidx reference.fa
python ./script/extract_seq.py \
    H2009_deepconsensus_nanomonsv_result.bed \
    reference.fa > nanomonsv_result_with_seq.txt
```
## Usage

```
nanomonsv_postprocess realignment \
    -i nanomonsv_result_with_seq.txt \
    -s nanomonsv.supporting_read.txt \
    -b bam_file \
    -d 98 \ # realignment identity threshold
    -l 180 \ # realignmment length threshold
    1>SV_group_info.txt 2>SV_pair_info.txt

nanomonsv_postprocess filt \
    -i SV_group_info.txt \
    -n nanomonsv.result.txt \
    -s nanomonsv.supporting_read.txt \
    -b bam_file \
    1>nanomonsv_new_results.txt 2>nanomonsv_postprocess_filter.log
```