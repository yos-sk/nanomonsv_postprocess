# Nanomonsv_postprocess
This tool removes redundant SV calls that arise when running nanomonsv with a diploid genome assembly, by realigning supporting reads around breakpoints and grouping calls that represent the same event on the two haplotypes.

It is primarily intended to run inside [PRCGAP](https://github.com/yos-sk/PRCGAP), where it is invoked as one step of the SV calling pipeline.

The input is the nanomonsv result file; see the [nanomonsv output format](https://github.com/friend1ws/nanomonsv#5-output-file-format) for column definitions.

## How to install
```
git clone https://github.com/yos-sk/nanomonsv_postprocess.git
cd nanomonsv_postprocess
cargo build --release
```

## Usage

1. Convert nanomonsv results to a bed file.
```
awk '{if (NR != 1) print $1 "\t" $2 - 1 "\t" $2 "\t" $8 "\n" $4 "\t" $5 - 1 "\t" $5 "\t" $8}' \
    nanomonsv.result.txt > nanomonsv_result.bed
```

2. Collect sequences around SV breakpoints with the `extract-seq` subcommand.
```
samtools faidx reference.fa
nanomonsv_postprocess extract-seq \
    -b nanomonsv_result.bed \
    -f reference.fa > nanomonsv_result_with_seq.txt
```

3. Realign supporting reads and group SVs that share the same breakpoint sequence.
```
nanomonsv_postprocess realignment \
    -i nanomonsv_result_with_seq.txt \
    -s nanomonsv.supporting_read.txt \
    -b bam_file \
    -d 98 \ # realignment identity threshold
    -l 180 \ # realignment length threshold
    1>SV_group_info.txt 2>SV_pair_info.txt
```

4. Filter the nanomonsv result using the grouping information.
```
nanomonsv_postprocess filt \
    -i SV_group_info.txt \
    -n nanomonsv.result.txt \
    -s nanomonsv.supporting_read.txt \
    -b bam_file \
    1>nanomonsv_new_results.txt 2>nanomonsv_postprocess_filter.log
```
