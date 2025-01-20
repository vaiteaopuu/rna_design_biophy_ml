#!/usr/bin/env bash

# USAGE: bash ./<this script>.sh <MSA> <STRUCT_MODE>
# STRUCT_MODE = 0 if no structure and 3 if structure
aln=$1

fam_name=$(basename $aln .seq)
struct_mode=$2
out_dir="out"
struct="dca_input/${fam_name}.str"

# Ajust these path for the benchmarks
# out_seq="${out_dir}/${fam_name}_${struct_mode}.seq"
# out_cl="${out_dir}/${fam_name}_${struct_mode}_only.seq"
out_dca="${out_dir}/${fam_name}_${struct_mode}.dca"

mkdir -p $out_dir

# Train the model
time ../mc_sec_c/bin/dca  -n 3000 -u 1000 -t 1 -l 0.05 -h 0.1 -m 0 -s ${struct_mode}\
     $aln $struct > ${out_dca} 2> ${out_dir}/${fam_name}.log

# Sample from the model
# time ../mc_sec_c/bin/dca  -n 1 -u 1000000 -t 1 -m 1 -s ${struct_mode} -f 1E-3 \
    #      $aln $struct -b ${out_dca} > $out_seq
