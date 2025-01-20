#!/usr/bin/env bash

aln=$1

clust=$(basename $aln .seq)
fam_name=$(basename $(dirname $aln))

struct_mode=$2

out_seq="./${fam_name}/${clust}${struct_id}_${struct_mode}.seq"
out_dca="./${fam_name}/${clust}_${struct_mode}.dca"
in_stat="./${fam_name}/in_${clust}.stat"
in_seq="./input/${fam_name}/${clust}.init"
out_stat="./${fam_name}/${clust}_${struct_id}_${struct_mode}.stat"

ref_seq=$(cat $in_seq)

target_id=$3
struct="input/${fam_name}/${target_id}.str"

if [ -n "$4" ]; then
    alpha="-w $4"
    out_seq="${out_seq}_w$4"
else
    alpha=""
fi

if [ -n "$5" ]; then
    phi="-p $5"
    out_seq="${out_seq}_p$5"
else
    phi="-p 1"
fi

if [[ $target_id == "subopt" ]]; then
    time ../mc_sec_c/bin/dca_prof -r "$ref_seq" -n 100 -u 10000 -t 1 -m 1 -s ${struct_mode} -f 1E-3 \
         $phi $alpha -b ${out_dca} $in_seq $struct | grep -v "#" | awk '{print $3}' > ${out_seq}_T${target_id}
else
    time ../mc_sec_c/bin/dca_prof -r "$ref_seq" -n 1 -u 1000000 -t 1 -m 1 -s ${struct_mode} -f 1E-3 \
         $phi $alpha -b ${out_dca} $in_seq $struct | grep -v "#" | awk '{print $3}' > ${out_seq}_T${target_id}
fi
