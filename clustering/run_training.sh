#!/usr/bin/env bash

aln=$1

clust=$(basename $aln .seq)
fam_name=$(basename $(dirname $aln))
struct_id=$2
struct="input/${fam_name}/${clust}.str"
target_struct="input/${fam_name}/clust_${struct_id}.str"

struct_mode=$3

out_seq="./${fam_name}/${clust}_${struct_id}_${struct_mode}.seq"
out_dca="./${fam_name}/${clust}_${struct_mode}.dca"
in_stat="./${fam_name}/in_${clust}.stat"
out_stat="./${fam_name}/${clust}_${struct_id}_${struct_mode}.stat"


if [ -n "$4" ]; then
    alpha="-a $4"
    out_seq="${out_seq}_a$4"
else
    alpha=""
fi

if [ -n "$5" ]; then
    phi="-p $5"
    out_seq="${out_seq}_p$5"
else
    phi="-p 1"
fi

mkdir -p ${fam_name}

time ../mc_sec_c/bin/dca -n 3000 -u 1000 -t 1 -l 0.01 -h 0.1 -m 0 -s ${struct_mode}\
     $aln $struct > ${out_dca} 2> ./${fam_name}/${clust}_${struct_mode}.log
