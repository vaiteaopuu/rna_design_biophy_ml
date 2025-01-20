#!/usr/bin/env bash

# USAGE: bash <this script> <MSA> <STRUCTURE> <INITIAL_SEQ> <OUT_DIR> <STRUCTURE_MODE> <PHI> <TEMPERATURE>
aln=$1
struct=$2


init_seq=$3

out_dir=$4

struct_mode=$5

phi=$6
temp=$7

fam_name=$(basename $aln .seq)

out_seq=${out_dir}/out_${struct_mode}.seq
out_cl=${out_dir}/out_cl_${fam_name}_${struct_mode}.seq
out_dca=${out_dir}/sys_${fam_name}_${struct_mode}.dca
out_stat=${out_dir}/out_${struct_mode}.stat
in_stat=${out_dir}/in.stat

ref_seq=$(cat ${init_seq})

# mkdir -p ${out_dir}

# running the training
time ../mc_sec_c/bin/dca_prof  -n 3000 -u 1000 -t 1 -l 0.01 -h 0.1 -m 0 -r ${ref_seq} -s ${struct_mode}\
     $aln $struct > ${out_dca}


# # One trajectory
# time ../mc_sec_c/bin/dca_prof  -n 1 -u 1000000 -t 1 -m 1 -r ${ref_seq} -s ${struct_mode} -f 1E-3 \
#          -b ${out_dca} $aln $struct | grep -v "#" | awk '{print $3}' > ${out_cl}

# # Run multiple trajectories
# for wei in {-1..4..1}; do
#     # time ../mc_sec_c/bin/dca_prof  -n 1 -u 1000000 -t 1 -m 1 -r ${ref_seq} -s ${struct_mode} -f 1E-3 \
#         #      -w $wei -b ${out_dca} $aln $struct | grep -v "#" | awk '{print $3}' > ${out_cl}_w${wei}_t${temp}_p$alpha
#     nohup time ../mc_sec_c/bin/dca_prof  -n 1 -u 2000000 -t ${temp} -m 1 -r ${ref_seq} -s ${struct_mode} -f 1E-3 -w $wei -p $phi $init_seq $struct -b ${out_dca} | grep -v "#" | awk '{print $3}' > ${out_cl}_w${wei}_p${phi}_t${temp}_23 &
# done
