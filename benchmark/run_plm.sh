#!/usr/bin/env bash

aln=$1
fam_name=$(basename $aln .seq)
out_dir="out"
out_coup=${out_dir}/${fam_name}.couplings
out_parms=${out_dir}/${fam_name}.parms

/home/vopuu/programs/plmc/bin/plmc -c $out_coup -o $out_parms -a AUGC- $aln
