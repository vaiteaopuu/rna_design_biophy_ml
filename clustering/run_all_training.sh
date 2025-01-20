faml=(RF00005_1 RF01852 RF00167 RF00168 RF00234 RF00380 RF01051 RF01734 RF01750 RF01763 RF01786 RF01831_1 RF01854 RF00028_1 RF01767 RF02682 RF00379 RF02683 RF00442_1 RF00080 RF00011_1 RF00050 RF00504 RF01689 RF01725 RF02001_2)

# TRAINING ####################################################################
for fam in ${faml[@]}; do
    for clust in {0,1}; do
        struct_mode=0
        nohup bash run_training.sh input/${fam}/clust_${clust}.seq $clust $struct_mode &
    done
done

for fam in ${faml[@]}; do
    for clust in {0,1}; do
        struct_mode=3
        nohup bash run_training.sh input/${fam}/clust_${clust}.seq $clust $struct_mode &
    done
done
