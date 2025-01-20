faml=(RF00005_1 RF01852 RF00167 RF00168 RF00234 RF00380 RF01051 RF01734 RF01750 RF01763 RF01786 RF01831_1 RF01854 RF00028_1 RF01767 RF02682 RF00379 RF02683 RF00442_1 RF00080 RF00011_1 RF00050 RF00504 RF01689 RF01725 RF02001_2)

# # normal sampling
# for fam in ${faml[@]}; do
#     for struct_mode in {0,3}; do
#         for clust in {0,1}; do
#             nohup bash run_sample.sh input/${fam}/clust_${clust}.seq $struct_mode clust_$clust &
#         done
#     done
# done

# # higher temperature
# struct_mode=0
# for fam in ${faml[@]}; do
#     for clust in {0,1}; do
#         nohup bash run_sample.sh input/${fam}/clust_${clust}.seq $struct_mode clust_$clust 1.3 &
#     done
# done

# # alpha mut with regular DCA
# struct_mode=0
# for fam in ${faml[@]}; do
#     for clust in {0,1}; do
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode clust_$clust -0.2 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode clust_$clust -0.5 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode clust_$clust -1 &
#     done
# done

# # alpha mut opposite target
# struct_mode=3
# for fam in ${faml[@]}; do
#     for clust in {0,1}; do
#         nohup bash run_sample_prof.sh input/${fam}/clust_$clust.seq $struct_mode clust_$clust -0.2 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_$clust.seq $struct_mode clust_$clust -0.5 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_$clust.seq $struct_mode clust_$clust -1 &
#     done
# done

# # alpha STRUCT opposite target
# struct_mode=3
# for fam in ${faml[@]}; do
#     nohup bash run_sample.sh input/${fam}/clust_0.seq $struct_mode clust_1 1 1 2 &
#     nohup bash run_sample.sh input/${fam}/clust_0.seq $struct_mode clust_1 1 1 3 &
#     nohup bash run_sample.sh input/${fam}/clust_0.seq $struct_mode clust_1 1 1 4 &
#     nohup bash run_sample.sh input/${fam}/clust_1.seq $struct_mode clust_0 1 1 2 &
#     nohup bash run_sample.sh input/${fam}/clust_1.seq $struct_mode clust_0 1 1 3 &
#     nohup bash run_sample.sh input/${fam}/clust_1.seq $struct_mode clust_0 1 1 4 &
# done

# # alpha STRUCT opposite target
# struct_mode=3
# for fam in ${faml[@]}; do
#     for clust in {0,1}; do
#         nohup bash run_sample.sh input/${fam}/clust_${clust}.seq $struct_mode clust_${clust} 1 1 2 &
#         nohup bash run_sample.sh input/${fam}/clust_${clust}.seq $struct_mode clust_${clust} 1 1 3 &
#     done
# done

# # alpha opposite target
# struct_mode=3
# for fam in ${faml[@]}; do
#     nohup bash run_sample_prof.sh input/${fam}/clust_0.seq $struct_mode clust_1 -0.2 &
#     nohup bash run_sample_prof.sh input/${fam}/clust_0.seq $struct_mode clust_1 -0.5 &
#     nohup bash run_sample_prof.sh input/${fam}/clust_0.seq $struct_mode clust_1 -1 &
#     nohup bash run_sample_prof.sh input/${fam}/clust_1.seq $struct_mode clust_0 -0.2 &
#     nohup bash run_sample_prof.sh input/${fam}/clust_1.seq $struct_mode clust_0 -0.5 &
#     nohup bash run_sample_prof.sh input/${fam}/clust_1.seq $struct_mode clust_0 -1 &
# done

# # supoptimal
# struct_mode=3
# for fam in ${faml[@]}; do
#     for clust in {0,1}; do
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode subopt -0.2 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode subopt -0.5 &
#         nohup bash run_sample_prof.sh input/${fam}/clust_${clust}.seq $struct_mode subopt -1 &
#     done
# done
