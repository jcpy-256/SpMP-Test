# ./make.sh
# full_time=$(date "+%Y-%m-%d %H-%M-%S")
# # yhbatch --output="./out/ilu_%x_${full_time}_%j.out" ./sub.sh "${full_time}"
# yhbatch --output="./out/trsv_%x_${full_time}_%j.out" ./sub_trsv.sh "${full_time}"

full_time=$(date "+%Y-%m-%d-%H-%M-%S")
nohup ./sub_trsv.sh "${full_time}" > "../out/trsv_SpMP_TRSV_${full_time}.out" 2>&1 &
