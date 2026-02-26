# ./make.sh
# full_time=$(date "+%Y-%m-%d %H-%M-%S")
# # yhbatch --output="./out/ilu_%x_${full_time}_%j.out" ./sub.sh "${full_time}"
# yhbatch --output="./out/ic_%x_${full_time}_%j.out" ./sub_ic.sh "${full_time}"


full_time=$(date "+%Y-%m-%d-%H-%M-%S")
nohup ./sub_ic.sh "${full_time}" > "../out/ic_SpMP_IC_${full_time}.out" 2>&1 &
