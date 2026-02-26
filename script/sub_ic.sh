#!/bin/sh
mtx_dir="$(cd "$(dirname "$0")"/.. && pwd)/data"     # data dir 
full_time="$1"

input_file="../matrixList.txt"

output_file="../log/SpMP_IC_${full_time}.log"
# 清空result.txt文件
> "$output_file"
echo $full_time


nthread=2
turns=100
export OMP_NUM_THREADS=$nthread
export OMP_PROC_BIND=close
export OMP_PLACES=cores
exec 3< "${input_file}"

while IFS= read -r matrix <&3; do
    # 构建矩阵的完整路径
    matrix_path="$mtx_dir/$matrix.mtx"
    echo $matrix_path
    if [ ! -f "$matrix_path" ]; then
        echo "Warning: Matrix file $matrix_path does not exist, skipping..." >> "$output_file"
        continue
    fi
    ../test/ic_test $matrix_path "../csv/SPMP_IC_${OMP_NUM_THREADS}_${full_time}.csv" $OMP_NUM_THREADS $turns >> $output_file
done


echo $(date +"%y-%m-%d %H-%M-%S")
exec 3<&-
