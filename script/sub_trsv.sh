#!/bin/sh
mtx_dir="$(cd "$(dirname "$0")"/.. && pwd)/data"     # data dir 
full_time="$1"

input_file="../matrixList.txt"

output_file="../log/SpMP_TRSV_${full_time}.log"
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
    # yhrun ./test/ilu_test $matrix "csv/SpMP_${OMP_NUM_THREADS}_${time3}.csv" $OMP_NUM_THREADS >> $output_file
    ../test/sptrsv_test "$matrix_path" "../csv/SPMP_TRSV_${OMP_NUM_THREADS}_${full_time}.csv" $OMP_NUM_THREADS $turns >> $output_file
    # >> "$output_file"
done

echo $(date +"%y-%m-%d %H-%M-%S")

# 关闭文件描述符 3
exec 3<&-
