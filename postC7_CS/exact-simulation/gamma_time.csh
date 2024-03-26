#!/bin/bash

# submit job
path=$(pwd)              # output from pwd as string
name=${path##*/}         # retain the part after last /
seq_file=./$1.dat
n_file=$(wc -l ${seq_file})
IFS=' '     # space is set as delimiter
read -ra num <<< $n_file
num=${num[0]}

reps=$(( $3 / num ))

echo [*] name : $1
echo [*] steps : $2
echo [*] n_max : $3
echo [*] n_sample : $4
echo [*] MAS : $5
echo [*] n_file : $num
echo [*] reps : $reps
echo [*] CS : $6

CS=$6

mkdir -p /cluster/scratch/chavezm/${name}/${CS}

# sbatch --time=24:00:00 --mem-per-cpu=512 --job-name=CS_POST \
#        --output=/cluster/scratch/chavezm/${name}/gamma.log \
# 	   --error=/cluster/scratch/chavezm/${name}/gamma.err \
#        --wrap="./OPTIM_externshape_time CC.sys\
#      0 -4500 0 0 0 \
# 	 ${CS} 0 0 0 0 0 \
#      ${CS} 0 0 0 0 0 \
#      5 $2 $seq_file $3 $4 $5 0 /cluster/scratch/chavezm/${name}/$CS/R${name}${seq: -2}_S3_$reps $num"



