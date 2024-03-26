#!/bin/bash

# load module
. /cluster/apps/local/env2lmod.sh
module load matlab/R2023b

path=$(pwd)       # output from pwd as string
name=${path##*/}  # retain the part after last /
directory=/cluster/scratch/chavezm/$name/
mkdir -p ${directory}

##########################################################################################
#
#                            first-order jobs
#
##########################################################################################
cores=4
run_time=24:00:00
memory=4096

jobname=O1t1
rep=2
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_C7_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

jobname=O1t2
rep=6
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_C7_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

jobname=O1t3
rep=10
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_C7_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\"" 

##########################################################################################
#
#                            second-order jobs
#
##########################################################################################

cores=30
run_time=120:00:00
memory=4096

jobname=C7O2t1
rep=2
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_C7_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=C7O2t2
rep=6
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_C7_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=C7O2t3
rep=10
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_C7_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\"" 
