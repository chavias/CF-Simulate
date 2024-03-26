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
#                            CHECK THE REPETITIONS!!!!
#
#                !!! scheme_R26_nocomp is half as long as scheme_R26  !!!
#
##########################################################################################



##########################################################################################
#
#                            first-order jobs
#
##########################################################################################
cores=4
run_time=04:00:00
memory=4096

jobname=N26O1t1CS
rep=2
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_nocomp_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

jobname=N26O1t2CS
rep=4
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_nocomp_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

jobname=N26O1t3CS
rep=12
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_nocomp_CS_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\"" 

# ##########################################################################################
# #
# #                            second-order jobs
# #
# ##########################################################################################

cores=20
run_time=120:00:00
memory=2048

jobname=N26O2t1CS
rep=2
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_nocomp_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=N26O2t2CS
rep=4
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_nocomp_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=N26O2t3CS
rep=12
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_nocomp_CS_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\"" 
