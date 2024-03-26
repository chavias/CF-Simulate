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
#                            second-order jobs
#
##########################################################################################

# maybe 1,2,3 repetitions would be best

cores=50
run_time=120:00:00
memory=2048

# jobname=R26O2r1
# rep=1
# sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
#        --output=${directory}$jobname.log \
#        --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_MAS_powder_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=R26O2r2
rep=2   #2
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_MAS_powder_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\""

jobname=R26O2r3
rep=4    #3
sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
       --output=${directory}$jobname.log \
       --wrap="matlab -nodisplay -singleCompThread -r \"O2=generate_R26_MAS_powder_O2(${rep},${cores}); save('${directory}${jobname}.mat','O2')\"" 


##########################################################################################
#
#                            first-order jobs
#
##########################################################################################

# maybe 1,2,3 repetitions would be best

# cores=4
# run_time=04:00:00
# memory=2048

# jobname=R26O1r1
# rep=1
# sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
#        --output=${directory}$jobname.log \
#        --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_MAS_powder_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

# jobname=R26O1t2
# rep=6     # 2
# sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
#        --output=${directory}$jobname.log \
#        --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_MAS_powder_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\""

# jobname=R26O1t3
# rep=10   # 3
# sbatch --ntasks=1 --cpus-per-task=$cores --time=$run_time --mem-per-cpu=$memory --job-name=$jobname \
#        --output=${directory}$jobname.log \
#        --wrap="matlab -nodisplay -singleCompThread -r \"O1=generate_R26_MAS_powder_O1(${rep},${cores}); save('${directory}${jobname}.mat','O1')\"" 
