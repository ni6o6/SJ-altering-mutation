#!/bin/bash

#######
#sh ./juncmut_bash/run_juncmut.sh <folder name under ./junction/>
#######

fo=$1

for pass in `\find ./junction/${fo} -name '*SJ.out.tab'`; do
IFS='/' read -r -a file <<< "${pass}"
IFS='.' read -r -a array <<< "${file[3]}"
pr=${array[0]}
echo ${pr}

./juncmut_bash/juncmut_env.py ${fo} 

./juncmut_bash/juncmut_juncutils.py ${pr} ${fo} --control_file ./reference/SJ_control_2_4.bed.gz

./juncmut_bash/juncmut_assadj.py ${pr} ${fo}

./juncmut_bash/juncmut_freq.py ${pr} ${fo} --read_num_thres 3 --freq_thres 0.05

./juncmut_bash/juncmut_spliceai.py ${pr} ${fo}

./juncmut_bash/juncmut_extractcanomut.py ${pr} ${fo}

./juncmut_bash/juncmut_annotgnomadsnp.py ${pr} ${fo}

done 
