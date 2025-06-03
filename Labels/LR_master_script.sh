#!/bin/bash

source /sci/labs/asafle/ofirse/icore-data/anaconda_files/etc/profile.d/conda.sh
conda activate /sci/labs/asafle/alexlevylab/icore-data/miniconda3/envs/lusi_ML_env

#run merge_unknown_manually_with_all_data.py
#run edit_manually_labeled_df.py

mkdir slurms

#run HvN (run,decide,merge)
#run model
job_ID1=($(sbatch --wrap="python3 run_logistic_regression.py HvN" -c15 --time=5:0:0 -o slurms/run_logistic_regression_HvN_%J.out | grep -oE "[[:digit:]]+"))

#decide labels from score
job_ID2=($(sbatch --dependency=afterok:$(echo $job_ID1) --wrap="python3 decide_labels_LR.py HvN" -o slurms/slurm_decide_labels_HvN_%J.out | grep -oE "[[:digit:]]+"))

#merge results into one table
job_ID3=($(sbatch --dependency=afterok:$(echo $job_ID2) --wrap="python3 merge_prediction_results.py HvN" -o slurms/merge_prediction_results_HvN_%J.out | grep -oE "[[:digit:]]+"))

#run AvP - with all data that is not predicted as -1 in host
#run model
job_ID4=($(sbatch  --dependency=afterok:$(echo $job_ID3) --wrap="python3 run_logistic_regression.py AvP" -c15 --time=5:0:0 -o slurms/run_logistic_regression_AvP_%J.out | grep -oE "[[:digit:]]+"))

#decide labels from score
job_ID5=($(sbatch --dependency=afterok:$(echo $job_ID4) --wrap="python3 decide_labels_LR.py AvP" -o slurms/slurm_decide_labels_AvP_%J.out | grep -oE "[[:digit:]]+"))

#merge results into one table
job_ID6=($(sbatch --dependency=afterok:$(echo $job_ID5) --wrap="python3 merge_prediction_results.py AvP" -o slurms/merge_prediction_results_AvP_%J.out | grep -oE "[[:digit:]]+"))

######################

#run on labels gut,oral,skin
job_ID7=($(sbatch --dependency=afterok:$(echo $job_ID6) --wrap="python3 run_logistic_regression.py GvOvS" -c15 --time=5:0:0 -o slurms/run_logistic_regression_GvOvS_%J.out | grep -oE "[[:digit:]]+"))

job_ID8=($(sbatch --dependency=afterok:$(echo $job_ID7) --wrap="python3 decide_labels_LR.py GvOvS" -o slurms/slurm_decide_labels_GvOvS_%J.out | grep -oE "[[:digit:]]+"))

job_ID9=($(sbatch --dependency=afterok:$(echo $job_ID8) --wrap="python3 merge_prediction_results.py GvOvS" -o slurms/merge_prediction_results_GvOvS_%J.out | grep -oE "[[:digit:]]+"))

#run on labels root/shoot
job_ID10=($(sbatch --dependency=afterok:$(echo $job_ID9) --wrap="python3 run_logistic_regression.py RvS" -c15 --time=5:0:0 -o slurms/run_logistic_regression_RvS_%J.out | grep -oE "[[:digit:]]+"))

job_ID11=($(sbatch --dependency=afterok:$(echo $job_ID10) --wrap="python3 decide_labels_LR.py RvS" -o slurms/slurm_decide_labels_RvS_%J.out | grep -oE "[[:digit:]]+"))

job_ID12=($(sbatch --dependency=afterok:$(echo $job_ID11) --wrap="python3 merge_prediction_results.py RvS" -o slurms/merge_prediction_results_RvS_%J.out | grep -oE "[[:digit:]]+"))

#run on labels teresterial aquatic
#job_ID13=($(sbatch --dependency=afterok:$(echo $job_ID12) --wrap="python3 run_logistic_regression.py TvA" -c15 --time=5:0:0 -o slurms/run_logistic_regression_TvA_%J.out | grep -oE "[[:digit:]]+"))

#job_ID14=($(sbatch --dependency=afterok:$(echo $job_ID13) --wrap="python3 decide_labels_LR.py TvA" -o slurms/slurm_decide_labels_TvA_%J.out | grep -oE "[[:digit:]]+"))

#job_ID15=($(sbatch --dependency=afterok:$(echo $job_ID14) --wrap="python3 merge_prediction_results.py TvA" -o slurms/merge_prediction_results_TvA_%J.out | grep -oE "[[:digit:]]+"))

#run on human mouse 
#changed jobID dependency
job_ID16=($(sbatch --dependency=afterok:$(echo $job_ID12) --wrap="python3 run_logistic_regression.py HvM" -c15 --time=5:0:0 -o slurms/run_logistic_regression_HvM_%J.out | grep -oE "[[:digit:]]+"))

job_ID17=($(sbatch --dependency=afterok:$(echo $job_ID16) --wrap="python3 decide_labels_LR.py HvM" -o slurms/slurm_decide_labels_HvM_%J.out | grep -oE "[[:digit:]]+"))

job_ID18=($(sbatch --dependency=afterok:$(echo $job_ID17) --wrap="python3 merge_prediction_results.py HvM" -o slurms/merge_prediction_results_HvM_%J.out | grep -oE "[[:digit:]]+"))

#run on arabidobsis poplar
job_ID19=($(sbatch --dependency=afterok:$(echo $job_ID18) --wrap="python3 run_logistic_regression.py ArvPo" -c15 --time=5:0:0 -o slurms/run_logistic_regression_ArvPo_%J.out | grep -oE "[[:digit:]]+"))

job_ID20=($(sbatch --dependency=afterok:$(echo $job_ID19) --wrap="python3 decide_labels_LR.py ArvPo" -o slurms/slurm_decide_labels_ArvPo_%J.out | grep -oE "[[:digit:]]+"))

job_ID21=($(sbatch --dependency=afterok:$(echo $job_ID20) --wrap="python3 merge_prediction_results.py ArvPo" -o slurms/merge_prediction_results_ArvPo_%J.out | grep -oE "[[:digit:]]+"))

#run on pathogen commensal
#job_ID22=($(sbatch --dependency=afterok:$(echo $job_ID21) --wrap="python3 run_logistic_regression.py PvC" -c15 --time=5:0:0 -o slurms/run_logistic_regression_PvC_%J.out | grep -oE "[[:digit:]]+"))

#job_ID23=($(sbatch --dependency=afterok:$(echo $job_ID22) --wrap="python3 decide_labels_LR.py PvC" -o slurms/slurm_decide_labels_PvC_%J.out | grep -oE "[[:digit:]]+"))

#job_ID24=($(sbatch --dependency=afterok:$(echo $job_ID23) --wrap="python3 merge_prediction_results.py PvC" -o slurms/merge_prediction_results_PvC_%J.out | grep -oE "[[:digit:]]+"))

##################################

#run on invertebrates
#job_ID25=($(sbatch --dependency=afterok:$(echo $job_ID24) --wrap="python3 run_logistic_regression.py Inv" -c15 --time=5:0:0 -o slurms/run_logistic_regression_Inv_%J.out | grep -oE "[[:digit:]]+"))

#job_ID26=($(sbatch --dependency=afterok:$(echo $job_ID25) --wrap="python3 decide_labels_LR.py Inv" -o slurms/slurm_decide_labels_Inv_%J.out | grep -oE "[[:digit:]]+"))

#job_ID27=($(sbatch --dependency=afterok:$(echo $job_ID26) --wrap="python3 merge_prediction_results.py Inv" -o slurms/merge_prediction_results_Inv_%J.out | grep -oE "[[:digit:]]+"))

#run on artropods
#job_ID28=($(sbatch --dependency=afterok:$(echo $job_ID27) --wrap="python3 run_logistic_regression.py Art" -c15 --time=5:0:0 -o slurms/run_logistic_regression_Art_%J.out | grep -oE "[[:digit:]]+"))

#job_ID29=($(sbatch --dependency=afterok:$(echo $job_ID28) --wrap="python3 decide_labels_LR.py Art" -o slurms/slurm_decide_labels_Art_%J.out | grep -oE "[[:digit:]]+"))

#job_ID30=($(sbatch --dependency=afterok:$(echo $job_ID29) --wrap="python3 merge_prediction_results.py Art" -o slurms/merge_prediction_results_Art_%J.out | grep -oE "[[:digit:]]+"))

#run on insects
#job_ID31=($(sbatch --dependency=afterok:$(echo $job_ID30) --wrap="python3 run_logistic_regression.py Ins" -c15 --time=5:0:0 -o slurms/run_logistic_regression_Ins_%J.out | grep -oE "[[:digit:]]+"))

#job_ID32=($(sbatch --dependency=afterok:$(echo $job_ID31) --wrap="python3 decide_labels_LR.py Ins" -o slurms/slurm_decide_labels_Ins_%J.out | grep -oE "[[:digit:]]+"))

#job_ID33=($(sbatch --dependency=afterok:$(echo $job_ID32) --wrap="python3 merge_prediction_results.py Ins" -o slurms/merge_prediction_results_Ins_%J.out | grep -oE "[[:digit:]]+"))

#run on hemolymph
#job_ID34=($(sbatch --dependency=afterok:$(echo $job_ID33) --wrap="python3 run_logistic_regression.py Hemo" -c15 --time=5:0:0 -o slurms/run_logistic_regression_Hemo_%J.out | grep -oE "[[:digit:]]+"))

#job_ID35=($(sbatch --dependency=afterok:$(echo $job_ID34) --wrap="python3 decide_labels_LR.py Hemo" -o slurms/slurm_decide_labels_Hemo_%J.out | grep -oE "[[:digit:]]+"))

#job_ID36=($(sbatch --dependency=afterok:$(echo $job_ID35) --wrap="python3 merge_prediction_results.py Hemo" -o slurms/merge_prediction_results_Hemo_%J.out | grep -oE "[[:digit:]]+"))
