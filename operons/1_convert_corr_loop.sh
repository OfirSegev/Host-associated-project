for i in scaffold_correlation_matrixs/scaf_correlation_matrix_clade_*
#for i in scaffold_correlation_matrixs/scaf_correlation_matrix_clade_Paraburkholderia*
do
    echo $i
    CLADE=$(echo $i | cut -d'_' -f 7 | cut -d'.' -f 1)
    sbatch --wrap="python3 plot_heatmap.py $i" -c2 --time=2:0:0 -o slurms/slurm_heatmap_"$CLADE"_%J.out
    sbatch --wrap="python3 convert_corrMatrix_to_mcl_input.py $i" -c2 --time=2:0:0 -o slurms/slurm_mcl_"$CLADE"_%J.out
done
