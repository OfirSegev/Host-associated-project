for i in $(ls -tr ./slurms/*)
do 
	number=$(echo $i | grep -oE "[[:digit:]]{6,9}")
       	echo -n $i'	' >> slurm_summary.tsv
       	seff $number | grep -o "State: [A-Za-z]\+" >> slurm_summary.tsv 
done
