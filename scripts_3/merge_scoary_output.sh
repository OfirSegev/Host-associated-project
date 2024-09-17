#!/bin/bash

clade_sample=$1
COLLAPSE=$2
prefix=$3
LABEL=$4

#arrange output files inside one folder
cd outputs_and_analysis
mkdir all_csvs

while IFS= read -r folder; do
	file_name=$(ls $folder/*csv | cut -d'/' -f 2)
	echo $file_name
	cp $folder/$file_name all_csvs
	cd all_csvs
	mv $file_name $folder.csv
	cd ..
done < <( ls -d $LABEL*/ | cut -f1 -d'/' )

# merge data to one table (remember scoary needs to finish)
cd all_csvs

python3 $prefix/combine_to_one_list.py $clade_sample $COLLAPSE $LABEL
