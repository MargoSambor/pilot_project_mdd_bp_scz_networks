#!/bin/bash

: '
for ii in "top_200_30_30_70_70" "top_200_35_35_75_75" "top_200_50_50_75_75" "top_200_ht_50_30_50_70" "top_200_ht_75_30_75_70"
	do
		results_folder="results/"$ii
		data_folder=$ii
		while read -r test
		do
			stringarray=($test)
			#for j in "TCSS" "resnikbma" "resnikmax" "gic"
			#for j in "resnikbma" "resnikmax"
			for j in "TCSS" "gic"
			do
				for k in "0.5" "0.7" "0.85"
				do
	        	               for i in "0.75" "0.8" "0.85" "0.9" "0.95"
        	        	       do
						mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/upregulated/fisher/proteome/$i
                        	        	mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/upregulated/cluster/$i
	                                        mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/downregulated/fisher/proteome/$i
		                                mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/downregulated/cluster/$i
					done
					echo bsub -q research -M 2048 -R "rusage[mem=2048]" -o log_ph.txt -e err_ph.txt python propagation_phospho.py $j $k $test $results_folder $data_folder
				done
			done
		done < ../data/top_weighted.txt
	done
'

results_folder="netextract_results/commonmind_scz_padj_0.05"
data_folder=""
while read -r test
	do
	stringarray=($test)
	for j in "TCSS" "gic"
		do
		for k in "0.5" "0.7" "0.85"
			do
        	               for i in "0.75" "0.8" "0.85" "0.9" "0.95"
       	        	       do
					mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/upregulated/fisher/proteome/$i
                       	        	mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/upregulated/cluster/$i
                                        mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/downregulated/fisher/proteome/$i
	                                mkdir -p ../$results_folder/$j/$k/${stringarray[0]}/downregulated/cluster/$i
				done
				bsub -q research -M 2048 -R "rusage[mem=2048]" -o log_ph.txt -e err_ph.txt python propagation_phospho.py $j $k $test $results_folder $data_folder
			done
		done
	done < /nfs/research/petsalaki/shared_folder/diffusion/data/TCGA.txt
