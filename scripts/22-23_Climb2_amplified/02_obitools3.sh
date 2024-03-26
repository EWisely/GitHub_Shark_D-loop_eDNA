#I could try the ecotag assignment step on my Dada2 sequences!  
obi import --fasta-input Climb2_Dada2_ASVs.fasta Climb2/Dada2_ASVs
obi ecotag -m 0.97 --taxonomy Climb2/taxonomy/my_tax -R Climb2/Climb2_refs_clean_97 Climb2/Dada2_ASVs Climb2/EMBL_assigned_ASVs
obi less Climb2/EMBL_assigned_ASVs
obi stats -c SCIENTIFIC_NAME Climb2/EMBL_assigned_ASVs
#2023-08-09 08:54:04,970 [stats : INFO ]  obi stats
  99.9 % |#################################################| ] remain : 00:00:00
SCIENTIFIC_NAME         	count	total
None                    	711	711	
Carcharhinus limbatus   	388	388	
Carcharhinus            	80	80	
2023-08-09 08:54:04,976 [stats : INFO ]  Done.

obi export --fasta-output Climb2/EMBL_assigned_ASVs -o EMBL_assigned_ASVs.fasta
cd ../usearch_analysis
#Imported fasta file made from Nexus traits file for the entire worldwide database of known blacktip haplotypes from NCBI and this study's blood samples.  
#filename: 23_bt_eDNA_haps.fasta

#Annotated (Dada2 ASVs ecotagged with EMBL) with that.
 usearch -annot ../obi3_Climb2/EMBL_assigned_ASVs.fasta -db 23_bt_eDNA_haps.fasta --fastaout 21GAL_Dada2_obi3_eDNAhaps_annot.fasta

less 21GAL_Dada2_obi3_eDNAhaps_annot.fasta | grep -c "Hap_5"
4
(base) eldridge@panthera:/fastone/Shark_popgen_21GAL/usearch_analysis$ less 21GAL_Dada2_obi3_eDNAhaps_annot.fasta | grep -c "Hap_15"
33
(base) eldridge@panthera:/fastone/Shark_popgen_21GAL/usearch_analysis$ less 21GAL_Dada2_obi3_eDNAhaps_annot.fasta | grep -c "Hap_23"
439

obi import --fasta-input ../usearch_analysis/21GAL_Dada2_obi3_eDNAhaps_annot.fasta Climb2/annotated_23bthaps_results_2021
obi export --tab-output Climb2/annotated_23bthaps_results_2021 >annotated_23bthaps_2021_results.tab

#now, try obi3 classification before annotation just like the 21 ASVs
cd obi3_Climb2
obi import --fasta-input ../22-23_obi3_Climb2/10_ASVs_Climb-22-23.fasta Climb2/22-23_Dada2_ASVs
obi ecotag -m 0.97 --taxonomy Climb2/taxonomy/my_tax -R Climb2/Climb2_refs_clean_97 Climb2/22-23_Dada2_ASVs Climb2/EMBL_assigned_ASVs_22-23

obi stats -c SCIENTIFIC_NAME Climb2/EMBL_assigned_ASVs_22-23
2024-03-19 17:56:58,207 [stats : INFO ]  obi stats
  99.8 % |#################################################\ ] remain : 00:00:00
SCIENTIFIC_NAME            	count	total
None                       	649	649	
Carcharhinus limbatus      	3	3	
Carcharhinus falciformis   	1	1	
2024-03-19 17:56:58,213 [stats : INFO ]  Done.

obi export --fasta-output Climb2/EMBL_assigned_ASVs_22-23 -o Climb2_22-23_ASVs_obi3_results.fasta
cd ../usearch_analysis/

usearch -annot ../obi3_Climb2/Climb2_22-23_ASVs_obi3_results.fasta -db 23_bt_eDNA_haps.fasta --fastaout Climb2-22-23_ASVs_obi3_annot_bt_eDNAhaps.fasta
obi import --fasta-input ../usearch_analysis/Climb2-22-23_ASVs_obi3_annot_bt_eDNAhaps.fasta Climb2/annotated_23bthaps_results_22-23
obi export --tab-output Climb2/annotated_23bthaps_results_22-23 >annotated_23bthaps_22-23_results.tab

