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
usearch -annot ../obi3_Climb2/EMBL_assigned_ASVs.fasta -knowndb Climb2_21_blood_haps.fasta -db Climb_all_haps_July2023.fasta -fastaout EMBL_ecotagged_Climb2_ASVs_annot_w_blood_haps.fasta
less EMBL_ecotagged_Climb2_ASVs_annot_w_blood_haps.fasta | grep -c "Wisely"
468
less EMBL_ecotagged_Climb2_ASVs_annot_w_blood_haps.fasta | grep -c "Wisely_C_limb_2_hap_2"
33
less EMBL_ecotagged_Climb2_ASVs_annot_w_blood_haps.fasta | grep -c "Wisely_C_limb_2_hap_1"
435
#all of the TRUE IDs were annotated with a Wisely blood hap when using them as the knowndb and all haps from NCBI as the db
#editme
