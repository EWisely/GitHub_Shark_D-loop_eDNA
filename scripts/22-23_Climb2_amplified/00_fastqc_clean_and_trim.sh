#Title: 00_fastqc_clean_and_trim.sh
#Description: take demultiplexed reads and clean all primers, adapters, and low-quality bases out of them before metabarcoding analysis
#Author: Eldridge Wisely
#Date: March 18, 2024

#Take arguments to specify primer set (have the ones I've already used pre-loaded or enter the primers and RC of them), and the adapter style (TruSeq of Nextera)


#eDNA-Climb(2) primers F:  CGCGTCAAGAATGCCAGTCC  R: TCCAAACCCGGGGTGAACCT

#Make a sample list for samples in their own folders in the data directory

#run from the 2024_Analysis_Shark_D-loop_eDNA/data/22-23 amplified then unamplified directories

#get trimmomatic TruSeq adapter file here: https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa
cd trimmomatic-adapter-files/
cp ~/Downloads/TruSeq3-PE-2.fa .
cd ..

cd 00_raw_demuxed
cp */*R1_001.fastq.gz .
cp */*R2_001.fastq.gz .
#rm -r *.ds*
ls *R1_001.fastq.gz |cut -d '_' -f 1,2 >../Sample_list.txt

#-G left side of R2 =reverse primer
#-a right side of R1 =RC of reverse primer
#-g left side of R1 = forward primer
#-A right side of R2 =RC of forward primer

 
mkdir ../01_cleaned_trimmed/cutadapt

while read line; do
cutadapt \
-q 20 \
--trim-n \
--minimum-length 10 \
-e 0.1 \
-A GGACTGGCATTCTTGACGCG \
-g CGCGTCAAGAATGCCAGTCC \
-a AGGTTCACCCCGGGTTTGGA \
-G TCCAAACCCGGGGTGAACCT \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-e .15 \
-g GAAGAGCACACGTCTGAACTCCAGTCACATGATCATGGTCACCCTATCTC \
-g GATCGGAAGAGCACACG \
-g CGGAAGAGCACACGTCTGAACTCCAG \
-a GATCGGAAGTGCGTCGTGTAGGGATAGAGAGTGCCCTTAATTGCTATCTT \
-a GATCGGAAGAGCGTCGAGTAGGGAAAGAAAGTGAGCATACGTGAAAAAGT \
-o "../01_cleaned_trimmed/cutadapt/""$line""_R1.cleaned.fastq.gz" \
-p "../01_cleaned_trimmed/cutadapt/""$line""_R2.cleaned.fastq.gz" \
"$line""_L001_R1_001.fastq.gz" \
"$line""_L001_R2_001.fastq.gz"; \
done < ../Sample_list.txt


cd ../01_cleaned_trimmed/cutadapt
mkdir fastqc
conda activate qc
fastqc * -o fastqc/
cd fastqc
multiqc *

#Use trimmomatic to polish


cd ../
mkdir trimmomatic
while read line; do
trimmomatic PE -trimlog trimmomatic/trimmomatic.log "$line""_R1.cleaned.fastq.gz" "$line""_R2.cleaned.fastq.gz" "trimmomatic/""$line""_final_R1.fastq.gz"  "trimmomatic/""$line""_cleaned.trimmed_R1.discards.fastq.gz" "trimmomatic/""$line""_final_R2.fastq.gz" "trimmomatic/""$line""_cleaned.trimmed_R2.discards.fastq.gz" -validatePairs ILLUMINACLIP:../../../trimmomatic-adapter-files/TruSeq3-PE-2.fa:2:30:10:2:true SLIDINGWINDOW:4:15 TRAILING:20;done < ../../Sample_list.txt

cd trimmomatic
mkdir fastqc
conda activate qc
fastqc * -o fastqc/
cd fastqc
multiqc *

cd ..
gunzip *_final_*gz

#inspect the multiqc files to make sure they're ready for import into Dada2
