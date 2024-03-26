#Eldridge Wisely
#run Dada2 on popgen eDNA primers

#for primer and adapter removed reads


### following this tutorial https://benjjneb.github.io/dada2/tutorial.html ###
#install dada2
## try http:// if https:// URLs are not supported
#if (!requireNamespace("BiocManager", quietly=TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
#and this package which can do lca https://rdrr.io/github/ammararuby/MButils/man/
#install.packages("remotes")
#remotes::install_github("ammararuby/MButils")

library(dada2); packageVersion("dada2")
help(package="dada2")
## the tutorial is based on version '1.16.0' but this is version 1.26.0
#setwd("/Users/Eldridge/Desktop/Marine Population Genetics with eDNA/Dada2_Climb2_analysis")
#path <- "/Users/Eldridge/Desktop/Marine Population Genetics with eDNA/Dada2_Climb2_analysis/Climb2_fastqs/cutadapt/2nd_cutadapt/annotated_and_concatenated/"
list.files(path)

library(plyr)
library(dplyr)
library(tidyr)
library(dada2)
library(ggplot2)


list.files(path)
#patternR1<-"_L001_R1_001"
#patternR2<-"_L001_R2_001"

patternR1<-"_R1.final.annotated"
patternR2<-"_R2.final.annotated"

dadapath<-"/Users/Eldridge/Desktop/Marine Population Genetics with eDNA/Dada2_Climb2_analysis"

#basicdada <- function(path, dadapath, lengthMin, lengthMax, patternR1, patternR2, locus){
# Chemin vers les R1-R2 sans primers
fnFs <- sort(list.files(path, pattern=paste0(patternR1, ".fastq"), full.names = TRUE))
fnRs <- sort(list.files(path, pattern=paste0(patternR2, ".fastq"), full.names = TRUE))


## the function only extracts files that end with the chosen pattern and
## they are extracted with their whole path

## then you can only keep the part of your files name you want :

#sample.names <- sapply(strsplit(basename(fnFs), "_L001_R1_001"), '[', 1)
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)
print("- Process sample files :")

print(fnFs)
print(fnRs)
#write.csv(fnFs,"fnFs.csv")
#write.csv(fnRs, "fnRs.csv")
write.csv(sample.names, "sample.names.dada2.csv")

# Chemins d'entree-sortie
filt_path <- file.path(dadapath, "filtered") # On mettra les reads trimmes dans filtered
filtFs <- file.path(filt_path, paste0(sample.names, paste0("_", "R1_filt.fastq")))
filtRs <- file.path(filt_path, paste0(sample.names, paste0("_", "R2_filt.fastq")))

print(filtFs)
print(filtRs)

#write.csv(filtFs, "filtFs.csv")
#write.csv(filtRs, "filtRs.csv")
# Examiner la qualite des reads
print("- Plotting quality profiles :")

plotQualityProfile(fnFs[1:10])
ggsave(paste0(dadapath, "/filtered", "R1_filt", ".png"), device="png")
print(paste0(dadapath, "/filtered", "R1_filt", ".png"))

plotQualityProfile(fnRs[1:10])
ggsave(paste0(dadapath, "/filtered", "R2_filt", ".png"), device="png")
print(paste0(dadapath, "/filtered", "R2_filt", ".png"))

# Trimming
print("## Filter and Trim ##")
#  if (locus=="18S-V9") {
#    trunc_len <- 110
#  } else {
#    if (locus=="COI") {
#      trunc_len <- 200
#    }  else {
#      trunc_len <- 220
#    }
#  }



#old command left here for reference of options
#  print(paste0("- Truncation length : ", trunc_len))
#out <- filterAndTrim(fwd=fnFs,filt=filtFs, rev=fnRs, filt.rev=filtRs, truncLen=c(220,180), maxN=0, maxEE=2, verbose=TRUE, multithread=TRUE) #removed trunLen=100 replaced with truncLen=c(220,180),removed matchIDs=TRUE

#expected insert length is 435, so shouldn't truncate forward and reverse reads below 440/2+20=240 minimum... can filter for length of total fragment later.
#matchIDs=TRUE could be messed up by the annotation step before putting them into this pipeline...

out <- filterAndTrim(fwd=fnFs,filt=filtFs, rev=fnRs, filt.rev=filtRs, maxN=0, maxEE=2, verbose=TRUE, multithread=TRUE)

out <- as.data.frame(out)
#out$sample <- sub("_L001_R1_001.fastq.gz", "", rownames(out))
out$sample<- sapply(strsplit(basename(rownames(out)), "_"), '[', 1)

out <- out[,c(3,1,2)]
saveRDS(out, paste0(dadapath, "/1_filterAndTrim_", patternR1, patternR2, ".rds"))
write.table(out, file=paste0(dadapath, "/1_filterAndTrim_", patternR1, patternR2, ".csv"), sep=";", col.names=T, row.names=F)

# Mettre a jour les donnees : enlever les echantillons vides
exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]
sample.names <- sapply(strsplit(basename(filtFs), "_"), '[', 1)


# Apprentissage
print("## Learning Errors ##")
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE, randomize=TRUE)
saveRDS(errF, paste0(dadapath, "/2_err", patternR1, ".rds"))
saveRDS(errR, paste0(dadapath, "/2_err", patternR2, ".rds"))
print(paste0(dadapath, "/2_err", patternR1, ".rds"))
print(paste0(dadapath, "/2_err", patternR2, ".rds"))

# Plot des erreurs
print("- Plotting errors :")
plotErrors(errF, nominalQ=TRUE)
ggsave(paste0(dadapath, "/2_err", patternR1, ".png"), device="png")
plotErrors(errR, nominalQ=TRUE)
ggsave(paste0(dadapath, "/2_err", patternR2, ".png"), device="png")
print(paste0(dadapath, "/2_err", patternR1, ".png"))
print(paste0(dadapath, "/2_err", patternR2, ".png"))

# Dereplication des sequences identiques. NB: en plus, dada2 fait un consensus (qualite moyenne) des sequences identiques. Le modele va s'en servir.
print("## Dereplicating identical sequences ## ")
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
saveRDS(derepFs, paste0(dadapath, "/3_derep", patternR1, "s.rds"))
saveRDS(derepRs, paste0(dadapath, "/3_derep", patternR2, "s.rds"))
print(paste0(dadapath, "/3_derep", patternR1, "s.rds"))
print(paste0(dadapath, "/3_derep", patternR2, "s.rds"))

# Infer the sequence variants in each sample
print("## ASV inference ## ")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, paste0(dadapath, "/4_dada", patternR1, "s.rds"))
saveRDS(dadaRs, paste0(dadapath, "/4_dada", patternR2, "s.rds"))
print(paste0(dadapath, "/4_dada", patternR1, "s.rds"))
print(paste0(dadapath, "/4_dada", patternR2, "s.rds"))

setwd("~/Desktop/Marine Population Genetics with eDNA/Dada2_Climb2_analysis")
save.image("Dada2_Climb2.RData")
#load("Dada2_Climb2.RData")

# Contigage des sequences obtenues
print("## Merging Sequences ## ")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, maxMismatch=0, minOverlap=10, verbose=TRUE)
saveRDS(mergers, paste0(dadapath, "/5_contigs", "Climb2_21", ".rds"))
print(paste0(dadapath, "/5_contigs", "Climb2_21", ".rds"))

#try again requiring more overlap

mergers1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, maxMismatch=0, minOverlap=20, verbose=TRUE)
saveRDS(mergers1, paste0(dadapath, "/5_contigs", "Climb2_21-moreoverlap", ".rds"))
print(paste0(dadapath, "/5_contigs", "Climb2_21-moreoverlap", ".rds"))


# Construire une table des sequences
print("## Making sequence table ## ") 
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, paste0(dadapath, "/6_seqtab", "Climb2_21", ".rds"))
print(paste0(dadapath, "/6_seqtab", "Climb2_21", ".rds"))

# Inspecter la distribution de taille des sequences uniques 
print("- Inspect ASV length distribution : ")
length_distribution = table(nchar(getSequences(seqtab)))
table(nchar(getSequences(seqtab)))


plot(length_distribution)
#target length without primers =395bp for C.limb2 according to Sanger data analysis. But 435 according to primer design in Geneious (re-check this). So set the min and max -15bp, +5bp?
#in the plot, most sequences (almost all, were over 400 and under 455)
lengthMin<-"380"
lengthMax<-"440"
saveRDS(length_distribution, paste0(dadapath, "/7_length_distribution", "Climb2_21", ".rds"))
print(paste0(dadapath, "/7_length_distribution", "Climb2_21", ".rds"))


# Enlever tous les ASV qui sont beaucoup plus grands/petits que la longueur attendue
print(paste0("- Removing ASV outside length range [", lengthMin, ":", lengthMax, "]"))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(lengthMin,lengthMax)]
saveRDS(seqtab2, paste0(dadapath, "/7_seqtab_trimmed", "Climb2_21", ".rds"))
print(paste0(dadapath, "/7_seqtab_trimmed", "Climb2_21", ".rds"))

# Enlever les chimeres, last dada step
print("## Removing chimeras ## ")
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab2.nochim, paste0(dadapath, "/7_nochimeras", "Climb2_21", ".rds"))
print(paste0(dadapath, "/7_nochimeras", "Climb2_21", ".rds"))



# Inspecter les chimeres
dim(seqtab2)
dim(seqtab2.nochim)
sum(seqtab2.nochim)/sum(seqtab2)

#99.5% non-chimeras in the trimmed seqtab2 

#Check des reads, table finale
getN <- function(x) sum(getUniques(x))
print("Creating summary table")
read.track <- cbind(sample=sample.names, 
                    dada_F=sapply(dadaFs,getN),
                    dada_R=sapply(dadaRs,getN),
                    merged_reads=sapply(mergers, getN), 
                    reads_before_chimera_removal=rowSums(seqtab2),
                    non_chimeric_reads=rowSums(seqtab2.nochim))
read.track <- merge(out, read.track, by="sample")
read.track$dada_F <- as.numeric(as.character(read.track$dada_F))
read.track$dada_R <- as.numeric(as.character(read.track$dada_R))
read.track$merged_reads <- as.numeric(as.character(read.track$merged_reads))
read.track$reads_before_chimera_removal <- as.numeric(as.character(read.track$reads_before_chimera_removal))
read.track$non_chimeric_reads <- as.numeric(as.character(read.track$non_chimeric_reads))
read.track$final_perc_reads_retained <- round(as.numeric(as.character(read.track$non_chimeric_reads))/as.numeric(as.character(read.track$reads.in))*100, 1) 

write.table(read.track, file=paste0(dadapath, "/8_readtrack", "Climb2_21", ".csv"), sep=";", col.names=T, row.names=F)
saveRDS(read.track, paste0(dadapath, "/8_readtrack", "Climb2_21", ".rds"))


##########################
#save final asv tables
print("- Saving final ASV table :")
asvname <-("Climb2_21")
saveRDS(seqtab2.nochim, paste0(dadapath, "/9_", asvname,"_asv.rds")) #used one containing potential chimeras due to the nature of my true haplotypes. #for chimera filtering, edit this to use seqtab2.nochim
print(paste0(dadapath, "/9_", asvname,"_asv.rds"))

#####copied from dada2_output_files.R from the french abyss project and edited
#Path for inputs (dada_dir or output_merge) and outputs
#ASV Table 
finaltab = readRDS("9_Climb2_21_asv.rds")
#I used my GoC MyBaits sequence target file for taxonomic assignment.  It's all GenBank mito sequences from Fish, Sharks,Cephalopods, and Marine Mammals found in the GoC according to Morzaria-Luna et al. as of download in December 2022.
rdpbank = ("/Users/Eldridge/Desktop/Marine Population Genetics with eDNA/MyBaits Design/baits_fasta/fasta_out/max_marine_mito.rdp.fasta")


print("## Assigning ASV taxonomy ## ")
taxa_max_marine_mito_baits <- assignTaxonomy(finaltab, rdpbank, minBoot=50, verbose=TRUE, outputBootstraps=TRUE, tryRC=TRUE, multithread=TRUE)
saveRDS(taxa_max_marine_mito_baits, paste0(dadapath, "/9_taxonomy_max_marine_mito_baits.rds"))

# Create output data
asv_seqs <- colnames(finaltab)
asv_headers <- vector(dim(finaltab)[2], mode="character")
for (i in 1:dim(finaltab)[2]) {
  asv_headers[i] <- paste(">Climb2_ASV", i, sep="_")
}
# ASV fasta files
print("Write ASV fasta file")
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste0(dadapath, "/10_ASVs.fasta"))
# Count table
print("Write count table to file")
asv_tab <- as.data.frame(t(finaltab))
rownames(asv_tab) <- sub(">", "", asv_headers)
asv_tab$ClusterID <- rownames(asv_tab)
asv_tav <- asv_tab[,c(3,1,2)]
ClusterID <- sub(">", "", asv_headers)
asv_tab <- cbind(ClusterID, asv_tab)
asv_tab <-cbind(asv_tab, asv_seqs)
write.table(asv_tab, paste0(dadapath, "/10_ASVs_counts.tsv"), sep="\t", quote=F, col.names=T, row.names=F)
# Taxa table
print("Write taxa table to file")
asv_tax <- as.data.frame(taxa_max_marine_mito_baits)#2
asv_tax <- cbind(ClusterID, asv_tax)
rownames(asv_tax) <- sub(">", "", asv_headers)
asv_tax$ClusterID <- rownames(asv_tax)
#asv_tax1 <- asv_tax[, c(13,1,2,3,4,5,6,7,8,9,10,11,12)]
write.table(asv_tax, paste0(dadapath, "/10_ASVs_taxonomy.tsv"),  sep="\t", quote=F, col.names=T, row.names=T)

save.image("After_Dada2_Climb2.RData")


####HANDOFF TO OBI3 for ecotag step######
###Then usearch for annotation with blood and all NCBI D-loop haps####

############ HANDOFF TO Metabar ################
metabar_reads_table <- as.data.frame(t(asv_tab))
metabar_reads_table <-
  rownames_to_column(metabar_reads_table, var = "id")
write_csv(metabar_reads_table,
          "data/21_Climb2/04_Usearch_annotated/21_Climb_dada2_metabar_reads_table.csv")


