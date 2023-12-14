#!/usr/bin/env Rscript
# nolint

###############
# Purpose: Pull out sequences of genes of interest from bacterial taxa related to mycobacterium to create a phylogenetic tree
# Contact: JReceveur jreceveur@som.umaryland.edu
# Last updated: 1 Oct 2023

#####
# Overview:
#####
# Using the package biomart, download the genome for a taxa of interest
# Run blastn against the downloaded genome and genes of interest
#   Uses the blast_nucleotide_to_nucleotide() function from metablastr instead of writing own parsing
# Takes the best hit (determined by bit score) and returns the sequence from that area of the subject strand
#   plus strand results are returned exactly while sequences matching the query sequences 'backwards' return the reverse complement
# Iterates through each gene stepwise
# After iterating through all genes of interest, the downloaded genome is removed to save space
#   NOTE: Only the first refseq accession ID is used currently 
#   NOTE: No filtering criteria is applied in this script
# Results are saved in the table Combined and by gene



#####
# Inputs: 
#####
# List of taxa on interest (using NCBI names)
# Reference sequences of genes of interest (.fna)
#   NOTE: All current references are from H37Rv

#####
# Outputs:
#####
# Combined table of all blast results with sequence info (.csv)
#   NOTE: genome-gene combinations for which there was no blast match are silent
# Individual tables (per gene basis)
# List of files for which no genome was found

#####
# Enhancements:
#####
#   Catch and retry on NCBI timeouts
#   Clean up script/paths to be more portable
#   Match formatting with expected input to tree code
#   Add filtering criteria
#   Write out table genome-gene combinations for which there was no match (currently printed in console)





######
# Package install
######
# *also need blast, see installation vignette  https://github.com/ropensci/biomartr

# BiocManager::install(
#   c(
#     "Biostrings",
#     "GenomicFeatures",
#     "GenomicRanges",
#     "Rsamtools",
#     "IRanges",
#     "rtracklayer")
# )
# 
# devtools::install_github("drostlab/metablastr", build_vignettes = TRUE, dependencies = TRUE)

# Need the path to blast on computer/instance


#####
# Package load
#####

# Full paths to everything, really only needs the path to blast (overkill here, clean up later or put blast path in function)
Sys.setenv(PATH="/Users/jreceveur/.rvm/gems/ruby-3.0.0/bin:/Users/jreceveur/.rvm/gems/ruby-3.0.0@global/bin:/Users/jreceveur/.rvm/rubies/ruby-3.0.0/bin:/Users/jreceveur/google-cloud-sdk/bin:/Users/jreceveur/opt/anaconda3/bin:/Users/jreceveur/opt/anaconda3/condabin:/Library/Frameworks/Python.framework/Versions/3.9/bin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin:/Users/jreceveur/.rvm/bin:/Library/Frameworks/Python.framework/Versions/3.9/bin:/Library/Frameworks/Python.framework/Versions/3.9/bin:/Users/jreceveur/.rvm/gems/ruby-3.0.0/bin:/Users/jreceveur/.rvm/gems/ruby-3.0.0@global/bin:/Users/jreceveur/.rvm/rubies/ruby-3.0.0/bin:/Library/Frameworks/Python.framework/Versions/3.9/bin:/usr/local/bin:/System/Cryptexes/App/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/local/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/bin:/var/run/com.apple.security.cryptexd/codex.system/bootstrap/usr/appleinternal/bin:/Users/jreceveur/.rvm/bin:/Users/jreceveur/Applications/quarto/bin:/Library/TeX/texbin:/usr/texbin:/Applications/RStudio.app/Contents/Resources/app/quarto/bin:/Applications/RStudio.app/Contents/Resources/app/bin/postback:/Applications/RStudio.app/Contents/Resources/app/bin/postback")



library(biomartr)
library(metablastr)
library(Biostrings)
library(R.utils)
library(dplyr)

print("Packages loaded")
options(timeout = 3000) # Increase timeout to help with NCBI database pulls


setwd("~/Documents/MycoSandel2023/Loop")

TaxaList<-read.csv("TaxaList.csv",header=T) # List of taxa to work on
getKingdomAssemblySummary(db = "refseq", skip_bacteria = FALSE)
# create an empty vector for taxa without a refseq reference
TaxaMissed<-NULL 
GenesMissed<-NULL
summary_table<-list()
if (!dir.exists("BlastOutputs")) {dir.create("BlastOutputs")}
tax_counter<-0 # to keep track of what species stopped at in case of errors
for (i in TaxaList$NCBI.Name){
  tax_counter<-tax_counter+1
  print(tax_counter)
  #getKingdomAssemblySummary(db = "refseq", skip_bacteria = FALSE)
  Available<-is.genome.available(organism = i, db = "refseq",skip_bacteria = "FALSE")
 
  if (Available ==T){ # Check if ref genome is available 
    
    if(dir.exists("genomes")){ # If genome directory is still around from previous loop, remove
      unlink("genomes", recursive=TRUE)
      dir.create("genomes")
    }
    else {
      dir.create("genomes")}
    suppressMessages(genome <- getGenome(db = "refseq", 
                                 organism = i,  
                                 path   = "genomes", 
                                 reference = FALSE))
    genome_seq<-read_genome(genome)
    CompressedGenome<-list.files("genomes",pattern="*.fna.gz")
    gunzip(paste0("genomes/",CompressedGenome),remove=T)
    ref_genome<-list.files("genomes",pattern="*.fna")
    ref_genome_path<-paste0("genomes/",ref_genome)
    GeneList<-list.files("~/Documents/MycoSandel2023/RefGenes") # directory contains fna files associated with genes of interest (all mtb h37)
    
    
    
    for (j in GeneList){ # Loop across all genes of interest
      suppressWarnings(rm(blast_test)) # remove blast results from previous loop
      try(
        # Suppressing messages here just because this function is extra verbose
        suppressMessages(blast_test <- blast_nucleotide_to_nucleotide(
        query   = paste0("~/Documents/MycoSandel2023/RefGenes/",j),
        subject = ref_genome_path,
        db.import  = FALSE,strand = 'both',
        output.path = "BlastOutputs/")))
      # no applicable method for 'collector_value' applied to an object 
      # of class seems to show up when no results are returned
      
      #####
      # Grabbing top sequence from reference and adding some metadata on the best_hit
      #####
      if (exists("blast_test")==T){
        blast_test_sort<-arrange(blast_test, -bit_score) # Sort by bit score
        best_hit<-blast_test[1,] # Take only highest score
        
        if(best_hit$s_start-best_hit$s_end >0){# If sequence matches to query 'backwards' take the reverse complement 
          best_hit$strand<-"minus"
          Fasta<-readDNAStringSet(ref_genome_path)
          Fasta2<-Fasta[grep(paste0(best_hit$subject_id,"*"), names(Fasta)), ] # Pull out the scaffold that the best blast hit came from
          
          gene_seq<-subseq(Fasta2, start = best_hit$s_end, end = best_hit$s_start)
          gene_seq_final<-reverseComplement(gene_seq) 
          start<-best_hit$s_end-(best_hit$q_len-best_hit$q_end)
          #start= best_hit$s_end-(difference between qlen and qalign)
          #end<-best_hit$s_start+best_hit$q_start
          end<-best_hit$s_start+best_hit$q_start
          
          
          if (start<0| end>best_hit$s_len){
            gene_seq_extended_final<-NA
          }
          else{gene_seq_extended<-subseq(Fasta2, start = start, end = end)
            gene_seq_extended_final<-reverseComplement(gene_seq_extended)}
          if (length(gene_seq_final)==1){
            best_hit$sequence<-as.character(gene_seq_final)
            best_hit$name<-names(gene_seq_final)
            best_hit$taxa<-i
            best_hit$gene<-j
            best_hit$gene_extended<-as.character(gene_seq_extended_final)
            summary_table[[paste0(j,":",i)]]<-best_hit}
          rm(best_hit)
          }
        else if(best_hit$s_start-best_hit$s_end <0){ # value below zero, matching to plus strand
          best_hit$strand<-"plus"
          Fasta<-readDNAStringSet(ref_genome_path)
          Fasta2<-Fasta[grep(paste0(best_hit$subject_id,"*"), names(Fasta)), ]
          gene_seq_final<-subseq(Fasta2, start = best_hit$s_start, end = best_hit$s_end)
          start<-best_hit$s_start-best_hit$q_start
          end<-best_hit$q_len-best_hit$q_end+best_hit$s_end
          if (start<0|end>best_hit$s_len){
            gene_seq_extended_final<-NA
          }
          else{gene_seq_extended_final<-subseq(Fasta2, start = start, end = best_hit$q_len-best_hit$q_end+best_hit$s_end)}
          if (length(gene_seq_final)==1){
            best_hit$sequence<-as.character(gene_seq_final)
            best_hit$name<-names(gene_seq_final)
            best_hit$taxa<-i
            best_hit$gene<-j
            best_hit$gene_extended<-as.character(gene_seq_extended_final)
            summary_table[[paste0(j,":",i)]]<-best_hit}
          rm(best_hit)
        }}
      else {
        print(paste0(i," has no BLAST match for ",j))
        Miss<-data.frame(taxa=i,gene=j)
        GenesMissed[[paste0(i,j)]]<-Miss
      }
    }
    #####
    # Cleanup
    #####
    system("rm -R genomes")
    system("rm -R BlastOutputs")
  }
else{TaxaMissed<-append(TaxaMissed,i)}}


#####
# Outputs
#####
Misses<-do.call(rbind,GenesMissed) 
Hits<-do.call(rbind,summary_table)
Combined<-full_join(Hits,Misses)

write.csv(Combined,"CombinedGeneList.csv")
write.csv(TaxaMissed,"TaxaMissingGenome.csv")
for (i in unique(Combined$gene)){
  subset_list<-subset(Combined,gene==i)
  write.csv(subset_list, paste0(gsub(".fna","",i),"GeneList.csv"))
}





#########
# Scratch
########

# scp /Users/jreceveur/Documents/MycoSandel2023/Sandel2023.R receveur@hpcc.msu.edu:/mnt/home/receveur/Documents/MycoPhylo

# conda activate ncbi_datasets
# 
# datasets download genome taxon Actinoalloteichus
# 
# 
# tar -xvf ncbi_dataset.zip
# 
# makeblastdb -in rpoC.fsa -dbtype nucl -parse_seqids
# 
# ncbi_dataset/data/GCA_000239155.2/GCA_000239155.2_ASM23915v2_genomic.fna
# 
# blastn -db ncbi_dataset/data/GCA_000239155.2/GCA_000239155.2_ASM23915v2_genomic.fna -query rpoC.fsa -out rpoCResults2.out
# 
# datasets download genome taxon 38288
# unzip ncbi_dataset
# 
# makeblastdb -in ncbi_dataset/data/GCF_000143825.1/GCF_000143825.1_ASM14382v1_genomic.fna  -dbtype nucl -parse_seqids
# 
# makeblastdb -in rpoB2.fsa  -dbtype nucl -parse_seqids
# 
# blastn -query rpoB2.fsa  -db ncbi_dataset/data/GCF_000143825.1/GCF_000143825.1_ASM14382v1_genomic.fna  -out rpoBResultsv3.out