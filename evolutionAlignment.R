#!/usr/bin/R
# Universitat Potsdam
# Author Gaurav Sablok
# date: 2024-3-17
evolutionAlignment <- function(fasta_file, path){
    library(ape)
    library(Biostrings)
    library(seqinr)
    library(msa)
    # a evolutionaryAlignment function that calculates the evolutionary
    # rates across the fasta file, given a fasta file, it will automatically
    # align the sequences and then will estimate the evolutionary rates on 
    # the clean alignments.To use this 
    # evolutionAlignment("/Users/gauravsablok/Desktop/CodeRelease/fasta_sample_datasets/read_check.fasta", 
    #               path = "/Users/gauravsablok")
    # it will also plot the evolutionary rates and it uses two methods ka/ks and dn/ds
    setwd(path)
    fasta_file <- readDNAStringSet(fasta_file)
    fasta_alignment <- msa(fasta_file)
    fasta_convert <- msaConvert(fasta_alignment, type = "ape::DNAbin")
    write.FASTA(fasta_convert, file = "fasta_convert.fasta")
    fasta_re_check <- read.alignment(paste(getwd(),
                            "fasta_convert.fasta", sep = "/"), format = "fasta")
    kaks_fasta <- kaks(fasta_re_check, rmgap = TRUE)
    dnds_fasta <- dnds(fasta_re_check)
    ka_histogram <- hist(kaks_fasta$ka, main = paste("Histogram of ka values"))
    ka_histogram <- hist(kaks_fasta$ka, main = paste("Histogram of ka values"))
    save.image(file = "ka_histogram", safe = TRUE)
    save.image(file = "ks_histogram", safe = TRUE)
}
