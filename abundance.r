##Loaded package
library(microeco)
library(magrittr)
set.seed(123)
library(ggplot2)
theme_set(theme_bw())

##extract protist relativate abundance across each taxon
setwd("F:/³ÌÐò/github/gut-protist/")
sample_info_18S<-read.delim('1264_sample_infor_18s',header = TRUE,sep = "\t", row.names = 1,stringsAsFactors = FALSE)
otu_table_18S<-read.delim('1264_OTU_table_18s',header = TRUE,sep = "\t", row.names = 1,stringsAsFactors = FALSE)
taxonomy_table_18S<-read.delim('taxonomy_table_18s',header = TRUE,sep = "\t", row.names = 1,stringsAsFactors = FALSE)
taxonomy_table_18S %<>% tidy_taxonomy
protist <- microtable$new(sample_table = sample_info_18S, otu_table = otu_table_18S, tax_table = taxonomy_table_18S)
dataset$tidy_dataset()
dataset$cal_abund()
dataset$save_abund(dirpath = "protist_abundance")

##extract prokaryotes relativate abundance across each taxon
sample_info_16S<-read.delim('1264_sample_infor_16s',header = TRUE,sep = "\t", row.names = 1,stringsAsFactors = FALSE)
otu_table_16S<-read.delim('1264_OTU_table_16s',header = TRUE,sep = "\t", row.names = 1,stringsAsFactors = FALSE)
taxonomy_table_16S<-read.csv('taxonomy_table_16s.csv',header = TRUE,sep = "\t", row.names = 1)
taxonomy_table_16S %<>% tidy_taxonomy
prokaryote <- microtable$new(sample_table = sample_info_16S, otu_table = otu_table_16S, tax_table = taxonomy_table_16S)
dataset$tidy_dataset()
dataset$cal_abund()
dataset$save_abund(dirpath = "prokaryote_abundance")

