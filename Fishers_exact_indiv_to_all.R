library(tidyverse)
library(patchwork)
library(gggenes)
library(ggfittext)

fet_wild <- read.delim ("./grenedalf.frequency.A2_Negative_combined.frequency.csv.neg_log_pval_new", sep=",", header=F)
fet_foco17<- read.delim("./grenedalf.frequency.Foco17Positive_merged_Foco17Negative_merged.frequency.csv.neg_log_pval", sep= ",", header=F)
fet_foco23<- read.delim("./grenedalf.frequency.Foco23Positive_Foco23Negative_merged.frequency.csv.neg_log_pval", sep= ",", header=F)
fet_25197<-read.delim("./grenedalf.frequency.25197-er-uninfected_25197-er-infected.frequency.csv.neg_log_pval", sep= ",", header=F)


colnames(fet_25197)<-(c("chrom", "pos", "ref", "alt", "inf_ref_count", "inf_alt_count", "neg_ref_count", "neg_alt_count", "log_p"))
colnames(fet_wild)<-(c("chrom", "pos", "ref", "alt", "inf_ref_count", "inf_alt_count", "neg_ref_count", "neg_alt_count", "log_p"))
colnames(fet_foco17)<-(c("chrom", "pos", "ref", "alt", "neg_ref_count", "neg_alt_count", "inf_ref_count", "inf_alt_count", "log_p"))
fet_foco17<-fet_foco17[c("chrom", "pos", "ref", "alt", "inf_ref_count", "inf_alt_count", "neg_ref_count", "neg_alt_count", "log_p")]
colnames(fet_foco23)<-(c("chrom", "pos", "ref", "alt", "neg_ref_count", "neg_alt_count", "inf_ref_count", "inf_alt_count", "log_p"))
fet_foco23<-fet_foco23[c("chrom", "pos", "ref", "alt", "inf_ref_count", "inf_alt_count", "neg_ref_count", "neg_alt_count", "log_p")]



fet_foco17$pool<-"Foco17"
fet_foco23$pool<- "Foco23"
fet_25197$pool<-"DGRP-517"
fet_wild$pool<-"Foco23 Wild"

df_combined <- rbind(fet_25197, fet_foco17, fet_foco23, fet_wild)

