#####make manhattan style plots with both individuals and with 25000bp windows
#####Input the files
fet_wild_i <- read.delim ("./A2_indiv.fet", sep="\t", header=F)
fet_foco17_1_i<- read.delim("./17_1_indiv.fet", sep= "\t", header=F)
fet_foco17_2_i<- read.delim("./17_2_indiv.fet", sep= "\t", header=F)
fet_foco23_1_i<- read.delim("./23_1_indiv.fet", sep= "\t", header=F)
fet_foco23_2_i<- read.delim("./23_2_indiv.fet", sep= "\t", header=F)
fet_foco23_n_i<- read.delim("./23_negative_i.fet", sep= "\t", header=F)
fet_foco17_n_i<- read.delim("./17_negative_i.fet", sep= "\t", header=F)
###Make all the input files the correct format
df_wild_i <- handle_input_fet(fet_wild_i, sample_set="wild")
df_17_1_i <- handle_input_fet(fet_foco17_1_i, sample_set="foco17_1")
df_17_2_i <- handle_input_fet(fet_foco17_2_i, sample_set="foco17_2")
df_23_1_i <- handle_input_fet(fet_foco23_1_i, sample_set="foco23_1")
df_23_2_i <- handle_input_fet(fet_foco23_2_i, sample_set="foco23_2")
df_17_n_i<- handle_input_fet(fet_foco17_n_i, sample_set="foco")
df_23_n_i<- handle_input_fet(fet_foco23_n_i, sample_set="foco")
#####Change the last column to a name that makes sense
names(df_17_1_i)[9] <- "p_val"
names(df_17_2_i)[9] <- "p_val"
names(df_23_1_i)[9] <- "p_val"
names(df_23_2_i)[9] <- "p_val"
names(df_wild_i)[9] <- "p_val"
names(df_23_n_i)[9] <- "p_val"
names(df_17_n_i)[9] <- "p_val"
######Add new column to say which pool it is
df_17_1_i$pool<- "Foco17 Lab 1"
df_17_2_i$pool<-"Foco17 Lab 2"
df_23_1_i$pool<- "Foco23 Lab 1"
df_23_2_i$pool<-"Foco23 Lab 2"
df_wild_i$pool<-"Foco23 Wild"
df_17_n_i$pool<- "Foco17 Negative"
df_23_n_i$pool<- "Foco23 Negative"

df_man<-rbind(df_wild_i, df_17_1_i, df_17_2_i, df_23_1_i, df_23_2_i)

#####The data imported above is the -log10 of the pval by default

##########Full genome plot of p-values################

df_man<- filter(df_man, p_val>8)
df_man<- filter(df_man, p_val>0)
df_man<- subset(df_man, chrom!="Y")


###rm(list=ls())
############################################################
#######p-value with significance lines
subset_i<-subset(df_man, p_val>80)

ggplot(filter(df_man)) +
  geom_point(aes(x=window/1e6, y=p_val, fill=pool), shape=21, size=2,stroke=NA, show.legend=FALSE) +
  scale_fill_manual(values = c("cyan", "cyan3", "coral", "coral3", "deeppink"))+
  geom_point(data=subset_i, aes(x=window/1e6, y=p_val), color="black")+
  theme_bw(base_size = 10) +
  ylab("-Log10(p)") +
  xlab("Position on chromosome (Mb)") +
  facet_grid(pool~chrom, scales = "free_x")+
  ggtitle("Individual p-value across the genome")+theme(panel.spacing.x = unit(0, "lines"))+
  theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())






####windowed p-value plot: should look the same but less
#------------------------------------------

#####Input the files
fet_wild <- read.delim ("./wild2.fet", sep="\t", header=F)
fet_foco17_1<- read.delim("./Foco17_1.fet", sep= "\t", header=F)
fet_foco17_2<- read.delim("./Foco17_2.fet", sep= "\t", header=F)
fet_foco23_1<- read.delim("./Foco23_1.fet", sep= "\t", header=F)
fet_foco23_2<- read.delim("./Foco23_2.fet", sep= "\t", header=F)
fet_foco23_n<- read.delim("./23_negative.fet", sep= "\t", header=F)
fet_foco17_n<- read.delim("./17_negative.fet", sep= "\t", header=F)

###Make all the input files the correct format
df_wild <- handle_input_fet(fet_wild, sample_set="wild")
df_17_1 <- handle_input_fet(fet_foco17_1, sample_set="foco17_1")
df_17_2 <- handle_input_fet(fet_foco17_2, sample_set="foco17_2")
df_23_1 <- handle_input_fet(fet_foco23_1, sample_set="foco23_1")
df_23_2 <- handle_input_fet(fet_foco23_2, sample_set="foco23_2")
df_17_n<- handle_input_fet(fet_foco17_n, sample_set="foco")
df_23_n<- handle_input_fet(fet_foco23_n, sample_set="foco")

#####Change the last column to a name that makes sense
names(df_17_1)[9] <- "p_val"
names(df_17_2)[9] <- "p_val"
names(df_23_1)[9] <- "p_val"
names(df_23_2)[9] <- "p_val"
names(df_wild)[9] <- "p_val"
names(df_23_n)[9] <- "p_val"
names(df_17_n)[9] <- "p_val"
####add column called 'pool'

df_17_1$pool<- "Foco17 Lab 1"
df_17_2$pool<-"Foco17 Lab 2"
df_23_1$pool<- "Foco23 Lab 1"
df_23_2$pool<-"Foco23 Lab 2"
df_wild$pool<-"Foco23 Wild"
df_17_n$pool<- "Foco17 Negative"
df_23_n$pool<- "Foco23 Negative"


df_combined <- rbind(df_17_1, df_17_2, df_23_1, df_23_2, df_wild, df_23_n, df_17_n)
df_no_neg<-rbind(df_17_1, df_17_2, df_23_1, df_23_2, df_wild)

#####The data imported above is the -log10 of the pval by default

##########Full genome plot of p-values################
df_combined<- filter(df_combined, geo_mean_pval>0)
df_combined<- subset(df_combined, chrom!="Y")

a<-1
b<-1.5

############################################################
#######p-value with significance lines
subset_window<-subset(df_combined, geo_mean_pval>1.5)

ggplot(filter(df_combined)) +
  geom_point(aes(x=window, y=geo_mean_pval, fill=pool), shape=21, size=2,stroke=NA) +scale_fill_manual(values = c("cyan", "cyan3", "grey", "coral", "coral3", "darkgrey", "deeppink"))+
  geom_point(data=subset_window, aes(x=window, y=geo_mean_pval), color="black")+
  theme_bw(base_size = 10) +
  ylab("") +
  xlab("Position in genome") +
  facet_grid(pool~chrom, scales = "free")+
  ggtitle("Windowed p-value across the genome")+theme(panel.spacing.x = unit(0, "lines"))+
  geom_hline(yintercept= a, linetype="dashed", color = 'orange')+geom_hline(yintercept= b, linetype="dashed", color = 'magenta')+
  theme(axis.text.x = element_text(angle = 90))

