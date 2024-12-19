
####Input FST files, separate into multiple columns
fst_foco17_1<-read.delim("./Foco17_1fst.csv", sep= ",", header=T)
fst_foco17_2<-read.delim("./Foco17_2fst.csv", sep= ",", header=T)
fst_foco23_1<-read.delim("./Foco23_1fst.csv", sep= ",", header=T)
fst_foco23_2<-read.delim("./Foco23_2fst.csv", sep= ",", header=T)
fst_a2<-read.delim("./A2fst.csv", sep= ",", header=T)
fst_foco23_n<-read.delim("./Foco23_Negatives_fst.csv", sep= ",", header=T)
fst_foco17_n<-read.delim("./Foco17_Negatives_fst.csv", sep= ",", header=T)
#####Change the last column to a name that makes sense
names(fst_foco17_1)[5] <- "fst"
names(fst_foco17_2)[5] <- "fst"
names(fst_foco23_1)[5] <- "fst"
names(fst_foco23_2)[5] <- "fst"
names(fst_a2)[5] <- "fst"
names(fst_foco17_n)[5]<- "fst"
names(fst_foco23_n)[5]<-"fst"

######Add new column to say which pool it is
fst_foco17_1$pool<- "Foco17 Lab 1"
fst_foco17_2$pool<-"Foco17 Lab 2"
fst_foco23_1$pool<- "Foco23 Lab 1"
fst_foco23_2$pool<-"Foco23 Lab 2"
fst_a2$pool<-"Foco23 Wild"
fst_foco17_n$pool<-"Foco17 Negative"
fst_foco23_n$pool<-"Foco23 Negative"

df_fst<-rbind(fst_a2, fst_foco17_1, fst_foco17_2, fst_foco23_1, fst_foco23_2, fst_foco17_n, fst_foco23_n)
df_fst_no_neg<-rbind(fst_a2, fst_foco17_1, fst_foco17_2, fst_foco23_1, fst_foco23_2)

##########Full genome plot of FST p-values################

#df_fst<- filter(df_fst, fst<0.01)
df_fst<- filter(df_fst, fst>0)
df_fst<- subset(df_fst, chrom!="Y")

df_fst_no_neg<- filter(df_fst_no_neg, fst>0)
df_fst_no_neg<- subset(df_fst_no_neg, chrom!="Y")

############################################################
#######FST p-value with significance lines based on typical FST 'categories' (not diverse, moderately diverse, extremely diverse)
subset_fst<-subset(df_fst, fst>.15)
subset_fst_no_neg<-subset(df_fst_no_neg, fst>.15)

ggplot(filter(df_fst_no_neg)) +
  geom_point(aes(x=end/1e6, y=fst, fill=pool), shape=21, size=2,stroke=NA, show.legend=FALSE) +
  scale_fill_manual(values = c("cyan", "cyan3", "coral", "coral3", "deeppink"))+
  geom_point(data=subset_fst_no_neg, aes(x=end/1e6, y=fst), color="black")+
  theme_bw(base_size = 10) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  ylab("FST of unique 25000bp windows") +
  xlab("Position on chromosome (Mb)") +
  facet_grid(pool~chrom, scales = "free_x")+
  ggtitle("FST across the genome")+theme(panel.spacing.x = unit(0, "lines"))



######new fst
library(tidyverse)



# Plot output of grenendalf fst command

# MDS 4/3/2024



# read in grenedalf output files

# one for wild (A2 vs. negative), one for foco datasets

fst_wild  <- read.delim("Wild_fst.csv", header=T, sep = ",")

fst_foco17_1  <- read.delim("Foco17_1_fst.csv", header=T, sep = ",")

fst_foco17_2  <- read.delim("Foco17_2_fst.csv", header=T, sep = ",")
fst_foco23_1  <- read.delim("Foco23_1_fst.csv", header=T, sep = ",")
fst_foco23_2  <- read.delim("Foco23_2_fst.csv", header=T, sep = ",")


# this function tidies input

# the grenedalf output columns provide detailed information

# see here for more info:

# https://github.com/lczech/grenedalf/wiki/Output

handle_fst_input <- function(input_df) {
  
  
  
  # for now, don't consider total value columns
  
  # could include these if wanted
  
  df <- input_df %>% select(!starts_with("total"))
  
  
  
  # handle "missing" columns:
  
  df_missing <- df %>%
    
    select(chrom, start, end, ends_with(".missing")) %>%
    
    pivot_longer(cols=-c(chrom, start, end),
                 
                 names_pattern = "(.*).missing",
                 
                 names_to = "p1_p2",
                 
                 values_to = "missing")
  
  
  
  # handle "numeric" columns:
  
  df_numeric <- df %>%
    
    select(chrom, start, end, ends_with(".numeric")) %>%
    
    pivot_longer(cols=-c(chrom, start, end),
                 
                 names_pattern = "(.*).numeric",
                 
                 names_to = "p1_p2",
                 
                 values_to = "numeric")
  
  
  
  # handle "passed" columns:
  
  df_passed <- df %>%
    
    select(chrom, start, end, ends_with(".passed")) %>%
    
    pivot_longer(cols=-c(chrom, start, end),
                 
                 names_pattern = "(.*).passed",
                 
                 names_to = "p1_p2",
                 
                 values_to = "passed")
  
  
  
  # handle "fst" columns:
  
  df_fst <- df %>%
    
    select(chrom, start, end, ends_with(".fst")) %>%
    
    pivot_longer(cols=-c(chrom, start, end),
                 
                 names_pattern = "(.*).fst",
                 
                 names_to = "p1_p2",
                 
                 values_to = "fst")
  
  
  
  # join everything together
  
  df <- left_join(df_missing, df_numeric)
  
  df <- left_join(df,     df_passed)
  
  df <- left_join(df,     df_fst)
  
  
  
  # pull out p1 & p2 separately
  
  p1_p2 <- str_match(df$p1_p2, "([A-Za-z0-9]+).([A-Za-z0-9]+)")
  
  
  
  df$p1 <- p1_p2[,2]
  
  df$p2 <- p1_p2[,3]
  
  
  
  # don't keep windows with NA Fst values
  
  df <- df %>% filter(!is.na(fst))
  
  
  
  # require that windows have >= 50 passed SNPs
  
  df <- df %>% filter(passed >= 50)
  
}



# tidy input dataframes

df_fst_wild <- handle_fst_input(fst_wild)
df_fst_foco17_1 <- handle_fst_input(fst_foco17_1)
df_fst_foco17_2 <- handle_fst_input(fst_foco17_2)
df_fst_foco23_1 <- handle_fst_input(fst_foco23_1)
df_fst_foco23_2 <- handle_fst_input(fst_foco23_2)

df_all<-rbind(df_fst_wild, df_fst_foco17_1,df_fst_foco17_2, df_fst_foco23_1, df_fst_foco23_2)

subset_fst<-subset(df_all, fst>.15)

ggplot(filter(df_all)) +
  geom_point(aes(x=end/1e6, y=fst, fill=p1), shape=21, size=2,stroke=NA, show.legend=FALSE) +
  scale_fill_manual(values = c( "deeppink", "cyan", "cyan3", "coral", "coral3"))+
  geom_point(data=subset_fst, aes(x=end/1e6, y=fst), color="black")+
  theme_bw(base_size = 10) +
  theme(plot.background = element_blank(),

        panel.grid.minor = element_blank())+
  ylab("FST of 25000bp windows") +
  xlab("Position on chromosome (Mb)") +
  facet_grid(p1~chrom, scales = "free_x")+
  ggtitle("FST across the genome")+theme(panel.spacing.x = unit(0, "lines"))
