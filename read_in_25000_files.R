##set window size
window_size <- 25000 
window_size_mb <- window_size / 1e6 

###helper code to make raw input files a better format
handle_input_fet <- function(df_to_handle, sample_set="wild") {
  
  # create column names (popoolation2 doesn't output column names)
  extra_col_names = paste0("X", seq(1,ncol(df_to_handle) - 5))
  colnames(df_to_handle) <- c("chrom", "window", "n_snps", "frac_sufficient_cov", "avg_min_coverage", extra_col_names) 
  
  df <- df_to_handle %>% pivot_longer(cols=starts_with("X"), names_to="unhelpful_name", values_to="p_p")  %>% select(-unhelpful_name)
  
  # chr2L	5000	1	0.002	59.3	1:2=2.29269731	1:3=2.91987693	1:4=9.43200886	1:5=33.37451954	1:6=0.94702188	1:7=1.51466859	1:8=0.17907030	2:3=0.50657442	2:4=7.84675129	2:5=70.85642166	2:6=1.11831056	2:7=0.53391643	2:8=4.34439480	3:4=5.34500299	3:5=43.77295891	3:6=1.88328402	3:7=1.24666426	3:8=5.12927387	4:5=5.51537397	4:6=9.84815294	4:7=9.47870196	4:8=14.57806963	5:6=59.53783747	5:7=75.64739283	5:8=83.47553973	6:7=0.38735115	6:8=1.09731604	7:8=2.61145212
  
  # parse metadata out of p-value value
  subparts <- str_match(df$p_p, "(\\d+):(\\d+)=([0-9\\.]+)")
  
  df$pop_1 <- as.character(subparts[,2])
  df$pop_2 <- as.character(subparts[,3])
  df$log_pval <- as.numeric(subparts[,4])
  
  # populations
  if (sample_set == "wild") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "A2",
                           "2" ~ "Neg")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "A2",
                           "2" ~ "Neg")
    
  } else if (sample_set == "foco") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Negative2",
                           "3" ~ "Foco17Positive",
                           "4" ~ "Foco23Negative1",
                           "5" ~ "Foco23Negative2",
                           "6" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Negative2",
                           "3" ~ "Foco17Positive",
                           "4" ~ "Foco23Negative1",
                           "5" ~ "Foco23Negative2",
                           "6" ~ "Foco23Positive")
  } else if (sample_set == "foco17_1") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative1",
                           "2" ~ "Foco17Positive")
    
  } else if (sample_set == "foco17_2") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco17Negative2",
                           "2" ~ "Foco17Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco17Negative2",
                           "2" ~ "Foco17Positive")
  } else if (sample_set == "foco23_1") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco23Negative1",
                           "2" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco23Negative1",
                           "2" ~ "Foco23Positive")
    
  } else if (sample_set == "foco23_2") {
    
    df$pop_1 <- case_match(df$pop_1,
                           "1" ~ "Foco23Negative2",
                           "2" ~ "Foco23Positive")
    
    df$pop_2 <- case_match(df$pop_2,
                           "1" ~ "Foco23Negative2",
                           "2" ~ "Foco23Positive")
  }
  
  df <- df %>% mutate(pop_1_2 = paste0(pop_1, "-", pop_2))
  
  # keep track of position in megabases
  df <- df %>% mutate(window_mb = window / 1e6)
  
  df <- df %>% mutate(window_start = window - window_size / 2,
                      window_end   = window + window_size / 2)
  
  # what does distribution of fraction sufficient coverage look like?
  # ggplot(df) + 
  # geom_histogram(aes(x=frac_sufficient_cov), bins=100)  +
  # theme_bw()
  
  # filter out windows with insufficient coverage
  df <- df %>% filter(frac_sufficient_cov > 0.05)
  
  # ggplot(df) + 
  # geom_histogram(aes(x=frac_sufficient_cov), bins=100)  +
  # theme_bw()
  
  # only keep big chromosomes
  df <- filter(df, str_detect(chrom, "[23X]"))
  
  # calculate geometric mean of p-values
  # which is the nth root of the product of p-values 
  # (the product is what popoolation2 reports by default)
  # n = the # of SNPs 
  df <- df %>% mutate(geo_mean_pval = log_pval ^ (1/n_snps),
                      mean_pval = log_pval / n_snps)
  
  # return the df
  df
}

#####Input the files
fet_wild <- read.delim ("./wild2.fet", sep="\t", header=F)
fet_foco17_1<- read.delim("./Foco17_1.fet", sep= "\t", header=F)
fet_foco17_2<- read.delim("./Foco17_2.fet", sep= "\t", header=F)
fet_foco23_1<- read.delim("./Foco23_1.fet", sep= "\t", header=F)
fet_foco23_2<- read.delim("./Foco23_2.fet", sep= "\t", header=F)
###Make all the input files the correct format
df_wild <- handle_input_fet(fet_wild, sample_set="wild")
df_17_1 <- handle_input_fet(fet_foco17_1, sample_set="foco17_1")
df_17_2 <- handle_input_fet(fet_foco17_2, sample_set="foco17_2")
df_23_1 <- handle_input_fet(fet_foco23_1, sample_set="foco23_1")
df_23_2 <- handle_input_fet(fet_foco23_2, sample_set="foco23_2")


####add column called 'pool'

df_17_1$pool<- "Foco17 Pool 1"
df_17_2$pool<-"Foco17 Pool 2"
df_23_1$pool<- "Foco23 Pool 1"
df_23_2$pool<-"Foco23 Pool 2"
df_wild$pool<-"Wild"

df_combined <- rbind(df_17_1, df_17_2, df_23_1, df_23_2, df_wild)
