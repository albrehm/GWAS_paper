####Input FST files, separate into multiple columns
foco17_5000_1<-read.delim("./Foco17_1_5000fst.csv", sep= ",", header=T)
foco17_5000_2<-read.delim("./Foco17_2_5000fst.csv", sep= ",", header=T)
foco23_5000_1<-read.delim("./Foco23_1_5000fst.csv", sep= ",", header=T)
foco23_5000_2<-read.delim("./Foco23_2_5000fst.csv", sep= ",", header=T)
a2_5000<-read.delim("./A2_1_5000fst.csv", sep= ",", header=T)

#####Change the last column to a name that makes sense
names(foco17_5000_1)[5] <- "p_val"
names(foco17_5000_2)[5] <- "p_val"
names(foco23_5000_1)[5] <- "p_val"
names(foco23_5000_2)[5] <- "p_val"
names(a2_5000)[5] <- "p_val"

######Add new column to say which pool it is
foco17_5000_1$pool<- "Foco17 Pool 1"
foco17_5000_2$pool<-"Foco17 Pool 2"
foco23_5000_1$pool<- "Foco23 Pool 1"
foco23_5000_2$pool<-"Foco23 Pool 2"
a2_5000$pool<-"Wild"

df_5000<-rbind(foco17_5000_1, foco17_5000_2, foco23_5000_1, foco23_5000_2, a2_5000)

df_5000<-na.omit(df_5000)

#####cutoff by pval#######
p_val_cutoffs_5000<- df_5000 %>%
  group_by(pool)%>%
  reframe(p_val=quantile(p_val,c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)), q=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) %>%
  filter(q==0.95)%>% 
  rename(log_pval_cutoff=p_val)

df_5000<- left_join(df_5000,p_val_cutoffs_5000)

above_cutoff_5000<- filter(df_5000, p_val > log_pval_cutoff)

above_cutoff_5000<-above_cutoff_5000%>%
  group_by(chrom, start)%>%
  mutate(n_groups= n())



df_5000 <- df_5000 %>% mutate(p_val = -log10(p_val))

ggplot(filter(df_5000)) +
  geom_point(aes(x=end, y=p_val, fill=pool), shape=21, size=2) +
  theme_bw(base_size = 10) +
  ylab("-log10 FST p-value") +
  xlab("position in genome") +
  facet_wrap(pool~chrom, scales="free")+
  ggtitle("-log10 FST p-value across the genome")+theme(panel.spacing.x = unit(0, "lines"))

