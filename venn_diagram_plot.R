install.packages("ggVennDiagram")
library(ggVennDiagram)
####total # SNPs in each pool
table(df_combined_group$pool)

####make pos and chrom a single variable instead of two columns
group_pos<-df_combined_group %>% 
  mutate(chrom_pos = paste(chrom, pos, sep= ":"))%>%
  distinct()

###split each pool into a list

venn_list<-group_pos %>% 
  group_by(pool) %>%
  summarise(positions = list(unique(chrom_pos)))%>%
  deframe()



ggVennDiagram(venn_list, label_alpha = 0, label= "count", set_size = 5)+scale_fill_distiller(palette = "Spectral")+
  scale_x_continuous(expand=expansion(mult=.2))+  guides(fill = "none", color = "none") +
  theme(legend.position = "none")
