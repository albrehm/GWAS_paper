library(tidyverse)
library(stringr)
# Plot output of grenendalf diversity command

# input data from wild (A2 vs negative) & DGRP (F2->F7) GWASes
div_wild <- read.delim("grenedalf.wild.diversity.csv", sep=",", header=T)
div_dgrp <- read.delim("grenedalf.dgrp.diversity.csv", sep=",", header=T)
div_foco <- read.delim("grenedalf.foco.diversity.csv", sep=",", header=T)


# Function to tidy the diversity data
handle_diversity_input <- function(df_to_handle) {
  df_long <- df_to_handle %>%
    pivot_longer(
      cols = -c(chrom, start, end),
      names_to = "sample_metric",
      values_to = "val"
    )
  
  # Split sample name from metric
  parsed <- str_match(df_long$sample_metric, "^(.+)\\.(.+)$")
  
  df_long <- df_long %>%
    mutate(
      sample = parsed[, 2],
      metric = parsed[, 3]
    ) %>%
    select(-sample_metric)
  
  # Spread metrics into columns
  df_wide <- df_long %>%
    pivot_wider(names_from = metric, values_from = val)
  
  # Drop "total" pseudo-sample
  df_wide <- df_wide %>% filter(sample != "total")
  
  # Make metric columns numeric
  metric_cols <- setdiff(names(df_wide), c("chrom", "start", "end", "sample"))
  df_wide[metric_cols] <- lapply(df_wide[metric_cols], function(x) as.numeric(as.character(x)))
  
  # Keep only big chromosomes
  df_wide <- df_wide %>% filter(str_detect(chrom, "^[23XY]"))
  
  return(df_wide)
}

# Apply to datasets
df_div_wild <- handle_diversity_input(div_wild) %>% mutate(dataset = "wild")
df_div_dgrp <- handle_diversity_input(div_dgrp) %>% mutate(dataset = "dgrp")
df_div_foco <- handle_diversity_input(div_foco) %>% mutate(dataset = "foco")

# Combine all
df_all <- bind_rows(df_div_wild, df_div_dgrp, df_div_foco)

# Define infected sample names (exact match, case-sensitive)
infected_samples <- c("X25197.er.infected", "Foco17Positive_merged", "Foco23Positive", "A2_redo")

# Tag infection status
df_all <- df_all %>%
  mutate(
    condition = if_else(sample %in% infected_samples, "infected", "uninfected")
  )




df_all$sample <- factor(
  df_all$sample,
  levels = c("A2_redo", "Negative_combined", "Foco17Negative_merged",
             "Foco17Positive_merged", "Foco23Positive", "Foco23Negative_merged",
             "X25197.er.infected", "X25197.er.uninfected"),
  labels = c("FoCo23 Wild Positive", "FoCo23 Wild Negative",
             "FoCo17 Negative", "FoCo17 Positive", "FoCo23 Lab Positive",
             "FoCo23 Lab Negative", "DGRP-517 Positive", "DGRP-517 Negative"))


# Plot infected vs uninfected (assuming 'sample' encodes that info)
ggplot(df_all, aes(x = sample, y = theta_pi, fill=sample)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.4) +
  facet_wrap(~chrom) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +
  labs(
    x = "",
    y = "Nucleotide diversity"
  )+
  scale_fill_manual(values = c("FoCo23 Wild Positive" = "deeppink3", "FoCo23 Wild Negative" = "deeppink", "FoCo17 Negative" = "blue", "DGRP-517 Positive" = "green3",
                                 "FoCo17 Positive" ="blue3", "DGRP-517 Negative" = "green", "FoCo23 Lab Positive" = "coral3", "FoCo23 Lab Negative"= "coral" ))

