library(readxl)
library(binom)
library(dplyr)
library(scales)
library(ggplot2)
library(patchwork)

###total population data for DGRP-517, Lab23, Wild23
GWAS<-(read_excel("./GWAS_pops.xlsx"))

outline <- c("No" = "orchid3", "Yes" = "black")

status_palette <- c("Positive" = "gray", "Negative" = "orangered2")

inf_stat<-ggplot(GWAS, aes(x = Status, y = dCt, group=Population)) +
  geom_jitter(
    aes(fill = Status, color = Used),
    shape = 21, size = 2.5, stroke = 1, width=.5
  ) +
  scale_fill_manual(values = status_palette) +
  scale_color_manual(values = outline) +
  facet_wrap(~Population, scales = "free_x") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12, face = "bold"))+
coord_cartesian(clip="off")+
  scale_y_log10()+ylab("Galbut virus RNA relative to RpL32 mRNA")+xlab("Infection Status")+
  labs(color= "Used in pool", fill= "Infection status")





###3Foco17 prevalence plot
foco1<-read_excel("./foco17_prevalence.xlsx")

ci <- binom.confint(foco1$Positive, foco1$Total, method = "wilson")


foco1$Lower_CI <- ci$lower
foco1$Upper_CI <- ci$upper
foco_over_time <- foco_over_time %>%
  left_join(foco1 %>% select(Year, Lower_CI, Upper_CI), by = "Year")


foco <- ggplot(data = foco1, aes(x = Year, y = Percentage)) +
  geom_line(color = "cyan4") +
  geom_point(color = "cyan4") +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, color = "cyan4") +
  scale_x_continuous(breaks = c(2017, 2018, 2019, 2021, 2022, 2023, 2024)) +
  ylab("Proportion Infected") +
  xlab("Year Screened") +
  theme_bw() +
  theme(legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"))+ facet_wrap(~All)

foco

foco/inf_stat+ plot_annotation(tag_levels = "A")&
  theme(legend.position = "bottom")

