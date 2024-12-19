######Foco23 lab and Foco23 wild plots
library(readxl)

lab_23<- read_excel("./23_lab.xlsx")
wild_23<-read_excel("./wild_23.xlsx")



lab<-ggplot(data=lab_23, aes(x= Title, y=dCt))+geom_jitter( width=0.3, height=0, aes(color=Status), show.legend=FALSE)+xlab("")+theme(axis.text.x = element_blank())+ggtitle("")+
  ylab("")+ facet_wrap(~Status, nrow=1)+scale_color_manual(values = c("coral", "coral3", "black"))+
  theme_bw()+xlab("")+scale_y_continuous(trans="log10", limits=c(.001,2000))+theme(axis.title.x=element_blank(),
                                                                                  axis.text.x=element_blank(),
                                                                                  axis.ticks.x=element_blank())+coord_cartesian(clip=FALSE)

wild<-ggplot(data=wild_23, aes(x= Title, y=dCt, color= Status))+geom_jitter(width=0.3, height=0, aes(color=Status), show.legend=FALSE)+xlab("")+theme(axis.text.x = element_blank())+ggtitle("")+
  scale_y_continuous(trans="log10", limits=c(.001,2000))+facet_wrap(~Status, nrow=1)+scale_color_manual(values=c("deeppink", "black"))+
  theme_bw()+ylab("Galbut virus RNA relative to RpL32 (log10)")+theme(axis.title.x=element_blank(),
                            axis.text.x=element_blank(),
                            axis.text=element_text(size=8),
                            axis.ticks.x=element_blank())+coord_cartesian(clip=FALSE)



######Foco17 plot of prevalence over time
library(scales)
foco_over_time <- read.delim ("./foco17_prop.txt", sep="\t", header=T)
foco_over_time<- foco_over_time[-c(3,4,5,6,7,8,9)]
foco_over_time<- foco_over_time[-c(6,7,8,9,10,11,12),]


foco<-ggplot(foco_over_time, aes(factor(Year), Infected)) +
  stat_summary(fun=mean, geom="bar", fill="grey70") +
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.2) +
  scale_y_continuous(labels=percent_format(), limits=c(0,1)) + xlab("Year")+ ylab("Percent infected")+
  theme_bw()+ ggtitle("")


foco<- ggplot(data=foco_over_time, aes(x=Year, y=Prop, color = "red"))+geom_line()+geom_point()+
  scale_x_continuous(breaks=c(2017,2019,2021,2022, 2023, 2024))+ylab("Proportion Infected")+ xlab("Year Screened")+
  theme(legend.position = "none")+ ggtitle("Proportion infected from 2017 wild colony reared in the lab over time")


foco

foco/ (wild+lab ) + plot_layout(heights=c(1.5, 1.6))


