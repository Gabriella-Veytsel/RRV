library(ggplot2)
library(dplyr)
library(tidytree)
library(reshape2)
library(tidyr)
library(lubridate)
library(ggpubr)
library(scales)
library(tidyverse)
library(cowplot)

#I adapted some of this code from Jiani's code on github: https://github.com/JianiC/ATL-flyway/tree/main/Analytical_Scripts/Markov_jumpshp.R
#Connecticut River
##################
jumps <- read.delim("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/relaxed_discrete_location/jumpTimes.txt",sep = '', header = TRUE)
state <- n_distinct(jumps$state) #Number of total state counts
jumps$from_to = paste(jumps$from,jumps$to, sep=".")
jumps <- jumps %>% mutate(time = 2022.1972602739727 - as.numeric(time))
jumps$year <- format(date_decimal(jumps$time), "%Y")

#From, To, Year
count<-jumps %>% group_by(from_to,year)%>% count()
count2<-cbind(count, read.table(text = as.character(count$from_to), sep = "."))

#From, To
count_total <- jumps %>% group_by(from_to) %>% count()
count_total <- count_total %>% separate_wider_delim(from_to, ".", names = c("From", "To"))
count_total <- count_total %>%
  mutate(F = case_when(
    From %in% c("USA:CT_East", "USA:ME") ~ "East",
    From %in% c("USA:CT_West", "USA:NY", "USA:VT") ~ "West"
  ))

count_total <- count_total %>%
  mutate(T = case_when(
    To %in% c("USA:CT_East", "USA:ME") ~ "East",
    To %in% c("USA:CT_West", "USA:NY", "USA:VT") ~ "West"
  ))

count_total$from_to = paste(count_total$F,count_total$T, sep=".")
river <- count_total %>%
  group_by(from_to) %>%
  summarize(total = sum(n))
river <- river %>% mutate(ave=total/state)

count_ave <- count_total %>% mutate(ave=n/state)
summary(count_ave$ave)
ggplot(count_ave, aes(x=To, y=From, fill= ave)) +
  theme_bw()+
  scale_fill_gradient2(low = 'wheat', mid = 'papayawhip', high = 'skyblue4',midpoint = 0.115, limits=c(0,6)) +
  xlab("Sink") +
  ylab("Source") +
  ggtitle("Average Markov Jump Counts") +
  geom_tile() +
  guides(fill = guide_colourbar(title = "Average Markov\nJump Counts"))
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/figures/heatmap_jumps_average_location.pdf", width=7, height=4)

#North America
##############
jumps <- read.delim("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/relaxed/jumpTimes.txt",sep = '', header = TRUE)
state <- n_distinct(jumps$state) #Number of total state counts
jumps$from_to = paste(jumps$from,jumps$to, sep=".")
jumps <- jumps %>% mutate(time = 2021.9698630136986 - as.numeric(time))
jumps$year <- format(date_decimal(jumps$time), "%Y")

#From, To
count_total <- jumps %>% group_by(from_to) %>% count()
count_total <- count_total %>% separate_wider_delim(from_to, ".", names = c("From", "To"))
count_total <- count_total %>% mutate(ave=n/state)
summary(count_total$ave)
ggplot(count_total, aes(x=To, y=From, fill= ave)) +
  theme_bw()+
  scale_fill_gradient2(low = 'wheat', mid = 'papayawhip', high = 'skyblue4',midpoint = 0.705586, limits=c(0,5)) +
  xlab("Sink") +
  ylab("Source") +
  ggtitle("Average Markov Jump Counts") +
  geom_tile() +
  guides(fill = guide_colourbar(title = "Average Markov\nJump Counts"))
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/heatmap_jumps_average_location.pdf", width=10, height=4)

#From, To, Year
count<-jumps %>% group_by(from_to,year)%>% count()
count2<-cbind(count, read.table(text = as.character(count$from_to), sep = "."))
count2 <- count2 %>% mutate(ave=n/state)

summary(count2$ave)
ggplot(count2, aes(x = year, y = V1, fill= ave)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
  ggtitle("") +
  #theme(panel.background = element_blank()) +
  #theme_bw(base_size=10)+
  scale_fill_gradient2(low = 'thistle', mid = 'papayawhip', high = 'purple4',midpoint = 0.0001644) +
  labs(fill="Average Markov\nJump Counts")+
  xlab("Year") +
  ylab("Source") +
  facet_grid(rows = vars(V2))
  #scale_x_continuous(expand=c(0, 0),position = "top")
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/heatmap_jumps_average_location_time.pdf", width=15, height=8)


#Heatmap 
# count2$year = as.numeric(count2$year)
# summary(count2$n)
# p<-ggplot(count2, aes(x = year, y = V1, fill= n)) + 
#   geom_tile() +
#   ggtitle("") +
#   theme_bw(base_size=10)+
#   scale_fill_gradient2(low = 'thistle', mid = 'papayawhip', high = 'blue',midpoint = 495, limits=c(0,20000)) +
#   labs(fill="Total Markov\nJump Counts")+
#   theme_bw(base_size=13) +
#   xlab("Year") +
#   ylab("Source") + 
#   facet_grid(rows = vars(V2)) +
#   scale_x_continuous(expand=c(0, 0),position = "top")
# p
# ggsave("heatmap_jumps_total.pdf", width=6, height=10)



# plot_grid(p,p1,labels = c('B', 'C'),  ncol = 2, nrow=1, rel_widths = c(1,1))
ggsave("D:/rabies/Results/State_Markov_BSSVS/jumphistory/totalandaveragejumps.pdf", width=17, height=13)
