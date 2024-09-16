library(ggplot2)
library(ape)
library(treeio)
library(ggtree)
library(phylotools)
library(phytools)
library(dplyr)
library(tidytree)
library(reshape2)
library(tidyr)
library(lubridate)
library(ggpubr)
library(scales)
library(tidyverse)
library(cowplot)

#River
####################################################################################################################################
tree <- read.beast("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/relaxed_discrete_location/combined.annot.tree") 
tree@phylo$edge.length[tree@phylo$edge.length==0] = tree@phylo$edge.length[tree@phylo$edge.length==0] + 1e-8 
tree@phylo <- di2multi(tree@phylo, tol=0)

tree_tibble <- tree %>% as_tibble()
tree_tibble$Location <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", tree_tibble$label)

table(tree_tibble$Location)

p <- ggtree(tree@phylo, mrsd="2022-03-14") + theme_tree2() 
p %<+% tree_tibble + 
  geom_tippoint(colour="black", shape=21, size = 3, aes(fill = Location)) +
  scale_fill_manual(values=c("purple", "deeppink", "navy", "dodgerblue", "gold3")) 
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/relaxed_discrete_location/tree.pdf", width=5, height=9)

tree@phylo$tip.label <- sub("^[^|]*\\|([^|]*)\\|.*$", "\\1", tree@phylo$tip.label)

p + geom_text(aes(label=node), hjust=-.3)
p + 
  geom_hilight(node=223, fill="pink") +
  geom_hilight(node=361, fill="gold") +
  geom_hilight(node=438, fill="lightblue")

tipclade <- tips(tree@phylo, 223) 
table(tipclade)

tipclade <- tips(tree@phylo, 361) 
table(tipclade)

tipclade <- tips(tree@phylo, 438) 
table(tipclade)

date_decimal(1978.299)
date_decimal(1974.51)
date_decimal(1981.51)

#Connecticut:Location
####################################################################################################################################
tree <- read.beast("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/strict_discrete_location/combined.annot.tree")
tree@phylo$edge.length[tree@phylo$edge.length==0] = tree@phylo$edge.length[tree@phylo$edge.length==0] + 1e-8 
tree@phylo <- di2multi(tree@phylo, tol=0)

tree_tibble <- tree %>% as_tibble()
tree_tibble$Location <- sub("^(?:[^|]*\\|){5}([^|]*)\\|.*$", "\\1", tree_tibble$label)

table(tree_tibble$Location)
p <- ggtree(tree@phylo, mrsd="2021-10-18") + theme_tree2() 
p1 <- p %<+% tree_tibble + 
  geom_tippoint(colour="black", shape=21, size = 4, aes(fill = Location)) +
  scale_fill_manual(values=c("red", "yellow", "blue","red3", "yellow2", 
                             "yellow3", "dodgerblue", "gold","royalblue", "navy", "lightblue")) 

tree@phylo$tip.label <- sub("^(?:[^|]*\\|){6}([^|]*).*", "\\1", tree@phylo$tip.label)
ggtree(tree@phylo, mrsd="2021-10-18") + theme_tree2() + geom_tiplab() #label tips by clade


# Clade: Node
# 2A: 83
# 2B: 105
# 2D: 79
# 2C: 75
# 1A: 125
# 1B: 140
# 1C: 120

p1 + 
  scale_x_continuous(limits = c(1970, 2030)) +
  #theme(plot.margin = margin(0, 0, 0, 0)) +
  geom_cladelab(node=83, label="2A", color="black", offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=105, label="2B", color="black",  offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=79, label="2D", color="black",  offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=75, label="2C", color="black",  offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=125, label="1A", color="black", offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=140, label="1B", color="black",  offset.text = -20, offset = -35, align = TRUE) +
  geom_cladelab(node=120, label="1C", color="black", offset.text = -20, offset = -35, align = TRUE) +
  geom_point2(aes(subset=(node==83)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==105)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==79)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==75)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==125)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==140)), shape=16, size=3, fill='black') +
  geom_point2(aes(subset=(node==120)), shape=16, size=3, fill='black') 
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/tree.pdf", width=7, height=9)

check <- tree@phylo %>% as.treedata() %>% as_tibble()
check <- tree_tibble %>% filter(node %in% c(83,105,79,75,125,140,120))
check <- check %>% mutate(latest = 2021.7945205479452)
check <- check %>%
  mutate(age = latest - as.numeric(height)) %>%
  mutate(lowerhpd =  sub(".*\\(([^,]+).*", "\\1", check$height_0.95_HPD)) %>%
  mutate(higherhpd = sub(".*, ([^)]*)", "\\1", check$height_0.95_HPD)) %>%
  mutate(higherhpd = sub("\\)$", "", higherhpd)) %>%
  mutate(lhpd = latest - as.numeric(lowerhpd)) %>%
  mutate(hhpd = latest - as.numeric(higherhpd)) %>%
  mutate(Age = date_decimal(age)) %>%
  mutate(LHPD = date_decimal(lhpd)) %>%
  mutate(HHPD = date_decimal(hhpd)) 


p + 
  geom_hilight(node=73, fill="pink") +
  geom_hilight(node=119, fill="gold")

tipclade <- tips(tree@phylo, 75) 
table(tipclade)

tipclade <- tips(tree@phylo, 119) %>% as.data.frame()
table(tipclade)

ggtree(tree) + 
  geom_cladelabel(node=17, label="Some random clade", color="red")

ontario_outbreak <- tips(tree@phylo, 577) %>% as.data.frame()
date_decimal(2021.7945205479452-22.6601510646775, tz = "UTC")

#Connecticut:Host
####################################################################################################################################
tree <- read.beast("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/strict_discrete_host/combined.annot.tree")
tree@phylo$edge.length[tree@phylo$edge.length==0] = tree@phylo$edge.length[tree@phylo$edge.length==0] + 1e-8 
tree@phylo <- di2multi(tree@phylo, tol=0)

tree_tibble <- tree %>% as_tibble()
tree_tibble$Host <- sub("^(?:[^|]*\\|){1}([^|]*)\\|.*$", "\\1", tree_tibble$label)

table(tree_tibble$Host)
p <- ggtree(tree@phylo, mrsd="2021-10-18") + theme_tree2() 
p %<+% tree_tibble + 
  geom_tippoint(colour="black", shape=21, size = 4, aes(fill = Host)) +
  scale_fill_manual(values=c("purple", "pink", "deeppink","lightblue",  "navy", "dodgerblue", "magenta", "gold3")) 
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/strict_discrete_host/tree.pdf", width=5, height=9)

#North America tree
####################################################################################################################################
tree <- read.beast("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/relaxed/combined.annot.tree")
tree@phylo$edge.length[tree@phylo$edge.length==0] = tree@phylo$edge.length[tree@phylo$edge.length==0] + 1e-8 
tree@phylo <- di2multi(tree@phylo, tol=0)

tree_tibble <- tree %>% as_tibble()
tree_tibble$Location <- sub("^[^|]*\\|[^|]*\\|([^|]*).*", "\\1", tree_tibble$label)
tree_tibble$Location[tree_tibble$Location %in% c("USA:MA", "USA:RI")] <- "USA:CT"

table(tree_tibble$Location)
p <- ggtree(tree@phylo, mrsd="2022-03-14") + theme_tree2() 
p %<+% tree_tibble + 
  geom_tippoint(colour="black", shape=21, size = 3, aes(fill = Location)) +
  scale_fill_manual(values=c("magenta", "pink", "deeppink", "dodgerblue","blue", "lightblue", "navy"))

p + geom_text(aes(label=node), hjust=-.3)
p + 
  geom_hilight(node=1035, fill="pink") +
  geom_hilight(node=916, fill="lightblue") +
  geom_hilight(node=559, fill="gold")

tree@phylo$tip.label <- sub("^[^|]*\\|[^|]*\\|([^|]*).*", "\\1", tree@phylo$tip.label)
tipclade <- tips(tree@phylo, 1035) %>% as.data.frame()
table(tipclade)

tipclade <- tips(tree@phylo, 916) %>% as.data.frame()
table(tipclade)

ontario_outbreak <- tips(tree@phylo, 577) %>% as.data.frame()
date_decimal(2013.2398, tz = "UTC")

ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/tree.pdf", width=5, height=9)

date_decimal(1978.425 )
date_decimal(1975.2893)
date_decimal(1981.648)

#North America tree - subsampled
####################################################################################################################################
tree <- read.beast("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/combined.annot.tree")
tree@phylo$edge.length[tree@phylo$edge.length==0] = tree@phylo$edge.length[tree@phylo$edge.length==0] + 1e-8 
tree@phylo <- di2multi(tree@phylo, tol=0)

tree_tibble <- tree %>% as_tibble()
tree_tibble$Location <- sub("^[^|]*\\|[^|]*\\|([^|]*).*", "\\1", tree_tibble$label)
tree_tibble$Location[tree_tibble$Location %in% c("USA:MA", "USA:RI")] <- "USA:CT"

table(tree_tibble$Location)
p <- ggtree(tree@phylo, mrsd="2021-12-21") + theme_tree2() 
z = p %<+% tree_tibble + 
  geom_tippoint(colour="black", shape=21, size = 3, aes(fill = Location)) +
  scale_fill_manual(values=c("magenta", "pink", "deeppink", "dodgerblue","blue", "lightblue", "navy")) 

p + geom_text(aes(label=node), hjust=-.3)
 + 
  geom_hilight(node=320, fill="pink") +
  geom_hilight(node=501, fill="lightblue") +
  geom_hilight(node=306, fill="gold")

tree@phylo$tip.label <- sub("^[^|]*\\|[^|]*\\|([^|]*).*", "\\1", tree@phylo$tip.label)
tipclade <- tips(tree@phylo, 320) %>% as.data.frame()
table(tipclade)

tipclade <- tips(tree@phylo, 501) %>% as.data.frame()
table(tipclade)

tipclade <- tips(tree@phylo, 306) %>% as.data.frame()
table(tipclade)

ontario_outbreak <- tips(tree@phylo, 324) %>% as.data.frame()
date_decimal(2021.9698630136986-8.6873, tz = "UTC")

date_decimal(1979.24, tz = "UTC")
date_decimal(1975.8343, tz = "UTC")
date_decimal(1982.6336, tz = "UTC")
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/tree.pdf", width=5, height=9)

