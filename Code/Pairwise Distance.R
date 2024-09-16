library(reshape2)
library(seqinr)
library(tidyverse)
library(ggplot2)
library(pegas)
library(cowplot)

dna <- read.dna("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/connecticut_rabies_sequences.align.trim.fasta", format = "fasta")
dna.d = 100*(dist.dna (dna, model='raw', pairwise.deletion = F, as.matrix=T)) 
dna.d.melt <- melt(dna.d)

summary(dna.d.melt$value)
dna.d.melt$valueCAT<- cut(dna.d.melt$value, breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2), right = FALSE)
dna.d.melt$name1 <- sub("\\|.*", "", dna.d.melt$Var1)
dna.d.melt$name2 <- sub("\\|.*", "", dna.d.melt$Var2)

# Define common color scale
color_scale <- scale_fill_manual(values = c("[0,0.25)" = '#f0f9e8',
                                            "[0.25,0.5)" = '#a8ddb5',
                                            "[0.5,0.75)" = '#238b45',
                                            "[0.75,1)" = '#4eb3d3',
                                            "[1,1.25)" =  '#2b8cbe',
                                            "[1.25,1.5)" =  '#08589e',
                                            "[1.5,1.75)" = '#c994c7',
                                            "[1.75,2.0)" = 'red'), name = 'Sequence\ndissimilarity, %')

# Function to create plots
create_plot <- function(data, title) {
  ggplot(data, aes(x = name1, y = name2)) +
    geom_tile(aes(fill = valueCAT), color = 'black') +
    color_scale +
    theme(axis.text.x = element_text(angle = 90, size = 6), axis.text.y = element_text(size = 6)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle(title)
}

# Subset data for each section
section1 <- dna.d.melt %>%
  filter(grepl("Section1", Var1)) %>%
  filter(grepl("Section1", Var2))

section2 <- dna.d.melt %>%
  filter(grepl("Section2", Var1)) %>%
  filter(grepl("Section2", Var2))

section3 <- dna.d.melt %>%
  filter(grepl("Section3", Var1)) %>%
  filter(grepl("Section3", Var2))

section1.2 <- dna.d.melt %>%
  filter(grepl("Section1", Var1)) %>%
  filter(grepl("Section2", Var2))

section2.3 <- dna.d.melt %>%
  filter(grepl("Section2", Var1)) %>%
  filter(grepl("Section3", Var2))

section1.3 <- dna.d.melt %>%
  filter(grepl("Section1", Var1)) %>%
  filter(grepl("Section3", Var2))

# Create plots
s1 <- create_plot(section1, "S1")
s2 <- create_plot(section2, "S2")
s3 <- create_plot(section3, "S3")
s1.2 <- create_plot(section1.2, "S1/S2")
s2.3 <- create_plot(section2.3, "S2/3")
s1.3 <- create_plot(section1.3, "S1/3")

# Combine plots into a grid
p <- plot_grid(s1, s2, s3, s1.2, s2.3, s1.3, ncol = 3, nrow = 2)

# Print the grid
print(p)
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/pairwisedistance_sections.pdf", width=20, height=7)

############################################################################
# Create plots
c1 <- create_plot(clade1, "Clade 1")
c2 <- create_plot(clade2, "Clade 2")
c1.2 <- create_plot(clade1.2, "Clade 1/2")

# Combine plots into a grid
p1 <- plot_grid(c1, c2, c1.2, ncol = 3, nrow = 1)

# Print the grid
print(p1)
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/pairwisedistance_clades.pdf", width=20, height=7)

summary(section1$value) #mean = 0.5
summary(section2$value) #mean = 1.1
summary(section3$value) #mean = 0.69

summary(section1.2$value) #mean = 1.3
summary(section2.3$value) #mean = 1.0
summary(section1.3$value) #mean = 1.5

summary(clade1$value) #mean = 0.5
summary(clade2$value) #mean = 0.77
summary(clade1.2$value) #mean = 1.5
summary(betweenClades$value) #mean = 1.50, range = 1.28, 1.80