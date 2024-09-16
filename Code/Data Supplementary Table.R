library(phylotools)
library(tidyverse)

northamerica <- read.fasta("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/prop/usa_canada_rabies_sequences_remove_subset_prop.aligned.trim.fasta") %>% select(seq.name)
ctriver <- read.fasta("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/analysis2.align.trim.fasta") %>% select(seq.name)
ct <- read.fasta("C:/Users/gev25289/Desktop/rabies/Fasta/connecticut/connecticut_rabies_sequences.align.trim.fasta") %>% select(seq.name)

head(ctriver$seq.name)
na <- northamerica %>% separate_wider_delim(seq.name, "|", names = c("Accession", "Host", "Location", "Collection Date"))
ctriv <- ctriver %>% separate_wider_delim(seq.name, "|", names = c("Accession", "Location", "Host", "Clade",  "Collection Date"))
ct <- ct %>% separate_wider_delim(seq.name, "|", names = c("Accession", "Host", "A", "Location", "X", "Y", "Clade", "Collection Date"))

ctriv <- ctriv %>% select("Accession", "Host", "Location", "Collection Date")
ct <- ct %>% select("Accession", "Host", "Location", "Collection Date")

ctriv$Location <- gsub("_West", "", ctriv$Location)
ctriv$Location <- gsub("_East", "", ctriv$Location)
ct$Location <- gsub(".West", "", ct$Location)
ct$Location <- gsub(".East", "", ct$Location)
ct$Location <- paste("USA:CT - ", ct$Location)

data <- bind_rows(ct, ctriv, na)
data$Host <- gsub("Raccoons", "Procyonlotor", data$Host)
data <- data %>%  distinct(Accession, .keep_all = TRUE)
write.csv(data, "C:/Users/gev25289/Desktop/rabies/data.csv", row.names = FALSE)
