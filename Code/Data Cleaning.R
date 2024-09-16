library(phylotools)
library(tidyverse)
library(tidygeocoder)
library(readxl)
library(ggtree)
library(phytools) #midpoint.root
library(geiger)
library(cowplot)

`%notin%` <- Negate(`%in%`)

#Note: the USA column of NCBI Virus is not trustworthy -- look at Geo_Location

###########
#NCBI VIRUS
###########

#Sequences
north_america_rabies_sequences <- read.fasta("C:/Users/gev25289/Desktop/rabies/Data/northamerica_rabies_sequences.fasta")
ct_rabies_sequences <- north_america_rabies_sequences %>%
  filter(grepl("CT|OR227628|OR227629", seq.name)) 

north_america_rabies_sequences$seq.name <- sub("\\.1.*", "", north_america_rabies_sequences$seq.name)
north_america_rabies_sequences$seq.name <- sub("\\.2.*", "", north_america_rabies_sequences$seq.name)
north_america_rabies_sequences$seq.name <- sub("\\.3.*", "", north_america_rabies_sequences$seq.name)

#Metadata
north_america_rabies_metadata <- read_csv("C:/Users/gev25289/Desktop/rabies/Data/northamerica_rabies_metadata.csv")
ct_rabies_metadata <- north_america_rabies_metadata %>%
  filter(USA == "CT"| Isolate == "19-4017" | Isolate == "18-3755") %>%
  mutate(Isolate = gsub("-","_", Isolate))
ct_rabies_metadata$Isolate <- ifelse(grepl("MN", ct_rabies_metadata$Accession), ct_rabies_metadata$Accession, ct_rabies_metadata$Isolate)

#Missing location info from https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0008113 S2 Table
nadindavis <- read_excel("C:/Users/gev25289/Downloads/pntd.0008113.s004.xlsx")
north_america_rabies_metadata <- north_america_rabies_metadata %>% 
  left_join(nadindavis, by = c("Accession" = "NCBI Accession Number"))
north_america_rabies_metadata$Geo_Location <- ifelse(north_america_rabies_metadata$Geo_Location == "Canada", paste("Canada", north_america_rabies_metadata$`STATE/ PROVINCE`, sep = ":"), north_america_rabies_metadata$Geo_Location)
north_america_rabies_metadata$Geo_Location <- ifelse(north_america_rabies_metadata$Geo_Location == "USA", paste("USA", north_america_rabies_metadata$`STATE/ PROVINCE`, sep = ":"), north_america_rabies_metadata$Geo_Location)

#Cleaning NCBI Virus data
north_america_rabies_fasta <- north_america_rabies_sequences %>% left_join(north_america_rabies_metadata, by = c("seq.name" = "Accession"))
north_america_rabies_fasta$short_Geo_Location <- sub("^([^:]+): .* ([A-Za-z]{2})$", "\\1: \\2", north_america_rabies_fasta$Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("New York", "NY", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Pennsylvania", "PA", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Alabama", "AL", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Alaska", "AK", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Hampden,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Tekax,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Hablekal,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Cancun,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Sinanche,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Washington,", "", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Ontario", "ON", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Newfoundland and Labrador", "NL", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub("Quebec", "QC", north_america_rabies_fasta$short_Geo_Location)
north_america_rabies_fasta$short_Geo_Location <- gsub(" ", "", north_america_rabies_fasta$short_Geo_Location)

north_america_rabies_fasta <- north_america_rabies_fasta %>% filter(Geo_Location %notin% c("Grenada", "Costa Rica"))
north_america_rabies_fasta <- north_america_rabies_fasta %>% filter(!is.na(Collection_Date))
north_america_rabies_fasta$short_Geo_Location <- gsub("Litchfield,|Fairfield,|Hartford,|NewHaven,|NewLondon,|Middlesex,|Tolland,|Windham,", "", north_america_rabies_fasta$short_Geo_Location)

table(north_america_rabies_fasta$short_Geo_Location) 

north_america_rabies_fasta %>% filter(grepl("USA", short_Geo_Location)) %>% summarise(n=n()) #618
north_america_rabies_fasta %>% filter(grepl("Canada", short_Geo_Location)) %>% summarise(n=n()) #527
north_america_rabies_fasta %>% filter(grepl("Mexico", short_Geo_Location)) %>% summarise(n=n()) #20
north_america_rabies_fasta %>% filter(grepl("Greenland", short_Geo_Location)) %>% summarise(n=n()) #14

#New sequence header
north_america_rabies_fasta$seq.name <- paste(north_america_rabies_fasta$seq.name, north_america_rabies_fasta$Host, 
                                             north_america_rabies_fasta$short_Geo_Location, 
                                             north_america_rabies_fasta$Collection_Date, sep = "|") 
north_america_rabies_fasta$seq.name <- gsub("NA", "", north_america_rabies_fasta$seq.name)
north_america_rabies_fasta$seq.name <- gsub(" ", "", north_america_rabies_fasta$seq.name)
north_america_rabies_fasta$city <- sub(".*: ([^,]+),.*", "\\1", north_america_rabies_fasta$Geo_Location)
north_america_rabies_fasta$city[north_america_rabies_fasta$city == "Lewistion"] <- "Lewiston"
north_america_rabies_fasta$state <- sub(".*(.{2})$", "\\1", north_america_rabies_fasta$Geo_Location)

check <- north_america_rabies_fasta %>% filter(!is.na(COUNTY)) 
table(check$address)

north_america_rabies_fasta$city <- ifelse(north_america_rabies_fasta$city %in% c("Canada", "USA", "USA: New York"), north_america_rabies_fasta$LOCATION, north_america_rabies_fasta$city)
north_america_rabies_fasta$state <- ifelse(north_america_rabies_fasta$state %in% c("SA", "rk"), north_america_rabies_fasta$`STATE/ PROVINCE`, north_america_rabies_fasta$state)

north_america_rabies_fasta$address <- paste(north_america_rabies_fasta$city, north_america_rabies_fasta$state, sep = ", ")
north_america_rabies_fasta %>% filter(state %in% c("ME", "NY")) %>% group_by(address) %>% summarize(n=n()) %>% print(n=100)

north_america_rabies_fasta <- north_america_rabies_fasta %>%
  select(c("seq.name", "seq.text", "Isolate", "Geo_Location", "short_Geo_Location", "Host", "Collection_Date", 
        "city", "state", "address", "COUNTY", "XCOORD", "YCOORD"))
north_america_rabies_fasta <- north_america_rabies_fasta %>%
  mutate(Accession = sub("\\|.*$", "", seq.name))
table(north_america_rabies_fasta$short_Geo_Location)
#Rest of US <= 18 sequences (relevant: MA 2, NJ 4, RI 1, PA 4, NH 2)

#####################################################
#DATASET: USA_CANADA_RABIES_SEQUENCES.FASTA (N = )
#####################################################

usa_canada_rabies_sequences <- north_america_rabies_fasta %>%
  filter(grepl("OR227628|OR227629", seq.name) |  grepl("NB|QC|NL|ON|NY|CT|VT|ME", short_Geo_Location))

table(usa_canada_rabies_sequences$short_Geo_Location)
#NB: 32, ON: 213, QC: 93, NL: 41, CT: 71 (inc. MA: 1, RI: 1), ME: 39, NY: 213, VT: 66

dat2fasta(usa_canada_rabies_sequences, "C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/usa_canada_rabies_sequences_new.fasta")

#Check dates
check <- usa_canada_rabies_sequences %>% filter(grepl("^\\d{4}$", usa_canada_rabies_sequences$Collection_Date))
table(check$short_Geo_Location)
usa_canada_rabies_sequences_completedate <- usa_canada_rabies_sequences %>% 
  filter(!grepl("^\\d{4}$", usa_canada_rabies_sequences$Collection_Date))

canada_details <- read_excel("C:/Users/gev25289/Downloads/pntd.0008113.s003.xlsx")
usa_canada_rabies_sequences <- usa_canada_rabies_sequences %>%
  left_join(canada_details, by = c("Accession" = "NCBI ACCESSION NUMBER"))

#Tree to subset USA Canada dataset
#tree <- read.tree("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/beforeRemoveDivergent/usa_canada_rabies_sequences_aligned_trim.fasta.treefile")
tree <- read.tree("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/afterSubsettingNY/usa_canada_rabies_sequences_remove_subset_new.aligned.trim.fasta.treefile")
tree_before <- ggtree(tree) + geom_treescale()
tree <- midpoint.root(tree) #splits tree to look at most divergent sequences
ggtree(tree) + geom_text(aes(label=node), hjust=-.3)
#tree_after <- ggtree(tree) + geom_hilight(node=1312, fill="pink") + geom_treescale() #check
tree_after <- ggtree(tree) + geom_hilight(node=1373, fill="pink") + geom_treescale() #check
p <- plot_grid(tree_before, tree_after, ncol=2, nrow=1, labels = c("A", "B"))
p
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/treeprep.pdf", width=10, height=7)
#tipclade <- tips(tree, 1312) %>% as.data.frame()
tipclade <- tips(tree, 1373) %>% as.data.frame()
tipclade$short_Geo_Location <-  sub("^[^|]*\\|[^|]*\\|([^|]*)\\|.*$", "\\1", tipclade$.)
tipclade$. <- gsub("_", ":", tipclade$.)
table(tipclade$short_Geo_Location) 
usa_canada_rabies_sequences_remove <- usa_canada_rabies_sequences %>%
  filter(seq.name %notin% tipclade$.) #590
table(usa_canada_rabies_sequences_remove$short_Geo_Location)

#Subset
usa_canada_rabies_sequences_remove <- usa_canada_rabies_sequences_remove %>% filter(short_Geo_Location != "USA:NJ")
usa_canada_rabies_sequences_remove_NY <- usa_canada_rabies_sequences_remove %>% filter(short_Geo_Location == "USA:NY")

set.seed(123456)
usa_canada_rabies_sequences_remove_NY <- usa_canada_rabies_sequences_remove_NY %>% slice_sample(n=50)
usa_canada_rabies_sequences_remove_rmNY <- usa_canada_rabies_sequences_remove %>% filter(short_Geo_Location != "USA:NY")
usa_canada_rabies_sequences_remove_new <- rbind(usa_canada_rabies_sequences_remove_NY, usa_canada_rabies_sequences_remove_rmNY) #404

usa_canada_rabies_sequences_remove_new.else <- usa_canada_rabies_sequences_remove_new %>% filter(!grepl("OR227628|OR227629", seq.name))
usa_canada_rabies_sequences_remove_new.else <- usa_canada_rabies_sequences_remove_new.else %>% filter(!grepl("CT", seq.name)) #333

set.seed(123456)
usa_canada_rabies_sequences_remove_new.ct <- usa_canada_rabies_sequences_remove_new %>% filter(grepl("OR227628|OR227629", seq.name) |  grepl("CT", short_Geo_Location))
usa_canada_rabies_sequences_remove_new.else <- usa_canada_rabies_sequences_remove_new.else %>% slice_sample(prop = 0.7)
usa_canada_rabies_sequences_remove_new <- rbind(usa_canada_rabies_sequences_remove_new.else, usa_canada_rabies_sequences_remove_new.ct) #304
table(usa_canada_rabies_sequences_remove_new$short_Geo_Location)
usa_canada_rabies_sequences_remove_new <- usa_canada_rabies_sequences_remove_new %>% mutate(year = year(as.Date(Collection_Date)))
usa_canada_rabies_sequences_remove_new$year <- ifelse(is.na(usa_canada_rabies_sequences_remove_new$year), usa_canada_rabies_sequences_remove_new$Collection_Date, usa_canada_rabies_sequences_remove_new$year)
table(usa_canada_rabies_sequences_remove_new$year)
table(usa_canada_rabies_sequences_remove_new$short_Geo_Location)

dat2fasta(usa_canada_rabies_sequences_remove_new, "C:/Users/gev25289/Desktop/rabies/fasta/northamerica/afterSubsettingNY/prop/usa_canada_rabies_sequences_remove_subset_prop.fasta")

###########################
#GEOCODE LAT/LONG FROM CITY
###########################
CTrivergeocode <- north_america_rabies_fasta %>%
  filter(state %in% c("NH", "VT", "MA", "ME", "NY", "NJ", "RI")) %>%
  geocode(address, method = 'osm', lat = latitude , long = longitude)
CTrivergeocode$latitude[CTrivergeocode$city == "Water"] <- NA
CTrivergeocode$longitude[CTrivergeocode$city == "Water"] <- NA

CTrivergeocode$latitude <- ifelse(!is.na(CTrivergeocode$YCOORD), CTrivergeocode$YCOORD, CTrivergeocode$latitude)
CTrivergeocode$longitude <- ifelse(!is.na(CTrivergeocode$XCOORD), CTrivergeocode$XCOORD, CTrivergeocode$longitude)
table(CTrivergeocode$state)  

CTrivergeocode_ny.me <- CTrivergeocode
CTrivergeocode <- CTrivergeocode %>% filter(state %notin% c("ME", "NY"))

##############################################
#METADATA FROM DR. LEE (LAT, LONG, ZIP, CLADE)
##############################################

ct_rabies_metadata_lee <- read_csv("C:/Users/gev25289/Desktop/rabies/Data/ct_rabies_metadata_lee.csv") %>%
  select(SampleID, CommonName, Zip, State, Latitude, Longitude, Clade) %>%
  mutate(SampleID = gsub("-","_", SampleID)) %>%
  mutate(SampleID = gsub("_0", "_", SampleID))

#join metadata by sampleID
ct_rabies_metadata_join <- ct_rabies_metadata %>% 
  left_join(ct_rabies_metadata_lee, by = c("Isolate" = "SampleID"))

#Dr. Lee uploaded these to NCBI Virus as complete genomes, but they're in "All Data", not in "Tableau" in his excel sheet
ct_rabies_metadata_join$Zip[ct_rabies_metadata_join$Isolate == "17_3553"] <- 06226
ct_rabies_metadata_join$Zip[ct_rabies_metadata_join$Isolate == "18_811"] <- 06355

ct_rabies_metadata_join$CommonName[ct_rabies_metadata_join$Host == "Procyon lotor"] <- "Raccoon"
ct_rabies_metadata_join$CommonName[ct_rabies_metadata_join$Host == "Mephitidae"] <- "Skunk"

ct_rabies_metadata_join <- ct_rabies_metadata_join %>%
  mutate(CommonName = as.factor(CommonName)) %>%
  mutate(Clade = as.factor(Clade))

ct_rabies_fasta <- ct_rabies_metadata_join %>% left_join(north_america_rabies_sequences, by = c("Accession" = "seq.name"))
ct_rabies_fasta$Geo_Location <- gsub(" ", "", ct_rabies_fasta$Geo_Location)
table(ct_rabies_fasta$Geo_Location)

#Get county and lat/long info based on zip code 
#With zip code and lat/long data, confirmed that the geo_location info can be incorrect, use county instead for CT
library(zipcodeR)
zip_code_db$county <- gsub(" County", "", zip_code_db$county)

sum(!is.na(ct_rabies_fasta$Latitude)) #Out of 69 cases, 47 have lat/long data
sum(!is.na(ct_rabies_fasta$Zip)) #Out of 69 cases, 67 have zip code data

ct_rabies_fasta$Zip <- ifelse(!is.na(ct_rabies_fasta$Zip), paste("0", ct_rabies_fasta$Zip, sep = ""), NA)
ct_rabies_fasta <- ct_rabies_fasta %>% left_join(zip_code_db, by = c("Zip" = "zipcode"))

ct_rabies_fasta$Latitude <- ifelse(is.na(ct_rabies_fasta$Latitude), ct_rabies_fasta$lat, ct_rabies_fasta$Latitude)
ct_rabies_fasta$Longitude <- ifelse(is.na(ct_rabies_fasta$Longitude), ct_rabies_fasta$lng, ct_rabies_fasta$Longitude)

#Missing zip code #1 - New London centroid lat/long 
library(housingData)
ct_rabies_fasta$county[ct_rabies_fasta$Accession == "MN418150"] <- "New London"
geoCounty %>% filter(county =="New London County") 
ct_rabies_fasta$Latitude[ct_rabies_fasta$Accession == "MN418150"] <- 41.48928
ct_rabies_fasta$Longitude[ct_rabies_fasta$Accession == "MN418150"] <- -72.09803	

#Missing zip code #2 Windsor Locks lat/long from Open Street Map
ct_rabies_fasta$county[ct_rabies_fasta$Accession == "PP447329"] <- "Hartford"
ct_rabies_fasta$Latitude[ct_rabies_fasta$Accession == "PP447329"] <- 41.9232
ct_rabies_fasta$Longitude[ct_rabies_fasta$Accession == "PP447329"] <- -72.6553	

levels(ct_rabies_fasta$Clade)
levels(ct_rabies_fasta$CommonName)

###########################
#CONNECTICUT RIVER ANALYSIS
###########################
ct_rabies_fasta <- ct_rabies_fasta %>%
  mutate(CT_River = case_when(
    county %in% c("Fairfield", "Litchfield", "New Haven") ~ "West",
    county %in% c("Tolland", "New London", "Windham") ~ "East",
    Isolate == "18_3755" ~ "East", #Massachusetts sample
    county == "Hartford" & Clade %in% c("2A", "2C") ~ "East",
    county == "Middlesex" & CommonName == "Skunk" ~ "East",
    Isolate == "19_4017" ~ "East", #Rhode Island sample
    Geo_Location == "USA:WindsorLocks,CT" ~ "West",
    county == "Hartford" & Clade %in% c("1B", "2B", "2D") ~ "West",
    county == "Middlesex" & CommonName != "Skunk" ~ "West"
  ))

ct_rabies_fasta <- ct_rabies_fasta %>%
  mutate(Counties_CT_River = case_when(
    county %in% c("Fairfield", "Litchfield", "New Haven") ~ ct_rabies_fasta$county,
    county %in% c("Tolland", "New London", "Windham") ~ ct_rabies_fasta$county,
    Isolate == "18_3755" ~ "Tolland", #Massachusetts sample
    Isolate == "19_4017" ~ "Windham", #Rhode Island sample
    county == "Hartford" & Clade %in% c("2A", "2C") ~ "Hartford.East",
    county == "Middlesex" & CommonName == "Skunk" ~ "Middlesex.East",
    county == "Hartford" & Clade %in% c("1B", "2B", "2D") ~ "Hartford.West",
    Geo_Location == "USA:WindsorLocks,CT" ~ "Hartford.West",
    county == "Middlesex" & CommonName != "Skunk" ~ "Middlesex.West"
  ))

#Connecticut and Housatonic River Analysis
ct_rabies_fasta <- ct_rabies_fasta %>%
  mutate(CT_Housa_River = case_when(
    county %in% c("Fairfield") ~ "Section1",
    Clade == "1C" ~ "Section1",
    county == "Litchfield" & Clade != "1C" ~ "Section2",
    county %in% c("New Haven") ~ "Section2",
    county == "Hartford" & Clade %in% c("1B", "2B", "2D") ~ "Section2",
    county == "Middlesex" & CommonName != "Skunk" ~ "Section2",
    Geo_Location == "USA:WindsorLocks,CT" ~ "Section2",
    county %in% c("Tolland", "New London", "Windham") ~ "Section3",
    Isolate == "18_3755" ~ "Section3", #Massachusetts sample
    county == "Hartford" & Clade %in% c("2A", "2C") ~ "Section3",
    county == "Middlesex" & CommonName == "Skunk" ~ "Section3",
    Isolate == "19_4017" ~ "Section3" #Rhode Island sample
  )) 

ct_rabies_fasta <- ct_rabies_fasta %>%
  mutate(Counties_CT_Housa_River = case_when(
    county %in% c("Fairfield") ~ ct_rabies_fasta$county,
    Clade == "1C" ~ "Litchfield.S1",
    county == "Litchfield" & Clade != "1C" ~ "Litchfield.S2",
    county %in% c("New Haven") ~ ct_rabies_fasta$county,
    county == "Hartford" & Clade %in% c("1B", "2B", "2D") ~ "Hartford.S2",
    Geo_Location == "USA:WindsorLocks,CT" ~ "Hartford.S2",
    county == "Middlesex" & CommonName != "Skunk" ~ "Middlesex.S2",
    county %in% c("Tolland", "New London", "Windham") ~ ct_rabies_fasta$county,
    Isolate == "18_3755" ~ "Tolland", #Massachusetts sample
    Isolate == "19_4017" ~ "Windham", #Rhode Island sample
    county == "Hartford" & Clade %in% c("2A", "2C") ~ "Hartford.S3",
    county == "Middlesex" & CommonName == "Skunk" ~ "Middlesex.S3" 
  ))

table(ct_rabies_fasta$CT_River) #2 traits 
table(ct_rabies_fasta$CT_Housa_River) #3 traits
table(ct_rabies_fasta$Counties_CT_River) #From 8 counties -> 10 traits
table(ct_rabies_fasta$Counties_CT_Housa_River) #From 8 counties -> 11 traits

ct_rabies_fasta$seq.name <- paste(ct_rabies_fasta$Accession, 
                                  ct_rabies_fasta$CommonName, 
                                  ct_rabies_fasta$CT_River,
                                  ct_rabies_fasta$Counties_CT_River,
                                  ct_rabies_fasta$CT_Housa_River,
                                  ct_rabies_fasta$Counties_CT_Housa_River, 
                                  ct_rabies_fasta$Clade, 
                                  ct_rabies_fasta$Collection_Date, 
                                  sep = "|") 
ct_rabies_fasta$seq.name <- gsub(" ", "", ct_rabies_fasta$seq.name)
dat2fasta(ct_rabies_fasta, "C:/Users/gev25289/Desktop/rabies/connecticut_rabies_sequences.fasta") #71

rabiesCTRiver_ny.me <- bind_rows(CTrivergeocode_ny.me, ct_rabies_fasta) %>%
  filter(Isolate %notin% c("18-3755", "19-4017")) #duplicates from combining datasets -> 394 sequences
rabiesCTRiver_ny.me$Accession <- ifelse(is.na(rabiesCTRiver_ny.me$Accession), sub("\\|.*", "", rabiesCTRiver_ny.me$seq.name), rabiesCTRiver_ny.me$Accession)
rabiesCTRiver_ny.me$latitude <- ifelse(is.na(rabiesCTRiver_ny.me$latitude), rabiesCTRiver_ny.me$Latitude, rabiesCTRiver_ny.me$latitude)
rabiesCTRiver_ny.me$longitude <- ifelse(is.na(rabiesCTRiver_ny.me$longitude), rabiesCTRiver_ny.me$Longitude, rabiesCTRiver_ny.me$longitude)

check <- rabiesCTRiver_ny.me %>% 
  filter(is.na(latitude))

submitters <- north_america_rabies_metadata %>% select(c("Accession", "Submitters"))
check <- check %>% left_join(submitters, by = "Accession")
table(check$Submitters.y)

rabiesCTRiver_ny.me <- rabiesCTRiver_ny.me %>%
  select(c("seq.name", "seq.text", "Accession", "Isolate", "Geo_Location", "short_Geo_Location", "Host", "CommonName", 
           "Collection_Date", "Zip", "city", "county", "state", "latitude", "longitude", "Clade", "CT_River", "Counties_CT_River"))

rabiesCTRiver_ny.me <- rabiesCTRiver_ny.me %>%
  mutate(CommonName = case_when(
    Host == "Bos taurus" ~ "Cow",
    Host == "Canidae" ~ "Canidae",
    Host == "Canis lupus familiaris" ~ "Dog",
    Host == "Cervidae" ~ "Deer", 
    Host == "Feliformia" ~ "Feliformia",
    Host == "Felis catus" ~ "Cat",
    Host == "Lontra canadensis" ~ "River Otter",
    Host == "Lynx rufus" ~ "Bobcat", 
    Host == "Marmota monax" ~ "Groundhog",
    Host == "Mephitidae" ~ "Skunk",
    Host == "Mephitis mephitis" ~ "Skunk",
    Host == "Procyon lotor" ~ "Raccoon",
    Host == "Urocyon cinereoargenteus" ~ "Gray Fox",
    Host == "Vulpes vulpes" ~ "Red Fox"
  ))
rabiesCTRiver_ny.me$CommonName[rabiesCTRiver_ny.me$seq.name == "OR227628|Raccoon|East|Tolland|Section3|Tolland|2018-08-15|2A"] <- "Raccoon"

rabiesCTRiver_ny.me <- rabiesCTRiver_ny.me %>%
  mutate(CT_River = case_when(
    short_Geo_Location == "USA:ME" ~ "East",
    short_Geo_Location == "USA:NY" ~ "West",
    short_Geo_Location == "USA:VT" ~ "West",
    short_Geo_Location == "USA:NH" ~ "East",
    short_Geo_Location == "USA:RI" ~ "East",
    short_Geo_Location == "USA:MA" ~ "East",
    .default = CT_River
  ))

rabiesCTRiver_ny.me$short_Geo_Location <- ifelse(is.na(rabiesCTRiver_ny.me$short_Geo_Location), 
                                                 paste("USA:CT", rabiesCTRiver_ny.me$Counties_CT_River, sep = "_"), 
                                                 rabiesCTRiver_ny.me$short_Geo_Location)
rabiesCTRiver_ny.me$seq.name <- paste(rabiesCTRiver_ny.me$Accession, 
                                      rabiesCTRiver_ny.me$short_Geo_Location,
                                      rabiesCTRiver_ny.me$CommonName, 
                                      rabiesCTRiver_ny.me$CT_River,
                                      rabiesCTRiver_ny.me$Clade,
                                      rabiesCTRiver_ny.me$Collection_Date, 
                                      sep = "|") 
rabiesCTRiver_ny.me$seq.name <- gsub(" ", "", rabiesCTRiver_ny.me$seq.name)
dat2fasta(rabiesCTRiver_ny.me, "C:/Users/gev25289/Desktop/rabies/Fasta/river/allsequences/mafft/connecticut_river_sequences_all.fasta") #394 -> 432

#Remove too few obs: MA, RI, and NH (except for CT included seqs)
rabiesCTRiver_ny.me %>% filter(state %in% c("MA", "RI", "NH")) %>% select(seq.name, short_Geo_Location, state)
rabiesCTRiver_ny.me.remove <- rabiesCTRiver_ny.me %>% filter(short_Geo_Location %notin% c("USA:MA", "USA:NH", "USA:RI")) 
rabiesCTRiver_ny.me.remove$seq.name <- gsub(" ", "", rabiesCTRiver_ny.me.remove$seq.name)
dat2fasta(rabiesCTRiver_ny.me.remove, "C:/Users/gev25289/Desktop/rabies/Fasta/river/allsequences/removeMA.NH.RI/mafft/#connecticut_river_sequences_all.remove.fasta") #389 -> 427

#Try a dataset that doesn't include maine and new york
rabiesCTRiver <- rabiesCTRiver_ny.me %>% filter(state %notin% c("ME", "NY", "rk")) 
rabiesCTRiver$seq.name <- gsub(" ", "", rabiesCTRiver$seq.name)
dat2fasta(rabiesCTRiver, "C:/Users/gev25289/Desktop/rabies/Fasta/river/removeME.NY/mafft/connecticut_river_sequences.fasta") #142

#Remove too few obs: MA, RI, and NH (except for CT included seqs)
rabiesCTRiver.remove <- rabiesCTRiver %>% filter(short_Geo_Location %notin% c("USA:MA", "USA:NH", "USA:RI"))
rabiesCTRiver.remove$seq.name <- gsub(" ", "", rabiesCTRiver.remove$seq.name)
dat2fasta(rabiesCTRiver.remove, "C:/Users/gev25289/Desktop/rabies/Fasta/river/removeME.NY/removeMA.NH.RI/mafft/connecticut_river_sequences_remove.fasta") #137

saveRDS(ct_rabies_fasta, "C:/Users/gev25289/Desktop/rabies/ct_rabies_fasta.rds")
saveRDS(usa_canada_rabies_sequences, "C:/Users/gev25289/Desktop/rabies/rds/usa_canada_rabies_sequences.rds")
saveRDS(rabiesCTRiver_ny.me.remove, "C:/Users/gev25289/Desktop/rabies/rds/rabiesCTRiver_ny.me.remove.rds")

#check what's not adding up between the NCBI Virus dataset and Dr. Lee's
which(ct_rabies_metadata_lee$SampleID %notin% ct_rabies_fasta$Isolate) #7
  #19_3631 is not publicly available, I'm leaving it out
