library(readxl)
library(grafify)

colors <- c("darkorange", "cadetblue1", "cornflowerblue", "lightcoral", "darkmagenta",
            "cadetblue", "plum", "darkgray", "firebrick", "darkslateblue",
            "goldenrod1", "deeppink", "dodgerblue4", "darkturquoise", "darkcyan", "royalblue",
            "darkorchid1", "violet", "slateblue", "deeppink4", "thistle", "chocolate3")

rabiescases <- read_excel("D:/rabies/Cases by Year.xlsx")
names <- colnames(rabiescases)[c(2:20)] %>% as.vector() 
rabiescasestidy <- rabiescases %>%
  pivot_longer(names, names_to = "Host", values_to = "cases")
ggplot(data = rabiescasestidy, aes(x=Year, y=cases, fill=Host)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = colors) +
  ylab("No. of Reported Rabies Cases")

ggsave("D:/rabies/cases.pdf", width=10, height=7)