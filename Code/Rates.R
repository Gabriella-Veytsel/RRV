library(tidyverse)
library(cowplot)

rates <- read.delim("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/relaxed_discrete_location/combined.rates.log" ,sep = '', header = TRUE)
library(ggplot2)

rates_west <- rates %>% select(
  Location.actualRates6, 
  Location.actualRates7,
  Location.actualRates10,
  Location.actualRates16,
  Location.actualRates17,
  Location.actualRates20,
  )

rates_east <- rates %>% select(
  Location.actualRates2,
  Location.actualRates12
)

rates_east_west <- rates %>% select(
  Location.actualRates1,
  Location.actualRates3,
  Location.actualRates4,
  Location.actualRates8,
  Location.actualRates9,
  Location.actualRates15
)

rates_west_east <- rates %>% select(
  Location.actualRates5,
  Location.actualRates11,
  Location.actualRates13,
  Location.actualRates14,
  Location.actualRates18,
  Location.actualRates19
)

#############################################################################
#PP >50%#####################################################################
#############################################################################
rates_west <- rates %>% select(
  Location.actualRates10,
  Location.actualRates16,
  Location.actualRates20,
)

rates_east <- rates %>% select(
  Location.actualRates2
)

rates_east_west <- rates %>% select(
  Location.actualRates1,
  Location.actualRates4,
  Location.actualRates9
)

rates_west_east <- rates %>% select(
  Location.actualRates13
)

#rates_east$WithinEast <- rowMeans(rates_east[,])
rates_west$WithinWest <- rowMeans(rates_west[,])
rates_east_west$East_to_West <- rowMeans(rates_east_west[,])
#rates_west_east$West_to_East <- rowMeans(rates_west_east[,])

density <- ggplot() + xlab("Mean transition rate/location") + #ggtitle("Density distribution of statistically supported mean transition rates between locations") +
  geom_density(aes(Location.actualRates2, fill = 'Within East'), alpha = .8, data = rates_east) +
  geom_density(aes(WithinWest, fill = 'Within West'), alpha = .8, data = rates_west) +
  geom_density(aes(East_to_West, fill = 'East to West'), alpha = .8, data = rates_east_west) +
  geom_density(aes(Location.actualRates13, fill = 'West to East'), alpha = .8, data = rates_west_east) +
  scale_fill_manual(name = "", values = c('royalblue',
                                          'cyan',
                                          'red',
                                          'yellow')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(text=element_text(size=12))
ggsave("density distribution.pdf", width=6, height=6)

mean(rates_east$Location.actualRates2) #0.73
mean(rates_west$WithinWest) #1.23
mean(rates_east_west$East_to_West) #1.24
mean(rates_west_east$Location.actualRates13) #0.38

within <- cbind(rates_east$Location.actualRates2, rates_west$WithinWest)
within <- as.data.frame(within)
scatter <- ggplot(within, aes(x=V1, y=V2)) + geom_point() + geom_abline(color="red") + #ggtitle("Statistically supported mean migration rates per MCMC step") +
  annotate("text", x = 3.8, y = 3.7, label = "x=y",  color="red", angle=45) +
  xlab("Within East") + ylab("Within West") + 
  theme(text=element_text(size=12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + theme(legend.position = "none")
ggsave("scatterplot.pdf", width=6, height=6)

p <- plot_grid(density, scatter, labels = c('B', 'C'),  ncol = 1)
plot_grid(ratesmap, p, labels = c('A'), ncol = 2, rel_widths = c(1.75, 1))
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/addresssamplebias/figures/density_and_scatterplot.pdf", width=20, height=10)
