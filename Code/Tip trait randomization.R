library(ggplot2)
library(ggpubr)  # For stat_cor
library(patchwork)  # For combining plots

#Connecticut River
###################################################
# Data for ctriverA
state_table_A <- c(31, 40, 28, 127, 72)
root_table_A <- c(0.0053, 0.0134, 0.0003, 0.9796, 0.0014)
sample_bias_A <- data.frame(state_table_A, root_table_A)

# Data for ctriverB
state_table_B <- c(31, 82, 50, 40, 72)
root_table_B <- c(0.0071, 0.0463, 0.4208, 0.5255, 0.0003)
sample_bias_B <- data.frame(state_table_B, root_table_B)

# Combine the x-values and y-values from both plots to determine common limits
x_limits <- range(c(state_table_A, state_table_B), na.rm = TRUE)
y_limits <- range(c(root_table_A, root_table_B), na.rm = TRUE)

# Create the plots
ctriverA <- ggplot(sample_bias_A, aes(x=state_table_A, y=root_table_A)) +
  stat_cor(method = "pearson", label.x = max(x_limits) * 0.7, label.y = max(y_limits) * 0.9) +  # Adjusted position
  geom_point() +
  xlab("number of taxa") +
  ylab("root state probability") +
  xlim(x_limits) +  # Set x-axis limits
  ylim(y_limits) +  # Set y-axis limits
  annotate("text", x = min(x_limits) + 0.1 * (max(x_limits) - min(x_limits)), 
         y = max(y_limits) * 0.95, label = "A", size = 6, fontface = "bold")  # Label A

ctriverB <- ggplot(sample_bias_B, aes(x=state_table_B, y=root_table_B)) +
  stat_cor(method = "pearson", label.x = max(x_limits) * 0.7, label.y = max(y_limits) * 0.9) +  # Adjusted position
  geom_point() +
  xlab("number of taxa") +
  ylab("root state probability") +
  xlim(x_limits) +  # Set x-axis limits
  ylim(y_limits) + # Set y-axis limits
  annotate("text", x = min(x_limits) + 0.1 * (max(x_limits) - min(x_limits)), 
           y = max(y_limits) * 0.95, label = "B", size = 6, fontface = "bold")  # Label B

# Combine the plots into one
combined_plot <- ctriverA + ctriverB + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/river/buffer/tiptraitrandomizationcombined.pdf", width=15, height=5)

library(ggplot2)
library(ggpubr)  # For stat_cor
library(patchwork)  # For combining plots

# Data for northamericaA
state_table_A <- c(52, 69, 39, 212, 32, 66, 84)
root_table_A <- c(0.0003, 0.0005, 0.0005, 0.99985, 0.0001, 0.0001, 0)
sample_bias_A <- data.frame(state_table_A, root_table_A)

# Data for northamericaB
state_table_B <- c(38, 71, 26, 62, 32, 24, 51)
root_table_B <- c(0.1216, 0.4633, 0.0004, 0.0002, 0.4121, 0.0020, 0.0005)
sample_bias_B <- data.frame(state_table_B, root_table_B)

# Combine the x-values and y-values from both plots to determine common limits
x_limits <- range(c(state_table_A, state_table_B), na.rm = TRUE)
y_limits <- range(c(root_table_A, root_table_B), na.rm = TRUE)

# Create the plots
northamericaA <- ggplot(sample_bias_A, aes(x=state_table_A, y=root_table_A)) +
  stat_cor(method = "pearson", label.x = max(x_limits) * 0.75, label.y = max(y_limits) * 0.9) +  # Adjusted position
  geom_point() +
  xlab("number of taxa") +
  ylab("root state probability") +
  xlim(x_limits) +  # Set x-axis limits
  ylim(y_limits) +  # Set y-axis limits
  annotate("text", x = min(x_limits) + 0.1 * (max(x_limits) - min(x_limits)), 
           y = max(y_limits) * 0.95, label = "A", size = 6, fontface = "bold")  # Label A

northamericaB <- ggplot(sample_bias_B, aes(x=state_table_B, y=root_table_B)) +
  stat_cor(method = "pearson", label.x = max(x_limits) * 0.75, label.y = max(y_limits) * 0.9) +  # Adjusted position
  geom_point() +
  xlab("number of taxa") +
  ylab("root state probability") +
  xlim(x_limits) +  # Set x-axis limits
  ylim(y_limits) +  # Set y-axis limits
  annotate("text", x = min(x_limits) + 0.1 * (max(x_limits) - min(x_limits)), 
           y = max(y_limits) * 0.95, label = "B", size = 6, fontface = "bold")  # Label B

# Combine the plots into one
combined_plot <- northamericaA + northamericaB + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)
ggsave("C:/Users/gev25289/Desktop/rabies/Fasta/northamerica/tiptraitrandomizationcombined.pdf", width=15, height=5)
