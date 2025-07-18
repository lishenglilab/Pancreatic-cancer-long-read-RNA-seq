# This R script generates density plots for read length and average quality scores
# from nanopore sequencing data. It uses ggplot2 for plotting.

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(scales)

data <- read.csv("nanopore_read_stats.csv")

# Filter the data to remove outliers in read length.
# This helps in creating a more focused and readable plot for read length distribution.
data_filtered <- data %>%
  group_by(sample) %>%
  mutate(
    lower = quantile(read_length, 0.005),  
    upper = quantile(read_length, 0.995)  
  ) %>%
  filter(read_length >= lower & read_length <= upper) %>%
  ungroup()


color_palette <- c(
  "PANC0203-1" = "#1f77b4", "PANC1005-1" = "#ff7f0e", "Panc0327-1" = "#2ca02c",
  "APSC-1-1" = "#d62728", "BXPC3-1" = "#9467bd", "Capan-1-1" = "#8c564b",
  "Capan-2-1" = "#e377c2", "Hup-T4-1" = "#7f7f7f", "Mia-Paca-2-1" = "#bcbd22",
  "SW1990-1" = "#17becf", "PANC0203-2" = "#aec7e8", "PANC1005-2" = "#ffbb78",
  "Panc0327-2" = "#98df8a", "APSC-1-2" = "#ff9896", "BXPC3-2" = "#c5b0d5",
  "Capan-1-2" = "#c49c94", "Capan-2-2" = "#f7b6d2", "Hup-T4-2" = "#c7c7c7",
  "Mia-Paca-2-2" = "#dbdb8d", "SW1990-2" = "#9edae5"
)





# Before filtering (using raw 'data'), calculate the 50th percentile (median) for each sample
median_values <- data %>%
  group_by(sample) %>%
  summarise(
    sample_median = quantile(read_length, 0.5, na.rm = TRUE)
  ) %>%
  ungroup()

# Calculate the mean of all sample medians
mean_of_medians <- mean(median_values$sample_median, na.rm = TRUE)


# Generate the plot
pdf("density_plot_length_mean_of_medians.pdf", width = 10, height = 6)
ggplot(data_filtered, aes(x = read_length, color = sample)) +
  geom_density(aes(y = after_stat(count)), size = 0.4) +
  # Add vertical line using the calculated 'mean_of_medians'
  geom_vline(
    xintercept = mean_of_medians,
    color = "black",
    linetype = "dashed",
    size = 0.6
  ) +
  annotate(
    "text",
    x = mean_of_medians,
    y = Inf,
    label = sprintf("Mean of Medians: %.1f", mean_of_medians),
    vjust = 2,
    color = "black"
  ) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(labels = label_comma()) +
  labs(
    title = "Read Length Distribution by Sample",
    x = "Read Length",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8)
  )
dev.off()



# Generate and save the density plot for average quality scores
pdf("density_plot_Q.pdf", width = 10, height = 6)
ggplot(data, aes(x = avg_quality, color = sample)) +
  geom_density(aes(y = after_stat(count)), size = 0.4) + 
  scale_color_manual(values = color_palette) +
  coord_cartesian(
    xlim = c(0, 40),         
    ylim = c(0, 1500000)      
  ) +
  scale_y_continuous(
    labels = label_comma(),  
    breaks = seq(0, 1500000, by = 500000)  
  ) +
  labs(
    title = "Quality Score Distribution by Sample",
    x = "Average Quality (Q)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()  
  )
dev.off()


