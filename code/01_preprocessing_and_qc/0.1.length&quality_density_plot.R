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
    sample_median = quantile(read_length, 0.5)
  ) %>%
  ungroup()

# Calculate the mean of all sample medians
mean_of_medians <- mean(median_values$sample_median)


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



# Calculate the mean of the median quality scores across all samples
mean_of_medians <- data %>%
  group_by(sample) %>%
  summarise(median_quality = median(avg_quality_correct), .groups = 'drop') %>%
  summarise(mean_val = mean(median_quality)) %>%
  pull(mean_val) # Extract the final calculated value from the data frame

# Print the result for verification
print(paste("The mean of medians for all samples is:", round(mean_of_medians, 2)))

# Create the plot and save it as a PDF
pdf("density_plot_Q_with_mean_median.pdf", width = 10, height = 6)

ggplot(data, aes(x = avg_quality_correct, color = sample)) +
  geom_density(aes(y = after_stat(count)), size = 0.4) +
  
  # --- New layer: Add a vertical line representing the mean of medians ---
  geom_vline(
    aes(xintercept = mean_of_medians, linetype = "Mean of Medians"), # Define in aes to show in the legend
    color = "firebrick", # Use a prominent color
    size = 0.8 # Make the line thicker
  ) +
  
  # Manually set sample colors (assuming color_palette is defined elsewhere)
  scale_color_manual(values = color_palette) +
  
  # --- New layer: Define the linetype for the vertical line and add it to the legend ---
  scale_linetype_manual(
    name = "Overall Metric", # Title in the legend
    values = c("Mean of Medians" = "dashed") # Set the linetype to dashed
  ) +
  
  # Set the coordinate axis limits
  coord_cartesian(
    xlim = c(0, 40),
    ylim = c(0, 1500000)
  ) +
  
  # Format the Y-axis
  scale_y_continuous(
    labels = label_comma(),
    breaks = seq(0, 1500000, by = 500000)
  ) +
  
  # Set titles and labels
  labs(
    title = "Quality Score Distribution by Sample",
    subtitle = paste("Overall Mean of Sample Medians is", round(mean_of_medians, 2)), # Add the mean to the subtitle
    x = "Average Quality (Q)",
    y = "Count",
    color = "Sample" # Rename the color legend title
  ) +
  
  # Apply and customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "firebrick", size = 12, hjust = 0.5), # Center and highlight the subtitle
    legend.position = "right",
    legend.box = "vertical", # Arrange legends vertically
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

# Close the PDF device
dev.off()

