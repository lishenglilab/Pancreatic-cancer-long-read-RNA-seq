# This script reads Kraken analysis results from an Excel file,
# processes the data, and generates a stacked bar plot to visualize
# the taxonomic composition (contamination status) of each sample.

library(readxl)
library(ggplot2)
library(tidyverse)

# Read the Excel file
# Assume the file is named "KRAKEN.xlsx" and the data is in "Sheet1"
data <- read_xlsx("KRAKEN.xlsx", sheet = "Sheet1")

# View the data
print(data)

# Clean and transform the data
data_clean <- data %>%
  rename(category = ...1) %>%
  mutate(across(-category, as.character))

# Convert to long format
data_long <- data_clean %>%
  pivot_longer(cols = -category, 
               names_to = "sample", 
               values_to = "value") %>%
  # Ensure value column is numeric
  mutate(value = as.numeric(value))


# Set factor levels to place Homo Sapien at the bottom
data_long$category <- factor(data_long$category, 
                             levels = c( "Mycoplasmatota", "Viruses", "others", "unclassified", "Homo Sapien"))

# Define color palette
colors <- c(
  "unclassified" = "#9467BD",
  "Homo Sapien" = "#1f77b4",
  "Mycoplasmatota" = "#ff7f0e",
  "Viruses" = "#2ca02c",
  "others" = "#d62728"
)

# Create PDF output
pdf("contamination_status_v2.pdf", width = 10, height = 6)
ggplot(data_long, aes(x = sample, y = value, fill = category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  labs(
    title = "Composition of Different Categories in Each Sample",
    x = "Sample",
    y = "Percentage (%)",
    fill = "Category"
  )

dev.off()
