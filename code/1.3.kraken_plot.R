# Install and load necessary packages
install.packages("readxl")
install.packages("ggplot2")
install.packages("tidyverse")

library(readxl)
library(ggplot2)
library(tidyverse)

# Read the Excel file
# Assume the file is named "KRAKEN.xlsx" and the data is in "Sheet1"
data <- read_xlsx("KRAKEN.xlsx", sheet = "Sheet1")

# View the data
print(data)

# Convert data from wide format to long format
# Note: The first column is named `...1`, not `Sample`
data_long <- data %>%
  pivot_longer(cols = -`...1`, names_to = "Sample", values_to = "Value") %>%
  rename(Category = `...1`)  # Rename `...1` column to `Category`

# Set factor levels for Category
data_long$Category <- factor(data_long$Category, 
                             levels = c("unclassified", "others", "Bacteria", "Viruses", "Homo"))

# Plot and save as PDF
pdf("contamination_status.pdf", width = 10, height = 6)
ggplot(data_long, aes(x = Sample, y = Value, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "The contamination status of each sample",
       x = "Sample",
       y = "Percentage (%)",
       fill = "Category") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "unclassified" = "#9467BD",  # purple
    "Homo" = "#1f77b4",         # blue
    "Bacteria" = "#ff7f0e",     # orange
    "Viruses" = "#2ca02c",      # green
    "others" = "#d62728"        # red
  ))

# Close the PDF device
dev.off()
