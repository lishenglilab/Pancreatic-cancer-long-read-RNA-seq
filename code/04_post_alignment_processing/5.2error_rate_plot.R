# This script visualizes alignment error rates from an Excel file.
# It creates a bar plot showing the error rate for each sample.

# --- 1. Load Libraries ---
library(readxl)
library(ggplot2)

# --- 2. Read and Prepare Data ---
data <- read_excel("error_rate.xlsx")
data$error_rate <- data$error_rate * 100

# --- 3. Generate and Save Plot ---
pdf("error_rate.pdf", width = 12, height = 8)
ggplot(data, aes(x = sample, y = error_rate)) +
  geom_col(fill = "steelblue") +
  labs(title = "Error Rate by Sample", 
       x = "Sample", 
       y = "Error Rate (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))
dev.off()

