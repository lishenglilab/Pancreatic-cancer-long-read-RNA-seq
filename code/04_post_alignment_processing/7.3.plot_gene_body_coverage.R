# 1. Load necessary packages
library(data.table)
library(ggplot2)

# 2. Define the color scheme for samples
color_scheme <- c(
  "PANC0203-1" = "#1f77b4", "PANC1005-1" = "#ff7f0e",
  "Panc0327-1" = "#2ca02c", "APSC-1-1" = "#d62728",
  "BXPC3-1" = "#9467bd", "Capan-1-1" = "#8c564b",
  "Capan-2-1" = "#e377c2", "HuP-T4-1" = "#7f7f7f",
  "Mia-Paca-2-1" = "#bcbd22", "SW1990-1" = "#17becf",
  "PANC0203-2" = "#aec7e8", "PANC1005-2" = "#ffbb78",
  "Panc0327-2" = "#98df8a", "APSC-1-2" = "#ff9896",
  "BXPC3-2" = "#c5b0d5", "Capan-1-2" = "#c49c94",
  "Capan-2-2" = "#f7b6d2", "HuP-T4-2" = "#c7c7c7",
  "Mia-Paca-2-2" = "#dbdb8d", "SW1990-2" = "#9edae5"
)

# 3. Generate a list of filenames for all samples
sample_names <- names(color_scheme)
file_names <- paste0(sample_names, "_metrics.txt")

# 4. Loop through, read, and process all sample files
# Use lapply to iterate over all filenames and store the results in a list
all_samples_list <- lapply(file_names, function(file) {
  # Check if the file exists; if not, skip it
  if (!file.exists(file)) {
    warning(paste("File not found, skipping:", file))
    return(NULL)
  }
  
  # Extract the sample name from the filename
  sample_name <- gsub("_metrics.txt$", "", file)
  
  # Locate the starting line of the data within the file
  lines <- readLines(file)
  hist_line <- grep("^## HISTOGRAM", lines)
  
  # Use fread to efficiently read the data
  hist_data <- fread(file, skip = hist_line, sep = "\t")
  
  # Filter and rename data, and add a new column to identify the sample
  processed_data <- hist_data[, .(
    sample = sample_name,  # Add a column for the sample name
    pos = as.numeric(normalized_position),
    cov = as.numeric(All_Reads.normalized_coverage)
  )]
  
  return(processed_data)
})

# 5. Combine all data.tables from the list into a single one
# rbindlist is a fast merging function provided by data.table
all_data <- rbindlist(all_samples_list)

# 6. Plot the data using ggplot2
# Added group and color mappings in aes()
p <- ggplot(all_data, aes(x = pos, y = cov, group = sample, color = sample)) +
  geom_line(linewidth = 0.8, alpha = 0.8) + # Use alpha for transparency to prevent lines from being too cluttered
  
  # Apply the custom color scheme using scale_color_manual
  scale_color_manual(values = color_scheme) +
  
  labs(title = "",
       x = "Normalized position",
       y = "Normalized coverage",
       color = "Sample") + # Change the legend title
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right" # The legend position can be adjusted as needed
  )

# 7. Save the plot
# Increase the width slightly to accommodate the legend
ggsave("All_Samples_Coverage.pdf", p, width = 8, height = 5, device = "pdf")

print("Plot generated: All_Samples_Coverage.pdf")