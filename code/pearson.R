# Load necessary package
library(pheatmap)

# Read the TSV file
tmptsv <- read_tsv("flair_quantify.tpm.tsv")


# Extract numeric data (exclude the first column which contains IDs)
numeric_data <- tmptsv[, -1]

# Calculate correlation matrices
pearson_cor <- cor(numeric_data, method = "pearson", use = "complete.obs")

# Set up color palette for heatmap (blue-white-red gradient)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Set PDF size and proportions, increasing height to ensure full visibility
pdf("pearson_correlation_heatmap.pdf", width = 9, height = 7)  # Adjusted height from 8 to 12

# Plot Pearson correlation heatmap
pheatmap(pearson_cor,
         color = color_palette,  # Assuming color palette is defined
         main = "Pearson Correlation",
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         margins = c(1, 2),        # Bottom and right margins
         fontsize_row = 7,         # Increase row label font size
         fontsize_col = 7,         # Increase column label font size
         angle_col = 45,           # Rotate column labels
         display_numbers = TRUE,   # Display numbers in cells
         number_format = "%.2f",   # Number formatting
         number_color = "black",   # Color of the displayed numbers
         fontsize_number = 6,      # Increase font size of numbers inside cells
         cellwidth = 20,           # Set cell width
         cellheight = 15)          # Set cell height

# Close the PDF device
dev.off()