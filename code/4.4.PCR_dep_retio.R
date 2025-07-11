data <- read.csv("reads_count.csv", header = TRUE, sep = ",")

colnames(data) <- c("sample", "total_reads", "reads_part", "percentage")

data$percentage_label <- paste0(round(data$percentage * 100, 2), "%")


library(ggplot2)

pdf("Simple_Barplot_Reads.pdf", width = 12, height = 8)  

ggplot(data, aes(x = sample, y = reads_part)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = percentage_label), 
            vjust = -0.5, size = 3, color = "black") +
  labs(title = "Bar Plot of Duplicates Reads with Percentage",
       x = "Sample",
       y = "Number of Reads") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

dev.off()
