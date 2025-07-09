library(readxl)
library(ggplot2)

data <- read_excel("error_rate.xlsx")

data$error_rate <- data$error_rate * 100

pdf("error_rate.pdf", width = 12, height = 8)
ggplot(data, aes(x = sample, y = error_rate)) +
  geom_col(fill = "steelblue") +
  labs(title = "Error Rate by Sample", 
       x = "Sample", 
       y = "Error Rate (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))
dev.off()

