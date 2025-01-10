# Load necessary libraries
library(dplyr)
library(ggplot2)


data <- read.csv("PupilBioTest_PMP_revA.csv")


# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Load the dataset
data <- read.csv("PupilBioTest_PMP_revA.csv")

# Step 2: Prepare the data
# Summing across all methylation patterns to get total coverage per CpG
data <- data %>%
  mutate(Total_Coverage = `X.000` + `X.001` + `X.010` + `X.011` + `X.100` + `X.101` + `X.110` + `X.111`)

# Step 3: Calculate Median and CV for each tissue
coverage_stats <- data %>%
  group_by(Tissue) %>%
  summarise(
    Median_Coverage = median(Total_Coverage),
    Mean_Coverage = mean(Total_Coverage),
    SD_Coverage = sd(Total_Coverage),
    CV_Coverage = (SD_Coverage / Mean_Coverage) * 100
  )

# Print the calculated statistics
print(coverage_stats)

# Step 4: Visualize Coverage Statistics
# Boxplot for Total Coverage by Tissue
ggplot(data, aes(x = Tissue, y = Total_Coverage, fill = Tissue)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Coverage Distribution by Tissue",
       x = "Tissue",
       y = "Total Coverage")

# Optional: Density plot
ggplot(data, aes(x = Total_Coverage, color = Tissue, fill = Tissue)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Density Plot of Coverage by Tissue",
       x = "Total Coverage",
       y = "Density")
