# Load necessary libraries
library(dplyr)
library(tidyr)



# Step 2: Data Transformation
# Convert the data to a long format where each methylation pattern is a separate row
data_long <- data %>%
  gather(key = "Methylation_Status", value = "Frequency", `X.000`:`X.111`) %>%
  mutate(Methylation_Status = as.factor(Methylation_Status))  # Convert Methylation_Status to factor

# Step 3: Data Preprocessing
# Ensure both tissue types have data for each CpG_Coordinates and Methylation_Status combination
data_long_filtered <- data_long %>%
  group_by(CpG_Coordinates, Methylation_Status) %>%
  filter(n_distinct(Tissue) == 2 & all(Frequency > 0))  # Ensure both tissues have data and no zero frequencies

# Check group sizes before performing the t-test
group_sizes <- data_long %>%
  group_by(CpG_Coordinates, Methylation_Status, Tissue) %>%
  summarise(count = n(), .groups = 'drop')


# Step 4: Perform T-test only on valid groups (those with both Tissue #1 and Tissue #2)
# First, filter groups where both tissues are present

valid_groups <- data_long_filtered %>%
  group_by(CpG_Coordinates, Methylation_Status) %>%
  filter(n_distinct(Tissue) == 2) %>%
  summarise(count_islet = sum(Tissue == "Islet"), count_cfdna = sum(Tissue == "cfDNA"), .groups = 'drop')

# View valid groups (optional)
print(valid_groups)

# Now, run the t-test only for valid groups
p_values <- data_long_filtered %>%
  group_by(CpG_Coordinates, Methylation_Status) %>%
  filter(CpG_Coordinates %in% valid_groups$CpG_Coordinates & Methylation_Status %in% valid_groups$Methylation_Status) %>%
  summarise(p_value = tryCatch(t.test(Frequency ~ Tissue)$p.value, error = function(e) NA), .groups = 'drop')

# View the result
print(p_values)


# Step 1: Define significance threshold for p-value (e.g., p < 0.05)
threshold_p_value <- 0.05

# Step 2: Filter PMPs with significant p-values (high specificity)
# Select only those PMPs where p-value < threshold
significant_pmp <- p_values %>%
  filter(p_value < threshold_p_value)

# View significant PMPs
print(significant_pmp)

# Step 3: Assign confidence score to each PMP
# We can assign confidence score as inverse log of p-value (or other method)
significant_pmp <- significant_pmp %>%
  mutate(confidence_score = -log10(p_value))

# View significant PMPs with confidence scores
print(significant_pmp)

# Step 4: Summarize the results
# For example, you can get the top 10 most confident PMPs
top_pmp <- significant_pmp %>%
  arrange(desc(confidence_score)) %>%
  head(10)

# View the top 10 PMPs
print(top_pmp)

write.csv(top_pmp,"top_pmp.csv")


# Step 1: Calculate the mean VRF for each CpG_Coordinates and Methylation_Status in both tissues
mean_vrf <- data_long_with_frequency %>%
  group_by(CpG_Coordinates, Methylation_Status, Tissue) %>%
  summarise(mean_vrf = mean(Frequency), .groups = 'drop')

# View the mean VRF for each PMP
print(mean_vrf)

# Step 2: Reshape the data for easier interpretation
mean_vrf_wide <- mean_vrf %>%
  pivot_wider(names_from = Tissue, values_from = mean_vrf, names_prefix = "mean_vrf_")

# View the wide format data with mean VRF for each tissue
print(mean_vrf_wide)





# Step 2: Filter the data for only the top 10 PMPs
top_pmp <- data_long_with_frequency %>%
  filter(CpG_Coordinates %in% top_pmp$CpG_Coordinates &
           Methylation_Status %in% top_pmp$Methylation_Status)

# Step 3: Calculate mean VRF for each CpG_Coordinates and Methylation_Status in both tissues
mean_vrf_top_10 <- top_pmp %>%
  group_by(CpG_Coordinates, Methylation_Status, Tissue) %>%
  summarise(mean_vrf = mean(Frequency), .groups = 'drop')

# Step 4: Reshape the results for easier interpretation
mean_vrf_top_10_wide <- mean_vrf_top_10 %>%
  pivot_wider(names_from = Tissue, values_from = mean_vrf, names_prefix = "mean_vrf_")

# View the result
print(mean_vrf_top_10_wide)

write.csv(mean_vrf_top_10_wide, "mean_vrf_top10.csv")

