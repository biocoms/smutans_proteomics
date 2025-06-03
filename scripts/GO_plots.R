# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# **Load data**
supernatant <- read.csv("data/supernatant_data.csv", stringsAsFactors = FALSE, na.strings = c(NA, "") )
intracellular <- read.csv("data/intracellular_data.csv", stringsAsFactors = FALSE, na.strings = c(NA, ""))
supernatant$biological_process[is.na(supernatant$biological_process)] <- "Unannotated"
supernatant$molecular_function[is.na(supernatant$molecular_function)] <- "Unannotated"
supernatant$cellular_component[is.na(supernatant$cellular_component)] <- "Unannotated"
intracellular$biological_process[is.na(intracellular$biological_process)] <- "Unannotated"
intracellular$molecular_function[is.na(intracellular$molecular_function)] <- "Unannotated"
intracellular$cellular_component[is.na(intracellular$cellular_component)] <- "Unannotated"

# Add labels to differentiate supernatant and intracellular
supernatant$Location <- "Supernatant"
intracellular$Location <- "Intracellular"

# **Combine both datasets**
data_combined <- rbind(supernatant, intracellular)

# **Convert GO categories into long format**
go_data <- data_combined %>%
  select(Entry, biological_process, molecular_function, cellular_component, Location) %>%
  pivot_longer(cols = c(biological_process, molecular_function, cellular_component), 
               names_to = "GO_Category", values_to = "GO_Term")

# **Split multiple GO terms separated by ";" into separate rows**
go_data <- go_data %>%
  separate_rows(GO_Term, sep = ";") %>%
  mutate(GO_Term = str_trim(GO_Term))  # Trim spaces

# **Count occurrences per GO term per location**
go_counts <- go_data %>%
  group_by(GO_Term, GO_Category, Location) %>%
  summarise(Count = n(), .groups = "drop")

# **Identify shared and unique GO terms**
shared_terms <- go_counts %>%
  group_by(GO_Term, GO_Category) %>%
  filter(n_distinct(Location) > 1)  # Keep separate rows

supernatant_unique <- go_counts %>%
  filter(Location == "Supernatant" & !GO_Term %in% shared_terms$GO_Term) %>%
  mutate(Type = "Supernatant Only")

intracellular_unique <- go_counts %>%
  filter(Location == "Intracellular" & !GO_Term %in% shared_terms$GO_Term) %>%
  mutate(Type = "Intracellular Only")

go_final <- bind_rows(supernatant_unique, intracellular_unique)

# **Structure shared data for stacking**
shared_data <- shared_terms %>%
  pivot_wider(names_from = Location, values_from = Count, values_fill = list(Count = 0)) %>%
  rename(Int_Count = Intracellular, Sup_Count = Supernatant) %>%
  pivot_longer(cols = c(Int_Count, Sup_Count), names_to = "Type", values_to = "Count") %>%
  mutate(Type = ifelse(Type == "Int_Count", "Intracellular", "Supernatant"))


# **Separate Data by GO Categories (BP, MF, CC)**
intracellular_bp <- go_final %>% filter(Type == "Intracellular Only", GO_Category == "biological_process")
intracellular_mf <- go_final %>% filter(Type == "Intracellular Only", GO_Category == "molecular_function")
intracellular_cc <- go_final %>% filter(Type == "Intracellular Only", GO_Category == "cellular_component")

supernatant_bp <- go_final %>% filter(Type == "Supernatant Only", GO_Category == "biological_process")
supernatant_mf <- go_final %>% filter(Type == "Supernatant Only", GO_Category == "molecular_function")
supernatant_cc <- go_final %>% filter(Type == "Supernatant Only", GO_Category == "cellular_component")

shared_bp <- shared_data %>% filter(GO_Category == "biological_process")
shared_mf <- shared_data %>% filter(GO_Category == "molecular_function")
shared_cc <- shared_data %>% filter(GO_Category == "cellular_component")


color_map <- c("Supernatant Only" = "#298c8c", "Intracellular Only" = "#f1a226", 
               "Supernatant" = "#298c8c", "Intracellular" = "#f1a226")

custom_theme <- theme_minimal(base_family = "Arial") +
  theme(
    text = element_text(family = "Arial", color = "black"),   # Serif font & black text
    plot.title = element_text(size = 16, face = "bold"),      # Bold plot titles
    axis.title = element_text(size = 14, face = "bold"),      # Bold axis titles
    axis.text = element_text(size = 14, margin = margin(r = 10)),  # Prevent text overlap
    axis.ticks = element_line(color = "black"),               # Black ticks
    panel.border = element_blank(),                           # Remove all borders first
    axis.line.x = element_line(color = "black", linewidth = 1),  # Bottom border only
    axis.line.y = element_line(color = "black", linewidth = 1),  # Left border only
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), # Remove grid lines
    panel.background = element_rect(fill = "transparent", color = NA), # Transparent background
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    legend.text = element_text(size = 12),  # Legend font size
    legend.title = element_text(size = 14, face = "bold")  # Legend title font size
  )

# **Plot Function**
plot_GO_category <- function(data, title, filename, hide_low_counts=FALSE, stack=FALSE) {
  if (nrow(data) == 0) {
    message(paste("Skipping empty plot:", title))
    return(NULL)
  }
  
  # Apply filter for display (only in intracellular plots)
  if (hide_low_counts) {
    data <- data %>% filter(Count > 3)
  }
  
  # Create the plot
  p <- ggplot(data, aes(x = reorder(GO_Term, Count), y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = ifelse(stack, "stack", "dodge")) +
    scale_fill_manual(values = color_map) +
    coord_flip() +
    labs(title = title, x = "GO Terms", y = "Protein Count", fill = "GO Term Origin") +
    custom_theme +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 4, color = "black")
  
  
  # Save plot as high-resolution PNG (DPI 600)
  ggsave(filename, p, width = 24, height = 15, dpi = 600)
  message(paste("Saved:", filename))
}

# **Save Plots**
plot_GO_category(intracellular_bp, "Intracellular Only - Biological Process", "results/go_plots/Int_BP.png", hide_low_counts=TRUE)
plot_GO_category(intracellular_mf, "Intracellular Only - Molecular Function", "results/go_plots/Int_MF.png", hide_low_counts=TRUE)
plot_GO_category(intracellular_cc, "Intracellular Only - Cellular Component", "results/go_plots/Int_CC.png", hide_low_counts=TRUE)

plot_GO_category(supernatant_bp, "Supernatant Only - Biological Process", "results/go_plots/Sup_BP.png")
plot_GO_category(supernatant_mf, "Supernatant Only - Molecular Function", "results/go_plots/Sup_MF.png")
plot_GO_category(supernatant_cc, "Supernatant Only - Cellular Component", "results/go_plots/Sup_CC.png")

plot_GO_category(shared_bp, "Shared - Biological Process", "results/go_plots/Shared_BP.png", stack=TRUE)
plot_GO_category(shared_mf, "Shared - Molecular Function", "results/go_plots/Shared_MF.png", stack=TRUE)
plot_GO_category(shared_cc, "Shared - Cellular Component", "results/go_plots/Shared_CC.png", stack=TRUE)


protein_location_mapping <- data_combined %>%
  group_by(Entry) %>%
  summarise(Location = case_when(
    any(Location == "Supernatant") & any(Location == "Intracellular") ~ "Shared",
    any(Location == "Supernatant") ~ "Supernatant",
    any(Location == "Intracellular") ~ "Intracellular",
    TRUE ~ "Uncategorized"
  ), .groups = "drop")


output_csv <- data_combined %>%
  select(Entry, biological_process, molecular_function, cellular_component) %>%
  left_join(protein_location_mapping, by = "Entry")  
write.csv(output_csv, "results/tables/GO_Annotations_Combined.csv", row.names = FALSE)

plot_GO_category_return <- function(data, title, hide_low_counts=FALSE, stack=FALSE) {
  if (nrow(data) == 0) return(NULL)
  if (hide_low_counts) data <- data %>% filter(Count > 3)
  
  p <- ggplot(data, aes(x = reorder(GO_Term, Count), y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = ifelse(stack, "stack", "dodge")) +
    scale_fill_manual(values = color_map) +
    coord_flip() +
    labs(title = title, x = "GO Terms", y = "Protein Count", fill = "GO Term Origin") +
    custom_theme +
    geom_text(aes(label = Count), position = position_stack(vjust = 0.5), size = 3, color = "black")
  
  return(p)
}
# Intracellular
pA <- plot_GO_category_return(intracellular_bp, "A. Intracellular Only - Biological Process", hide_low_counts=TRUE)
pB <- plot_GO_category_return(intracellular_cc, "B. Intracellular Only - Cellular Component", hide_low_counts=TRUE)
pC <- plot_GO_category_return(intracellular_mf, "C. Intracellular Only - Molecular Function", hide_low_counts=TRUE)

# Supernatant
pD <- plot_GO_category_return(supernatant_bp, "D. Supernatant Only - Biological Process")
pE <- plot_GO_category_return(supernatant_mf, "E. Supernatant Only - Molecular Function")

# Shared (stacked)
pF <- plot_GO_category_return(shared_bp, "F. Shared - Biological Process", stack=TRUE)
pG <- plot_GO_category_return(shared_cc, "G. Shared - Cellular Component", stack=TRUE)
pH <- plot_GO_category_return(shared_mf, "H. Shared - Molecular Function", stack=TRUE)

combined_plot <- (
  (pA + pB + pC) /
    (pD + pE + pF) /
    (pG + pH)
) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

print(combined_plot)

ggsave("results/plots/Combined_GO_Figure.png", combined_plot, width = 40, height = 48, dpi = 600, limitsize = FALSE)
