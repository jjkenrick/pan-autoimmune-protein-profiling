source('fig_functions.R')

####   themes, colors, abbreviations   ####

disease_colors <- c(
  "Healthy" = "#4B5446",
  "Infectious" = "#E5E7EB",
  "Myositis" = "#B7CBE4",
  "Rheumatoid arthritis" = "#F1B6DA",
  "Sjögrens syndrome" = "#9C755F",
  "Systemic lupus erythematosus" = "#69B193",
  "Systemic sclerosis" = "#734F96"
)

disease_colors_abbrev <- c(
  "HC" = "#4B5446",
  "IC" = "#E5E7EB",
  "IIM" = "#B7CBE4",
  "RA" = "#F1B6DA",
  "SjD" = "#9C755F",
  "SLE" = "#69B193",
  "SSc" = "#734F96"
)

disease_abbreviations <- c(
  "Healthy" = "HC",
  "Infectious" = "IC",
  "Myositis" = "IIM",
  "Rheumatoid arthritis" = "RA",
  "Sjögrens syndrome" = "SjD",
  "Systemic lupus erythematosus" = "SLE",
  "Systemic sclerosis" = "SSc"
)

disease_order <- c(
  "Healthy",
  "Infectious",
  "Myositis",
  "Rheumatoid arthritis",
  "Sjögrens syndrome",
  "Systemic lupus erythematosus",
  "Systemic sclerosis"
)

theme_main <- theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    panel.spacing = unit(0.2, "lines"), 
                    panel.background=element_rect(fill="white"),
                    panel.border = element_blank(),
                    plot.title = element_text(face = "bold",
                                              size = rel(1), hjust = 0.5),
                    plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.ticks = element_line(),
                    axis.ticks.length = unit(.25, "cm"),
                    axis.line = element_line(size = 0.5),
                    axis.text = element_text(size = rel(2.0), color = 'black'), # Increased from 1.5 to 2.0
                    legend.key = element_blank(),
                    legend.position = "right",
                    legend.text = element_text(size=rel(1.8)),
                    legend.key.size= unit(0.7, "cm"),
                    legend.title = element_text(size=rel(1)),
                    plot.margin=unit(c(10,5,5,5),"mm"),
                    strip.background=element_rect(colour="grey90",fill="grey90"),
                    strip.text = element_text(face="bold"))

theme_simple <- 
  theme_main +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1,
                                   size = rel(1)), # Added size specification to maintain larger text
        strip.background = element_rect(color="white", fill="white"))


### Figure 1 ####

disease_meta <- create_disease_metadata(
  disease_order,
  disease_abbreviations, 
  disease_colors
)
fig1 <- create_disease_plots(long_data, disease_meta)
fig1a <- fig1$age_dist + fig1$sample_counts

# Print both volcano plots
vol1 <- create_labeled_volcano_plot(DA_results$group_comparisons$autoimmune_vs_healthy$p_0.01,'Autoimmune vs Healthy',
                                    pval_threshold = 0.01,
                                    n_labels_up = 10, 
                                    fc_threshold = 0.75,
                                    n_labels_down = 5,
                                    show.legend = F)
vol2 <- create_labeled_volcano_plot(DA_results$group_comparisons$autoimmune_vs_infections$p_0.01,'Autoimmune vs Infections',
                                    pval_threshold = 0.01,
                                    n_labels_up = 10, 
                                    fc_threshold = 0.25,
                                    n_labels_down = 5,
                                    show.legend = T)

vol1+vol2


# gene set enrichment
filtered_ai_healthy <- DA_results$group_comparisons$autoimmune_vs_healthy$p_0.01$results |> 
  filter(adj.P.Val < 0.01 & logFC > 0.75) |> 
  pull(protein)

go_results_AI_healthy <- do_jk_gsea(
  DE_genes = filtered_ai_healthy,
  adjust_method = "BH",
  cut_off = 0.05,
  go_or_kegg = "go_bp"
)

go_results_AI_healthy |> cnetplot()
barplot(go_results_AI_healthy, showCategory=25) 

# same for infectious
filtered_ai_infectious <- DA_results$group_comparisons$autoimmune_vs_infections$p_0.01$results |> 
  filter(adj.P.Val < 0.05 & logFC > 0.25) |> 
  pull(protein)

go_results_AI_infectious <- do_jk_gsea(
  DE_genes = filtered_ai_infectious,
  adjust_method = "BH",
  cut_off = 0.05,
  go_or_kegg = "go_bp"
)





###########       Figure 2        ############
process_summary <- function(da_results, disease_name) {
  # Get the DA results table from the new format
  da_table <- da_results
  
  # Count up and down regulated genes using the new 'regulation' column
  up_prots <- sum(da_table$regulation == "Higher", na.rm = TRUE)
  down_prots <- sum(da_table$regulation == "Lower", na.rm = TRUE)
  
  tibble(
    direction = c("Higher", "Lower"),
    count = c(up_prots, down_prots),
    Disease = disease_name
  )
}

# Function to get significant proteins for each disease
get_sig_proteins <- function(results, p_cutoff = 0.01) {
  results |> 
    filter(adj.P.Val < p_cutoff) |> 
    pull(protein)
}

# Create summaries using the new DA results format
all_summaries <- bind_rows(
  process_summary(DA_results$individual_diseases$Myositis$vs_other_auto$p_0.01$results, "Myositis"),
  process_summary(DA_results$individual_diseases$`Systemic sclerosis`$vs_other_auto$p_0.01$results, "Systemic sclerosis"),
  process_summary(DA_results$individual_diseases$`Rheumatoid arthritis`$vs_other_auto$p_0.01$results, "Rheumatoid arthritis"),
  process_summary(DA_results$individual_diseases$`Sjögrens syndrome`$vs_other_auto$p_0.01$results, "Sjögrens syndrome"),
  process_summary(DA_results$individual_diseases$`Systemic lupus erythematosus`$vs_other_auto$p_0.01$results, "Systemic lupus erythematosus")
) |>
  mutate(Disease_abbrev = disease_abbreviations[Disease])

# Create the plot
### barplot ###
fig2a <- ggplot(all_summaries, aes(x = Disease_abbrev, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c("Lower" = "#3B4992", "Higher" = "#EE0000"),
    name = "Direction"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(angle = 45, hjust = 1, size = 16),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.text = element_text(size = 10)
  ) +
  labs(
    x = "Disease",
    y = "Number of Differentially Abundant Proteins"
  )

sig_proteins <- list(
  IIM = get_sig_proteins(DA_results$individual_diseases$Myositis$vs_other_auto$p_0.01$results),
  SLE = get_sig_proteins(DA_results$individual_diseases$`Systemic lupus erythematosus`$vs_other_auto$p_0.01$results),
  SjD = get_sig_proteins(DA_results$individual_diseases$`Sjögrens syndrome`$vs_other_auto$p_0.01$results),
  RA = get_sig_proteins(DA_results$individual_diseases$`Rheumatoid arthritis`$vs_other_auto$p_0.01$results),
  SSc = get_sig_proteins(DA_results$individual_diseases$`Systemic sclerosis`$vs_other_auto$p_0.01$results)
)



# Create binary matrix for UpSet plot
all_proteins <- unique(unlist(sig_proteins))
binary_matrix <- matrix(0, nrow = length(all_proteins), 
                        ncol = length(sig_proteins))
colnames(binary_matrix) <- names(sig_proteins)
rownames(binary_matrix) <- all_proteins

for(disease in names(sig_proteins)) {
  binary_matrix[sig_proteins[[disease]], disease] <- 1
}

# Convert to data frame for UpSetR
binary_df <- as.data.frame(binary_matrix)

# Create UpSet plot
fig2b <- upset(binary_df,
               nsets = 5,
               order.by = "freq",
               main.bar.color = "navy",
               sets.bar.color = "darkred",
               point.size = 3,
               line.size = 1,
               mainbar.y.label = "Protein Intersections",
               sets.x.label = "Proteins per Disease",
               text.scale = 2,
               matrix.color = "darkblue")


vol_SSc <- create_enhanced_autoimmune_volcano(
  disease_results = list(
    vs_other_auto = DA_results$individual_diseases$`Systemic sclerosis`$vs_other_auto$p_0.01,
    vs_healthy = DA_results$individual_diseases$`Systemic sclerosis`$vs_healthy$p_0.01,
    vs_infections = DA_results$individual_diseases$`Systemic sclerosis`$vs_infections$p_0.01
  ),
  disease_name = "Systemic sclerosis",
  custom_proteins = c("KLK4", "MCAM", "CCN3","CD1C","KIT","CA6","FABP9","KLK13","CST6","COL4A1","ITGAV","ROBO2","IGFBP1","ROBO2","HSPG2","CDH6","LAMA4"),
  label_mode = "both",
  max_labels = 35
)
#### Figure 3 ####
panel_sle <- c('MFAP5','TNFSF11','PCSK9') 
panel_iim <- c('HSPB6','NOS1','CA3')
panel_ssc <- c('CCN3','KLK4','MCAM')
panel_sjd <- c('IRAG2','CRACR2A','TRAF2')
panel_ra <- c('CRTAC1','COMP','IL6')
p_iim <- do_many_boxplots_AI(
  data_long = long_data, 
  protein_list = panel_iim, 
  colorpal = disease_colors,
  ncol = 3,
  yourtitle = "IIM", 
  disease_order = disease_order,
  display_type = 'abbrev',
  disease_abbreviations = disease_abbreviations,
  highlight_disease_abbrev = "IIM",
  show_x_axis = FALSE,  # Hide x-axis
  point_size = 0.5      # Smaller points
)

p_ra <- do_many_boxplots_AI(
  data_long = long_data, 
  protein_list = panel_ra, 
  colorpal = disease_colors,
  ncol = 3,
  yourtitle = "RA", 
  disease_order = disease_order,
  display_type = 'abbrev',
  disease_abbreviations = disease_abbreviations,
  highlight_disease_abbrev = "RA",
  show_x_axis = FALSE,  # Hide x-axis
  point_size = 0.5      # Smaller points
)

p_sjd <- do_many_boxplots_AI(
  data_long = long_data, 
  protein_list = panel_sjd, 
  colorpal = disease_colors,
  ncol = 3,
  yourtitle = "SjD", 
  disease_order = disease_order,
  display_type = 'abbrev',
  disease_abbreviations = disease_abbreviations,
  highlight_disease_abbrev = "SjD",
  show_x_axis = FALSE,  # Hide x-axis
  point_size = 0.5      # Smaller points
)

p_sle <- do_many_boxplots_AI(
  data_long = long_data, 
  protein_list = panel_sle, 
  colorpal = disease_colors,
  ncol = 3,
  yourtitle = "SLE", 
  disease_order = disease_order,
  display_type = 'abbrev',
  disease_abbreviations = disease_abbreviations,
  highlight_disease_abbrev = "SLE",
  show_x_axis = FALSE,  # Hide x-axis
  point_size = 0.5      # Smaller points
)

p_ssc <- do_many_boxplots_AI(
  data_long = long_data, 
  protein_list = panel_ssc, 
  colorpal = disease_colors,
  ncol = 3,
  yourtitle = "SSc", 
  disease_order = disease_order,
  display_type = 'abbrev',
  disease_abbreviations = disease_abbreviations,
  highlight_disease_abbrev = "SSc",
  show_x_axis = TRUE,   # Show x-axis only on the bottom plot
  point_size = 0.5      # Smaller points
)

# Combine plots and save with increased height
combined_plot <- p_iim / p_ra / p_sjd / p_sle / p_ssc

#### Figure 4 ####
auc_summary <- ML_results$rocs %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    x = 0.75,  # Add fixed x coordinate for ggplot
    y = 0.25   # Add fixed y coordinate for ggplot
  ) |> 
  dplyr::ungroup()

fig4a <- ggplot(ML_results$rocs, 
                aes(x = 1 - Specificity, 
                    y = Sensitivity, 
                    color = fold,
                    group = interaction(fold, Group))) +
  geom_line(alpha = 0.7, size = 0.8) +  # Increased line thickness and opacity
  facet_wrap(~Group, labeller = labeller(Group = disease_abbreviations)) +  # Use abbreviated labels
  geom_text(data = auc_summary,
            aes(x = x, y = y, 
                label = sprintf("Mean AUC: %.3f", mean_auc)),
            color = "black",
            inherit.aes = FALSE) +
  theme_minimal() +
  scale_color_viridis_d(option = "turbo") +  # Using a more vibrant color palette
  labs(title = "ROC Curves by Disease and Fold",
       color = "Fold") +
  theme(
    panel.grid.minor = element_blank(),  # Remove minor grid lines for cleaner look
    strip.text = element_text(face = "bold"),  # Bold facet labels
    legend.position = "bottom"  # Move legend to bottom
  )

conf_matrix <- ML_results$predictions %>%
  group_by(fold) %>%
  count(Disease, .pred_class) %>%
  group_by(Disease, .pred_class) %>%
  summarise(n = round(mean(n), 1)) %>%
  ungroup() %>%
  mutate(Disease = recode(Disease, !!!disease_abbreviations),
         .pred_class = recode(.pred_class, !!!disease_abbreviations)) %>%
  # Add the following lines:
  complete(Disease, .pred_class, fill = list(n = 0))
fig4b <- ggplot(conf_matrix, 
                aes(x = Disease, 
                    y = .pred_class, 
                    fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), 
            color = ifelse(conf_matrix$n > max(conf_matrix$n)/2, "white", "black"),size = 5) +
  scale_fill_gradient(low = "white", high = "navy") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 8)) +
  labs(title = "Average Confusion Matrix Across Folds",
       x = "True Disease",
       y = "Predicted Disease",
       fill = "Count")+
  theme_main

fig4b <- ggplot(conf_matrix,
                aes(x = Disease,
                    y = .pred_class,
                    fill = n)) +
  geom_tile() +
  geom_text(aes(label = n),
            color = ifelse(conf_matrix$n > max(conf_matrix$n)/2, "white", "black"),
            size = 6) + # Set a specific absolute size for the numbers inside the tiles
  scale_fill_gradient(low = "white", high = "navy") +
  theme_minimal() + # Or theme_simple if you prefer its base settings
  labs(title = "Average Confusion Matrix Across Folds",
       x = "True Disease",
       y = "Predicted Disease",
       fill = "Count") +
  theme_main + # Apply your main theme first
  theme(
    # Override specific elements from theme_main with absolute sizes for this plot
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # Larger plot title
    axis.title.x = element_text(size = 16, face = "bold"),           # Larger x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),           # Larger y-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),    # Larger x-axis labels
    axis.text.y = element_text(size = 14),                           # Larger y-axis labels
    legend.title = element_text(size = 14, face = "bold"),           # Larger legend title
    legend.text = element_text(size = 12)                            # Larger legend text
  )


# 3. Important Proteins Bar Plot
important_proteins <- ML_results$protein_importances %>%
  mutate(class = factor(class, levels = disease_order)) %>%
  group_by(term, class) %>%
  summarise(
    mean_importance = mean(abs(estimate)),
    sd_importance = sd(abs(estimate)),
    se_importance = sd_importance / sqrt(n()),
    ci_lower = mean_importance - 1.96 * se_importance,
    ci_upper = mean_importance + 1.96 * se_importance,
    .groups = 'drop'
  ) %>%
  group_by(class) %>%
  slice_max(order_by = mean_importance, n = 10) %>%
  mutate(class = recode(class, !!!disease_abbreviations))

fig4c <- ggplot(important_proteins,
                aes(x = reorder(term, mean_importance), 
                    y = mean_importance,
                    fill = class)) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2) +
  facet_wrap(~class, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = disease_colors_abbrev) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 8)
  ) +
  labs(title = "Top 10 Important Proteins by Disease",
       x = "Protein",
       y = "Mean Absolute Importance ± 95% CI")+
  theme(
    # Override specific elements from theme_main with absolute sizes for this plot
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5), # Larger plot title
    axis.title.x = element_text(size = 16, face = "bold"),           # Larger x-axis title
    axis.title.y = element_text(size = 16, face = "bold"),           # Larger y-axis title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),    # Larger x-axis labels
    axis.text.y = element_text(size = 14),                           # Larger y-axis labels
    legend.title = element_text(size = 14, face = "bold"),           # Larger legend title
    legend.text = element_text(size = 14)                            # Larger legend text
  )

fig4c <- ggplot(important_proteins,
                aes(x = reorder(term, mean_importance),
                    y = mean_importance,
                    fill = class)) +
  geom_col() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.2) +
  facet_wrap(~class, scales = "free_y") +
  coord_flip() +
  scale_fill_manual(values = disease_colors_abbrev) +
  theme_minimal() + # Starting with theme_minimal
  # Merge your theme customizations into one block for clarity and proper overriding
  theme(
    # Facet label text (disease abbreviations)
    strip.text = element_text(size = 14, face = "bold"), # Increased size and added bold
    
    # Axis text (protein names on y-axis, importance on x-axis)
    axis.text.y = element_text(size = 12), # Increased size for protein names
    axis.text.x = element_text(size = 12), # Added or adjusted size for importance numbers
    
    # Axis titles
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    # Plot title
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    
    # Legend titles and text (if a legend is shown, though fill is faceted)
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  labs(title = "Top 10 Important Proteins by Disease",
       x = "Protein",
       y = "Mean Absolute Importance ± 95% CI")
# Removed the duplicate theme() block at the end.
# If you want to explicitly apply theme_main, it should go BEFORE these specific overrides.
# Example: theme_minimal() + theme_main + theme( ... overrides ... )
# However, theme_minimal() already sets a good base, and the explicit sizes below override it.

# Get the top 10 proteins for each disease class
top_proteins_by_class <- important_proteins %>%
  group_by(class) %>%
  top_n(10, mean_importance) |> 
  ungroup()

# Find proteins that appear in multiple classes
overlapping_proteins <- top_proteins_by_class %>%
  count(term) %>%
  filter(n > 1) %>%
  arrange(desc(n))

print(overlapping_proteins)
# 4. Disease Overlap UpSet Plot
# First, create binary matrix for disease presence
important_proteins_by_disease <- ML_results$protein_importances %>%
  group_by(term, class) %>%
  summarise(
    mean_importance = mean(abs(estimate)),
    .groups = 'drop'
  ) %>%
  filter(mean_importance > 0.0) %>%  # Small threshold to filter out noise
  mutate(class = recode(class, !!!disease_abbreviations))

# Create binary matrix for protein presence across diseases
upset_matrix_ml <- important_proteins_by_disease %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = term,
    names_from = class,
    values_from = present,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(upset_matrix_ml) <- upset_matrix_ml$term
upset_matrix_ml <- upset_matrix_ml[,-1]  # Remove the term column

# Create upset plot
#pdf("upset_plot.pdf", width = 12, height = 8)  # Made plot larger to accommodate more intersections
fig4d <- upset(upset_matrix_ml, 
               nsets = ncol(upset_matrix_ml),
               order.by = "freq",
               main.bar.color = "navy",
               sets.bar.color = "darkred",
               # text.scale = 3,
               text.scale = c(2, 2, 2, 1.8, 2, 2.5),
               point.size = 2.5,
               line.size = 1,
               mainbar.y.label = "Number of Shared Proteins",
               sets.x.label = "Proteins per Disease",
               nintersects = NA)  # Show all intersections

## UMAP 
AI_wide <- wide_data |>  filter(Disease %in% AUTOIMMUNE_DISEASES)
prepared_data <- prepare_pca_data(AI_wide)
pca_results <- run_pca(prepared_data)
loadings <- process_pc_loadings(pca_results, n_components = 30)
pc_scores <- pca_results$scores %>% select(matches("^PC\\d+$")) %>% as.matrix()

print(paste("PC scores dimensions (all proteins):", dim(pc_scores)[1], "x", dim(pc_scores)[2]))
## scatter and scree plot code in pca_umap.R script
pca_df <- data.frame(
  PC1 = pca_results$scores$PC1,
  PC2 = pca_results$scores$PC2,
  Disease = pca_results$scores$Disease
)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Disease)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Plot All Proteins",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()

umap_coords <- uwot::umap(pc_scores, 
                          n_neighbors = 15, 
                          min_dist = 0.1, 
                          n_components = 2,
                          metric = "euclidean")

# Create the plotting dataframe
umap_plot_data <- data.frame(
  UMAP1 = umap_coords[,1],
  UMAP2 = umap_coords[,2],
  Disease = pca_results$scores$Disease
)

umap_all <- ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, color = Disease)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = disease_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Disease"
  ) +
  theme_simple+
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )
top50_glm_prots <- important_proteins$term
## GLM filtering ##
glm_filtered <- wide_data %>% filter(Disease %in% AUTOIMMUNE_DISEASES) %>% select(DAid, Disease, Sex, Age, all_of(top50_glm_prots))
glm_prepared_data <- prepare_pca_data(glm_filtered)
glm_pca_results <- run_pca(glm_prepared_data)
glm_pc_scores <- glm_pca_results$scores %>% select(matches("^PC\\d+$")) %>% as.matrix()
glm_pca_df <- data.frame(
  PC1 = glm_pca_results$scores$PC1,
  PC2 = glm_pca_results$scores$PC2,
  Disease = glm_pca_results$scores$Disease
)

ggplot(glm_pca_df, aes(x = PC1, y = PC2, color = Disease)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Plot GLM important proteins",
    x = "Principal Component 1",
    y = "Principal Component 2"
  ) +
  theme_minimal()
print(paste("PC scores dimensions (GLM proteins):", dim(glm_pc_scores)[1], "x", dim(glm_pc_scores)[2]))

glm_umap_coords <- uwot::umap(glm_pc_scores, 
                              n_neighbors = 15, 
                              min_dist = 0.1, 
                              n_components = 2,
                              metric = "euclidean")

# Create the plotting dataframe
glm_plot_data <- data.frame(
  UMAP1 = glm_umap_coords[,1],
  UMAP2 = glm_umap_coords[,2],
  Disease = glm_pca_results$scores$Disease
)
umap_glm <- ggplot(glm_plot_data, aes(x = UMAP1, y = UMAP2, color = Disease)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_manual(values = disease_colors) +
  labs(
    x = "UMAP1",
    y = "UMAP2",
    color = "Disease"
  ) +
  theme_simple+
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    legend.position = c(1, 1),
    legend.justification = c("right", "top"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

both_umaps <- umap_all+ umap_glm


#### Supplementary Figures ####
#plot_method_comparison(comparison_results)

vol_IIM <- create_enhanced_autoimmune_volcano(
  disease_results = list(
    vs_other_auto = DA_results$individual_diseases$Myositis$vs_other_auto$p_0.01,
    vs_healthy = DA_results$individual_diseases$Myositis$vs_healthy$p_0.01,
    vs_infections = DA_results$individual_diseases$Myositis$vs_infections$p_0.01
  ),
  disease_name = "Myositis",
 # custom_proteins = c("KLK4", "MCAM", "CCN3","CD1C","KIT","CA6","FABP9","KLK13","CST6","COL4A1","ITGAV","ROBO2","IGFBP1","ROBO2","HSPG2","CDH6","LAMA4"),
  label_mode = "both",
  max_labels = 35
)
vol_SLE <- create_enhanced_autoimmune_volcano(
  disease_results = list(
    vs_other_auto = DA_results$individual_diseases$`Systemic lupus erythematosus`$vs_other_auto$p_0.01,
    vs_healthy = DA_results$individual_diseases$`Systemic lupus erythematosus`$vs_healthy$p_0.01,
    vs_infections = DA_results$individual_diseases$`Systemic lupus erythematosus`$vs_infections$p_0.01
  ),
  disease_name = "Systemic lupus erythematosus",
 # custom_proteins = c("KLK4", "MCAM", "CCN3","CD1C","KIT","CA6","FABP9","KLK13","CST6","COL4A1","ITGAV","ROBO2","IGFBP1","ROBO2","HSPG2","CDH6","LAMA4"),
  label_mode = "both",
  max_labels = 35
)
  
vol_RA <- create_enhanced_autoimmune_volcano(
  disease_results = list(
    vs_other_auto = DA_results$individual_diseases$`Rheumatoid arthritis`$vs_other_auto$p_0.01,
    vs_healthy = DA_results$individual_diseases$`Rheumatoid arthritis`$vs_healthy$p_0.01,
    vs_infections = DA_results$individual_diseases$`Rheumatoid arthritis`$vs_infections$p_0.01
  ),
  disease_name = "Rheumatoid arthritis",
  custom_proteins = c("CRTAC1","COMP","IL6","ADGRG2","CES3","MTPN","ITGAM","CD1C","WASF3","TNFSF14","S100A12"),
  label_mode = "both",
  max_labels = 25
)

vol_SjD <- create_enhanced_autoimmune_volcano(
  disease_results = list(
    vs_other_auto = DA_results$individual_diseases$`Sjögrens syndrome`$vs_other_auto$p_0.01,
    vs_healthy = DA_results$individual_diseases$`Sjögrens syndrome`$vs_healthy$p_0.01,
    vs_infections = DA_results$individual_diseases$`Sjögrens syndrome`$vs_infections$p_0.01
  ),
  disease_name = "Sjögren's syndrome",
  custom_proteins = c("CNDP1",'IRAG2','CRACR2A','ITM2A','TRAF2','UXS1','STX16','DNMBP','CA4','CPA2','CBLN4','IQGAP2','ICAM4','CNTN1','PRDX1','ADAMTS15'),
  label_mode = "both",
  max_labels = 25
)

sjdprots_vsuto <- DA_results$individual_diseases$`Sjögrens syndrome`$vs_other_auto$p_0.01$results |> filter(logFC>0.25) |> select(protein)
## overlap proteins


# Function to analyze overlapping significant proteins across diseases
analyze_significant_protein_overlap <- function(DA_results, p_threshold = 0.01, logfc_threshold = 0.25) {
  diseases <- names(DA_results$individual_diseases)
  
  # Initialize lists to store significant proteins for each disease
  up_regulated <- list()
  down_regulated <- list()
  all_significant <- list()
  
  # Extract significant proteins for each disease
  for (disease in diseases) {
    results <- DA_results$individual_diseases[[disease]]$vs_healthy$p_0.01$results
    
    # Get up-regulated proteins (significant P-value and logFC > threshold)
    up <- results$protein[results$adj.P.Val < p_threshold & results$logFC > logfc_threshold]
    up_regulated[[disease]] <- up
    
    # Get down-regulated proteins (significant P-value and logFC < -threshold)
    down <- results$protein[results$adj.P.Val < p_threshold & results$logFC < -logfc_threshold]
    down_regulated[[disease]] <- down
    
    # All significant proteins
    all_significant[[disease]] <- c(up, down)
  }
  
  # Display counts
  protein_counts <- data.frame(
    Disease = diseases,
    Up_regulated = sapply(up_regulated, length),
    Down_regulated = sapply(down_regulated, length),
    Total_significant = sapply(all_significant, length)
  )
  
  print("Number of significant proteins per disease:")
  print(protein_counts)
  
  # Find overlaps in up-regulated proteins across diseases
  cat("\nOverlaps in UP-regulated proteins across diseases:\n")
  for (i in 1:(length(diseases)-1)) {
    for (j in (i+1):length(diseases)) {
      disease1 <- diseases[i]
      disease2 <- diseases[j]
      
      overlap <- intersect(up_regulated[[disease1]], up_regulated[[disease2]])
      
      if (length(overlap) > 0) {
        cat(sprintf("\n%s AND %s: %d proteins\n", disease1, disease2, length(overlap)))
        cat("Examples: ", paste(head(overlap, 5), collapse=", "), 
            if(length(overlap) > 5) "..." else "", "\n")
      }
    }
  }
  
  # Find overlaps in down-regulated proteins across diseases
  cat("\nOverlaps in DOWN-regulated proteins across diseases:\n")
  for (i in 1:(length(diseases)-1)) {
    for (j in (i+1):length(diseases)) {
      disease1 <- diseases[i]
      disease2 <- diseases[j]
      
      overlap <- intersect(down_regulated[[disease1]], down_regulated[[disease2]])
      
      if (length(overlap) > 0) {
        cat(sprintf("\n%s AND %s: %d proteins\n", disease1, disease2, length(overlap)))
        cat("Examples: ", paste(head(overlap, 5), collapse=", "), 
            if(length(overlap) > 5) "..." else "", "\n")
      }
    }
  }
  
  # Let's also find proteins that appear across all diseases
  common_up <- Reduce(intersect, up_regulated)
  common_down <- Reduce(intersect, down_regulated)
  
  cat("\nProteins UP-regulated across ALL diseases:\n")
  if (length(common_up) > 0) {
    print(common_up)
  } else {
    cat("None\n")
  }
  
  cat("\nProteins DOWN-regulated across ALL diseases:\n")
  if (length(common_down) > 0) {
    print(common_down)
  } else {
    cat("None\n")
  }
  
  # Create a summary of overlaps for visualization
  overlap_matrix_up <- matrix(0, nrow = length(diseases), ncol = length(diseases))
  overlap_matrix_down <- matrix(0, nrow = length(diseases), ncol = length(diseases))
  
  for (i in 1:length(diseases)) {
    for (j in 1:length(diseases)) {
      if (i != j) {
        disease1 <- diseases[i]
        disease2 <- diseases[j]
        
        overlap_up <- intersect(up_regulated[[disease1]], up_regulated[[disease2]])
        overlap_down <- intersect(down_regulated[[disease1]], down_regulated[[disease2]])
        
        overlap_matrix_up[i, j] <- length(overlap_up)
        overlap_matrix_down[i, j] <- length(overlap_down)
      }
    }
  }
  
  rownames(overlap_matrix_up) <- colnames(overlap_matrix_up) <- diseases
  rownames(overlap_matrix_down) <- colnames(overlap_matrix_down) <- diseases
  
  return(list(
    up_regulated = up_regulated,
    down_regulated = down_regulated,
    overlap_matrix_up = overlap_matrix_up,
    overlap_matrix_down = overlap_matrix_down,
    common_up = common_up,
    common_down = common_down
  ))
}

# Run the analysis
overlap_results <- analyze_significant_protein_overlap(DA_results)

# To visualize the overlaps, you could use a heatmap
if (requireNamespace("pheatmap", quietly = TRUE)) {
  cat("\nHeatmap of UP-regulated protein overlaps:\n")
  pheatmap::pheatmap(overlap_results$overlap_matrix_up, 
                     display_numbers = TRUE, 
                     main = "Number of proteins with higher levels than healthy")
  
  cat("\nHeatmap of DOWN-regulated protein overlaps:\n")
  pheatmap::pheatmap(overlap_results$overlap_matrix_down,
                     display_numbers = TRUE,
                     main = "Number of proteins with lower levels than healthy")
} else {
  cat("\nInstall the 'pheatmap' package to visualize the overlap heatmaps\n")
}

svg('supplementary_data/hm1.svg')
hm1
dev.off()
hm1 <- pheatmap::pheatmap(overlap_results$overlap_matrix_up, 
                                                display_numbers = TRUE, 
                                                main = "Number of proteins with higher levels than healthy")
hm2 <-   pheatmap::pheatmap(overlap_results$overlap_matrix_down,
                            display_numbers = TRUE,
                            main = "Number of proteins with lower levels than healthy")
svg('supplementary_data/hm2.svg')
hm2
dev.off()
# Check if certain proteins specifically appears in multiple diseases
check_specific_protein <- function(protein_name, overlap_results) {
  diseases <- names(overlap_results$up_regulated)
  
  cat(sprintf("\nChecking for '%s' across diseases:\n", protein_name))
  
  for (disease in diseases) {
    in_up <- protein_name %in% overlap_results$up_regulated[[disease]]
    in_down <- protein_name %in% overlap_results$down_regulated[[disease]]
    
    if (in_up) {
      cat(sprintf("  %s: PRESENT (UP-regulated)\n", disease))
    } else if (in_down) {
      cat(sprintf("  %s: PRESENT (DOWN-regulated)\n", disease))
    } else {
      cat(sprintf("  %s: NOT found\n", disease))
    }
  }
}

required_packages <- c("ggplot2", "UpSetR", "VennDiagram", "ComplexHeatmap", "eulerr", "cowplot")
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Package '%s' is not installed. Install it with: install.packages('%s')\n", pkg, pkg))
  }
}

# Enhanced visualization function
visualize_protein_overlaps <- function(overlap_results) {
  diseases <- names(overlap_results$up_regulated)
  
  # 1. UpSet plot for up-regulated proteins
  if(requireNamespace("UpSetR", quietly = TRUE)) {
    # Prepare data for UpSetR
    protein_sets_up <- list()
    for(disease in diseases) {
      protein_sets_up[[disease]] <- overlap_results$up_regulated[[disease]]
    }
    
    # Create UpSet plot
    upset_plot_up <- UpSetR::upset(UpSetR::fromList(protein_sets_up), 
                                   order.by = "freq", 
                                   nsets = length(diseases),
                                   mainbar.y.label = "Protein Intersections", 
                                   sets.x.label = "Proteins per Disease",
                                   text.scale = 1.2,
                                   point.size = 3.5,
                                   line.size = 1)
    print(upset_plot_up)
    
    # Same for down-regulated proteins
    protein_sets_down <- list()
    for(disease in diseases) {
      protein_sets_down[[disease]] <- overlap_results$down_regulated[[disease]]
    }
    
    upset_plot_down <- UpSetR::upset(UpSetR::fromList(protein_sets_down), 
                                     order.by = "freq", 
                                     nsets = length(diseases),
                                     mainbar.y.label = "Protein Intersections", 
                                     sets.x.label = "Proteins per Disease",
                                     text.scale = 1.2,
                                     point.size = 3.5,
                                     line.size = 1)
    print(upset_plot_down)
  } else {
    cat("Install 'UpSetR' package for intersection plots\n")
  }
  
  # 2. Venn diagram (if 5 or fewer diseases)
  if(length(diseases) <= 5 && requireNamespace("VennDiagram", quietly = TRUE)) {
    # Set up colors
    colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
    
    # Create Venn diagram for up-regulated proteins
    venn_up <- VennDiagram::venn.diagram(
      x = overlap_results$up_regulated,
      category.names = diseases,
      filename = NULL,
      output = TRUE,
      fill = colors[1:length(diseases)],
      alpha = 0.5,
      lwd = 2,
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans"
    )
    grid::grid.draw(venn_up)
    
    # Create Venn diagram for down-regulated proteins
    venn_down <- VennDiagram::venn.diagram(
      x = overlap_results$down_regulated,
      category.names = diseases,
      filename = NULL,
      output = TRUE,
      fill = colors[1:length(diseases)],
      alpha = 0.5,
      lwd = 2,
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans"
    )
    grid::grid.draw(venn_down)
  } else if(length(diseases) > 5) {
    cat("Venn diagrams are limited to 5 or fewer sets\n")
  } else {
    cat("Install 'VennDiagram' package for Venn diagrams\n")
  }
  
  # 3. Euler diagram (alternative to Venn)
  if(requireNamespace("eulerr", quietly = TRUE)) {
    # Create Euler diagram for up-regulated proteins
    euler_up <- eulerr::euler(overlap_results$up_regulated)
    plot(euler_up, quantities = TRUE, main = "Up-regulated Proteins")
    
    # Create Euler diagram for down-regulated proteins
    euler_down <- eulerr::euler(overlap_results$down_regulated)
    plot(euler_down, quantities = TRUE, main = "Down-regulated Proteins")
  } else {
    cat("Install 'eulerr' package for Euler diagrams\n")
  }
  
  # 4. Heatmap of common proteins
  if(requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
    # Get all unique proteins
    all_up_proteins <- unique(unlist(overlap_results$up_regulated))
    all_down_proteins <- unique(unlist(overlap_results$down_regulated))
    
    # Create presence/absence matrices
    matrix_up <- matrix(0, nrow = length(all_up_proteins), ncol = length(diseases))
    rownames(matrix_up) <- all_up_proteins
    colnames(matrix_up) <- diseases
    
    matrix_down <- matrix(0, nrow = length(all_down_proteins), ncol = length(diseases))
    rownames(matrix_down) <- all_down_proteins
    colnames(matrix_down) <- diseases
    
    # Fill matrices
    for(i in 1:length(diseases)) {
      disease <- diseases[i]
      matrix_up[overlap_results$up_regulated[[disease]], i] <- 1
      matrix_down[overlap_results$down_regulated[[disease]], i] <- 1
    }
    
    # Sort rows by sum (proteins present in most diseases first)
    matrix_up <- matrix_up[order(rowSums(matrix_up), decreasing = TRUE), ]
    matrix_down <- matrix_down[order(rowSums(matrix_down), decreasing = TRUE), ]
    
    # Filter to show only proteins present in at least 2 diseases
    matrix_up_filtered <- matrix_up[rowSums(matrix_up) >= 2, ]
    matrix_down_filtered <- matrix_down[rowSums(matrix_down) >= 2, ]
    
    # Create heatmaps
    hm_up <- ComplexHeatmap::Heatmap(matrix_up_filtered, 
                                     name = "Presence",
                                     column_title = "Up-regulated Proteins Present in Multiple Diseases",
                                     row_title = "Proteins",
                                     col = c("white", "#FF0000"),
                                     show_row_names = TRUE,
                                     cluster_rows = FALSE,
                                     cluster_columns = FALSE)
    
    hm_down <- ComplexHeatmap::Heatmap(matrix_down_filtered, 
                                       name = "Presence",
                                       column_title = "Down-regulated Proteins Present in Multiple Diseases",
                                       row_title = "Proteins",
                                       col = c("white", "#0000FF"),
                                       show_row_names = TRUE,
                                       cluster_rows = FALSE,
                                       cluster_columns = FALSE)
    
    ComplexHeatmap::draw(hm_up)
    ComplexHeatmap::draw(hm_down)
  } else {
    cat("Install 'ComplexHeatmap' package for heatmaps\n")
  }
  
  # 5. Bar chart of protein overlap counts
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    # Create data frame for plotting
    count_data <- data.frame(
      Disease_Pair = character(),
      Direction = character(),
      Count = integer(),
      stringsAsFactors = FALSE
    )
    
    # Calculate pairwise overlaps
    for(i in 1:(length(diseases)-1)) {
      for(j in (i+1):length(diseases)) {
        d1 <- diseases[i]
        d2 <- diseases[j]
        
        up_overlap <- length(intersect(overlap_results$up_regulated[[d1]], 
                                       overlap_results$up_regulated[[d2]]))
        
        down_overlap <- length(intersect(overlap_results$down_regulated[[d1]], 
                                         overlap_results$down_regulated[[d2]]))
        
        count_data <- rbind(count_data, 
                            data.frame(Disease_Pair = paste(d1, "-", d2),
                                       Direction = "Up-regulated",
                                       Count = up_overlap))
        
        count_data <- rbind(count_data, 
                            data.frame(Disease_Pair = paste(d1, "-", d2),
                                       Direction = "Down-regulated",
                                       Count = down_overlap))
      }
    }
    
    # Create bar plot
    g <- ggplot2::ggplot(count_data, ggplot2::aes(x = Disease_Pair, y = Count, fill = Direction)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_manual(values = c("Up-regulated" = "#FF6666", "Down-regulated" = "#6666FF")) +
      ggplot2::labs(title = "Number of Overlapping Proteins Between Disease Pairs",
                    x = "Disease Pair", y = "Number of Proteins")
    
    print(g)
  }
  
  # 6. Network diagram of protein-disease relationships
  if(requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE)) {
    # Create edge lists
    edges_up <- data.frame(Protein = character(), Disease = character(), Regulation = character())
    edges_down <- data.frame(Protein = character(), Disease = character(), Regulation = character())
    
    for(disease in diseases) {
      for(protein in overlap_results$up_regulated[[disease]]) {
        edges_up <- rbind(edges_up, data.frame(Protein = protein, 
                                               Disease = disease, 
                                               Regulation = "Up"))
      }
      
      for(protein in overlap_results$down_regulated[[disease]]) {
        edges_down <- rbind(edges_down, data.frame(Protein = protein, 
                                                   Disease = disease, 
                                                   Regulation = "Down"))
      }
    }
    
    # Filter to include only proteins that appear in multiple diseases
    protein_counts_up <- table(edges_up$Protein)
    shared_proteins_up <- names(protein_counts_up[protein_counts_up > 1])
    
    protein_counts_down <- table(edges_down$Protein)
    shared_proteins_down <- names(protein_counts_down[protein_counts_down > 1])
    
    edges_up_filtered <- edges_up[edges_up$Protein %in% shared_proteins_up, ]
    edges_down_filtered <- edges_down[edges_down$Protein %in% shared_proteins_down, ]
    
    # Create network plots
    if(length(shared_proteins_up) > 0) {
      g_up <- igraph::graph_from_data_frame(edges_up_filtered, directed = FALSE)
      
      # Set node types
      V(g_up)$type <- ifelse(V(g_up)$name %in% diseases, "Disease", "Protein")
      
      # Plot
      plot(g_up, 
           vertex.color = ifelse(V(g_up)$type == "Disease", "lightblue", "salmon"),
           vertex.size = ifelse(V(g_up)$type == "Disease", 15, 10),
           vertex.label.cex = 0.8,
           main = "Network of Up-regulated Proteins Shared Between Diseases")
    }
    
    if(length(shared_proteins_down) > 0) {
      g_down <- igraph::graph_from_data_frame(edges_down_filtered, directed = FALSE)
      
      # Set node types
      V(g_down)$type <- ifelse(V(g_down)$name %in% diseases, "Disease", "Protein")
      
      # Plot
      plot(g_down, 
           vertex.color = ifelse(V(g_down)$type == "Disease", "lightblue", "lightskyblue"),
           vertex.size = ifelse(V(g_down)$type == "Disease", 15, 10),
           vertex.label.cex = 0.8,
           main = "Network of Down-regulated Proteins Shared Between Diseases")
    }
  } else {
    cat("Install 'igraph' package for network diagrams\n")
  }
}

# Run the visualization function with our existing overlap results
supp_plots <- visualize_protein_overlaps(overlap_results)

pdf('supplementary_data/go_cnet_healthy.pdf')
go_results_AI_healthy |> cnetplot()
dev.off()

pdf('supplementary_data/go_bar_healthy.pdf')
barplot(go_results_AI_healthy) 
dev.off()

pdf('supplementary_data/go_dot_healthy.pdf')
dotplot(go_results_AI_healthy) 
dev.off()

pdf('supplementary_data/go_cnet_IC.pdf')
go_results_AI_infectious |> cnetplot()
dev.off()

pdf('supplementary_data/go_bar_IC.pdf')
barplot(go_results_AI_infectious) 
dev.off()

pdf('supplementary_data/go_dot_IC.pdf')
dotplot(go_results_AI_infectious) 
dev.off()

pdf('supplementary_data/go_infectious.pdf')
g4/g5+g6
dev.off()
svg('supplementary_data/venn_down.svg')
# Check if the VennDiagram package is installed and number of diseases is appropriate
if (requireNamespace("VennDiagram", quietly = TRUE) && length(diseases) <= 5) {
    # 1. Generate the Venn diagram object
    venn_down <- VennDiagram::venn.diagram(
      x = overlap_results$down_regulated,
      category.names = diseases,
      filename = NULL,
      output = TRUE,
      fill = colors[1:length(diseases)],
      alpha = 0.5,
      lwd = 2,
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans"
    )
    
    # 2. **CRITICAL CHANGE:** Start a new, blank plotting page
    grid::grid.newpage() 
    
    # 3. Draw the Venn diagram on the new, blank page
    grid::grid.draw(venn_down)
} else if(length(diseases) > 5) {
    cat("Venn diagrams are limited to 5 or fewer sets\n")
} else {
    cat("Install 'VennDiagram' package for Venn diagrams\n")
}
dev.off()

svg('supplementary_data/vol_sjd.svg',width = 16,height = 12)
vol_SjD
dev.off()
svg('supplementary_data/vol_sle.svg',width = 16,height = 12)
vol_SLE
dev.off()
svg('supplementary_data/vol_ra.svg',width = 16,height = 12)
vol_RA
dev.off()
