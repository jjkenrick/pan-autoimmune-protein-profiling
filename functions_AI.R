# Function script includes differential abundance analysis, GLMnet lasso multiclassification, comparison of both methods, and other functions
# Load required libraries
library(limma)
library(MatchIt)
library(tibble)
library(purrr)
library(ggbeeswarm) 
library(ggridges)
library(patchwork)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)
library(ggupset)
library(UpSetR)
library(data.table)
library(ggrepel)
library(tidyverse)
library(tidymodels)
library(themis)
library(multiROC)
library(pheatmap)
library(gprofiler2)
library(tidyverse)
library(stats)
library(cluster)
library(factoextra)
library(ggalluvial)
library(ComplexHeatmap)
library(viridis)
library(gridExtra)
library(uwot)
select <- dplyr::select # need to specify otherwise error unable to find inherited method for signature 'x=tbl_df'


# ===== Differential abundance =====
#' Core limma differential abundance analysis
#' @param protein_data Matrix with proteins as rows, samples as columns
#' @param design_matrix Design matrix for the analysis
#' @param contrast_string Contrast to test (e.g., "GroupA-GroupB")
#' @param p_threshold P-value threshold for significance
#' @param fc_threshold Log2 fold change threshold
#' @return List containing limma results
run_limma_analysis <- function(protein_data, design_matrix, contrast_string, 
                               p_threshold = 0.01, fc_threshold = 0.25) {
  
  # Fit linear model
  fit <- lmFit(protein_data, design = design_matrix)
  
  # Create and apply contrasts
  contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = design_matrix)
  contrast_fit <- contrasts.fit(fit, contrast_matrix)
  
  # Empirical Bayes
  ebayes_fit <- eBayes(contrast_fit)
  
  # Extract results
  results <- topTable(ebayes_fit, n = Inf, adjust.method = 'BH', confint = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("protein") %>%
    mutate(
      regulation = case_when(
        logFC >= fc_threshold & adj.P.Val < p_threshold ~ "Higher",
        logFC <= -fc_threshold & adj.P.Val < p_threshold ~ "Lower",
        TRUE ~ "Not significant"
      ),
      significance = adj.P.Val < p_threshold
    )
  
  # Summary statistics
  summary_stats <- summary(decideTests(ebayes_fit, p.value = p_threshold, lfc = fc_threshold))
  
  list(
    results = results,
    summary = summary_stats,
    design_matrix = design_matrix,
    contrast_matrix = contrast_matrix,
    parameters = list(
      p_threshold = p_threshold,
      fc_threshold = fc_threshold,
      contrast = contrast_string
    )
  )
}

# ===== DATA PREPARATION FUNCTIONS =====

#' Prepare data for differential abundance analysis (no matching)
#' @param wide_data Wide format data with samples as rows
#' @param case_group Name of case group
#' @param control_group Name of control group
#' @return Prepared data for analysis
prepare_comparison_data <- function(wide_data, case_group, control_group) {
  
  # Handle special group combinations
  if (case_group == "Autoimmune" && control_group == "Infection") {
    prepared_data <- wide_data %>%
      mutate(
        Group = case_when(
          Disease %in% AUTOIMMUNE_DISEASES ~ "Autoimmune",
          Disease %in% INFECTIONS ~ "Infection",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Group))
  } else if (case_group == "Autoimmune" && control_group == "Healthy") {
    prepared_data <- wide_data %>%
      mutate(
        Group = case_when(
          Disease %in% AUTOIMMUNE_DISEASES ~ "Autoimmune",
          Disease == "Healthy" ~ "Healthy",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Group))
  } else if (control_group == "Infections") {
    # Individual disease vs all infections grouped together
    prepared_data <- wide_data %>%
      mutate(
        Group = case_when(
          Disease == case_group ~ case_group,
          Disease %in% INFECTIONS ~ "Infections",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Group))
  } else if (case_group == "Infections") {
    # All infections grouped together vs another group
    prepared_data <- wide_data %>%
      mutate(
        Group = case_when(
          Disease %in% INFECTIONS ~ "Infections",
          Disease == control_group ~ control_group,
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Group))
  } else {
    # Standard disease vs disease comparison
    prepared_data <- wide_data %>%
      filter(Disease %in% c(case_group, control_group)) %>%
      mutate(Group = Disease)
  }
  
  # Clean metadata and make group names R-friendly
  prepared_data <- prepared_data %>%
    filter(!is.na(Sex), !is.na(Age), !is.na(Group)) %>%
    mutate(
      Sex = factor(Sex, levels = c("F", "M")),
      Age = as.numeric(Age),
      Group = factor(Group, levels = c(control_group, case_group)),
      Group_clean = factor(gsub(" ", "_", Group), 
                           levels = gsub(" ", "_", c(control_group, case_group)))
    )
  
  return(prepared_data)
}

# ===== MAIN COMPARISON FUNCTION =====

#' Run differential abundance comparison between two groups
#' @param wide_data Wide format data frame
#' @param case_group Case group name
#' @param control_group Control group name
#' @param p_thresholds Vector of p-value thresholds to test
#' @param fc_threshold Log2 fold change threshold
#' @return List containing results for each p-value threshold
run_da_comparison <- function(wide_data, case_group, control_group, 
                              p_thresholds = 0.01, fc_threshold = 0.25) {

  # Prepare data (no matching)
  prepared_data <- prepare_comparison_data(wide_data, case_group, control_group)
  
  if (nrow(prepared_data) == 0) {
    warning("No samples available for comparison")
    return(NULL)
  }
  # Create design matrix using clean group names
  design_matrix <- model.matrix(~ 0 + Group_clean + Age + Sex, data = prepared_data)
  colnames(design_matrix) <- gsub("Group_clean", "", colnames(design_matrix))
  
  # Prepare NPX data
  metadata_cols <- c("DAid", "Disease", "Group", "Group_clean", "Sex", "Age","Class") 
  protein_data <- prepared_data %>%
    select(-any_of(metadata_cols)) %>%
    as.matrix() %>%
    t() # Transpose: proteins as rows, samples as columns
  
  # Create contrast string using clean names
  case_group_clean <- gsub(" ", "_", case_group)
  control_group_clean <- gsub(" ", "_", control_group)
  contrast_string <- paste(case_group_clean, "-", control_group_clean)
  
  # Run analysis for each p-value threshold
  results_list <- map(p_thresholds, function(p_thresh) {
    
    limma_result <- run_limma_analysis(
      protein_data = protein_data,
      design_matrix = design_matrix,
      contrast_string = contrast_string,
      p_threshold = p_thresh,
      fc_threshold = fc_threshold
    )
    
    # Add comparison metadata
    limma_result$comparison_info <- list(
      case_group = case_group,
      control_group = control_group,
      n_case = sum(prepared_data$Group == case_group),
      n_control = sum(prepared_data$Group == control_group),
      matched = FALSE
    )
    
    return(limma_result)
  })
  
  names(results_list) <- paste0("p_", p_thresholds)
  
  return(results_list)
}

# ===== INDIVIDUAL DISEASE ANALYSIS =====

#' Run differential abundance analysis for individual autoimmune diseases
#' @param wide_data Wide format data frame
#' @param disease_name Name of the autoimmune disease
#' @param p_thresholds Vector of p-value thresholds to test
#' @param fc_threshold Log2 fold change threshold
#' @return List containing results for each comparison and p-value threshold
run_individual_disease_analysis <- function(wide_data, disease_name, 
                                            p_thresholds = 0.01, 
                                            fc_threshold = 0.25) {
  
  # Check if disease exists in data
  if (!disease_name %in% unique(wide_data$Disease)) {
    warning(paste("Disease", disease_name, "not found in data"))
    return(NULL)
  }
  
  # Run comparisons
  results <- list()
  
  # Disease vs Healthy
  results$vs_healthy <- run_da_comparison(
    wide_data, disease_name, "Healthy", p_thresholds, fc_threshold
  )
  
  # Disease vs Infections
  results$vs_infections <- run_da_comparison(
    wide_data, disease_name, "Infections", p_thresholds, fc_threshold
  )
  
  # Disease vs Other Autoimmune
  other_auto_data <- wide_data %>%
    filter(Disease %in% AUTOIMMUNE_DISEASES)  |> 
    mutate(Disease = ifelse(Disease == disease_name, disease_name, "Other_Auto"))
  
  results$vs_other_auto <- run_da_comparison(
    other_auto_data, disease_name, "Other_Auto", p_thresholds, fc_threshold
  )
  
  return(results)
}

# ===== COMPREHENSIVE ANALYSIS FUNCTION =====
# calls run_da_comparison function and pastes in the names of control groups 
#' Run comprehensive differential abundance analysis
#' @param wide_data Wide format data frame
#' @param p_thresholds Vector of p-value thresholds to test
#' @param fc_threshold Log2 fold change threshold
#' @return List containing all analysis results
run_comprehensive_da_analysis <- function(wide_data, p_thresholds = 0.01, 
                                          fc_threshold = 0.25) {
  
  results <- list()
  
  # Group-level comparisons
  results$group_comparisons <- list(
    autoimmune_vs_healthy = run_da_comparison(
      wide_data, "Autoimmune", "Healthy", p_thresholds, fc_threshold
    ),
    autoimmune_vs_infections = run_da_comparison(
      wide_data, "Autoimmune", "Infection", p_thresholds, fc_threshold
    )
  )
  
  # Individual disease comparisons
  available_diseases <- intersect(AUTOIMMUNE_DISEASES, unique(wide_data$Disease))
  
  results$individual_diseases <- map(available_diseases, function(disease) {
    run_individual_disease_analysis(wide_data, disease, p_thresholds, fc_threshold)
  })
  names(results$individual_diseases) <- available_diseases
  
  # Analysis summary
  results$summary <- sig_prot_stats(results)
  
  return(results)
}

# ===== UTILITY FUNCTIONS =====

# renamed sig prot stats
#' Create analysis summary
#' @param results Complete results list
#' @return Summary data frame
sig_prot_stats <- function(results) {
  
  summary_data <- list()
  
  # Group comparisons summary
  for (comp_name in names(results$group_comparisons)) {
    comp_results <- results$group_comparisons[[comp_name]]
    if (!is.null(comp_results)) {
      for (p_thresh in names(comp_results)) {
        result_df <- comp_results[[p_thresh]]$results
        
        sig_up <- result_df  |> 
          filter(regulation == "Higher") |> 
          nrow()
        
        sig_down <- result_df |> 
          filter(regulation == "Lower") |> 
          nrow()
        
        total_sig <- sig_up + sig_down
        
        summary_data[[paste(comp_name, p_thresh, sep = "_")]] <- data.frame(
          comparison = comp_name,
          p_threshold = gsub("p_", "", p_thresh),
          n_significant = total_sig,
          sig_up = sig_up,
          sig_down = sig_down,
          type = "group"
        )
      }
    }
  }
  
  # Individual disease summary
  for (disease in names(results$individual_diseases)) {
    disease_results <- results$individual_diseases[[disease]]
    if (!is.null(disease_results)) {
      for (comp_type in names(disease_results)) {
        comp_results <- disease_results[[comp_type]]
        if (!is.null(comp_results)) {
          for (p_thresh in names(comp_results)) {
            result_df <- comp_results[[p_thresh]]$results
            
            sig_up <- result_df  |> 
              filter(regulation == "Higher") |> 
              nrow()
            
            sig_down <- result_df |> 
              filter(regulation == "Lower") |> 
              nrow()
            
            total_sig <- sig_up + sig_down
            
            summary_data[[paste(disease, comp_type, p_thresh, sep = "_")]] <- data.frame(
              comparison = paste(disease, comp_type, sep = "_"),
              p_threshold = gsub("p_", "", p_thresh),
              n_significant = total_sig,
              sig_up = sig_up,
              sig_down = sig_down,
              type = "individual"
            )
          }
        }
      }
    }
  }
  
  # Combine all summaries
  do.call(rbind, summary_data) |> 
    arrange(comparison, p_threshold)
}

# ===== Machine Learning =====

# Function to run multiclassification using glmnet lasso regression
# hyperparameter estimation is for inner folds
lasso_multiclass_nested_cv <- function(wide_data, 
                                       n_outer_folds = 5, 
                                       n_inner_folds = 10) {
  
  # Outer CV splits
  set.seed(213)
  outer_folds <- vfold_cv(wide_data, v = n_outer_folds, strata = Disease)
  
  # Loop over outer folds to fit models (do it 5 times and save in outer results)
  outer_results <- map(c(1:n_outer_folds), function(i) {
    
    # outer folds contains list of the 5 splits, so it will go 1-5
    outer_split <- outer_folds$splits[[i]]
    train_data <- training(outer_split)
    test_data <- testing(outer_split)
    
    # Inner CV for hyperparameter tuning - now this is with 10
    inner_folds <- vfold_cv(train_data, v = n_inner_folds, strata = Disease)
    
    # Define recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = train_data) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> # normalize values
      step_nzv(all_numeric()) |> # remove near-zero variance predictors
      step_impute_knn(all_numeric()) |> # imputation
      step_downsample(Disease) # downsamples myositis in the training set
    
    # Model specification
    glmnet_lasso_specs <- 
      multinom_reg(penalty = tune(), mixture = 1) |> 
      set_mode("classification") |> 
      set_engine("glmnet")
    
    # Workflow
    glmnet_wflow <- 
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs)
    
    # Grid search for inner CV
    glmnet_grid <- grid_regular(penalty(), levels = 10) # grid is at 10 bc otherwise it's a bit slow, this is tuning alpha penalty (regularization pararameters)
    
    # tune with the inner folds
    inner_res <- 
      tune_grid(
        glmnet_wflow,
        resamples = inner_folds,
        grid = glmnet_grid,
        metrics = metric_set(roc_auc),
        control = control_grid(save_pred = TRUE)
      )
    
    # Select best hyperparameter from inner CV
    best_params <- select_best(inner_res, metric = "roc_auc")
    final_wflow <- finalize_workflow(glmnet_wflow, best_params)
    
    # Fit final model on the outer training data
    final_fit <- fit(final_wflow, data = train_data)
    
    # Evaluate on outer test data
    test_predictions <- 
      predict(final_fit, test_data, type = "prob") |> # generate predictions in probability 
      bind_cols(predict(final_fit, test_data)) |> # bind to test results
      bind_cols(test_data |> select(DAid, Disease)) |> 
      mutate(Disease = factor(Disease))
    
    roc_multiclass_results <- 
      generate_roc_multiclass(test_predictions) # make ROC plots from test predictions
    
    # Extract important proteins
    protein_importance <- 
      final_fit |> # take final model
      extract_fit_parsnip() |> 
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(desc(abs(estimate))) # arrange in descending order
    
    # return list with predictions, roc curves, and important proteins
    list(
      predictions = test_predictions,
      roc = roc_multiclass_results,
      important_proteins = protein_importance
    )
  })
  
  
  # Aggregate results across outer folds
  all_predictions <- map_dfr(outer_results, "predictions", .id = "fold")
  all_roc <- map_dfr(outer_results, "roc", .id = "fold")
  all_importances <- map_dfr(outer_results, "important_proteins", .id = "fold")
  
  list(
    predictions = all_predictions,
    rocs = all_roc,
    protein_importances = all_importances
  )
}

# Function to generate ROC curve for multiclass classification using the multiROC package
# this one is used within the big function, about getting the predictions in the df and will help with roc curves
generate_roc_multiclass <- function(predictions) {
  print("Input predictions structure:")
  print(str(predictions))
  
  dat <- predictions |> 
    select(DAid, Disease) |>  
    mutate(value = 1) |>  
    tidyr::pivot_wider(names_from = Disease, 
                       values_from = value, 
                       values_fill = 0)
  print("After pivot_wider:")
  print(names(dat))
  
  # Create true labels dataframe with _true suffix
  true_dat <- 
    dat |> 
    rename_with(~paste0(., "_true"), -DAid)
  
  # Create probability predictions dataframe
  dat_prob <- 
    predictions |> 
    rename_with(~stringr::str_remove(., ".pred_")) |> 
    select(-class, -Disease)
  
  # Add _pred_glmnet suffix to probability columns
  prob_data <- 
    dat_prob |>  
    rename_with(~paste0(., "_pred_glmnet"), -DAid)
  
  # Join the dataframes
  final_df <- 
    true_dat |>  
    left_join(prob_data, by = "DAid") |>  
    select(-DAid) |> 
    as.data.frame()
  
  # Generate ROC curves
  roc_res <- multiROC::multi_roc(final_df, force_diag = TRUE)
  
  # Prepare plot data
  plot_roc_df <- multiROC::plot_roc_data(roc_res)
  
  # Format results
  roc_dat <- 
    plot_roc_df |> 
    filter(!Group %in% c("Macro", "Micro")) |> 
    mutate(Performance = paste0(Group, ": ", round(AUC, 4))) |> 
    arrange(desc(AUC))
  
  return(roc_dat)
}

# takes output of generate_roc_multiclass and plots it
plot_roc_multiclass <- function(roc_dat) {
  roc_dat |> 
    arrange(Performance) |> 
    ggplot(aes(x = 1 - Specificity, y = Sensitivity), color = Group) +
    geom_path(size = 1, show.legend = F) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                 color = "grey",
                 linetype = 'dotdash',
                 show.legend = F) +
    geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_x_continuous(breaks = c(0, 1)) +
    facet_wrap(~Group, nrow = 5) 
  
}


# ===== DE & ML integration =====
####                              ####
####          GSEA                ####
####                              ####

# define background
bg_genes <- readRDS('path/to/gene_sets/explore_bg_genes.rds')

do_jk_gsea <- function(DE_genes, adjust_method, cut_off, go_or_kegg) {
  # read in background genes from gene_sets folder
  # if KEGG read in kegg_medicus
  if (go_or_kegg == 'kegg') {
    t2g <- readRDS('path/to/gene_sets/KEGG_MEDICUS_hs.v2023.2.rds')
    # if GO read in go_bp otherwise give error message 
  } else if (go_or_kegg == 'go_bp') {
    t2g <- readRDS('path/to/gene_sets//GO_BP_hs.v2023.2.rds')
  } else {
    stop('Background gene set not found')
  }
  
  res <- enricher(DE_genes,universe = bg_genes, pvalueCutoff = cut_off, pAdjustMethod = adjust_method,TERM2GENE = t2g)
  
  
  return(res)
  
}
# Function that takes differential abundance results and machine learning results along with a p-value threshold and compares the number
# of unique proteins found by each method
# ranks overlapping proteins by p.val or by importance 

# Updated for new DA_results format
# Combined function that performs both comparison and detailed analysis
# note: this saves the importance for each fold of the machine learning, meaning there will be 5 diff importance scores
# for each protein if it was selected in each fold
analyze_da_ml_methods <- function(da_results, ml_results, wide_data, pval_threshold = 0.01) {

  # Initialize results list
  analysis_results <- list()
  
  # Get diseases from the new results structure
  diseases <- names(da_results$individual_diseases)
  
  # Process each disease
  for (disease_name in diseases) {
    
    # Extract DA significant proteins
    da_proteins <- da_results$individual_diseases[[disease_name]]$vs_other_auto$p_0.01$results |> 
      filter(adj.P.Val < pval_threshold)  |> 
      pull(protein)
    
    # Extract ML significant proteins
    ml_proteins <- ml_results$protein_importances  |> 
      filter(class == disease_name)  |> 
      filter(abs(estimate) > 0) |> 
      pull(term) |> 
      unique()
    
    # Calculate overlaps
    overlapping_proteins <- intersect(da_proteins, ml_proteins)
    da_only_proteins <- setdiff(da_proteins, ml_proteins)
    ml_only_proteins <- setdiff(ml_proteins, da_proteins)
    
    # Initialize disease results
    disease_results <- list()
    
    # Store basic comparison results
    disease_results$comparison <- list(
      overlap = list(
        proteins = overlapping_proteins,
        count = length(overlapping_proteins)
      ),
      da_only = list(
        proteins = da_only_proteins,
        count = length(da_only_proteins)
      ),
      ml_only = list(
        proteins = ml_only_proteins,
        count = length(ml_only_proteins)
      ),
      summary = data.frame(
        category = c("Both Methods", "DA Only", "ML Only"),
        count = c(
          length(overlapping_proteins),
          length(da_only_proteins),
          length(ml_only_proteins)
        )
      )
    )
    
    if (length(overlapping_proteins) > 0) {
      
      # DA stats for overlapping proteins
      da_stats <- da_results$individual_diseases[[disease_name]]$vs_other_auto$p_0.01$results |> 
        filter(protein %in% overlapping_proteins) |> 
        arrange(adj.P.Val)
      
      # ML stats for overlapping proteins
      ml_stats <- ml_results$protein_importances |> 
        filter(class == disease_name, term %in% overlapping_proteins)  |> 
        arrange(desc(abs(estimate)))
      
      # Expression patterns
      # changed logFC to > 0.25 < instead of 0 on Sep 2
      expression_patterns <- da_stats |> 
        arrange(desc(abs(logFC))) |> 
        mutate(
          regulation = case_when(
            logFC > 0.25 & adj.P.Val < 0.01 ~ "Higher levels",
            logFC < 0.25 & adj.P.Val < 0.01 ~ "Lower levels",
            TRUE ~ "Not significant"
          )
        )
      
      # ML importance ranking
      importance_ranking <- data.frame(
        protein = ml_stats$term,
        importance = ml_stats$estimate
      )  |> 
        arrange(desc(abs(importance)))
      
      # Expression data and correlation (if wide_data is available)
      correlation_long <- NULL
      if (!missing(wide_data) && exists("wide_data") && nrow(wide_data) > 0) {
        tryCatch({
          expression_data <- wide_data %>%
            select(all_of(c("Disease", overlapping_proteins))) %>%
            filter(Disease == disease_name)
          
          if (nrow(expression_data) > 0) {
            correlation_matrix <- expression_data %>%
              select(-Disease) %>%
              cor(method = "spearman")
            
            correlation_long <- correlation_matrix %>%
              as.data.frame() %>%
              mutate(protein1 = rownames(.)) %>%
              pivot_longer(-protein1, 
                           names_to = "protein2", 
                           values_to = "correlation") %>%
              filter(protein1 < protein2) %>%
              arrange(desc(abs(correlation)))
          }
        }, error = function(e) {
          cat("Warning: Could not compute correlations for", disease_name, ":", e$message, "\n")
        })
      }
      
      # ML vs DA comparison
      ml_da_comparison <- da_stats %>%
        select(protein, logFC, adj.P.Val) %>%
        inner_join(importance_ranking, by = "protein") %>%
        mutate(
          abs_logFC = abs(logFC),
          abs_importance = abs(importance)
        )
      
      # Store detailed analysis results
      disease_results$detailed_analysis <- list(
        expression_patterns = expression_patterns,
        importance_ranking = importance_ranking,
        protein_correlations = correlation_long,
        ml_da_comparison = ml_da_comparison,
        summary_stats = list(
          n_proteins = length(overlapping_proteins),
          n_upregulated = sum(expression_patterns$regulation == "Higher levels"),
          n_downregulated = sum(expression_patterns$regulation == "Lower levels"),
          median_importance = median(abs(importance_ranking$importance)),
          median_correlation = if(!is.null(correlation_long) && nrow(correlation_long) > 0) {
            median(abs(correlation_long$correlation))
          } else {
            NA
          }
        )
      )
      
    } else {
      cat("No overlapping proteins found for", disease_name, "- skipping detailed analysis\n")
    }
    
    # Store disease results
    analysis_results[[disease_name]] <- disease_results

  }
    return(analysis_results)
}



# Create UMAP Plot
create_umap_plot <- function(pc_scores, disease_labels, title) {
  #set.seed(42)
  umap_coords <- uwot::umap(pc_scores, n_neighbors = 15, min_dist = 0.1, n_components = 2, metric = "euclidean")
  umap_plot_data <- data.frame(UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2], Disease = disease_labels)
  disease_colors <- c(
    "Myositis" = "#1F77B4", "Rheumatoid arthritis" = "#D62728", "Sjögrens syndrome" = "#FF7F0E",
    "Systemic lupus erythematosus" = "#2CA02C", "Systemic sclerosis" = "#9467BD"
  )
  ggplot(umap_plot_data, aes(x = UMAP1, y = UMAP2, color = Disease)) +
    geom_point(alpha = 0.8, size = 2.5) +
    scale_color_manual(values = disease_colors) +
    labs(title = title, x = "UMAP1", y = "UMAP2", color = "Disease") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12), axis.text = element_text(size = 10),
      legend.title = element_text(size = 12), legend.text = element_text(size = 10),
      panel.grid.minor = element_blank(), legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}


# # Run the combined analysis
# combined_results <- analyze_da_ml_methods(DA_results_0708, ML_results, wide_data, pval_threshold = 0.01)
# 
# # Extract comparison results for plotting
# comparison_results <- extract_comparison_results(combined_results)
# plot_method_comparison(comparison_results)
# 
# # Extract detailed analysis for specific plots
# detailed_results <- extract_detailed_analysis(combined_results)
# plot_expression_patterns(detailed_results, "Myositis")
# plot_importance_vs_foldchange(detailed_results, "Myositis")
# 
# # Access specific results
# myositis_pathway_results <- combined_results$Myositis$pathway_analysis$go_bp



