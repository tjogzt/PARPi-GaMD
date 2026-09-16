#!/usr/bin/env Rscript
# 08-ml_pipeline.R — ML analysis for PARP1 trapping mechanism
# 
# Purpose:
#   Classify PARP1 inhibitors into Type II (neutral allosteric) vs Type III 
#   (pro-release) using GaMD-derived PMF features. Validate with LOOCV.
#   Blind-predict AZD5305 (saruparib) Type classification.
#
# Usage:
#   Rscript code/08-ml_pipeline.R
#
# Input:
#   results/analysis/*.xvg — PMF curves from GaMD reweighting
#
# Output:
#   results/ml/ — CV metrics, SHAP values, blind prediction

# ==============================================================================
# 0. Setup
# ==============================================================================
library(xgboost)
library(lightgbm)
library(shapviz)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

dir.create("results/ml", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Read PMF data
# ==============================================================================
read_xvg <- function(path) {
  lines <- readLines(path)
  data_lines <- lines[!grepl("^[@#]", lines)]
  data_lines <- data_lines[data_lines != ""]
  read.table(text = data_lines, header = FALSE)
}

# ==============================================================================
# 2. Extract features from a single PMF curve
# ==============================================================================
extract_pmf_features <- function(pmf_df) {
  colnames(pmf_df) <- c("RC", "PMF")
  pmf_df$PMF <- pmf_df$PMF - min(pmf_df$PMF)
  
  rc  <- pmf_df$RC
  pmf <- pmf_df$PMF
  
  min_idx     <- which.min(pmf)
  min_rc      <- rc[min_idx]
  well_depth  <- max(pmf) - min(pmf)
  
  # Basin width (FWHM)
  half_max    <- max(pmf) / 2
  above_half  <- which(pmf <= half_max)
  basin_width <- if (length(above_half) > 1) {
    rc[tail(above_half, 1)] - rc[above_half[1]]
  } else {
    NA_real_
  }
  
  # Barrier height
  left_max  <- if (min_idx > 1) max(pmf[1:min_idx]) else 0
  right_max <- if (min_idx < length(pmf)) max(pmf[min_idx:length(pmf)]) else 0
  barrier   <- max(left_max, right_max) - pmf[min_idx]
  
  # Number of metastable states (plain local-minimum count).
  # Note: 12-descriptive_analysis.R / 28-pca_features.R wrap this in a
  # "pmf <= pmf + 0.5" filter that is tautological and never fires — the same quantity
  # results. See data_manifest.md (n_states semantics).
  d1 <- diff(pmf)
  minima_idx <- which(diff(sign(d1)) == 2) + 1
  n_states   <- length(minima_idx)
  
  # Area under PMF curve
  auc <- sum(diff(rc) * (pmf[-1] + pmf[-length(pmf)]) / 2)
  
  # Skewness
  pmf_left  <- rc < min_rc
  pmf_right <- rc > min_rc
  skew <- if (sum(pmf_left) > 0 & sum(pmf_right) > 0) {
    mean(pmf[pmf_left]) - mean(pmf[pmf_right])
  } else {
    NA_real_
  }
  
  # Asymmetry ratio
  left_width  <- if (min_idx > 1) min_rc - rc[1] else NA_real_
  right_width <- if (min_idx < length(rc)) rc[length(rc)] - min_rc else NA_real_
  asym_ratio  <- if (!is.na(left_width) & !is.na(right_width) & right_width > 0) {
    left_width / right_width
  } else {
    NA_real_
  }
  
  data.frame(
    min_position   = min_rc,
    well_depth     = well_depth,
    basin_width    = basin_width,
    barrier        = barrier,
    n_states       = n_states,
    auc            = auc,
    skew           = skew,
    asym_ratio     = asym_ratio,
    row.names      = NULL
  )
}

# ==============================================================================
# 3. Parse PMF filenames (NO REGEX — pure string operations)
# ==============================================================================
parse_filename <- function(path) {
  base <- basename(path)
  
  # sys1: sys1_{LIGAND}_pmf_c{CV}.xvg
  if (grepl("^sys1_", base)) {
    stripped <- tools::file_path_sans_ext(base)
    parts <- strsplit(stripped, "_", fixed = TRUE)[[1]]
    pmf_pos <- which(parts == "pmf")
    cv_pos  <- which(grepl("^c[0-9]+$", parts))
    if (length(pmf_pos) > 0 && length(cv_pos) > 0) {
      cv_num <- as.integer(gsub("c", "", parts[cv_pos], fixed = TRUE))
      ligand <- paste(parts[2:(pmf_pos - 1)], collapse = "_")
      return(list(system = "sys1", ligand = ligand,
                  cv = paste0("CV", cv_num),
                  cumulant = NA_integer_, path = path))
    }
  }
  
  # sys2: pmf-c{C}_sys2_{LIGAND}_CV{N}_cv.dat.xvg
  if (grepl("^pmf-c", base)) {
    stripped <- tools::file_path_sans_ext(base)
    # pmf-c{CUMULANT}-sys2_{LIGAND}_CV{N}_cv.dat
    after_prefix <- gsub("^pmf-c", "", stripped)
    cum_str <- strsplit(after_prefix, "-", fixed = TRUE)[[1]][1]
    cum <- as.integer(cum_str)
    
    # Extract ligand: between "sys2_" and "_CV"
    lig_start <- regexpr("sys2_", stripped, fixed = TRUE)
    cv_start  <- regexpr("_CV", stripped, fixed = TRUE)
    if (lig_start > 0 && cv_start > 0) {
      lig_start <- lig_start + 5  # length of "sys2_"
      ligand <- substr(stripped, lig_start, cv_start - 1)
    } else {
      return(NULL)
    }
    
    # Extract CV number: _CV{N}_cv.dat
    cv_part <- gsub(".*_CV", "", stripped)
    cv_num  <- as.integer(strsplit(cv_part, "_", fixed = TRUE)[[1]][1])
    
    if (!is.na(cum) && !is.na(cv_num)) {
      return(list(system = "sys2", ligand = ligand,
                  cv = paste0("CV", cv_num),
                  cumulant = cum, path = path))
    }
  }
  
  NULL
}

# ==============================================================================
# 4. Build feature matrix
# ==============================================================================
known_labels <- list(
  "sys2" = list(
    "APO"         = "APO",
    "talazoparib" = "Type_II"
  ),
  "sys1" = list(
    "APO"         = "APO",
    "AZD5305"     = "unknown",
    "olaparib"    = "Type_II",
    "talazoparib" = "Type_II",
    "veliparib"   = "Type_III",
    "niraparib"   = "Type_II",
    "rucaparib"   = "Type_II"
  )
)

# Scan PMF files
analysis_dir <- "results/analysis"
all_files    <- list.files(analysis_dir, pattern = "pmf.*[.]xvg$", full.names = TRUE)

# Parse all files
pmf_meta_list <- lapply(all_files, function(f) {
  info <- parse_filename(f)
  if (is.null(info)) return(NULL)
  as.data.frame(info, stringsAsFactors = FALSE)
})
pmf_meta <- do.call(rbind, pmf_meta_list)

if (is.null(pmf_meta) || nrow(pmf_meta) == 0) {
  stop("No PMF files found in results/analysis/")
}

cat(sprintf("Found %d PMF files from %d unique systems\n", 
            nrow(pmf_meta), length(unique(pmf_meta$system))))

# Extract features row by row, collect in list then rbind with consistent names
all_feature_names <- c("min_position", "well_depth", "basin_width", 
                       "barrier", "n_states", "auc", "skew", "asym_ratio")

feature_rows <- do.call(rbind, lapply(seq_len(nrow(pmf_meta)), function(i) {
  row      <- pmf_meta[i, ]
  pmf_data <- tryCatch(read_xvg(row$path), error = function(e) NULL)
  if (is.null(pmf_data)) return(NULL)
  
  feats <- extract_pmf_features(pmf_data)
  
  # Build consistent row: system, ligand, cv, cumulant + all feature columns
  out <- data.frame(
    system   = row$system,
    ligand   = row$ligand,
    cv       = row$cv,
    cumulant = if (is.na(row$cumulant)) 0L else row$cumulant,
    stringsAsFactors = FALSE
  )
  # Add feature columns, fill missing with NA
  for (fn in names(feats)) {
    out[[fn]] <- feats[[fn]]
  }
  out
}))

cat(sprintf("Extracted features from %d PMF files\n", nrow(feature_rows)))

# Select best cumulant per (system, ligand, CV)
feature_rows <- feature_rows %>%
  group_by(system, ligand, cv) %>%
  slice_max(cumulant, n = 1) %>%
  ungroup()

# Pivot to wide: each CV's columns get CV prefix
feat_cols <- setdiff(names(feature_rows), c("system", "ligand", "cv", "cumulant"))

features_wide <- feature_rows %>%
  select(system, ligand, cv, all_of(feat_cols)) %>%
  pivot_longer(
    cols      = all_of(feat_cols),
    names_to  = "metric",
    values_to = "value"
  ) %>%
  mutate(col_name = paste0(cv, "_", metric)) %>%
  select(-cv, -metric) %>%
  pivot_wider(
    id_cols     = c(system, ligand),
    names_from  = col_name,
    values_from = value,
    values_fn   = first
  )

cat(sprintf("Feature matrix: %d rows x %d cols\n", 
            nrow(features_wide), ncol(features_wide)))

# --- Assign labels ---
features_wide$label_type <- NA_character_

for (sys in names(known_labels)) {
  for (lig in names(known_labels[[sys]])) {
    lbl <- known_labels[[sys]][[lig]]
    if (lbl != "APO") {
      features_wide$label_type[features_wide$system == sys & 
                               features_wide$ligand == lig] <- lbl
    }
  }
}

# Filter to single system for clean analysis (S1 and S2 use different CV definitions)
# Change to "sys2" or NULL (all) as needed
system_filter <- "sys1"
if (!is.null(system_filter)) {
  features_wide <- features_wide %>% filter(system == system_filter)
}

features_known   <- features_wide %>% filter(!is.na(label_type) & label_type != "unknown")
features_unknown <- features_wide %>% filter(ligand == "AZD5305" | label_type == "unknown")

cat(sprintf("Known samples:  %d\n", nrow(features_known)))
cat(sprintf("Unknown samples: %d\n", nrow(features_unknown)))

# ==============================================================================
# 5. Prepare ML data
# ==============================================================================
prepare_ml_data <- function(df) {
  feature_cols <- grep("^CV[0-9]+_", names(df), value = TRUE, perl = TRUE)
  if (length(feature_cols) == 0) {
    feature_cols <- names(df)[sapply(df, is.numeric)]
  }
  
  x <- as.matrix(df[, feature_cols, drop = FALSE])
  
  for (j in seq_len(ncol(x))) {
    col_na <- is.na(x[, j])
    if (any(col_na)) {
      x[col_na, j] <- median(x[, j], na.rm = TRUE)
    }
  }
  
  y <- ifelse(df$label_type == "Type_II", 1, 0)
  
  list(
    x             = x,
    y             = y,
    feature_names = feature_cols,
    ligand_names  = df$ligand
  )
}

# ==============================================================================
# 6. XGBoost LOOCV
# ==============================================================================
run_xgboost_loocv <- function(ml_data, n_folds = NULL, seed = 49) {
  n <- length(ml_data$y)
  if (is.null(n_folds)) {
    n_folds <- n
    cv_mode <- "LOOCV"
  } else {
    cv_mode <- sprintf("%d-fold", n_folds)
  }
  
  set.seed(seed)
  
  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    max_depth        = min(3, n - 1),
    eta              = 0.1,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    nthread          = 2,
    seed             = 49
  )
  
  nrounds <- 50
  
  if (n_folds == n) {
    fold_ids <- seq_len(n)
  } else {
    fold_ids <- sample(rep(seq_len(n_folds), length.out = n))
  }
  
  cv_results <- data.frame(
    fold      = integer(),
    ligand    = character(),
    actual    = integer(),
    predicted = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (fold in unique(fold_ids)) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    dtrain <- xgb.DMatrix(ml_data$x[train_idx, , drop = FALSE], 
                          label = ml_data$y[train_idx])
    dtest  <- xgb.DMatrix(ml_data$x[test_idx, , drop = FALSE], 
                          label = ml_data$y[test_idx])
    
    model <- xgb.train(
      params   = params,
      data     = dtrain,
      nrounds  = nrounds,
      verbose  = 0
    )
    
    pred <- predict(model, dtest)
    
    cv_results <- rbind(cv_results, data.frame(
      fold      = fold,
      ligand    = ml_data$ligand_names[test_idx],
      actual    = ml_data$y[test_idx],
      predicted = pred,
      stringsAsFactors = FALSE
    ))
  }
  
  cv_results$pred_class <- ifelse(cv_results$predicted >= 0.5, 1, 0)
  accuracy <- mean(cv_results$pred_class == cv_results$actual)
  
  unique_classes <- unique(cv_results$actual)
  auc_val <- if (length(unique_classes) >= 2) {
    pos <- cv_results$predicted[cv_results$actual == 1]
    neg <- cv_results$predicted[cv_results$actual == 0]
    if (length(pos) > 0 & length(neg) > 0) {
      mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }
  
  list(
    cv_results = cv_results,
    accuracy   = accuracy,
    auc        = auc_val,
    cv_mode    = cv_mode,
    n_folds    = n_folds
  )
}

# ==============================================================================
# 7. LightGBM LOOCV
# ==============================================================================
run_lgb_loocv <- function(ml_data, n_folds = NULL, seed = 49) {
  n <- length(ml_data$y)
  if (is.null(n_folds)) {
    n_folds <- n
    cv_mode <- "LOOCV"
  } else {
    cv_mode <- sprintf("%d-fold", n_folds)
  }
  
  set.seed(seed)
  
  params <- list(
    objective      = "binary",
    metric         = "binary_logloss",
    max_depth      = min(3, n - 1),
    learning_rate  = 0.1,
    num_leaves     = max(2, min(7, 2^(min(3, n - 1)) - 1)),
    min_data_in_leaf = 1,
    verbose        = -1,
    num_threads    = 2,
    seed           = 49
  )
  
  nrounds <- 50
  
  if (n_folds == n) {
    fold_ids <- seq_len(n)
  } else {
    fold_ids <- sample(rep(seq_len(n_folds), length.out = n))
  }
  
  cv_results <- data.frame(
    fold      = integer(),
    ligand    = character(),
    actual    = integer(),
    predicted = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (fold in unique(fold_ids)) {
    test_idx  <- which(fold_ids == fold)
    train_idx <- setdiff(seq_len(n), test_idx)
    
    dtrain <- lgb.Dataset(
      data  = ml_data$x[train_idx, , drop = FALSE],
      label = ml_data$y[train_idx]
    )
    
    model <- lgb.train(
      params  = params,
      data    = dtrain,
      nrounds = nrounds,
      verbose = -1
    )
    
    pred <- predict(model, ml_data$x[test_idx, , drop = FALSE])
    
    cv_results <- rbind(cv_results, data.frame(
      fold      = fold,
      ligand    = ml_data$ligand_names[test_idx],
      actual    = ml_data$y[test_idx],
      predicted = pred,
      stringsAsFactors = FALSE
    ))
  }
  
  cv_results$pred_class <- ifelse(cv_results$predicted >= 0.5, 1, 0)
  accuracy <- mean(cv_results$pred_class == cv_results$actual)
  
  unique_classes <- unique(cv_results$actual)
  auc_val <- if (length(unique_classes) >= 2) {
    pos <- cv_results$predicted[cv_results$actual == 1]
    neg <- cv_results$predicted[cv_results$actual == 0]
    if (length(pos) > 0 & length(neg) > 0) {
      mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
    } else {
      NA_real_
    }
  } else {
    NA_real_
  }
  
  list(
    cv_results = cv_results,
    accuracy   = accuracy,
    auc        = auc_val,
    cv_mode    = cv_mode,
    n_folds    = n_folds
  )
}

# ==============================================================================
# 8. Full model for SHAP + blind prediction
# ==============================================================================
train_full_xgboost <- function(ml_data, seed = 49) {
  set.seed(seed)
  n <- length(ml_data$y)
  
  params <- list(
    objective        = "binary:logistic",
    eval_metric      = "logloss",
    max_depth        = min(3, n),
    eta              = 0.1,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    nthread          = 2,
    seed             = 49
  )
  
  dtrain <- xgb.DMatrix(ml_data$x, label = ml_data$y)
  xgb.train(params = params, data = dtrain, nrounds = 50, verbose = 0)
}

train_full_lightgbm <- function(ml_data, seed = 49) {
  set.seed(seed)
  n <- length(ml_data$y)
  
  params <- list(
    objective      = "binary",
    metric         = "binary_logloss",
    max_depth      = min(3, n),
    learning_rate  = 0.1,
    num_leaves     = min(7, 2^(min(3, n)) - 1),
    verbose        = -1,
    num_threads    = 2,
    seed           = 49
  )
  
  dtrain <- lgb.Dataset(data = ml_data$x, label = ml_data$y)
  lgb.train(params = params, data = dtrain, nrounds = 50, verbose = -1)
}

# ==============================================================================
# 9. SHAP analysis
# ==============================================================================
compute_shap_xgb <- function(model, ml_data) {
  shap_contrib <- predict(model, ml_data$x, predcontrib = TRUE)
  shap_values  <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
  bias         <- shap_contrib[1, ncol(shap_contrib)]
  
  importance <- data.frame(
    feature       = colnames(shap_values),
    mean_abs_shap = colMeans(abs(shap_values)),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(mean_abs_shap))
  
  shp <- shapviz(shap_values, X = ml_data$x, X_names = colnames(shap_values))
  
  list(shap_values = shap_values, bias = bias, importance = importance,
       shapviz_obj = shp, feature_names = colnames(shap_values))
}

# ==============================================================================
# 10. Blind prediction
# ==============================================================================
blind_predict <- function(model, features_unknown, ml_data, model_type = "xgboost") {
  feature_cols <- ml_data$feature_names
  x_unknown <- as.matrix(features_unknown[, feature_cols, drop = FALSE])
  
  for (j in seq_len(ncol(x_unknown))) {
    col_na <- is.na(x_unknown[, j])
    if (any(col_na)) {
      x_unknown[col_na, j] <- median(x_unknown[, j], na.rm = TRUE)
    }
  }
  
  if (model_type == "xgboost") {
    pred_prob <- predict(model, xgb.DMatrix(x_unknown))
  } else {
    pred_prob <- predict(model, x_unknown)
  }
  
  pred_class <- ifelse(pred_prob >= 0.5, "Type_II", "Type_III")
  
  data.frame(
    ligand       = features_unknown$ligand,
    pred_prob_II = pred_prob,
    pred_class   = pred_class,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 11. MAIN EXECUTION
# ==============================================================================
cat("\n=========== PARP1 Trapping ML Analysis ===========\n")

if (nrow(features_known) < 2) {
  cat("\n*** SMOKE TEST MODE: < 2 known samples — code validation only ***\n")
  cat("Classification results are not meaningful with n < 2.\n\n")
  
  cat("Feature summary:\n")
  print(summary(features_wide))
  
  # Save feature matrix even in smoke test
  write.csv(features_wide, "results/ml/feature_matrix.csv", row.names = FALSE)
  
  cat("\nResults saved to results/ml/\n")
  cat("=========== DONE (smoke test) ===========\n")
  quit(save = "no", status = 0)
}

ml_data <- prepare_ml_data(features_known)

cat(sprintf("Feature matrix: %d samples x %d features\n", 
            nrow(ml_data$x), ncol(ml_data$x)))
cat(sprintf("Class balance: %d Type_II, %d Type_III\n",
            sum(ml_data$y == 1), sum(ml_data$y == 0)))
cat("Features:\n")
cat(paste("  ", ml_data$feature_names, collapse = "\n"), "\n")

# XGBoost LOOCV
cat("\n--- XGBoost LOOCV ---\n")
xgb_cv_results <- run_xgboost_loocv(ml_data, n_folds = NULL)
cat(sprintf("CV mode:  %s\n", xgb_cv_results$cv_mode))
cat(sprintf("Accuracy: %.3f\n", xgb_cv_results$accuracy))
if (!is.na(xgb_cv_results$auc)) {
  cat(sprintf("AUC:      %.3f\n", xgb_cv_results$auc))
}
print(xgb_cv_results$cv_results, row.names = FALSE)

# LightGBM LOOCV
cat("\n--- LightGBM LOOCV ---\n")
lgb_cv_results <- run_lgb_loocv(ml_data, n_folds = NULL)
cat(sprintf("CV mode:  %s\n", lgb_cv_results$cv_mode))
cat(sprintf("Accuracy: %.3f\n", lgb_cv_results$accuracy))
if (!is.na(lgb_cv_results$auc)) {
  cat(sprintf("AUC:      %.3f\n", lgb_cv_results$auc))
}

# Full model
cat("\n--- Full Model ---\n")
xgb_full <- train_full_xgboost(ml_data)
lgb_full <- train_full_lightgbm(ml_data)

# SHAP
cat("\n--- SHAP Feature Importance (XGBoost) ---\n")
xgb_shap <- compute_shap_xgb(xgb_full, ml_data)
print(xgb_shap$importance, row.names = FALSE)

# Blind prediction for AZD5305
if (nrow(features_unknown) > 0) {
  cat("\n--- Blind Prediction: AZD5305 ---\n")
  
  xgb_pred <- blind_predict(xgb_full, features_unknown, ml_data, "xgboost")
  lgb_pred <- blind_predict(lgb_full, features_unknown, ml_data, "lightgbm")
  
  cat("XGBoost  - P(Type_II):", sprintf("%.3f", xgb_pred$pred_prob_II), 
      "->", xgb_pred$pred_class, "\n")
  cat("LightGBM - P(Type_II):", sprintf("%.3f", lgb_pred$pred_prob_II),
      "->", lgb_pred$pred_class, "\n")
}

# ==============================================================================
# 12. Save results
# ==============================================================================
write.csv(xgb_cv_results$cv_results, "results/ml/xgboost_loocv.csv", row.names = FALSE)
write.csv(lgb_cv_results$cv_results, "results/ml/lightgbm_loocv.csv", row.names = FALSE)
write.csv(xgb_shap$importance, "results/ml/xgboost_shap_importance.csv", row.names = FALSE)
write.csv(features_wide, "results/ml/feature_matrix.csv", row.names = FALSE)
saveRDS(xgb_full, "results/ml/xgboost_full_model.rds")
saveRDS(lgb_full, "results/ml/lightgbm_full_model.rds")

if (nrow(features_unknown) > 0) {
  blind_results <- rbind(
    cbind(model = "XGBoost", xgb_pred),
    cbind(model = "LightGBM", lgb_pred)
  )
  write.csv(blind_results, "results/ml/blind_prediction_AZD5305.csv", row.names = FALSE)
}

cat("\nResults saved to results/ml/\n")
cat("=========== DONE ===========\n")
