## =============================================================================
## Comorbidity Network Simulation from Multi-Omics Data
## Compatible with R >= 3.6
##
## Simulates three disease conditions (e.g. Epilepsy, Type 2 Diabetes, CVD),
## integrates genomics / transcriptomics / proteomics layers, builds a
## comorbidity network, trains a multi-label risk predictor and visualises
## the most influential features.
##
## Packages required (CRAN):
##   igraph, ggplot2, ggraph, tidygraph, reshape2, RColorBrewer,
##   randomForest, caret, corrplot, gridExtra, scales, viridis
##
## Install once with:
##   install.packages(c("igraph","ggplot2","ggraph","tidygraph","reshape2",
##                      "RColorBrewer","randomForest","caret","corrplot",
##                      "gridExtra","scales","viridis"))
## =============================================================================

# ------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------

set.seed(42)

suppressPackageStartupMessages({
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(reshape2)
  library(RColorBrewer)
  library(randomForest)
  library(caret)
  library(corrplot)
  library(gridExtra)
  library(scales)
  library(viridis)
})

# Palette used throughout
DISEASE_COLOURS <- c(
  "Epilepsy"       = "#E63946",
  "Type2Diabetes"  = "#2A9D8F",
  "CVD"            = "#F4A261"
)

OMICS_COLOURS <- c(
  "Genomics"       = "#457B9D",
  "Transcriptomics"= "#A8DADC",
  "Proteomics"     = "#F1FAEE"
)


# ------------------------------------------------------------------
# 1. Simulate multi-omics data
# ------------------------------------------------------------------

#' simulate_multiomics
#'
#' Generates a list of three omics matrices (samples × features) with
#' disease-specific signal embedded in feature subsets.
#'
#' @param n_samples   Integer. Total number of subjects.
#' @param n_genomic   Integer. Number of SNP / genetic-score features.
#' @param n_transcript Integer. Number of gene-expression features.
#' @param n_proteomic Integer. Number of protein-abundance features.
#' @param disease_prev Named numeric. Prevalence of each disease (0–1).
#'
#' @return A named list:
#'   \item{genomics}      data.frame – genomic features + disease labels
#'   \item{transcriptomics} data.frame – expression features + disease labels
#'   \item{proteomics}    data.frame – protein features + disease labels
#'   \item{labels}        data.frame – binary disease indicator per subject
#'   \item{subject_ids}   character vector

simulate_multiomics <- function(
    n_samples    = 500,
    n_genomic    = 50,
    n_transcript = 80,
    n_proteomic  = 40,
    disease_prev = c(Epilepsy = 0.15, Type2Diabetes = 0.25, CVD = 0.20)
) {

  diseases    <- names(disease_prev)
  subject_ids <- paste0("S", seq_len(n_samples))

  # ---- Disease labels (allow comorbidity / multi-label) ----
  labels <- as.data.frame(
    lapply(disease_prev, function(p) rbinom(n_samples, 1, p))
  )
  rownames(labels) <- subject_ids

  # ---- Helper: add disease-specific mean shift to a feature block ----
  add_signal <- function(mat, labels, effect_size = 0.8, frac_causal = 0.3) {
    n_feat   <- ncol(mat)
    n_causal <- max(1L, floor(n_feat * frac_causal))
    for (d in seq_len(ncol(labels))) {
      causal_idx <- sample(n_feat, n_causal)
      mat[labels[[d]] == 1, causal_idx] <-
        mat[labels[[d]] == 1, causal_idx] + effect_size
    }
    mat
  }

  # ---- Genomics: minor-allele-count proxy (0/1/2 scale, then normalise) ----
  geno_mat  <- matrix(
    sample(0:2, n_samples * n_genomic, replace = TRUE,
           prob = c(0.6, 0.3, 0.1)),
    nrow = n_samples, ncol = n_genomic
  )
  geno_mat  <- add_signal(geno_mat, labels, effect_size = 0.5, frac_causal = 0.25)
  geno_df   <- as.data.frame(geno_mat)
  colnames(geno_df) <- paste0("SNP_", seq_len(n_genomic))

  # ---- Transcriptomics: log-normalised counts proxy ----
  trans_mat <- matrix(
    rnorm(n_samples * n_transcript, mean = 5, sd = 1.5),
    nrow = n_samples, ncol = n_transcript
  )
  trans_mat <- add_signal(trans_mat, labels, effect_size = 1.0, frac_causal = 0.30)
  trans_df  <- as.data.frame(trans_mat)
  colnames(trans_df) <- paste0("GENE_", seq_len(n_transcript))

  # ---- Proteomics: log-abundance proxy ----
  prot_mat  <- matrix(
    rnorm(n_samples * n_proteomic, mean = 3, sd = 0.8),
    nrow = n_samples, ncol = n_proteomic
  )
  prot_mat  <- add_signal(prot_mat, labels, effect_size = 0.9, frac_causal = 0.35)
  prot_df   <- as.data.frame(prot_mat)
  colnames(prot_df) <- paste0("PROT_", seq_len(n_proteomic))

  # ---- Attach subject IDs ----
  for (df_name in c("geno_df", "trans_df", "prot_df")) {
    assign(df_name, `rownames<-`(get(df_name), subject_ids))
  }

  list(
    genomics         = cbind(geno_df,  labels),
    transcriptomics  = cbind(trans_df, labels),
    proteomics       = cbind(prot_df,  labels),
    labels           = labels,
    subject_ids      = subject_ids
  )
}


# ------------------------------------------------------------------
# 2. Build comorbidity network
# ------------------------------------------------------------------

#' build_comorbidity_network
#'
#' Constructs a weighted, undirected comorbidity network where nodes are
#' diseases and edges are weighted by three complementary measures:
#'   (a) phenotypic co-occurrence (phi coefficient)
#'   (b) shared genomic signal (feature correlation)
#'   (c) shared transcriptomic signal
#'
#' @param omics_list  List returned by \code{simulate_multiomics}.
#' @param phi_thresh  Numeric. Minimum phi coefficient to retain an edge.
#'
#' @return An igraph object with edge attributes \code{phi}, \code{genomic_sim},
#'         \code{transcript_sim}, and composite \code{weight}.

build_comorbidity_network <- function(omics_list, phi_thresh = 0.05) {

  labels   <- omics_list$labels
  diseases <- colnames(labels)
  n_dis    <- length(diseases)

  # (a) Phenotypic co-occurrence: phi coefficient for each disease pair
  phi_matrix <- matrix(0, n_dis, n_dis, dimnames = list(diseases, diseases))
  for (i in seq_len(n_dis - 1)) {
    for (j in (i + 1):n_dis) {
      tbl <- table(labels[[i]], labels[[j]])
      if (all(dim(tbl) == c(2, 2))) {
        # Coerce to numeric immediately: table() returns integers, and
        # products like n11*n00 with n~600 exceed .Machine$integer.max
        # (~2.1e9), producing NA and breaking the denom > 0 guard.
        n11 <- as.numeric(tbl[2, 2]); n10 <- as.numeric(tbl[2, 1])
        n01 <- as.numeric(tbl[1, 2]); n00 <- as.numeric(tbl[1, 1])
        denom <- sqrt((n11 + n10) * (n01 + n00) *
                        (n11 + n01) * (n10 + n00))
        phi   <- if (is.finite(denom) && denom > 0) (n11 * n00 - n10 * n01) / denom else 0
      } else { phi <- 0 }
      phi_matrix[i, j] <- phi_matrix[j, i] <- abs(phi)
    }
  }

  # (b) Shared omics signal: per-disease mean feature vector, then cosine sim
  cosine_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)) + 1e-9)

  get_mean_profiles <- function(df, labels) {
    feat_cols <- setdiff(colnames(df), colnames(labels))
    lapply(colnames(labels), function(d) {
      colMeans(df[df[[d]] == 1, feat_cols, drop = FALSE])
    })
  }

  geno_profiles  <- get_mean_profiles(omics_list$genomics,       labels)
  trans_profiles <- get_mean_profiles(omics_list$transcriptomics, labels)

  genomic_sim_mat  <- matrix(0, n_dis, n_dis, dimnames = list(diseases, diseases))
  transcript_sim_mat <- genomic_sim_mat

  for (i in seq_len(n_dis - 1)) {
    for (j in (i + 1):n_dis) {
      gs <- cosine_sim(geno_profiles[[i]],  geno_profiles[[j]])
      ts <- cosine_sim(trans_profiles[[i]], trans_profiles[[j]])
      genomic_sim_mat[i, j]    <- genomic_sim_mat[j, i]    <- gs
      transcript_sim_mat[i, j] <- transcript_sim_mat[j, i] <- ts
    }
  }

  # Composite weight (equally weighted average of the three measures)
  weight_matrix <- (phi_matrix + genomic_sim_mat + transcript_sim_mat) / 3

  # Build igraph from upper triangle
  edges <- data.frame()
  for (i in seq_len(n_dis - 1)) {
    for (j in (i + 1):n_dis) {
      if (phi_matrix[i, j] >= phi_thresh) {
        edges <- rbind(edges, data.frame(
          from            = diseases[i],
          to              = diseases[j],
          phi             = round(phi_matrix[i, j], 4),
          genomic_sim     = round(genomic_sim_mat[i, j], 4),
          transcript_sim  = round(transcript_sim_mat[i, j], 4),
          weight          = round(weight_matrix[i, j], 4),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Node attributes: prevalence, mean comorbidity degree
  node_df <- data.frame(
    name       = diseases,
    prevalence = colMeans(labels),
    stringsAsFactors = FALSE
  )

  g <- graph_from_data_frame(edges, directed = FALSE, vertices = node_df)
  g
}


# ------------------------------------------------------------------
# 3. Risk prediction + feature importance
# ------------------------------------------------------------------

#' train_risk_model
#'
#' Trains a Random Forest for each disease using the integrated omics feature
#' matrix and returns out-of-bag importance scores.
#'
#' @param omics_list   List from \code{simulate_multiomics}.
#' @param ntree        Number of trees.
#' @param importance   Logical. Compute variable importance?
#'
#' @return Named list (one entry per disease):
#'   \item{model}      randomForest object
#'   \item{importance} data.frame with MeanDecreaseGini and omics layer tag

train_risk_models <- function(omics_list, ntree = 200, importance = TRUE) {

  labels <- omics_list$labels

  # Integrated feature matrix (genomics + transcriptomics + proteomics)
  feat_cols_g <- grep("^SNP_",  colnames(omics_list$genomics),         value = TRUE)
  feat_cols_t <- grep("^GENE_", colnames(omics_list$transcriptomics),  value = TRUE)
  feat_cols_p <- grep("^PROT_", colnames(omics_list$proteomics),       value = TRUE)

  X <- cbind(
    omics_list$genomics[,         feat_cols_g],
    omics_list$transcriptomics[,  feat_cols_t],
    omics_list$proteomics[,       feat_cols_p]
  )

  models <- lapply(colnames(labels), function(disease) {
    y   <- factor(labels[[disease]], levels = c(0, 1),
                  labels = c("Control", "Case"))
    rf  <- randomForest(
      x          = X,
      y          = y,
      ntree      = ntree,
      importance = importance,
      do.trace   = FALSE
    )

    imp_df <- data.frame(
      feature    = rownames(importance(rf)),
      MeanDecGini = importance(rf)[, "MeanDecreaseGini"],
      layer      = dplyr::case_when(
        grepl("^SNP_",  rownames(importance(rf))) ~ "Genomics",
        grepl("^GENE_", rownames(importance(rf))) ~ "Transcriptomics",
        TRUE                                      ~ "Proteomics"
      ),
      stringsAsFactors = FALSE
    )
    imp_df <- imp_df[order(imp_df$MeanDecGini, decreasing = TRUE), ]

    list(model = rf, importance = imp_df, disease = disease)
  })
  names(models) <- colnames(labels)
  models
}


# ------------------------------------------------------------------
# 4. Visualisation functions
# ------------------------------------------------------------------

## 4a. Network plot ---------------------------------------------------

#' plot_comorbidity_network
#'
#' Draws the comorbidity network with node size ~ prevalence and
#' edge width / opacity ~ composite weight.
#'
#' @param g      igraph object from \code{build_comorbidity_network}.
#' @param title  Plot title string.
#' @return A ggplot object.

plot_comorbidity_network <- function(g, title = "Multi-Omics Comorbidity Network") {

  tg <- as_tbl_graph(g)

  ggraph(tg, layout = "circle") +
    geom_edge_link(
      aes(width = weight, alpha = weight, colour = phi),
      show.legend = TRUE
    ) +
    scale_edge_width(range  = c(1, 6),   name = "Composite\nWeight") +
    scale_edge_alpha(range  = c(0.4, 1), guide = "none") +
    scale_edge_colour_gradient(
      low  = "#A8DADC", high = "#E63946",
      name = "Phenotypic\nCo-occurrence\n(phi)"
    ) +
    geom_node_point(
      aes(size = prevalence, fill = name),
      shape = 21, stroke = 1.5, colour = "white"
    ) +
    scale_size(range = c(8, 20), name = "Prevalence") +
    scale_fill_manual(values = DISEASE_COLOURS, name = "Disease") +
    geom_node_label(
      aes(label = name),
      repel       = FALSE,
      fontface    = "bold",
      size        = 4.5,
      label.size  = 0,
      fill        = NA,
      colour      = "white",
      nudge_y     = -0.15
    ) +
    labs(title = title,
         subtitle = "Node size ∝ prevalence  |  Edge width & colour ∝ comorbidity strength") +
    theme_graph(base_family = "sans") +
    theme(
      plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, colour = "grey50"),
      legend.position = "right"
    )
}


## 4b. Feature importance (top-N bar chart) ---------------------------

#' plot_feature_importance
#'
#' Bar chart of top-N most important features for one disease,
#' coloured by omics layer.
#'
#' @param model_result  Single-disease list from \code{train_risk_models}.
#' @param top_n         Number of features to display.
#' @return A ggplot object.

plot_feature_importance <- function(model_result, top_n = 20) {

  imp  <- model_result$importance
  top  <- head(imp, top_n)
  top$feature <- factor(top$feature, levels = rev(top$feature))

  layer_cols <- c(
    "Genomics"        = "#457B9D",
    "Transcriptomics" = "#2A9D8F",
    "Proteomics"      = "#F4A261"
  )

  ggplot(top, aes(x = feature, y = MeanDecGini, fill = layer)) +
    geom_col(width = 0.7, colour = "white", size = 0.3) +
    coord_flip() +
    scale_fill_manual(values = layer_cols, name = "Omics Layer") +
    labs(
      title    = paste("Top", top_n, "Risk Features –", model_result$disease),
      subtitle = "Metric: Mean Decrease in Gini (Random Forest OOB)",
      x        = NULL,
      y        = "Mean Decrease Gini"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey50", size = 9),
      panel.grid.major.y = element_blank(),
      legend.position    = "bottom"
    )
}


## 4c. Cross-disease feature importance heatmap -----------------------

#' plot_importance_heatmap
#'
#' Heatmap of the top-N features (union across all diseases) versus each
#' disease, scaled per column.
#'
#' @param models  Named list from \code{train_risk_models}.
#' @param top_n   Features per disease to include in the union set.
#' @return A ggplot object.

plot_importance_heatmap <- function(models, top_n = 15) {

  # Collect union of top-N features across all diseases
  top_feats <- unique(unlist(lapply(models, function(m) head(m$importance$feature, top_n))))

  # Build a wide matrix
  imp_wide <- do.call(cbind, lapply(models, function(m) {
    idx <- match(top_feats, m$importance$feature)
    vals <- m$importance$MeanDecGini[idx]
    vals[is.na(vals)] <- 0
    vals
  }))
  rownames(imp_wide) <- top_feats
  colnames(imp_wide) <- names(models)

  # Scale each disease column to [0,1]
  imp_scaled <- apply(imp_wide, 2, rescale)

  melt_df <- reshape2::melt(imp_scaled)
  colnames(melt_df) <- c("Feature", "Disease", "Importance")

  # Add layer annotation
  melt_df$Layer <- dplyr::case_when(
    grepl("^SNP_",  melt_df$Feature) ~ "Genomics",
    grepl("^GENE_", melt_df$Feature) ~ "Transcriptomics",
    TRUE                             ~ "Proteomics"
  )

  ggplot(melt_df, aes(x = Disease, y = Feature, fill = Importance)) +
    geom_tile(colour = "white", size = 0.5) +
    scale_fill_viridis_c(option = "plasma", name = "Scaled\nImportance") +
    facet_grid(Layer ~ ., scales = "free_y", space = "free_y") +
    labs(
      title    = "Cross-Disease Multi-Omics Feature Importance",
      subtitle = "Importance scaled 0–1 within each disease column",
      x        = NULL, y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x       = element_text(angle = 30, hjust = 1, face = "bold"),
      axis.text.y       = element_text(size = 7),
      strip.text.y      = element_text(angle = 0, face = "bold"),
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(colour = "grey50", size = 9),
      panel.grid        = element_blank(),
      legend.position   = "right"
    )
}


## 4d. Omics-layer contribution pie / donut per disease ---------------

#' plot_layer_contribution
#'
#' Donut chart showing what fraction of the top-K importance comes from
#' each omics layer, per disease.
#'
#' @param models  Named list from \code{train_risk_models}.
#' @param top_k   Number of top features to consider.
#' @return A ggplot object (faceted).

plot_layer_contribution <- function(models, top_k = 30) {

  layer_cols <- c(
    "Genomics"        = "#457B9D",
    "Transcriptomics" = "#2A9D8F",
    "Proteomics"      = "#F4A261"
  )

  agg <- do.call(rbind, lapply(names(models), function(d) {
    imp  <- head(models[[d]]$importance, top_k)
    tot  <- sum(imp$MeanDecGini)
    agg2 <- aggregate(MeanDecGini ~ layer, data = imp, FUN = sum)
    agg2$frac    <- agg2$MeanDecGini / tot
    agg2$disease <- d
    agg2
  }))

  ggplot(agg, aes(x = 2, y = frac, fill = layer)) +
    geom_col(colour = "white", size = 0.5) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    facet_wrap(~ disease) +
    scale_fill_manual(values = layer_cols, name = "Omics Layer") +
    labs(
      title    = "Omics Layer Contribution to Risk Prediction",
      subtitle = paste("Fraction of top-", top_k, " feature importance", sep = "")
    ) +
    theme_void() +
    theme(
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle   = element_text(colour = "grey50", size = 9, hjust = 0.5),
      strip.text      = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
}


## 4e. Comorbidity co-occurrence matrix heatmap ----------------------

#' plot_cooccurrence_heatmap
#'
#' Symmetric heatmap of pairwise phi co-occurrence coefficients.
#'
#' @param labels  Binary disease-label data.frame.
#' @return A ggplot object.

plot_cooccurrence_heatmap <- function(labels) {

  diseases <- colnames(labels)
  n        <- length(diseases)
  phi_mat  <- matrix(NA_real_, n, n, dimnames = list(diseases, diseases))
  diag(phi_mat) <- 1

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      tbl  <- table(labels[[i]], labels[[j]])
      if (all(dim(tbl) == c(2, 2))) {
        chisq_val <- chisq.test(tbl, correct = FALSE)$statistic
        # Coerce to numeric before multiplication to prevent integer overflow
        phi <- sqrt(chisq_val / sum(tbl)) * sign(
          as.numeric(tbl[2,2]) * as.numeric(tbl[1,1]) -
          as.numeric(tbl[2,1]) * as.numeric(tbl[1,2]))
      } else { phi <- 0 }
      phi_mat[i, j] <- phi_mat[j, i] <- round(phi, 4)
    }
  }

  melt_df <- reshape2::melt(phi_mat)
  colnames(melt_df) <- c("Disease1", "Disease2", "Phi")

  ggplot(melt_df, aes(x = Disease1, y = Disease2, fill = Phi)) +
    geom_tile(colour = "white", size = 1) +
    geom_text(aes(label = sprintf("%.3f", Phi)), size = 4.5, fontface = "bold") +
    scale_fill_gradient2(
      low  = "#457B9D", mid = "white", high = "#E63946",
      midpoint = 0, limits = c(-1, 1),
      name = "Phi\nCoefficient"
    ) +
    labs(
      title    = "Phenotypic Comorbidity Co-occurrence Matrix",
      subtitle = "Phi coefficient (chi-square-based)",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title  = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey50", size = 9),
      axis.text   = element_text(face = "bold"),
      panel.grid  = element_blank()
    )
}


# ------------------------------------------------------------------
# 5. Main execution
# ------------------------------------------------------------------

cat("\n=== Multi-Omics Comorbidity Network Pipeline ===\n\n")

## 5a. Simulate data
cat("Step 1: Simulating multi-omics data...\n")
omics <- simulate_multiomics(
  n_samples    = 600,
  n_genomic    = 60,
  n_transcript = 100,
  n_proteomic  = 50
)
cat(sprintf("  Subjects: %d  |  Genomic features: %d  |  Transcript: %d  |  Proteomic: %d\n",
            nrow(omics$labels),
            sum(grepl("^SNP_",  colnames(omics$genomics))),
            sum(grepl("^GENE_", colnames(omics$transcriptomics))),
            sum(grepl("^PROT_", colnames(omics$proteomics)))))
cat("  Disease prevalence:\n")
print(round(colMeans(omics$labels), 3))

## 5b. Build comorbidity network
cat("\nStep 2: Building comorbidity network...\n")
comorb_net <- build_comorbidity_network(omics, phi_thresh = 0.0)
cat(sprintf("  Nodes: %d  |  Edges: %d\n", vcount(comorb_net), ecount(comorb_net)))
cat("  Edge weights:\n")
print(as.data.frame(igraph::as_data_frame(comorb_net, what = "edges")))

## 5c. Train risk models
cat("\nStep 3: Training Random Forest risk models (may take ~30s)...\n")
risk_models <- train_risk_models(omics, ntree = 300)
for (d in names(risk_models)) {
  oob_err <- risk_models[[d]]$model$err.rate[300, "OOB"]
  cat(sprintf("  %-16s  OOB error: %.3f\n", d, oob_err))
}

## 5d. Produce all plots
cat("\nStep 4: Generating visualisations...\n")

p_network  <- plot_comorbidity_network(comorb_net)
p_imp_epi  <- plot_feature_importance(risk_models[["Epilepsy"]],      top_n = 20)
p_imp_t2d  <- plot_feature_importance(risk_models[["Type2Diabetes"]], top_n = 20)
p_imp_cvd  <- plot_feature_importance(risk_models[["CVD"]],           top_n = 20)
p_heatmap  <- plot_importance_heatmap(risk_models, top_n = 15)
p_donut    <- plot_layer_contribution(risk_models, top_k = 30)
p_phi_mat  <- plot_cooccurrence_heatmap(omics$labels)

## 5e. Render (use pdf() / png() in batch, or just print() interactively)
cat("  Rendering plots...\n")

pdf("comorbidity_multiomics_results.pdf", width = 12, height = 8, useDingbats = FALSE)

print(p_network)
print(p_phi_mat)

grid.arrange(p_imp_epi, p_imp_t2d, p_imp_cvd, ncol = 3,
             top = "Per-Disease Top Feature Importance (Random Forest)")

print(p_heatmap)
print(p_donut)

dev.off()

cat("\n=== Pipeline complete.  Output: comorbidity_multiomics_results.pdf ===\n")


# ------------------------------------------------------------------
# 6. Utility helpers (stand-alone calls)
# ------------------------------------------------------------------

#' get_top_features_for_disease
#'
#' Quick accessor: return top-N features for a named disease.
#'
#' @param models   Named list from \code{train_risk_models}.
#' @param disease  Character. One of the disease names.
#' @param n        Number of features to return.
#' @return data.frame with columns feature, MeanDecGini, layer.
#'
#' @examples
#' get_top_features_for_disease(risk_models, "Epilepsy", n = 10)

get_top_features_for_disease <- function(models, disease, n = 10) {
  if (!disease %in% names(models))
    stop("Disease not found. Available: ", paste(names(models), collapse = ", "))
  head(models[[disease]]$importance, n)
}


#' get_shared_features
#'
#' Returns features that appear in the top-N lists of ALL specified diseases —
#' a proxy for pan-disease risk features.
#'
#' @param models    Named list from \code{train_risk_models}.
#' @param diseases  Character vector of disease names (default: all).
#' @param top_n     Top-N list per disease to consider.
#' @return Character vector of shared feature names.
#'
#' @examples
#' get_shared_features(risk_models, top_n = 25)

get_shared_features <- function(models,
                                diseases = names(models),
                                top_n    = 20) {
  feat_lists <- lapply(diseases, function(d) head(models[[d]]$importance$feature, top_n))
  Reduce(intersect, feat_lists)
}


# Example calls (printed to console):
cat("\n--- Example utility calls ---\n")
cat("Top 10 features for Epilepsy:\n")
print(get_top_features_for_disease(risk_models, "Epilepsy", n = 10))

shared <- get_shared_features(risk_models, top_n = 30)
cat(sprintf("\nShared top-30 features across ALL three diseases (%d found):\n", length(shared)))
if (length(shared) > 0) print(shared) else cat("  (none – increase top_n)\n")


# ==================================================================
# Section 7: Graph Convolutional Network (GCN) for Individual Risk
# ==================================================================
#
# Architecture
# ────────────
#   Each SUBJECT is a graph node.  Edges connect omics-similar subjects
#   via mutual k-NN (intersection): edge (i,j) exists only when both
#   i nominates j AND j nominates i among their top-k neighbours.
#
#   2-layer GCN:
#     H1    = ReLU( A_hat @ X  @ W1 + b1 )        [n × hidden]
#     Logit = A_hat @ H1 @ W2 + b2                 [n × n_diseases]
#     Ŷ     = sigmoid(Logit)                        multi-label output
#
#   A_hat = D^{-1/2}(A + I)D^{-1/2}  (symmetric normalisation + self-loops)
#
# Training
# ────────
#   Binary cross-entropy loss, Adam optimiser, analytical backprop.
#   No external dependencies — pure base-R matrix algebra.
#
# Compatibility: R >= 3.6  (no torch / keras required)
# ==================================================================


# ------------------------------------------------------------------
# 7a.  Patient Similarity Graph
# ------------------------------------------------------------------

#' build_patient_graph
#'
#' Constructs a patient-similarity network from a feature matrix using
#' cosine similarity and one of two symmetrisation strategies:
#'
#'   "mutual" (default, recommended):
#'     Edge (i,j) is retained only when BOTH i nominates j AND j nominates i
#'     among their respective top-k neighbours.  Formally:
#'       A = pmin(A_directed, t(A_directed))
#'     Effect: sparser graph, lower hub degree, edges carry stronger
#'     reciprocal similarity signal.  Isolated nodes (patients with no
#'     reciprocated neighbour) are rescued by connecting them to their
#'     single most similar subject.
#'
#'   "union":
#'     Edge (i,j) is retained when EITHER i nominates j OR j nominates i.
#'     Formally: A = pmax(A_directed, t(A_directed))
#'     Effect: denser graph (mean degree approaches 2k in the limit),
#'     inflated hub degree, no isolated nodes by construction.
#'
#' Both strategies produce a symmetric (undirected) adjacency matrix
#' compatible with GCN normalisation.
#'
#' @param X          Numeric matrix [n x d].  Rows = subjects.
#' @param k          Integer. Number of nearest neighbours nominated per
#'                   subject before symmetrisation.
#' @param method     Character. "mutual" (intersection) or "union".
#' @param self_loops Logical. Add self-loops before returning (handled
#'                   separately in normalize_adjacency; leave FALSE).
#' @param verbose    Logical. Print graph diagnostics after construction.
#' @return Symmetric binary integer adjacency matrix [n x n].

build_patient_graph <- function(X,
                                k          = 15,
                                method     = c("mutual", "union"),
                                self_loops = FALSE,
                                verbose    = TRUE) {

  method <- match.arg(method)
  n      <- nrow(X)

  # ---- Row-normalise for cosine similarity ----
  norms  <- sqrt(rowSums(X^2)) + 1e-9
  X_norm <- X / norms

  # ---- Full cosine similarity matrix (n x n) ----
  S <- X_norm %*% t(X_norm)

  # ---- Directed nomination matrix ----
  # A_dir[i,j] = 1  iff  j is among i's top-k neighbours
  A_dir <- matrix(0L, n, n, dimnames = list(rownames(X), rownames(X)))
  for (i in seq_len(n)) {
    s_i       <- S[i, ]
    s_i[i]    <- -Inf                         # exclude self
    top_k_idx <- order(s_i, decreasing = TRUE)[seq_len(k)]
    A_dir[i, top_k_idx] <- 1L
  }

  # ---- Symmetrise according to chosen method ----
  if (method == "mutual") {
    # Intersection: both endpoints must have nominated each other
    A <- pmin(A_dir, t(A_dir))

    # Safety net: rescue isolated nodes (no reciprocated neighbour)
    # by connecting each to its single most-similar subject
    isolated <- which(rowSums(A) == 0L)
    if (length(isolated) > 0) {
      for (i in isolated) {
        s_i    <- S[i, ]; s_i[i] <- -Inf
        best_j <- which.max(s_i)
        A[i, best_j] <- A[best_j, i] <- 1L
      }
    }
  } else {
    # Union: either endpoint's nomination is sufficient
    A        <- pmax(A_dir, t(A_dir))
    isolated <- integer(0)          # union never produces isolated nodes
  }

  if (self_loops) diag(A) <- 1L

  # ---- Graph diagnostics ----
  if (verbose) {
    degrees <- rowSums(A)
    cat(sprintf(
      "  [PSN %s-kNN]  nodes=%d  edges=%d  degree: mean=%.1f  min=%d  max=%d  isolated_rescued=%d\n",
      method, n, sum(A) / 2L,
      mean(degrees), min(degrees), max(degrees), length(isolated)
    ))
  }

  A
}


# ------------------------------------------------------------------
# 7b.  GCN Primitives
# ------------------------------------------------------------------

#' normalize_adjacency
#'
#' Computes the symmetrically normalised adjacency with self-loops:
#'   A_hat = D^{-1/2} (A + I) D^{-1/2}
#'
#' @param A Square adjacency matrix (0/1 or weighted).
#' @return Dense numeric matrix of same dimensions.

normalize_adjacency <- function(A) {
  n       <- nrow(A)
  A_tilde <- A + diag(n)                          # add self-loops
  deg     <- rowSums(A_tilde)
  D_inv_sqrt <- diag(1 / sqrt(deg + 1e-9))
  D_inv_sqrt %*% A_tilde %*% D_inv_sqrt
}


# Activation helpers
.relu        <- function(x)  pmax(x, 0)
.relu_deriv  <- function(x)  (x > 0) * 1.0
.sigmoid     <- function(x)  1 / (1 + exp(-clamp(x, -30, 30)))
clamp        <- function(x, lo, hi) pmax(pmin(x, hi), lo)

# Inverted dropout (pass-through when rate = 0)
.dropout <- function(H, rate, training) {
  if (!training || rate <= 0) return(list(out = H, mask = NULL))
  mask <- matrix(
    rbinom(length(H), 1, 1 - rate) / (1 - rate),
    nrow = nrow(H), ncol = ncol(H)
  )
  list(out = H * mask, mask = mask)
}

# Binary cross-entropy (multi-label, averaged over subjects)
bce_loss_multilabel <- function(Y_hat, Y, eps = 1e-7) {
  -mean(Y * log(Y_hat + eps) + (1 - Y) * log(1 - Y_hat + eps))
}

# Xavier/Glorot initialisation
.xavier <- function(fan_in, fan_out) {
  limit <- sqrt(6 / (fan_in + fan_out))
  matrix(runif(fan_in * fan_out, -limit, limit), fan_in, fan_out)
}


# ------------------------------------------------------------------
# 7c.  Forward Pass
# ------------------------------------------------------------------

#' gcn_forward
#'
#' Two-layer GCN forward pass.
#'
#' @param A_hat    Normalised adjacency [n × n].
#' @param X        Feature matrix [n × d].
#' @param params   Named list: W1, b1, W2, b2.
#' @param drop_rate Dropout rate on H1 (0 = disabled).
#' @param training  Logical. Apply dropout only during training.
#' @return Named list: Y_hat, H1, H1_drop, logit, drop_mask.

gcn_forward <- function(A_hat, X, params, drop_rate = 0.1, training = TRUE) {

  # Layer 1
  AX    <- A_hat %*% X                            # n × d
  pre1  <- AX %*% params$W1 +
           matrix(params$b1, nrow(X), length(params$b1), byrow = TRUE)
  H1    <- .relu(pre1)                             # n × hidden

  # Dropout on H1
  dp    <- .dropout(H1, drop_rate, training)
  H1d   <- dp$out

  # Layer 2
  AH1   <- A_hat %*% H1d                          # n × hidden
  logit <- AH1 %*% params$W2 +
           matrix(params$b2, nrow(X), length(params$b2), byrow = TRUE)
  Y_hat <- .sigmoid(logit)                         # n × n_diseases

  list(Y_hat    = Y_hat,
       H1       = H1,
       H1_drop  = H1d,
       pre1     = pre1,
       logit    = logit,
       AX       = AX,
       AH1      = AH1,
       drop_mask= dp$mask)
}


# ------------------------------------------------------------------
# 7d.  Analytical Backward Pass
# ------------------------------------------------------------------

#' gcn_backward
#'
#' Computes analytical gradients for W1, b1, W2, b2 via chain rule.
#' Derivation (combined sigmoid + BCE shortcut):
#'   dL/d(logit) = (Ŷ - Y) / n
#'
#' @param A_hat   Normalised adjacency.
#' @param cache   List returned by \code{gcn_forward}.
#' @param params  Current parameter list.
#' @param Y       Binary label matrix [n × n_diseases].
#' @param drop_rate Dropout rate (needed to scale gradients).
#' @return Named list of gradients: dW1, db1, dW2, db2.

gcn_backward <- function(A_hat, cache, params, Y, drop_rate = 0.1) {

  n     <- nrow(Y)
  Y_hat <- cache$Y_hat

  # --- Output layer ---
  # Combined sigmoid + BCE gradient: dL/d(logit) = (Ŷ - Y) / n
  d_logit <- (Y_hat - Y) / n                      # n × n_diseases

  dW2     <- t(cache$AH1) %*% d_logit             # hidden × n_diseases
  db2     <- colSums(d_logit)

  # --- Propagate through A_hat (layer 2 aggregation) ---
  d_H1d   <- A_hat %*% d_logit %*% t(params$W2)   # n × hidden

  # --- Undo dropout ---
  if (!is.null(cache$drop_mask)) {
    d_H1d <- d_H1d * cache$drop_mask
  }

  # --- ReLU derivative ---
  d_pre1  <- d_H1d * .relu_deriv(cache$pre1)      # n × hidden

  # --- Layer 1 ---
  dW1     <- t(cache$AX) %*% d_pre1               # d × hidden
  db1     <- colSums(d_pre1)

  list(dW1 = dW1, db1 = db1, dW2 = dW2, db2 = db2)
}


# ------------------------------------------------------------------
# 7e.  Adam Optimiser
# ------------------------------------------------------------------

#' adam_update
#'
#' One Adam step with optional L2 weight decay.
#'
#' @param params     Named list of current parameters.
#' @param grads      Named list of gradients (same names).
#' @param state      Adam moment estimates (m, v, t); initialised internally.
#' @param lr         Learning rate.
#' @param beta1,beta2 Adam moment decay rates.
#' @param epsilon    Numerical stability constant.
#' @param lambda     L2 weight-decay coefficient.
#' @return List: updated \code{params} and \code{state}.

adam_update <- function(params, grads, state,
                        lr      = 1e-3,
                        beta1   = 0.9,
                        beta2   = 0.999,
                        epsilon = 1e-8,
                        lambda  = 1e-4) {

  state$t <- state$t + 1L
  t       <- state$t

  for (nm in names(grads)) {
    g <- grads[[nm]]

    # L2 regularisation (weight decay) – skip biases
    if (!grepl("^b", nm)) g <- g + lambda * params[[nm]]

    # Moment updates
    state$m[[nm]] <- beta1 * state$m[[nm]] + (1 - beta1) * g
    state$v[[nm]] <- beta2 * state$v[[nm]] + (1 - beta2) * g^2

    # Bias correction
    m_hat <- state$m[[nm]] / (1 - beta1^t)
    v_hat <- state$v[[nm]] / (1 - beta2^t)

    params[[nm]] <- params[[nm]] - lr * m_hat / (sqrt(v_hat) + epsilon)
  }

  list(params = params, state = state)
}


# ------------------------------------------------------------------
# 7f.  Training Loop
# ------------------------------------------------------------------

#' train_gcn
#'
#' Builds the patient similarity graph, initialises a 2-layer GCN and
#' trains it with Adam on multi-label binary cross-entropy loss.
#'
#' @param omics_list  List from \code{simulate_multiomics}.
#' @param k_neighbors Integer. kNN neighbours for graph construction.
#' @param hidden_dim  Integer. GCN hidden layer width.
#' @param epochs      Integer. Training epochs.
#' @param lr          Numeric. Adam learning rate.
#' @param drop_rate   Numeric. Dropout on hidden layer (0 = disabled).
#' @param lambda      Numeric. L2 weight decay.
#' @param val_frac    Numeric. Fraction of subjects held out for validation.
#' @param verbose     Integer. Print loss every \code{verbose} epochs.
#'
#' @return A named list (the "GCN model object"):
#'   \item{params}      Final weight matrices.
#'   \item{A_hat}       Normalised adjacency.
#'   \item{A_raw}       Raw binary adjacency.
#'   \item{X}           Feature matrix.
#'   \item{Y}           Label matrix.
#'   \item{subject_ids} Character vector.
#'   \item{feature_names} Character vector.
#'   \item{disease_names} Character vector.
#'   \item{history}     data.frame of train/val loss per epoch.
#'   \item{val_idx}     Validation subject indices.
#'   \item{train_idx}   Training subject indices.

train_gcn <- function(
    omics_list,
    k_neighbors = 15,
    graph_method = c("mutual", "union"),
    hidden_dim  = 64,
    epochs      = 300,
    lr          = 5e-3,
    drop_rate   = 0.1,
    lambda      = 1e-4,
    val_frac    = 0.15,
    verbose     = 50
) {

  graph_method <- match.arg(graph_method)

  # --- Feature matrix ---
  feat_g <- grep("^SNP_",  colnames(omics_list$genomics),        value = TRUE)
  feat_t <- grep("^GENE_", colnames(omics_list$transcriptomics), value = TRUE)
  feat_p <- grep("^PROT_", colnames(omics_list$proteomics),      value = TRUE)

  X_raw <- cbind(
    omics_list$genomics[,        feat_g],
    omics_list$transcriptomics[, feat_t],
    omics_list$proteomics[,      feat_p]
  )

  # Z-score normalise each feature column
  col_means <- colMeans(X_raw)
  col_sds   <- apply(X_raw, 2, sd) + 1e-9
  X <- scale(X_raw, center = col_means, scale = col_sds)

  Y <- as.matrix(omics_list$labels)
  n <- nrow(X)
  d <- ncol(X)
  n_out <- ncol(Y)

  # --- Train / val split ---
  val_n     <- max(1L, floor(n * val_frac))
  val_idx   <- sample(n, val_n)
  train_idx <- setdiff(seq_len(n), val_idx)

  # --- Graph ---
  cat("  Building patient similarity graph (k =", k_neighbors,
      "| method =", graph_method, ")...\n")
  A_raw <- build_patient_graph(X, k = k_neighbors,
                               method = graph_method, verbose = TRUE)
  A_hat <- normalize_adjacency(A_raw)
  cat(sprintf("  Normalised adjacency ready.  Density: %.4f\n",
              mean(A_raw[upper.tri(A_raw)])))

  # --- Parameter initialisation ---
  params <- list(
    W1 = .xavier(d,       hidden_dim),
    b1 = rep(0, hidden_dim),
    W2 = .xavier(hidden_dim, n_out),
    b2 = rep(0, n_out)
  )

  # --- Adam state ---
  adam_state <- list(
    t = 0L,
    m = lapply(params, function(p) p * 0),
    v = lapply(params, function(p) p * 0)
  )

  history <- data.frame(epoch = integer(), train_loss = numeric(),
                        val_loss = numeric())

  cat("  Training GCN...\n")
  for (ep in seq_len(epochs)) {

    # Forward (training mode)
    cache <- gcn_forward(A_hat, X, params, drop_rate, training = TRUE)
    tr_loss <- bce_loss_multilabel(cache$Y_hat[train_idx, , drop = FALSE],
                                   Y[train_idx, , drop = FALSE])

    # Backward
    grads <- gcn_backward(A_hat, cache, params, Y, drop_rate)

    # Adam step
    upd    <- adam_update(params, grads, adam_state,
                          lr = lr, lambda = lambda)
    params     <- upd$params
    adam_state <- upd$state

    # Validation loss (no dropout)
    val_cache  <- gcn_forward(A_hat, X, params, 0, training = FALSE)
    val_loss   <- bce_loss_multilabel(val_cache$Y_hat[val_idx, , drop = FALSE],
                                     Y[val_idx, , drop = FALSE])

    history <- rbind(history,
                     data.frame(epoch = ep, train_loss = tr_loss,
                                val_loss = val_loss))

    if (verbose > 0 && (ep %% verbose == 0 || ep == 1)) {
      cat(sprintf("    Epoch %4d  train_loss=%.4f  val_loss=%.4f\n",
                  ep, tr_loss, val_loss))
    }
  }

  # Final predictions (inference mode)
  final_cache <- gcn_forward(A_hat, X, params, 0, training = FALSE)

  # gcn_forward returns Y_hat from raw matrix algebra which drops dimnames.
  # Restore row (subject) and column (disease) names so that downstream
  # column-name indexing (e.g. probs[, "Epilepsy"]) works correctly.
  Y_hat_named <- final_cache$Y_hat
  rownames(Y_hat_named) <- omics_list$subject_ids
  colnames(Y_hat_named) <- colnames(Y)

  list(
    params        = params,
    A_hat         = A_hat,
    A_raw         = A_raw,
    X             = X,
    X_raw         = X_raw,
    Y             = Y,
    Y_hat         = Y_hat_named,
    subject_ids   = omics_list$subject_ids,
    feature_names = c(feat_g, feat_t, feat_p),
    disease_names = colnames(Y),
    col_means     = col_means,
    col_sds       = col_sds,
    history       = history,
    val_idx       = val_idx,
    train_idx     = train_idx
  )
}


# ------------------------------------------------------------------
# 7g.  Individual Risk Prediction
# ------------------------------------------------------------------

#' predict_individual_risk_gcn
#'
#' Returns predicted disease-risk probabilities for one or more subjects
#' identified by their subject IDs.
#'
#' @param gcn_model  Object returned by \code{train_gcn}.
#' @param subject_ids Character vector of subject IDs to query.
#'                    Defaults to all subjects.
#' @return data.frame: subject_id, one column per disease (probability),
#'         plus true-label columns (prefixed "true_").

predict_individual_risk_gcn <- function(gcn_model,
                                        subject_ids = NULL) {

  if (is.null(subject_ids)) subject_ids <- gcn_model$subject_ids

  idx <- match(subject_ids, gcn_model$subject_ids)
  bad <- subject_ids[is.na(idx)]
  if (length(bad) > 0)
    stop("Unknown subject IDs: ", paste(bad, collapse = ", "))

  probs  <- gcn_model$Y_hat[idx, , drop = FALSE]
  true_Y <- gcn_model$Y[idx, , drop = FALSE]

  # Use integer column positions rather than string names as a
  # belt-and-braces guard: matrix algebra can silently drop dimnames,
  # and probs[, "Epilepsy"] throws 'no dimnames attribute' in that case.
  df <- data.frame(subject_id = subject_ids, stringsAsFactors = FALSE)
  for (k in seq_along(gcn_model$disease_names)) {
    d <- gcn_model$disease_names[k]
    df[[paste0("risk_", d)]] <- round(probs[, k], 4)
  }
  for (k in seq_along(gcn_model$disease_names)) {
    d <- gcn_model$disease_names[k]
    df[[paste0("true_", d)]] <- true_Y[, k]
  }
  df
}


# ------------------------------------------------------------------
# 7h.  Perturbation-based Feature Saliency (per individual)
# ------------------------------------------------------------------

#' compute_gcn_saliency
#'
#' Measures feature importance for a specific subject by zeroing out
#' each feature (one at a time) in the subject's row and recording the
#' absolute change in predicted risk for the target disease.
#'
#' @param gcn_model  Object from \code{train_gcn}.
#' @param subject_id Character. Subject to explain.
#' @param disease    Character. Disease to explain (default: all, returns list).
#' @param n_top      Integer. Number of top features to return per disease.
#' @return data.frame: feature, saliency, layer (for one disease) or
#'         named list of data.frames (all diseases).

compute_gcn_saliency <- function(gcn_model,
                                 subject_id,
                                 disease = NULL,
                                 n_top   = 20) {

  idx <- match(subject_id, gcn_model$subject_ids)
  if (is.na(idx)) stop("Unknown subject ID: ", subject_id)

  X_orig  <- gcn_model$X
  A_hat   <- gcn_model$A_hat
  params  <- gcn_model$params
  d_names <- gcn_model$disease_names
  f_names <- gcn_model$feature_names
  n_feat  <- length(f_names)

  # Baseline prediction
  # gcn_forward drops dimnames; assign them so vector indexing by name is safe
  base_cache <- gcn_forward(A_hat, X_orig, params, 0, training = FALSE)
  base_prob  <- base_cache$Y_hat[idx, ]
  names(base_prob) <- gcn_model$disease_names

  diseases_to_explain <- if (is.null(disease)) d_names else disease

  results <- lapply(diseases_to_explain, function(dis) {

    d_col <- which(d_names == dis)

    saliency_vals <- vapply(seq_len(n_feat), function(j) {
      X_perturb        <- X_orig
      X_perturb[idx, j] <- 0          # ablate feature j for subject i
      cache_p <- gcn_forward(A_hat, X_perturb, params, 0, training = FALSE)
      abs(cache_p$Y_hat[idx, d_col] - base_prob[d_col])
    }, numeric(1))

    sal_df <- data.frame(
      feature  = f_names,
      saliency = saliency_vals,
      layer    = dplyr::case_when(
        grepl("^SNP_",  f_names) ~ "Genomics",
        grepl("^GENE_", f_names) ~ "Transcriptomics",
        TRUE                     ~ "Proteomics"
      ),
      baseline_risk = round(base_prob[d_col], 4),
      disease  = dis,
      stringsAsFactors = FALSE
    )
    sal_df <- sal_df[order(sal_df$saliency, decreasing = TRUE), ]
    head(sal_df, n_top)
  })

  names(results) <- diseases_to_explain
  if (length(results) == 1) results[[1]] else results
}


# ==================================================================
# Section 8: GCN Visualisation Functions
# ==================================================================

LAYER_COLS <- c(
  "Genomics"        = "#457B9D",
  "Transcriptomics" = "#2A9D8F",
  "Proteomics"      = "#F4A261"
)


## 8a. Training curve --------------------------------------------------

#' plot_gcn_training_curve
#'
#' Line chart of training and validation BCE loss over epochs.
#'
#' @param gcn_model Object from \code{train_gcn}.
#' @return ggplot object.

plot_gcn_training_curve <- function(gcn_model) {

  hist <- gcn_model$history
  melt_hist <- reshape2::melt(hist, id.vars = "epoch",
                              variable.name = "Split",
                              value.name    = "Loss")
  melt_hist$Split <- ifelse(melt_hist$Split == "train_loss", "Train", "Validation")

  ggplot(melt_hist, aes(x = epoch, y = Loss, colour = Split, linetype = Split)) +
    geom_line(size = 1) +
    scale_colour_manual(values = c(Train = "#457B9D", Validation = "#E63946")) +
    scale_linetype_manual(values = c(Train = "solid", Validation = "dashed")) +
    labs(
      title    = "GCN Training Curve",
      subtitle = "Multi-label binary cross-entropy loss",
      x        = "Epoch",
      y        = "BCE Loss",
      colour   = NULL, linetype = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey50", size = 9),
      legend.position = "top"
    )
}


## 8b. Patient graph embedding (2D force-directed) --------------------

#' plot_patient_graph
#'
#' Visualises the patient similarity graph as a 2D force-directed layout.
#' Nodes are coloured by a selected disease's true label, sized by GCN
#' predicted risk.  Optionally highlights a specific subject.
#'
#' @param gcn_model        Object from \code{train_gcn}.
#' @param color_disease    Character. Disease for node colour.
#' @param highlight_id     Character or NULL. Subject ID to highlight.
#' @param max_nodes        Integer. Subsample graph for speed (igraph layout).
#' @return ggplot object.

plot_patient_graph <- function(gcn_model,
                               color_disease = NULL,
                               highlight_id  = NULL,
                               max_nodes     = 300) {

  if (is.null(color_disease)) color_disease <- gcn_model$disease_names[1]
  d_col <- which(gcn_model$disease_names == color_disease)

  n     <- length(gcn_model$subject_ids)
  idx   <- if (n > max_nodes) sort(sample(n, max_nodes)) else seq_len(n)
  sids  <- gcn_model$subject_ids[idx]

  A_sub <- gcn_model$A_raw[idx, idx]
  g_sub <- graph_from_adjacency_matrix(A_sub, mode = "undirected",
                                       diag = FALSE)

  # Layout
  layout_mat <- layout_with_fr(g_sub, niter = 500)
  # layout_with_fr() returns a plain matrix with no rownames.
  # Assign vertex names so character indexing via edge_list$from/to works.
  rownames(layout_mat) <- V(g_sub)$name

  node_df <- data.frame(
    subject_id   = sids,
    x            = layout_mat[, 1],
    y            = layout_mat[, 2],
    true_label   = factor(gcn_model$Y[idx, d_col]),
    pred_risk    = gcn_model$Y_hat[idx, d_col],
    highlighted  = if (!is.null(highlight_id)) sids == highlight_id else FALSE,
    stringsAsFactors = FALSE
  )

  # Edges
  # Convert vertex names to integer row positions in layout_mat as a
  # belt-and-braces guard: match() is safe whether vertices are named or not.
  edge_list  <- as_data_frame(g_sub, what = "edges")
  vnames     <- rownames(layout_mat)          # vertex name -> row index map
  from_idx   <- match(edge_list$from, vnames)
  to_idx     <- match(edge_list$to,   vnames)
  edge_df    <- data.frame(
    x    = layout_mat[from_idx, 1],
    y    = layout_mat[from_idx, 2],
    xend = layout_mat[to_idx,   1],
    yend = layout_mat[to_idx,   2]
  )

  p <- ggplot() +
    geom_segment(data = edge_df,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey80", size = 0.2, alpha = 0.5) +
    geom_point(data = node_df,
               aes(x = x, y = y,
                   colour = true_label,
                   size   = pred_risk,
                   alpha  = pred_risk),
               shape = 16) +
    scale_colour_manual(
      values = c("0" = "#A8DADC", "1" = "#E63946"),
      labels = c("0" = "Control", "1" = "Case"),
      name   = paste(color_disease, "\nTrue Label")
    ) +
    scale_size(range = c(1.5, 6), name = "GCN\nPred. Risk") +
    scale_alpha(range = c(0.4, 1), guide = "none")

  if (!is.null(highlight_id) && highlight_id %in% sids) {
    hi <- node_df[node_df$subject_id == highlight_id, ]
    p  <- p +
      geom_point(data = hi,
                 aes(x = x, y = y), size = 10,
                 shape = 1, colour = "#F4A261", stroke = 2) +
      geom_text(data = hi,
                aes(x = x, y = y, label = subject_id),
                vjust = -1.3, fontface = "bold",
                colour = "#F4A261", size = 3.5)
  }

  p +
    labs(
      title    = "Patient Similarity Network (GCN)",
      subtitle = paste0("Force-directed layout  |  colour = ", color_disease,
                        " true label  |  size = GCN predicted risk"),
      x = NULL, y = NULL
    ) +
    theme_void() +
    theme(
      plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(colour = "grey50", size = 9,  hjust = 0.5),
      legend.position = "right"
    )
}


## 8c. Per-individual GCN saliency bar chart --------------------------

#' plot_gnn_saliency_individual
#'
#' Horizontal bar chart showing the features that most influenced
#' the GCN's risk prediction for a specific subject and disease.
#'
#' @param gcn_model  Object from \code{train_gcn}.
#' @param subject_id Character.
#' @param disease    Character.
#' @param top_n      Number of features to display.
#' @return ggplot object.

plot_gnn_saliency_individual <- function(gcn_model,
                                         subject_id,
                                         disease,
                                         top_n = 20) {

  sal <- compute_gcn_saliency(gcn_model, subject_id, disease, n_top = top_n)
  if (is.list(sal)) sal <- sal[[disease]]   # handle list return

  sal$feature <- factor(sal$feature, levels = rev(sal$feature))

  base_risk <- unique(sal$baseline_risk)
  true_label <- gcn_model$Y[match(subject_id, gcn_model$subject_ids),
                             disease]

  ggplot(sal, aes(x = feature, y = saliency, fill = layer)) +
    geom_col(width = 0.72, colour = "white", size = 0.3) +
    coord_flip() +
    scale_fill_manual(values = LAYER_COLS, name = "Omics Layer") +
    labs(
      title    = paste("GCN Feature Saliency –", subject_id),
      subtitle = sprintf(
        "Disease: %s  |  Predicted risk: %.3f  |  True label: %s",
        disease, base_risk, ifelse(true_label == 1, "Case", "Control")
      ),
      x = NULL,
      y = "Perturbation-based Saliency (|ΔRisk|)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(colour = "grey50", size = 9),
      panel.grid.major.y = element_blank(),
      legend.position    = "bottom"
    )
}


## 8d. Individual multi-disease risk profile (GCN vs RF) --------------

#' plot_individual_risk_profile
#'
#' Grouped bar chart comparing GCN and Random Forest predicted risk
#' probabilities for a single subject across all three diseases.
#'
#' @param gcn_model   Object from \code{train_gcn}.
#' @param rf_models   Named list from \code{train_risk_models}.
#' @param subject_id  Character.
#' @return ggplot object.

plot_individual_risk_profile <- function(gcn_model, rf_models, subject_id) {

  idx    <- match(subject_id, gcn_model$subject_ids)
  d_names <- gcn_model$disease_names

  # GCN predicted risks
  gcn_probs <- setNames(gcn_model$Y_hat[idx, ], d_names)

  # RF predicted risks (OOB probability for subject, using predict)
  rf_probs <- sapply(d_names, function(d) {
    feat_names <- gcn_model$feature_names
    X_row  <- gcn_model$X_raw[idx, feat_names, drop = FALSE]
    pred   <- predict(rf_models[[d]]$model, newdata = X_row, type = "prob")
    pred[1, "Case"]
  })

  # True labels
  true_Y <- gcn_model$Y[idx, ]

  df <- data.frame(
    disease   = rep(d_names, 2),
    model     = c(rep("GCN", length(d_names)), rep("Random Forest", length(d_names))),
    risk      = c(gcn_probs, rf_probs),
    true_label = rep(true_Y, 2),
    stringsAsFactors = FALSE
  )

  # Add jitter-friendly true-label indicator
  true_df <- data.frame(
    disease    = d_names,
    true_label = as.integer(true_Y),
    stringsAsFactors = FALSE
  )

  ggplot(df, aes(x = disease, y = risk, fill = model)) +
    geom_col(position = position_dodge(width = 0.7),
             width = 0.6, colour = "white") +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               colour = "grey40", size = 0.5) +
    geom_point(data = true_df,
               aes(x = disease, y = true_label + 0.02, shape = NULL,
                   colour = factor(true_label)),
               inherit.aes = FALSE, size = 5,
               shape = 23, fill = "white", stroke = 2) +
    scale_fill_manual(
      values = c("GCN" = "#E63946", "Random Forest" = "#457B9D"),
      name = "Model"
    ) +
    scale_colour_manual(
      values = c("0" = "#A8DADC", "1" = "#E63946"),
      labels = c("0" = "Control", "1" = "Case"),
      name = "True Label"
    ) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    labs(
      title    = paste("Individual Risk Profile –", subject_id),
      subtitle = "Diamond = true label (white = Control, red = Case)  |  dashed = 0.5 decision boundary",
      x        = NULL,
      y        = "Predicted Risk Probability"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(colour = "grey50", size = 9),
      legend.position = "top"
    )
}


## 8e. Population risk landscape (GCN) --------------------------------

#' plot_risk_landscape
#'
#' Faceted density ridgeline showing GCN predicted-risk distribution
#' stratified by true case/control status, per disease.
#' Falls back to density plots if ggridges is unavailable.
#'
#' @param gcn_model Object from \code{train_gcn}.
#' @return ggplot object.

plot_risk_landscape <- function(gcn_model) {

  d_names <- gcn_model$disease_names
  rows <- do.call(rbind, lapply(seq_along(d_names), function(k) {
    d <- d_names[k]
    data.frame(
      disease    = d,
      risk       = gcn_model$Y_hat[, k],
      true_label = factor(gcn_model$Y[, k], labels = c("Control", "Case")),
      stringsAsFactors = FALSE
    )
  }))

  ggplot(rows, aes(x = risk, fill = true_label, colour = true_label)) +
    geom_density(alpha = 0.45, size = 0.8, adjust = 1.2) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               colour = "grey30", size = 0.5) +
    facet_wrap(~ disease, ncol = 1, scales = "free_y") +
    scale_fill_manual(
      values = c("Control" = "#A8DADC", "Case" = "#E63946"),
      name = "True Label"
    ) +
    scale_colour_manual(
      values = c("Control" = "#2A9D8F", "Case" = "#C1121F"),
      guide  = "none"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title    = "GCN Population Risk Landscape",
      subtitle = "Predicted risk distribution stratified by true case/control",
      x        = "GCN Predicted Risk",
      y        = "Density"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 13),
      plot.subtitle   = element_text(colour = "grey50", size = 9),
      strip.text      = element_text(face = "bold", size = 11),
      legend.position = "top"
    )
}


## 8f. Cross-subject saliency heatmap ---------------------------------

#' plot_saliency_heatmap_subjects
#'
#' Heatmap of perturbation saliency for a set of subjects × top features
#' for a given disease.  Reveals which features are consistently important
#' versus subject-specific drivers.
#'
#' @param gcn_model   Object from \code{train_gcn}.
#' @param subject_ids Character vector of subjects to include.
#' @param disease     Character. Target disease.
#' @param top_n       Number of top features (union across subjects).
#' @return ggplot object.

plot_saliency_heatmap_subjects <- function(gcn_model,
                                           subject_ids,
                                           disease,
                                           top_n = 20) {

  # Compute saliency for each subject
  all_sal <- lapply(subject_ids, function(sid) {
    sal <- compute_gcn_saliency(gcn_model, sid, disease, n_top = top_n)
    if (is.list(sal)) sal <- sal[[disease]]
    sal$subject_id <- sid
    sal
  })

  sal_df <- do.call(rbind, all_sal)

  # Union of top features
  top_feats <- unique(sal_df$feature[sal_df$saliency > 0])
  top_feats <- head(top_feats, top_n)

  wide <- reshape2::dcast(
    sal_df[sal_df$feature %in% top_feats, ],
    subject_id ~ feature,
    value.var = "saliency",
    fill = 0
  )
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$subject_id

  # Scale by row (subject-level normalisation)
  mat_scaled <- t(apply(mat, 1, rescale))

  melt_df <- reshape2::melt(mat_scaled)
  colnames(melt_df) <- c("subject_id", "feature", "saliency")

  # Add layer info
  melt_df$layer <- dplyr::case_when(
    grepl("^SNP_",  melt_df$feature) ~ "Genomics",
    grepl("^GENE_", melt_df$feature) ~ "Transcriptomics",
    TRUE                             ~ "Proteomics"
  )

  ggplot(melt_df, aes(x = subject_id, y = feature, fill = saliency)) +
    geom_tile(colour = "white", size = 0.4) +
    scale_fill_viridis_c(option = "inferno", name = "Scaled\nSaliency") +
    facet_grid(layer ~ ., scales = "free_y", space = "free_y") +
    labs(
      title    = paste("Per-Subject GCN Saliency –", disease),
      subtitle = "Rows = subjects  |  Columns = features  |  scaled per subject",
      x        = "Subject ID", y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x       = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y       = element_text(size  = 7),
      strip.text.y      = element_text(angle = 0, face = "bold"),
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(colour = "grey50", size = 9),
      panel.grid        = element_blank(),
      legend.position   = "right"
    )
}


# ==================================================================
# Section 8b. Utility: Compare mutual vs union graph topology
# ==================================================================

#' compare_graph_methods
#'
#' Builds both mutual and union kNN graphs on the same feature matrix
#' and returns a side-by-side summary of their topological properties.
#' Useful for choosing k and method before committing to a full training run.
#'
#' Properties reported per method:
#'   n_edges       Total undirected edges
#'   density       Fraction of possible edges present
#'   mean_degree   Average node degree
#'   max_degree    Highest node degree (hub indicator)
#'   min_degree    Lowest node degree (isolation risk)
#'   n_isolated    Nodes with degree 0 BEFORE the rescue step (mutual only)
#'   gini_degree   Gini coefficient of degree distribution (0=uniform, 1=hub)
#'
#' @param X  Numeric matrix [n x d] — the same scaled feature matrix
#'            passed to train_gcn.
#' @param k  Integer. Number of neighbours to nominate per subject.
#' @return   data.frame with one row per method.
#'
#' @examples
#' compare_graph_methods(gcn_model$X, k = 15)

compare_graph_methods <- function(X, k = 15) {

  .gini <- function(x) {
    x   <- sort(x)
    n   <- length(x)
    2 * sum(x * seq_len(n)) / (n * sum(x)) - (n + 1) / n
  }

  results <- lapply(c("mutual", "union"), function(m) {
    A   <- build_patient_graph(X, k = k, method = m, verbose = FALSE)
    deg <- rowSums(A)
    n   <- nrow(X)

    # For mutual: count nodes with degree 0 BEFORE rescue
    # Re-run directed step to get pre-rescue count
    norms  <- sqrt(rowSums(X^2)) + 1e-9
    X_norm <- X / norms
    S      <- X_norm %*% t(X_norm)
    A_dir  <- matrix(0L, n, n)
    for (i in seq_len(n)) {
      s_i <- S[i, ]; s_i[i] <- -Inf
      top_k_idx <- order(s_i, decreasing = TRUE)[seq_len(k)]
      A_dir[i, top_k_idx] <- 1L
    }
    A_pre <- if (m == "mutual") pmin(A_dir, t(A_dir)) else pmax(A_dir, t(A_dir))
    n_isolated_pre <- sum(rowSums(A_pre) == 0L)

    data.frame(
      method       = m,
      n_edges      = sum(A) / 2L,
      density      = round(mean(A[upper.tri(A)]), 5),
      mean_degree  = round(mean(deg), 2),
      median_degree= median(deg),
      max_degree   = max(deg),
      min_degree   = min(deg),
      n_isolated_pre_rescue = n_isolated_pre,
      gini_degree  = round(.gini(deg), 4),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  cat("\n--- Graph topology comparison (k =", k, ") ---\n")
  print(out, row.names = FALSE)
  invisible(out)
}

# ==================================================================

cat("\n\n=== GCN Risk Prediction Pipeline ===\n\n")

## 9a. Train
cat("Step 5: Training Graph Convolutional Network...\n")

# Optional topology audit — compare mutual vs union before committing
# (uses the same scaled X built inside train_gcn; replicate scaling here)
{
  feat_g <- grep("^SNP_",  colnames(omics$genomics),        value = TRUE)
  feat_t <- grep("^GENE_", colnames(omics$transcriptomics), value = TRUE)
  feat_p <- grep("^PROT_", colnames(omics$proteomics),      value = TRUE)
  X_audit <- scale(cbind(omics$genomics[, feat_g],
                         omics$transcriptomics[, feat_t],
                         omics$proteomics[, feat_p]))
  compare_graph_methods(X_audit, k = 15)
}
gcn_model <- train_gcn(
  omics_list   = omics,
  k_neighbors  = 15,
  graph_method = "mutual",      # intersection kNN: edge only if both endpoints agree
  hidden_dim   = 64,
  epochs       = 300,
  lr           = 5e-3,
  drop_rate    = 0.10,
  lambda       = 1e-4,
  val_frac     = 0.15,
  verbose      = 50
)

## 9b. Evaluate
cat("\nStep 6: GCN evaluation summary...\n")
for (d in gcn_model$disease_names) {
  d_col   <- which(gcn_model$disease_names == d)
  pred_cls <- as.integer(gcn_model$Y_hat[, d_col] >= 0.5)
  true_cls <- gcn_model$Y[, d_col]
  acc      <- mean(pred_cls == true_cls)
  # AUC (rank-based, no extra packages)
  pos  <- gcn_model$Y_hat[true_cls == 1, d_col]
  neg  <- gcn_model$Y_hat[true_cls == 0, d_col]
  auc  <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
  cat(sprintf("  %-16s  Accuracy: %.3f  AUC: %.3f\n", d, acc, auc))
}

## 9c. Example predictions for 5 specific subjects
cat("\nStep 7: Individual risk predictions (GCN)...\n")
sample_ids <- head(gcn_model$subject_ids, 5)
preds <- predict_individual_risk_gcn(gcn_model, sample_ids)
print(preds)

## 9d. Generate GCN plots
cat("\nStep 8: Generating GCN visualisations...\n")

p_train_curve <- plot_gcn_training_curve(gcn_model)
p_patient_net <- plot_patient_graph(gcn_model,
                                    color_disease = "Epilepsy",
                                    highlight_id  = sample_ids[1])

# Saliency for the first subject across each disease
p_sal_epi <- plot_gnn_saliency_individual(gcn_model, sample_ids[1], "Epilepsy",   top_n = 20)
p_sal_t2d <- plot_gnn_saliency_individual(gcn_model, sample_ids[1], "Type2Diabetes", top_n = 20)
p_sal_cvd <- plot_gnn_saliency_individual(gcn_model, sample_ids[1], "CVD",        top_n = 20)

p_risk_profile  <- plot_individual_risk_profile(gcn_model, risk_models, sample_ids[1])
p_risk_landscape <- plot_risk_landscape(gcn_model)

# Cross-subject saliency heatmap (first 10 subjects, Epilepsy)
p_sal_heatmap <- plot_saliency_heatmap_subjects(
  gcn_model, head(gcn_model$subject_ids, 12), "Epilepsy", top_n = 20
)

## 9e. Write GCN plots to PDF
pdf("comorbidity_multiomics_GCN.pdf", width = 12, height = 8, useDingbats = FALSE)

print(p_train_curve)
print(p_patient_net)

grid.arrange(p_sal_epi, p_sal_t2d, p_sal_cvd, ncol = 3,
             top = paste("GCN Perturbation Saliency –", sample_ids[1]))

print(p_risk_profile)
print(p_risk_landscape)
print(p_sal_heatmap)

dev.off()

cat("\n=== GCN pipeline complete.  Output: comorbidity_multiomics_GCN.pdf ===\n")


# ------------------------------------------------------------------
# Utility: batch saliency for a named subject + disease
# ------------------------------------------------------------------

#' explain_subject
#'
#' One-shot explanation for a subject: prints predicted risks from both
#' GCN and RF, then returns the saliency table for every disease.
#'
#' @param subject_id Character.
#' @param gcn_model  Object from \code{train_gcn}.
#' @param rf_models  Named list from \code{train_risk_models}.
#' @return Invisibly, named list of saliency data.frames.
#'
#' @examples
#' explain_subject("S1", gcn_model, risk_models)

explain_subject <- function(subject_id, gcn_model, rf_models) {

  cat(sprintf("\n=== Explanation for subject: %s ===\n", subject_id))

  preds <- predict_individual_risk_gcn(gcn_model, subject_id)
  cat("GCN predicted risks:\n")
  print(preds)

  sal_list <- compute_gcn_saliency(gcn_model, subject_id, disease = NULL, n_top = 10)

  for (d in names(sal_list)) {
    cat(sprintf("\nTop features driving %s risk:\n", d))
    print(sal_list[[d]][, c("feature", "layer", "saliency", "baseline_risk")])
  }

  invisible(sal_list)
}

# Quick demo
explain_subject(sample_ids[1], gcn_model, risk_models)
