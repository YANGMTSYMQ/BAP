library(neuroCombat)
library(preprocessCore)
library(ggplot2)
library(dplyr)

meta <- read.csv("Data/Meta_info.csv", row.names = 1)
exp  <- read.csv("Data/DIS_AGE_Combat_exp2", row.names = 1)

rownames(meta)<-meta$Sample
common_samples <- intersect(rownames(exp), meta$Sample)
exp  <- exp[common_samples, , drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
meta$Sample <- rownames(meta)
stopifnot(identical(rownames(exp), meta$Sample))
exp <- exp[, colSums(is.na(exp)) == 0, drop = FALSE]
exp_qn <- exp

for (ds in unique(meta$dataset)) {
  
  idx <- meta$dataset == ds
  
  x <- as.matrix(exp[idx, , drop = FALSE])
  
  x_qn <- preprocessCore::normalize.quantiles(t(x))
  rownames(x_qn) <- colnames(x)
  colnames(x_qn) <- rownames(x)
  
  exp_qn[idx, ] <- t(x_qn)
}
keep_gene <- apply(exp_qn, 2, sd) > 0
exp_qn <- exp_qn[, keep_gene, drop = FALSE]

meta$DIS_HEA <- factor(meta$DIS_HEA)
table(meta$Set)
internal_batch <- "Training set"
tuning_batch   <- "External tuning set"
idx_internal <- meta$Set == internal_batch
idx_tuning   <- meta$Set == tuning_batch
idx_external <- grepl("ETES", meta$Batch)
meta_internal <- meta[idx_internal, , drop = FALSE]
exp_internal  <- exp_qn[meta_internal$Sample, , drop = FALSE]
mod_internal <- model.matrix(
  ~ AGE,
  data = meta_internal
)
combat_internal <- neuroCombat(
  dat        = t(exp_internal),
  batch      = meta_internal$Batch,
  mod        = mod_internal,
  parametric = TRUE,
  verbose    = TRUE
)
exp_internal_cb <- t(combat_internal$dat.combat)
meta_tuning <- meta[idx_tuning, , drop = FALSE]
exp_tuning  <- exp_qn[meta_tuning$Sample, , drop = FALSE]
meta_external <- meta[idx_external, , drop = FALSE]
exp_external  <- exp_qn[meta_external$Sample, , drop = FALSE]
mod_external <- model.matrix(
  ~ AGE,
  data = meta_external
)

combat_external <- neuroCombat(
  dat        = t(exp_external),
  batch      = meta_external$Batch,
  mod        = mod_external,
  parametric = TRUE,
  verbose    = TRUE
)

exp_external_cb <- t(combat_external$dat.combat)


expr_combat <- rbind(
  exp_internal_cb,
  exp_tuning,
  exp_external_cb
)

expr_combat <- expr_combat[meta$Sample, , drop = FALSE]





stopifnot(identical(rownames(expr_combat), rownames(meta)))


pca_res <- prcomp(
  expr_combat,
  center = TRUE,
  scale. = TRUE
)

pca_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

pc1_var <- round(pca_var[1] * 100, 2)
pc2_var <- round(pca_var[2] * 100, 2)
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1    = pca_res$x[, 1],
  PC2    = pca_res$x[, 2],
  meta[rownames(pca_res$x), , drop = FALSE],
  check.names = FALSE
)

p_set <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Set)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA after Quantile Normalization and ComBat",
    subtitle = "Colored by Set",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Set"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_set)

p_batch <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Batch)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA after Quantile Normalization and ComBat",
    subtitle = "Colored by Batch",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Batch"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_batch)

p_dataset <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = dataset)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA after Quantile Normalization and ComBat",
    subtitle = "Colored by Dataset",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Dataset"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_dataset)


p_dis <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = DIS_HEA)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    title = "PCA after Quantile Normalization and ComBat",
    subtitle = "Colored by Disease Status",
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "DIS_HEA"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(p_dis)


