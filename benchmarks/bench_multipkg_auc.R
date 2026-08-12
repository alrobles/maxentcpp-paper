#!/usr/bin/env Rscript
# Consistent AUC comparison: same auc_fun on all four packages.
suppressMessages(library(maxentcpp))
suppressMessages(library(terra))
suppressMessages(library(maxnet))
suppressMessages(library(dismo))
suppressMessages(library(predicts))
suppressMessages(library(rJava))
.jinit()

occ_df <- maxentcpp::example_occ_df
s <- readRDS(system.file("extdata", "stack_1_12_crop.rds", package = "maxentcpp"))
r <- terra::unwrap(s)
g1 <- maxent_grid_from_terra(r[[1]], name = "bio1")
g2 <- maxent_grid_from_terra(r[[2]], name = "bio12")
info <- maxent_grid_info(g1)
dim_obj <- maxent_dimension(info$nrows, info$ncols, info$xll, info$yll, info$cellsize)
occ <- maxent_read_occurrences(occ_df, dim_obj, lon_col = "long", lat_col = "lat")
bg <- maxent_background_indices(g1, n = 10000L, seed = 42L)
all_rows <- c(bg$rows, occ$rows); all_cols <- c(bg$cols, occ$cols)
n_total <- length(all_rows)
sample_indices <- seq(length(bg$rows), n_total - 1L)
env_vals <- lapply(list(g1, g2), function(g) grid_get_values_batch(g, as.integer(all_rows), as.integer(all_cols)))
names(env_vals) <- feature_names <- c("bio1", "bio12")
valid_idx <- maxent_complete_cases(env_vals, nodata_value = -9999)
env_vals <- lapply(env_vals, function(v) v[valid_idx])
old_to_new <- rep(NA_integer_, n_total); old_to_new[valid_idx] <- seq_along(valid_idx) - 1L
sample_indices <- stats::na.omit(old_to_new[sample_indices + 1L])
n_total <- length(valid_idx)

m1 <- as.matrix(r[[1]]); m2 <- as.matrix(r[[2]])
cells <- which(!is.na(m1) & !is.na(m2))
env_all <- data.frame(bio1 = m1[cells], bio12 = m2[cells])
occ_cell <- terra::cellFromRowCol(r, occ$rows, occ$cols)
pa <- rep(0, nrow(env_all)); pa[match(occ_cell, cells)] <- 1

auc_fun <- function(pred, lab) {
    rk <- rank(pred); n1 <- sum(lab==1); n0 <- sum(lab==0)
    (sum(rk[lab==1]) - n1*(n1+1)/2) / (n1*n0)
}

# maxentcpp: predictions on the same 2371-cell frame for a fair AUC comparison
features <- maxent_generate_features(env_vals, types = c("linear","quadratic","hinge"))
fs <- maxent_featured_space(n_total, as.integer(sample_indices), features)
fr <- maxent_fit(fs, max_iter = 500L, convergence = 1e-5)
proj <- maxent_project_cloglog(fs, list(g1, g2), feature_names)
cpp_full <- as.numeric(t(maxent_grid_to_matrix(proj)))
cpp_on_cells <- cpp_full[cells]
cat(sprintf("maxentcpp: AUC=%.4f (same frame, cloglog)\n", auc_fun(cpp_on_cells, pa)))

mn <- maxnet(pa, env_all, maxnet.formula(pa, env_all, classes = "lqh"))
mn_pred <- predict(mn, env_all, type = "cloglog")
cat(sprintf("maxnet:   AUC=%.4f\n", auc_fun(mn_pred, pa)))

dm <- dismo::maxent(env_all, pa, args = c("maximumiterations=500","linear=true",
    "quadratic=true","product=false","threshold=false","hinge=true","betamultiplier=1"))
dm_pred <- dismo::predict(dm, env_all)
cat(sprintf("dismo:    AUC=%.4f (same-frame cloglog predict)\n", auc_fun(dm_pred, pa)))
cat(sprintf("dismo:    Training.AUC (self-report)=%.4f, iterations=%d\n",
    dm@results["Training.AUC", 1], dm@results["Iterations", 1]))

pm <- predicts::MaxEnt(env_all, pa, args = c("maximumiterations=500","linear=true",
    "quadratic=true","product=false","threshold=false","hinge=true","betamultiplier=1"))
pm_pred <- predicts::predict(pm, env_all)
cat(sprintf("predicts: AUC=%.4f (same-frame cloglog predict)\n", auc_fun(pm_pred, pa)))
cat(sprintf("predicts: Training.AUC (self-report)=%.4f, iterations=%d\n",
    pm@results["Training.AUC", 1], pm@results["Iterations", 1]))
cat("DONE\n")
