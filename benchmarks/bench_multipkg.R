#!/usr/bin/env Rscript
# Multi-package benchmark: maxentcpp vs dismo vs predicts vs maxnet
# on the bundled Abeillia abeillei dataset (73 presences, 2,371 bg cells).
# Fresh-state per fit, warm-up excluded. Measured on Gamma 12 Aug 2026.
suppressMessages(library(maxentcpp))
suppressMessages(library(terra))
suppressMessages(library(maxnet))
suppressMessages(library(dismo))
suppressMessages(library(predicts))
suppressMessages(library(rJava))

.jinit()
cat("Java:", .jcall("java/lang/System", "S", "getProperty", "java.version"), "\n")

occ_df <- maxentcpp::example_occ_df
s <- readRDS(system.file("extdata", "stack_1_12_crop.rds", package = "maxentcpp"))
r <- terra::unwrap(s)
g1 <- maxent_grid_from_terra(r[[1]], name = "bio1")
g2 <- maxent_grid_from_terra(r[[2]], name = "bio12")
feature_names <- c("bio1", "bio12")

# ---- Prepare data (same pipeline as maxent_run) ----
info <- maxent_grid_info(g1)
dim_obj <- maxent_dimension(info$nrows, info$ncols, info$xll, info$yll, info$cellsize)
occ <- maxent_read_occurrences(occ_df, dim_obj, lon_col = "long", lat_col = "lat")
bg <- maxent_background_indices(g1, n = 10000L, seed = 42L)
all_rows <- c(bg$rows, occ$rows); all_cols <- c(bg$cols, occ$cols)
n_total <- length(all_rows)
sample_indices <- seq(length(bg$rows), n_total - 1L)
env_vals <- lapply(list(g1, g2), function(g) grid_get_values_batch(g, as.integer(all_rows), as.integer(all_cols)))
names(env_vals) <- feature_names
valid_idx <- maxent_complete_cases(env_vals, nodata_value = -9999)
env_vals <- lapply(env_vals, function(v) v[valid_idx])
old_to_new <- rep(NA_integer_, n_total); old_to_new[valid_idx] <- seq_along(valid_idx) - 1L
sample_indices <- stats::na.omit(old_to_new[sample_indices + 1L])
n_total <- length(valid_idx)

# ---- maxentcpp warm-up ----
features <- maxent_generate_features(env_vals, types = c("linear","quadratic","hinge"))
fs <- maxent_featured_space(n_total, as.integer(sample_indices), features)
invisible(maxent_fit(fs, max_iter = 500L, convergence = 1e-5))

# ---- 1. maxentcpp (fresh features each run) ----
cat("=== maxentcpp Sequential (30 fresh runs) ===\n")
times_cpp <- numeric(30); iters_cpp <- numeric(30); auc_cpp <- numeric(30)
for (i in 1:30) {
    features <- maxent_generate_features(env_vals, types = c("linear","quadratic","hinge"))
    fs <- maxent_featured_space(n_total, as.integer(sample_indices), features)
    t0 <- system.time(fr <- maxent_fit(fs, max_iter = 500L, convergence = 1e-5))[["elapsed"]]
    times_cpp[i] <- t0; iters_cpp[i] <- fr$iterations
    # training AUC from the model
    fmat <- do.call(cbind, lapply(features, function(f) vapply(seq_len(n_total), function(i) maxent_feature_eval(f, i), numeric(1))))
    preds <- maxent_predict_model(fs, fmat)
    lab <- rep(0, n_total); lab[sample_indices + 1L] <- 1
    rk <- rank(preds); n1 <- sum(lab==1); n0 <- sum(lab==0)
    auc_cpp[i] <- (sum(rk[lab==1]) - n1*(n1+1)/2) / (n1*n0)
}
cat(sprintf("  median: %.1f ms | iters: %d | AUC: %.4f\n",
    median(times_cpp)*1000, median(iters_cpp), median(auc_cpp)))

# ---- 2. maxnet ----
m1 <- as.matrix(r[[1]]); m2 <- as.matrix(r[[2]])
cells <- which(!is.na(m1) & !is.na(m2))
env_all <- data.frame(bio1 = m1[cells], bio12 = m2[cells])
# map the 73 occurrence rows/cols (1-based grid indices) to raster cells,
# then to positions within the no-NA cell frame used by maxnet
occ_cell <- terra::cellFromRowCol(r, occ$rows, occ$cols)
pa <- rep(0, nrow(env_all)); pa[match(occ_cell, cells)] <- 1
cat("presences in maxnet frame:", sum(pa), "\n")
cat("=== maxnet (10 runs) ===\n")
times_mn <- numeric(10)
for (i in 1:10) {
    t0 <- system.time({ mn <- maxnet(pa, env_all, maxnet.formula(pa, env_all, classes = "lqh")) })[["elapsed"]]
    times_mn[i] <- t0
}
cat(sprintf("  median: %.1f ms\n", median(times_mn)*1000))

# ---- 3. dismo (Java maxent.jar, fresh jar invocation per run) ----
cat("=== dismo maxent (10 runs) ===\n")
times_dismo <- numeric(10)
for (i in 1:10) {
    t0 <- system.time({
        dm <- dismo::maxent(env_all, pa, args = c("maximumiterations=500",
                           "linear=true", "quadratic=true", "product=false",
                           "threshold=false", "hinge=true", "betamultiplier=1"))
    })[["elapsed"]]
    times_dismo[i] <- t0
}
cat(sprintf("  median: %.1f ms\n", median(times_dismo)*1000))

# ---- 4. predicts (Java maxent.jar via rJava, S4 MaxEnt) ----
cat("=== predicts MaxEnt (10 runs) ===\n")
times_pred <- numeric(10)
for (i in 1:10) {
    t0 <- system.time({
        pm <- predicts::MaxEnt(env_all, pa,
                               args = c("maximumiterations=500",
                                        "linear=true", "quadratic=true", "product=false",
                                        "threshold=false", "hinge=true", "betamultiplier=1"))
    })[["elapsed"]]
    times_pred[i] <- t0
}
cat(sprintf("  median: %.1f ms\n", median(times_pred)*1000))

# ---- Summary ----
cat("\n=== SUMMARY (median) ===\n")
res <- data.frame(
    package = c("maxentcpp", "maxnet", "dismo", "predicts"),
    median_ms = c(median(times_cpp)*1000, median(times_mn)*1000,
                  median(times_dismo)*1000, median(times_pred)*1000),
    auc = c(median(auc_cpp), NA, NA, NA),
    iters = c(median(iters_cpp), NA, NA, NA)
)
print(res)
cat(sprintf("\nSpeedup maxentcpp vs maxnet: %.0fx | vs dismo: %.0fx | vs predicts: %.0fx\n",
    res$median_ms[2]/res$median_ms[1], res$median_ms[3]/res$median_ms[1], res$median_ms[4]/res$median_ms[1]))
saveRDS(list(times_cpp=times_cpp, times_mn=times_mn, times_dismo=times_dismo,
             times_pred=times_pred, auc_cpp=auc_cpp, iters_cpp=iters_cpp),
        "/tmp/multi_pkg_bench.rds")
cat("SAVED /tmp/multi_pkg_bench.rds\n")
cat("DONE\n")
