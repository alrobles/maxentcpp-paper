suppressMessages(library(maxentcpp))
set.seed(42)
n <- 2371
env <- data.frame(
  bio1 = runif(n, 10, 30),
  bio12 = runif(n, 500, 2500),
  landcover = factor(sample(c("forest","grass","urban","wetland","agri"), n, replace=TRUE))
)
occ <- sample(n, 100)

tm <- function(expr) { s <- system.time(expr); s[["elapsed"]]*1000 }

# 1. baseline fit
f1 <- maxent_generate_features(env[,1:2], types="lqh")
t0 <- tm({ fs <- maxent_featured_space(n, occ, f1); m <- maxent_fit(fs) })
cat(sprintf("baseline fit lqh: %.0f ms\n", t0))

# 2. with categorical
f2 <- maxent_generate_features(env, types="lqh", categorical="landcover")
t1 <- tm({ fs2 <- maxent_featured_space(n, occ, f2); m2 <- maxent_fit(fs2) })
cat(sprintf("fit WITH categorical: %.0f ms\n", t1))

# 3. jackknife (3 vars incl categorical)
t2 <- tm({ jk <- maxent_jackknife(env, occ, n, types=c("linear","quadratic","hinge"), categorical="landcover") })
cat(sprintf("jackknife (3 vars): %.0f ms\n", t2))

# 4. cross-validation k=5
t3 <- tm({ cv <- maxent_cross_validate(env, occ, n, k=5L, types=c("linear","quadratic","hinge"), categorical="landcover") })
cat(sprintf("cross-validate k=5: %.0f ms\n", t3))

# 5. replicates bootstrap n=5
t4 <- tm({ rp <- maxent_replicate(env, occ, n, n_replicates=5L, replicate_type="bootstrap", types=c("linear","quadratic","hinge"), categorical="landcover") })
cat(sprintf("replicate bootstrap n=5: %.0f ms\n", t4))

cat("DONE\n")
