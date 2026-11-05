# ------------------------------------------------------------------
# Reproduction script for:
#   An automated CT-based radiomics pipeline for functional outcome
#   prediction in intracerebral hematoma: a feasibility study
#
# Reproduces: ICC filtering, correlation filtering, LASSO feature
# selection, Rad-score construction, hierarchical logistic models
# (Models 1-3), DeLong comparisons, and bootstrap optimism correction.
#
# Input : ICH_radiomics_deidentified.xlsx (same folder)
# Tested: R 4.5.1, readxl, irr, glmnet 4.1-10, pROC 1.19.0.1
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readxl); library(irr); library(glmnet); library(pROC)
})

set.seed(3407)
f <- "ICH_radiomics_deidentified.xlsx"

clin <- read_excel(f, "clinical")
gt   <- read_excel(f, "ground_truth")
d1   <- read_excel(f, "doctor1")
d2   <- read_excel(f, "doctor2")
d3   <- read_excel(f, "doctor3")
dl   <- read_excel(f, "DL_predict")

stopifnot(identical(clin$ID, gt$ID), identical(gt$ID, dl$ID),
          identical(gt$ID, d1$ID), identical(gt$ID, d2$ID), identical(gt$ID, d3$ID))
y     <- clin$outcome_poor_mRS4to6
feats <- setdiff(names(gt), "ID")

## 1. Inter-rater reproducibility: ICC(2,1) across three neurosurgeons ----
icc <- sapply(feats, function(v) {
  m <- cbind(d1[[v]], d2[[v]], d3[[v]])
  irr::icc(m, model = "twoway", type = "agreement", unit = "single")$value
})
keep_icc <- feats[icc > 0.75]
cat("Features with ICC > 0.75:", length(keep_icc), "of", length(feats), "\n")   # 78

## 2. Redundancy filtering: remove any feature with R^2 > 0.9 vs another ----
X   <- as.matrix(gt[, keep_icc])
R2  <- cor(X)^2; diag(R2) <- 0
keep <- keep_icc[apply(R2, 1, max) <= 0.9]
cat("Features after R^2 filtering:", length(keep), "\n")                      # 23

## 3. LASSO logistic regression (10-fold CV, lambda at minimum deviance) ----
Xk <- as.matrix(gt[, keep])
cvfit <- cv.glmnet(Xk, y, family = "binomial", alpha = 1,
                   standardize = TRUE, nfolds = 10)
co  <- coef(cvfit, s = "lambda.min")
sel <- rownames(co)[-1][as.vector(co[-1]) != 0]
cat("lambda.min =", round(cvfit$lambda.min, 3), "\nSelected features:\n"); print(co[c("(Intercept)", sel), , drop = FALSE])

## 4. Rad-score = 100 * (b0 + sum(bi * xi)); same coefficients applied to DL features ----
rad <- function(df) as.vector(100 * (co["(Intercept)", 1] + as.matrix(df[, sel]) %*% co[sel, 1]))
clin$Rad_GT <- rad(gt)
clin$Rad_DL <- rad(dl)

roc_gt <- roc(y, clin$Rad_GT, quiet = TRUE)
roc_dl <- roc(y, clin$Rad_DL, quiet = TRUE)
cat(sprintf("\nRad-score AUC: physician-defined %.3f, DL-based %.3f; DeLong p = %.3f\n",
            auc(roc_gt), auc(roc_dl), roc.test(roc_gt, roc_dl)$p.value))

## 5. Hierarchical models ----
m1 <- glm(y ~ age + sex_male + BMI + GCS, data = clin, family = binomial)
m2 <- update(m1, . ~ . + DL_volume_mL)
m3 <- update(m2, . ~ . + Rad_DL)

perf <- function(m, lab) {
  p <- fitted(m); r <- roc(y, p, quiet = TRUE); ci <- ci.auc(r)
  pred <- as.integer(p >= 0.5); tp <- sum(pred & y); tn <- sum(!pred & !y)
  fp <- sum(pred & !y); fn <- sum(!pred & y)
  data.frame(Model = lab, AUC = auc(r), CI_low = ci[1], CI_high = ci[3],
             ACC = (tp + tn) / length(y), F1 = 2 * tp / (2 * tp + fp + fn),
             SEN = tp / (tp + fn), SPE = tn / (tn + fp),
             PPV = tp / (tp + fp), NPV = tn / (tn + fn))
}
tab3 <- rbind(perf(m1, "Model 1"), perf(m2, "Model 2"), perf(m3, "Model 3"))
print(tab3, digits = 3)                                                       # Table 3

r1 <- roc(y, fitted(m1), quiet = TRUE); r2 <- roc(y, fitted(m2), quiet = TRUE); r3 <- roc(y, fitted(m3), quiet = TRUE)
cat(sprintf("DeLong: M2 vs M1 p = %.3f; M3 vs M1 p = %.3f; M3 vs M2 p = %.3f\n",
            roc.test(r2, r1)$p.value, roc.test(r3, r1)$p.value, roc.test(r3, r2)$p.value))

cat("\nModel 3 coefficients (logit scale) and odds ratios:\n")             # Suppl. Table S4
print(cbind(beta = coef(m3), OR = exp(coef(m3)), exp(confint.default(m3)),
            p = summary(m3)$coefficients[, 4]), digits = 3)

## 6. Bootstrap optimism correction of the full pipeline (B = 1000) ----
B <- 1000
pipeline <- function(idx) {
  g <- gt[idx, ]; dd <- dl[idx, ]; cc <- clin[idx, ]; yy <- y[idx]
  a <- cbind(d1[idx, ], d2[idx, ], d3[idx, ])
  ic <- sapply(feats, function(v) irr::icc(cbind(d1[[v]][idx], d2[[v]][idx], d3[[v]][idx]),
                                           model = "twoway", type = "agreement", unit = "single")$value)
  k1 <- feats[ic > 0.75]; if (length(k1) < 2) return(NULL)
  r2 <- cor(as.matrix(g[, k1]))^2; diag(r2) <- 0
  k2 <- k1[apply(r2, 1, max) <= 0.9]; if (length(k2) < 2) return(NULL)
  cv <- cv.glmnet(as.matrix(g[, k2]), yy, family = "binomial", alpha = 1, standardize = TRUE, nfolds = 10)
  cb <- coef(cv, s = "lambda.min"); s <- rownames(cb)[-1][as.vector(cb[-1]) != 0]
  if (length(s) == 0) return(NULL)
  rs <- function(df) as.vector(100 * (cb["(Intercept)", 1] + as.matrix(df[, s]) %*% cb[s, 1]))
  cc$Rad_DL <- rs(dd)
  list(m3 = glm(yy ~ age + sex_male + BMI + GCS + DL_volume_mL + Rad_DL, data = cc, family = binomial),
       rs = rs, sel = s)
}
app  <- auc(r3)
opt  <- numeric(0); selfreq <- setNames(integer(length(feats)), feats)
for (b in seq_len(B)) {
  idx <- sample(nrow(clin), replace = TRUE)
  fit <- tryCatch(pipeline(idx), error = function(e) NULL); if (is.null(fit)) next
  selfreq[fit$sel] <- selfreq[fit$sel] + 1
  clin_t <- clin; clin_t$Rad_DL <- fit$rs(dl)
  auc_boot <- auc(roc(y[idx], fitted(fit$m3), quiet = TRUE))
  auc_test <- auc(roc(y, predict(fit$m3, newdata = clin_t, type = "response"), quiet = TRUE))
  opt <- c(opt, auc_boot - auc_test)
}
cat(sprintf("\nModel 3: apparent AUC %.3f, optimism %.3f, corrected AUC %.3f (%d valid resamples)\n",
            app, mean(opt), app - mean(opt), length(opt)))                    # Suppl. Table S5
cat("Selection frequency (%) of the three Rad-score features:\n")
print(round(100 * selfreq[sel] / length(opt), 1))                              # Suppl. Table S6
