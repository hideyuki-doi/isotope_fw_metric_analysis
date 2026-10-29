# Cluster-bootstrap mediation (a x b) with sampling effort (log10 individuals)
# as a covariate in both the a-path and b-path models. Resamples whole studies
# (Source_DOI). Replicates run in parallel. Output: mediation_bootstrap_nsample.csv
suppressMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(parallel)})
N_BOOT <- as.integer(Sys.getenv("ISOFRESH_N_BOOT", "999")); N_CORES <- as.integer(Sys.getenv("ISOFRESH_CORES", "9")); SEED <- 42
RAW <- Sys.getenv("ISOFRESH_RAW", "data")
DATA <- Sys.getenv("ISOFRESH_DATA", ".")  # folder holding IsoFresh_Metrics_v2.RData
dir.create("manuscript_figs", showWarnings = FALSE)
load(file.path(DATA, "IsoFresh_Metrics_v2.RData"))
rd <- function(f) read_delim(file.path(RAW, f), delim = ";", show_col_types = FALSE)
ns <- bind_rows(rd("IsoFresh_Isotope_Fish.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N),
                rd("IsoFresh_Isotope_Other.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N)) %>%
  filter(!is.na(d13C), !is.na(d15N)) %>% group_by(FoodWeb_ID) %>%
  summarize(N_Individuals = ifelse(all(is.na(N_Sample)), NA_real_, sum(N_Sample, na.rm = TRUE)), .groups = "drop")
d <- final_df %>% left_join(ns, by = "FoodWeb_ID") %>% filter(Include_S10) %>%
  mutate(Log_Richness = log10(Richness), Log_N_Ind = log10(N_Individuals))
z <- function(x) as.numeric(scale(x))
METRICS <- c("Connectance", "z_NODF", "Mean_Strength", "Mean_Closeness", "Modularity", "z_Q")
PREDS <- c("Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality"); MEDS <- c("Log_Richness", "Prey_Niche_Area")

fit <- function(dat, y, x, extra) {
  tryCatch(suppressMessages(suppressWarnings(lmer(as.formula(paste(y, "~", paste(c(x, extra, "N_z", "(1 | Region/Source_DOI)"), collapse = " + "))), dat))), error = function(e) NULL)
}
cf <- function(m, t) { if (is.null(m)) return(NA_real_); s <- summary(m)$coefficients; if (t %in% rownames(s)) s[t, "Estimate"] else NA_real_ }
ab <- function(dat) cf(fit(dat, "med_z", "p_z", NULL), "p_z") * cf(fit(dat, "m_z", "p_z", "med_z"), "med_z")

boot_one <- function(metric, pred, med) {
  dd <- d %>% filter(!is.na(.data[[metric]]), !is.na(.data[[pred]]), !is.na(.data[[med]]), !is.na(Log_N_Ind), !is.na(Region), !is.na(Source_DOI))
  if (nrow(dd) < 30 || n_distinct(dd$Source_DOI) < 10) return(NULL)
  dd$m_z <- z(dd[[metric]]); dd$p_z <- z(dd[[pred]]); dd$med_z <- z(dd[[med]]); dd$N_z <- z(dd$Log_N_Ind)
  pa <- cf(fit(dd, "med_z", "p_z", NULL), "p_z"); pb <- cf(fit(dd, "m_z", "p_z", "med_z"), "med_z")
  ids <- unique(dd$Source_DOI); by_study <- split(dd, dd$Source_DOI)
  bv <- unlist(mclapply(seq_len(N_BOOT), function(i) {
    set.seed(SEED + i); s <- sample(ids, length(ids), replace = TRUE)
    bd <- bind_rows(lapply(seq_along(s), function(j) { r <- by_study[[s[j]]]; r$Source_DOI <- paste0(s[j], "__rep", j); r }))
    ab(bd)
  }, mc.cores = N_CORES))
  bv <- bv[!is.na(bv)]; if (length(bv) < N_BOOT * 0.5) return(NULL)
  ci <- quantile(bv, c(0.025, 0.975))
  tibble(metric = metric, predictor = pred, mediator = med, n = nrow(dd), path_a = pa, path_b = pb, indirect_effect = pa * pb,
         ci_lwr = ci[[1]], ci_upr = ci[[2]], n_boot_ok = length(bv), sig = (ci[[1]] > 0 & ci[[2]] > 0) | (ci[[1]] < 0 & ci[[2]] < 0))
}
res <- list(); t0 <- Sys.time()
for (m in METRICS) for (p in PREDS) for (md in MEDS) {
  r <- tryCatch(boot_one(m, p, md), error = function(e) NULL)
  if (!is.null(r)) res[[length(res) + 1]] <- r
  cat(sprintf("done %s ~ %s via %s (%.1f min)\n", m, p, md, as.numeric(Sys.time() - t0, units = "mins")))
  write_csv(bind_rows(res), if (N_BOOT < 100) "mediation_bootstrap_nsample_TEST.csv" else "mediation_bootstrap_nsample.csv")
}
cat("ALL DONE\n")
