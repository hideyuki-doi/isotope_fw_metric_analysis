# piecewise SEM with null-corrected modularity (z_Q) and raw Q as extra outcomes,
# sampling effort as covariate in every equation. Tests whether the direct
# climate -> z_Q path offsets the richness-mediated (negative) path.
suppressMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(piecewiseSEM)})
RAW <- Sys.getenv("ISOFRESH_RAW", "data"); DATA <- Sys.getenv("ISOFRESH_DATA", ".")
load(file.path(DATA, "IsoFresh_Metrics_v2.RData"))
rd <- function(f) read_delim(file.path(RAW, f), delim = ";", show_col_types = FALSE)
ns <- bind_rows(rd("IsoFresh_Isotope_Fish.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N), rd("IsoFresh_Isotope_Other.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N)) %>%
  filter(!is.na(d13C), !is.na(d15N)) %>% group_by(FoodWeb_ID) %>% summarize(N_Individuals = ifelse(all(is.na(N_Sample)), NA_real_, sum(N_Sample, na.rm = TRUE)), .groups = "drop")
d <- final_df %>% left_join(ns, by = "FoodWeb_ID") %>% filter(Include_S10) %>% mutate(Log_Richness = log10(Richness), Log_N_Ind = log10(N_Individuals))
z <- function(x) as.numeric(scale(x))
DRV <- list(Latitude = "Abs_Latitude", Temperature = "Mean_Annual_Temp", Seasonality = "Temp_Seasonality")
rows <- list()
for (nm in names(DRV)) for (adj in c(FALSE, TRUE)) {
  col <- DRV[[nm]]
  dat <- d %>% filter(!is.na(Connectance), !is.na(z_NODF), !is.na(Mean_Strength), !is.na(z_Q), !is.na(Modularity), !is.na(.data[[col]]), !is.na(Log_Richness),
                      !is.na(Prey_Niche_Area), !is.na(Log_N_Ind), !is.na(Region), !is.na(Source_DOI)) %>%
    mutate(Exog_z = z(.data[[col]]), Richness_z = z(Log_Richness), PreyArea_z = z(Prey_Niche_Area), zQ_z = z(z_Q), Q_z = z(Modularity),
           Connectance_z = z(Connectance), zNODF_z = z(z_NODF), Strength_z = z(Mean_Strength), N_z = z(Log_N_Ind))
  nn <- if (adj) " + N_z" else ""; f <- function(y, rhs) as.formula(paste0(y, " ~ ", rhs, nn, " + (1 | Region/Source_DOI)"))
  ms <- list(lmer(f("Richness_z", "Exog_z"), dat), lmer(f("PreyArea_z", "Exog_z + Richness_z"), dat),
             lmer(f("zQ_z", "Exog_z + Richness_z + PreyArea_z"), dat), lmer(f("Q_z", "Exog_z + Richness_z + PreyArea_z"), dat))
  cf <- summary(do.call(psem, c(ms, list(data = as.data.frame(dat)))))$coefficients; names(cf)[names(cf) == ""] <- "Sig"
  g <- function(r, p, col = "Std.Estimate") { x <- cf %>% filter(Response == r, Predictor == p); if (nrow(x)) x[[col]][1] else NA_real_ }
  a <- g("Richness_z", "Exog_z"); ap <- g("PreyArea_z", "Exog_z"); ar <- g("PreyArea_z", "Richness_z")
  for (o in c("zQ_z", "Q_z")) {
    br <- g(o, "Richness_z"); bp <- g(o, "PreyArea_z"); dir <- g(o, "Exog_z")
    rows[[length(rows) + 1]] <- tibble(driver = nm, adjusted = adj, outcome = o, n = nrow(dat), direct = dir, direct_p = g(o, "Exog_z", "P.Value"),
      richness_to_outcome = br, richness_to_outcome_p = g(o, "Richness_z", "P.Value"), via_richness = a * br, via_prey = ap * bp + a * ar * bp,
      total = dir + a * br + ap * bp + a * ar * bp)
  }
}
res <- bind_rows(rows); write_csv(res, "sem_zq_traced_effects.csv")
print(as.data.frame(res %>% filter(adjusted) %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
cat("\n(unadjusted)\n"); print(as.data.frame(res %>% filter(!adjusted, outcome == "zQ_z") %>% mutate(across(where(is.numeric), ~round(.x, 3)))))
