# Theory-structured piecewise SEM for connectance: climate -> richness (diversity
# pathway), -> prey niche area (resource-spread pathway), -> fish trophic height
# (generalism/omnivory pathway, Petchey diet breadth), and direct. Sampling effort
# in every equation. Prey area -> trophic height path added (d-sep test showed it was missing);
# paths running through both prey area and trophic height are counted under prey area. Uses theory_results/theory_analysis_df.rds from theory_tests.R
suppressMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(piecewiseSEM)})
d <- readRDS("theory_results/theory_analysis_df.rds"); z <- function(x) as.numeric(scale(x))
DRV <- c(Latitude = "Abs_Latitude", Temperature = "Mean_Annual_Temp", Seasonality = "Temp_Seasonality")
rows <- list(); coefs <- list()
for (nm in names(DRV)) {
  col <- DRV[[nm]]
  dat <- d %>% filter(!is.na(Connectance), !is.na(.data[[col]]), !is.na(Log_S), !is.na(Prey_Niche_Area), !is.na(Fish_Trophic_Height),
                      !is.na(Log_N_Ind), !is.na(Region), !is.na(Source_DOI)) %>%
    mutate(X = z(.data[[col]]), Rich = z(Log_S), Prey = z(Prey_Niche_Area), TH = z(Fish_Trophic_Height), Conn = z(Connectance), N = z(Log_N_Ind))
  f <- function(s) as.formula(paste(s, "+ N + (1 | Region/Source_DOI)"))
  ms <- list(lmer(f("Rich ~ X"), dat), lmer(f("Prey ~ X + Rich"), dat), lmer(f("TH ~ X + Rich + Prey"), dat), lmer(f("Conn ~ X + Rich + Prey + TH"), dat))
  sm <- summary(do.call(psem, c(ms, list(data = as.data.frame(dat)))))
  cf <- sm$coefficients; names(cf)[names(cf) == ""] <- "Sig"; coefs[[nm]] <- cf %>% mutate(driver = nm)
  g <- function(r, p, k = "Std.Estimate") { x <- cf %>% filter(Response == r, Predictor == p); if (nrow(x)) x[[k]][1] else NA_real_ }
  a <- g("Rich", "X"); bR <- g("Conn", "Rich"); aP <- g("Prey", "X"); rP <- g("Prey", "Rich"); bP <- g("Conn", "Prey")
  aT <- g("TH", "X"); rT <- g("TH", "Rich"); pT <- g("TH", "Prey"); bT <- g("Conn", "TH"); dir <- g("Conn", "X")
  rows[[nm]] <- tibble(driver = nm, n = nrow(dat), direct = dir, direct_p = g("Conn", "X", "P.Value"),
    via_richness = a * bR, via_prey_area = (aP + a * rP) * (bP + pT * bT), via_fish_trophic_height = (aT + a * rT) * bT,
    climate_to_TH = aT, climate_to_TH_p = g("TH", "X", "P.Value"), TH_to_conn = bT, TH_to_conn_p = g("Conn", "TH", "P.Value"),
    fisherC = sm$Cstat$Fisher.C, fisherC_p = sm$Cstat$P.Value) %>%
    mutate(total = direct + via_richness + via_prey_area + via_fish_trophic_height)
}
res <- bind_rows(rows); write_csv(res, "theory_results/T6_theory_sem_connectance.csv"); write_csv(bind_rows(coefs), "theory_results/T6_theory_sem_coefficients.csv")
print(as.data.frame(res %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
