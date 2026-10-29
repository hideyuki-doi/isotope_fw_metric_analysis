# piecewiseSEM re-run with sampling effort (log10 total individuals measured)
# added as a covariate to EVERY equation. Base and adjusted models are fitted
# on the identical sample so differences are due to the covariate only.
suppressMessages({library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest); library(piecewiseSEM); library(ggplot2)})
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
DRV <- list(Latitude = "Abs_Latitude", Temperature = "Mean_Annual_Temp", Seasonality = "Temp_Seasonality")

run_sem <- function(exog_col, adjust, dat) {
  nn <- if (adjust) " + N_z" else ""
  f <- function(y, rhs) as.formula(paste0(y, " ~ ", rhs, nn, " + (1 | Region/Source_DOI)"))
  ms <- list(lmer(f("Richness_z", "Exog_z"), dat), lmer(f("PreyArea_z", "Exog_z + Richness_z"), dat),
             lmer(f("Connectance_z", "Exog_z + Richness_z + PreyArea_z"), dat),
             lmer(f("zNODF_z", "Exog_z + Richness_z + PreyArea_z"), dat),
             lmer(f("Strength_z", "Exog_z + Richness_z + PreyArea_z"), dat))
  sm <- summary(do.call(psem, c(ms, list(data = as.data.frame(dat)))))
  cf <- sm$coefficients; names(cf)[names(cf) == ""] <- "Sig"
  g <- function(r, p) { x <- cf %>% filter(Response == r, Predictor == p); if (nrow(x)) x$Std.Estimate[1] else NA_real_ }
  gp <- function(r, p) { x <- cf %>% filter(Response == r, Predictor == p); if (nrow(x)) x$P.Value[1] else NA_real_ }
  tr <- lapply(c("Connectance_z", "zNODF_z", "Strength_z"), function(o) {
    a <- g("Richness_z", "Exog_z"); ap <- g("PreyArea_z", "Exog_z"); ar <- g("PreyArea_z", "Richness_z")
    br <- g(o, "Richness_z"); bp <- g(o, "PreyArea_z"); dir <- g(o, "Exog_z")
    tibble(outcome = o, direct = dir, direct_p = gp(o, "Exog_z"), via_richness = a * br, via_prey = ap * bp + a * ar * bp,
           total = dir + a * br + ap * bp + a * ar * bp,
           effort_on_outcome = if (adjust) g(o, "N_z") else NA_real_)
  })
  list(effects = bind_rows(tr), coefs = cf, effort_on_richness = if (adjust) g("Richness_z", "N_z") else NA_real_,
       effort_on_prey = if (adjust) g("PreyArea_z", "N_z") else NA_real_, n = nrow(dat))
}

out <- list(); coefs_out <- list()
for (nm in names(DRV)) {
  col <- DRV[[nm]]
  dat <- d %>% filter(!is.na(Connectance), !is.na(z_NODF), !is.na(Mean_Strength), !is.na(.data[[col]]), !is.na(Log_Richness),
                      !is.na(Prey_Niche_Area), !is.na(Log_N_Ind), !is.na(Region), !is.na(Source_DOI)) %>%
    mutate(Exog_z = z(.data[[col]]), Richness_z = z(Log_Richness), PreyArea_z = z(Prey_Niche_Area), Connectance_z = z(Connectance),
           zNODF_z = z(z_NODF), Strength_z = z(Mean_Strength), N_z = z(Log_N_Ind))
  cat(sprintf("\n[%s] n = %d (identical sample for both models)\n", nm, nrow(dat)))
  for (adj in c(FALSE, TRUE)) {
    r <- run_sem(col, adj, dat)
    out[[length(out) + 1]] <- r$effects %>% mutate(driver = nm, model = ifelse(adj, "+ sampling effort", "base"), n = r$n,
                                                   effort_on_richness = r$effort_on_richness, effort_on_prey = r$effort_on_prey)
    coefs_out[[length(coefs_out) + 1]] <- r$coefs %>% mutate(driver = nm, model = ifelse(adj, "+ sampling effort", "base"))
  }
}
eff <- bind_rows(out); write_csv(eff, "sem_nsample_traced_effects.csv"); write_csv(bind_rows(coefs_out), "sem_nsample_coefficients.csv")
print(as.data.frame(eff %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

# figure: direct / via richness / via prey for base vs adjusted
LAB <- c(Connectance_z = "Connectance", zNODF_z = "Nestedness (z-NODF)", Strength_z = "Mean link strength")
pl <- eff %>% mutate(sgn = ifelse(driver == "Temperature", -1, 1)) %>%
  pivot_longer(c(direct, via_richness, via_prey), names_to = "path", values_to = "v") %>%
  mutate(v = v * sgn, outcome = factor(LAB[outcome], levels = LAB), driver = factor(driver, levels = names(DRV)),
         path = factor(path, levels = c("direct", "via_richness", "via_prey"), labels = c("Direct", "Via taxon richness", "Via prey niche area")),
         model = factor(model, levels = c("base", "+ sampling effort")))
fig <- ggplot(pl, aes(v, outcome, fill = model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("grey60", "#2a78d6"), name = NULL) +
  facet_grid(path ~ driver) +
  labs(x = "Standardised path effect of moving poleward", y = NULL,
       title = "SEM paths with and without sampling effort as covariate") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
ggsave("manuscript_figs/FigS6_sem_sampling_effort.png", fig, width = 180, height = 120, units = "mm", dpi = 300, bg = "white")
ggsave("manuscript_figs/FigS6_sem_sampling_effort.pdf", fig, width = 180, height = 120, units = "mm")
