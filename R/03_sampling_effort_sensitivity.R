# Sensitivity: does sampling effort (number of individuals isotope-measured per
# food web) change the climate effects?  Reads N_Sample from the raw IsoFresh
# CSVs; reuses final_df from IsoFresh_Metrics_v2.RData.
suppressMessages({library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest); library(ggplot2); library(patchwork)})
RAW <- Sys.getenv("ISOFRESH_RAW", "data")
DATA <- Sys.getenv("ISOFRESH_DATA", ".")  # folder holding IsoFresh_Metrics_v2.RData
dir.create("manuscript_figs", showWarnings = FALSE)
load(file.path(DATA, "IsoFresh_Metrics_v2.RData"))
rd <- function(f) read_delim(file.path(RAW, f), delim = ";", show_col_types = FALSE)
ns <- bind_rows(rd("IsoFresh_Isotope_Fish.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N),
                rd("IsoFresh_Isotope_Other.csv") %>% select(FoodWeb_ID, N_Sample, d13C, d15N)) %>%
  filter(!is.na(d13C), !is.na(d15N)) %>%
  group_by(FoodWeb_ID) %>%
  summarize(N_Individuals = ifelse(all(is.na(N_Sample)), NA_real_, sum(N_Sample, na.rm = TRUE)),
            N_NA_frac = mean(is.na(N_Sample)), N_Records = n(), .groups = "drop")
d <- final_df %>% left_join(ns, by = "FoodWeb_ID") %>% filter(Include_S10) %>%
  mutate(Log_N_Ind = log10(N_Individuals), Log_Richness = log10(Richness))
cat("\nN_Individuals summary (S>=10 webs):\n"); print(summary(d$N_Individuals))
cat("share of webs with missing N_Sample:", mean(is.na(d$N_Individuals)), "\n")
cat("cor(Log_N_Ind, Log_Richness)  =", round(cor(d$Log_N_Ind, d$Log_Richness, use = "complete"), 2), "\n")
cat("cor(Log_N_Ind, Abs_Latitude)  =", round(cor(d$Log_N_Ind, d$Abs_Latitude, use = "complete"), 2), "\n")
cat("cor(Log_N_Ind, Connectance)   =", round(cor(d$Log_N_Ind, d$Connectance, use = "complete"), 2), "\n")

z <- function(x) as.numeric(scale(x))
METRICS <- c("Connectance", "z_NODF", "Mean_Strength", "Mean_Closeness", "Modularity", "z_Q")
DRIVERS <- c(Abs_Latitude = "Latitude", Mean_Annual_Temp = "Temperature", Temp_Seasonality = "Seasonality")
fit1 <- function(dat, y, x, adj) {
  dd <- dat %>% filter(!is.na(.data[[y]]), !is.na(.data[[x]]), !is.na(Log_N_Ind), !is.na(Region), !is.na(Source_DOI))
  if (nrow(dd) < 30) return(NULL)
  dd$yz <- z(dd[[y]]); dd$xz <- z(dd[[x]]); dd$nz <- z(dd$Log_N_Ind)
  f <- if (adj) yz ~ xz + nz + (1 | Region/Source_DOI) else yz ~ xz + (1 | Region/Source_DOI)
  m <- tryCatch(suppressMessages(suppressWarnings(lmer(f, dd))), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  s <- summary(m)$coefficients
  tibble(est = s["xz", "Estimate"], se = s["xz", "Std. Error"], p = s["xz", "Pr(>|t|)"], n = nrow(dd))
}
rows <- list()
for (eco in c("All", "Lentic", "Lotic")) {
  dat <- if (eco == "All") d else d %>% filter(Ecosystem_Grp == eco)
  for (y in METRICS) for (x in names(DRIVERS)) for (adj in c(FALSE, TRUE)) {
    r <- fit1(dat, y, x, adj)
    if (!is.null(r)) rows[[length(rows) + 1]] <- bind_cols(tibble(ecosystem = eco, metric = y, driver = DRIVERS[[x]], model = ifelse(adj, "+ log N samples", "base")), r)
  }
}
res <- bind_rows(rows) %>% group_by(model) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
write_csv(res, "nsample_sensitivity.csv")

# subset: only webs with >= median individuals measured
med <- median(d$N_Individuals, na.rm = TRUE)
rows2 <- list()
for (y in METRICS) for (x in names(DRIVERS)) {
  r <- fit1(d %>% filter(N_Individuals >= med), y, x, FALSE)
  if (!is.null(r)) rows2[[length(rows2) + 1]] <- bind_cols(tibble(ecosystem = "All", metric = y, driver = DRIVERS[[x]], model = paste0("only N >= ", med)), r)
}
res2 <- bind_rows(rows2) %>% mutate(q = p.adjust(p, "BH"))
write_csv(res2, "nsample_sensitivity_subset.csv")

cmp <- res %>% filter(ecosystem == "All") %>% select(metric, driver, model, est, q) %>%
  pivot_wider(names_from = model, values_from = c(est, q)) %>% mutate(across(where(is.numeric), ~round(.x, 3)))
cat("\n=== Pooled ecosystems: base vs +log N individuals ===\n"); print(as.data.frame(cmp))
cat("\n=== Subset N >= median (", med, ") ===\n"); print(as.data.frame(res2 %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

# figure
LAB <- c(Connectance = "Connectance", z_NODF = "Nestedness (z-NODF)", Mean_Strength = "Mean link strength", Mean_Closeness = "Mean closeness", Modularity = "Modularity (Q)", z_Q = "Standardised modularity (z-Q)")
pl <- bind_rows(res %>% filter(ecosystem == "All"), res2) %>%
  mutate(sgn = ifelse(driver == "Temperature", -1, 1), est = est * sgn, lo = est - 1.96 * se, hi = est + 1.96 * se,
         metric_f = factor(LAB[metric], levels = rev(LAB)), driver = factor(driver, levels = c("Latitude", "Temperature", "Seasonality")),
         model = factor(model, levels = c("base", "+ log N samples", paste0("only N >= ", med))), sig = q < 0.05)
fig <- ggplot(pl, aes(est, metric_f, colour = model, shape = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, position = position_dodge(0.7)) +
  geom_point(size = 2, fill = "white", position = position_dodge(0.7)) +
  scale_colour_manual(values = c("grey20", "#2a78d6", "#eb6834"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = "BH-FDR") +
  scale_x_continuous(breaks = c(-0.4, 0, 0.4)) + facet_wrap(~driver) +
  labs(x = "Standardised effect of moving poleward (95% CI)", y = NULL, title = "Sampling effort (no. of samples measured) does not change the climate effects") +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
dir.create("manuscript_figs", showWarnings = FALSE)
ggsave("manuscript_figs/FigS5_sampling_effort.png", fig, width = 180, height = 90, units = "mm", dpi = 300, bg = "white")
ggsave("manuscript_figs/FigS5_sampling_effort.pdf", fig, width = 180, height = 90, units = "mm")
