# ==========================================================================
# IsoFresh: manuscript figures (reads the CSV/RData outputs of
# network_analysis_v2_FULL_PIPELINE.R; no heavy computation here)
#
# Run from the folder that holds the pipeline outputs:
#   Rscript make_manuscript_figures.R
# Inputs : model_battery_results.csv, piecewiseSEM_traced_effects_*.csv,
#          piecewiseSEM_coefficients_latitude.csv, mediation_bootstrap.csv,
#          seasonality_interaction.csv, bodysize_temperature_interaction.csv,
#          IsoFresh_Metrics_v2.RData, IsoFresh_Metrics_v2_thresholdsweep.RData
# Outputs: manuscript_figs/Fig1..Fig4 and FigS1..FigS4 (PDF + PNG)
#
# Design rules
#  * Storyline order: pattern (Fig1) -> climate effects (Fig2) -> why: direct
#    vs diversity-mediated (Fig3, Fig4). Everything else is supplementary.
#  * All three climate drivers are expressed in ONE direction -- "toward the
#    poles": higher latitude, COOLER temperature (sign of temperature flipped),
#    higher seasonality -- so effects are directly comparable.
#  * Colours: Lentic blue / Lotic orange; diverging blue-grey-red for signed
#    effects. Filled symbol / star = BH-FDR q < 0.05.
# ==========================================================================
suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  library(lme4); library(lmerTest); library(maps); library(readr); library(scales)
})

OUT <- "manuscript_figs"; dir.create(OUT, showWarnings = FALSE)
COL_LENTIC <- "#2a78d6"; COL_LOTIC <- "#eb6834"; COL_ALL <- "grey20"
DIV_NEG <- "#e34948"; DIV_MID <- "#f4f3f0"; DIV_POS <- "#1c5cab"
COL_DIRECT <- "#eb6834"; COL_RICH <- "#2a78d6"; COL_PREY <- "#9cc4f1"

theme_ms <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.key.size = unit(3.5, "mm"),
        plot.title = element_text(face = "bold", size = 10),
        plot.tag = element_text(face = "bold", size = 11),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3))
theme_set(theme_ms)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
  message("saved ", name)
}
stars <- function(q) dplyr::case_when(is.na(q) ~ "", q < 0.001 ~ "***", q < 0.01 ~ "**", q < 0.05 ~ "*", TRUE ~ "")

# ---- data ----
load("IsoFresh_Metrics_v2.RData")
if (file.exists("IsoFresh_Omnivory_v2.RData")) { load("IsoFresh_Omnivory_v2.RData"); final_df <- left_join(final_df, omnivory_v2_df, by = "FoodWeb_ID") }
d <- final_df %>% filter(Include_S10) %>% mutate(Log_Richness = log10(Richness))
battery <- read_csv("model_battery_results.csv", show_col_types = FALSE)

DRIVER_LAB <- c(Abs_Latitude = "Latitude", Mean_Annual_Temp = "Temperature", Temp_Seasonality = "Seasonality")
DRIVER_AXIS <- c(Latitude = "Latitude\n(per +1 SD)", Temperature = "Cooler temperature\n(per -1 SD)", Seasonality = "Temp. seasonality\n(per +1 SD)")
flip <- function(driver) ifelse(driver == "Temperature", -1, 1)  # all drivers point poleward

# ==========================================================================
# Fig 1  The pattern
# ==========================================================================
lmm_band <- function(dat, x, y, n = 80) {
  dd <- dat %>% filter(!is.na(.data[[x]]), !is.na(.data[[y]]), !is.na(Region), !is.na(Source_DOI))
  m <- lmer(as.formula(paste(y, "~", x, "+ (1 | Region/Source_DOI)")), data = dd)
  xr <- range(dd[[x]]); nd <- tibble(x = seq(xr[1], xr[2], length.out = n))
  X <- cbind(1, nd$x); V <- as.matrix(vcov(m))
  nd$fit <- as.numeric(X %*% fixef(m)); se <- sqrt(rowSums((X %*% V) * X))
  nd$lwr <- nd$fit - 1.96 * se; nd$upr <- nd$fit + 1.96 * se
  list(band = nd, p = summary(m)$coefficients[2, "Pr(>|t|)"], n = nrow(dd))
}
scatter_panel <- function(dat, x, y, xlab, ylab, log_x = FALSE) {
  lb <- lmm_band(dat, x, y)
  ptxt <- ifelse(lb$p < 0.001, "p < 0.001", sprintf("p = %.3f", lb$p))
  p <- ggplot(dat, aes(.data[[x]], .data[[y]])) +
    geom_point(aes(colour = Ecosystem_Grp), alpha = 0.35, size = 1, stroke = 0) +
    geom_ribbon(data = lb$band, aes(x = x, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "grey30", alpha = 0.18) +
    geom_line(data = lb$band, aes(x = x, y = fit), inherit.aes = FALSE, colour = "grey10", linewidth = 0.8) +
    scale_colour_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = NULL) +
    annotate("text", x = Inf, y = Inf, label = sprintf("n = %d, %s", lb$n, ptxt), hjust = 1.05, vjust = 1.6, size = 2.8) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    labs(x = xlab, y = ylab)
  if (log_x) p <- p + scale_x_continuous(labels = function(v) round(10^v))
  p
}
world <- map_data("world")
p1a <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group), fill = "grey90", colour = NA) +
  geom_point(data = d, aes(Longitude, Latitude, colour = Ecosystem_Grp), size = 0.9, alpha = 0.7, stroke = 0) +
  scale_colour_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = NULL) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
  coord_quickmap(ylim = c(-50, 80), expand = FALSE) + theme_void(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10), plot.tag = element_text(face = "bold", size = 11)) +
  labs(title = sprintf("%d freshwater food webs (S >= 10)", nrow(d)))
p1b <- scatter_panel(d, "Abs_Latitude", "Connectance", "Absolute latitude (deg)", "Connectance")
p1c <- scatter_panel(d, "Abs_Latitude", "Log_Richness", "Absolute latitude (deg)", "Number of taxa", log_x = FALSE) +
  scale_y_continuous(labels = function(v) round(10^v))
p1d <- scatter_panel(d, "Log_Richness", "Connectance", "Number of taxa", "Connectance", log_x = TRUE)
# one shared legend (from the map); scatter panels carry none
no_leg <- theme(legend.position = "none")
fig1 <- (p1a) / ((p1b + no_leg) | (p1c + no_leg) | (p1d + no_leg)) + plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "a")
save_fig(fig1, "Fig1_pattern", 180, 165)

# ==========================================================================
# Fig 2  Climate effects on network structure (3 drivers x 6 metrics)
# ==========================================================================
METRIC_LAB <- c(Connectance = "Connectance", z_NODF = "Nestedness (z-NODF)", Mean_Strength = "Mean link strength",
                Mean_Closeness = "Mean closeness", Modularity = "Modularity (Q)", z_Q = "Standardised modularity (z-Q)")
METRIC_GRP <- c(Connectance = "Density & nestedness", z_NODF = "Density & nestedness", Mean_Strength = "Density & nestedness",
                Mean_Closeness = "Density & nestedness", Modularity = "Modularity", z_Q = "Modularity")
f2 <- battery %>% filter(population == "S10", predictor %in% names(DRIVER_LAB), metric %in% names(METRIC_LAB)) %>%
  mutate(driver = factor(DRIVER_LAB[predictor], levels = DRIVER_LAB), sgn = flip(as.character(driver)),
         est = estimate * sgn, lwr = est - 1.96 * std.error, upr = est + 1.96 * std.error,
         metric_f = factor(METRIC_LAB[metric], levels = rev(METRIC_LAB)), grp = factor(METRIC_GRP[metric], levels = c("Density & nestedness", "Modularity")),
         eco = factor(ecosystem, levels = c("All", "Lentic", "Lotic")), sig = q.value < 0.05)
fig2 <- ggplot(f2, aes(est, metric_f, colour = eco, shape = sig, alpha = eco)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = lwr, xmax = upr), orientation = "y", width = 0, linewidth = 0.5, position = position_dodge(width = 0.7)) +
  geom_point(aes(size = eco), fill = "white", position = position_dodge(width = 0.7)) +
  scale_colour_manual(values = c(All = COL_ALL, Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = "Ecosystem") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = "BH-FDR") +
  scale_size_manual(values = c(All = 2.6, Lentic = 1.8, Lotic = 1.8), guide = "none") +
  scale_alpha_manual(values = c(All = 1, Lentic = 0.85, Lotic = 0.85), guide = "none") +
  facet_grid(grp ~ driver, scales = "free_y", space = "free_y") +
  scale_x_continuous(breaks = c(-0.4, 0, 0.4)) +
  labs(x = "Standardised effect of moving poleward (95% CI)\n(higher latitude | cooler temperature | higher seasonality)", y = NULL,
       title = "Climate effects on food-web structure", subtitle = "S >= 10 food webs; mixed models with Region/Source-study random effects") +
  guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2))
save_fig(fig2, "Fig2_climate_effects", 180, 105)

# ==========================================================================
# Fig 3  Direct vs diversity-mediated effects
# piecewise SEM with sampling effort (log10 individuals measured) as a covariate
# in every equation; produced by sem_nsample.R (identical sample, n = 474)
# ==========================================================================
tr <- read_csv("sem_nsample_traced_effects.csv", show_col_types = FALSE) %>% filter(model == "+ sampling effort")
OUT_LAB <- c(Connectance_z = "Connectance", zNODF_z = "Nestedness (z-NODF)", Strength_z = "Mean link strength")
f3 <- tr %>% mutate(sgn = flip(driver),
         `Direct effect` = direct * sgn, `Via taxon richness` = via_richness * sgn, `Via prey niche area` = via_prey * sgn,
         total = total * sgn, outcome_f = factor(OUT_LAB[outcome], levels = rev(OUT_LAB)),
         driver = factor(driver, levels = c("Latitude", "Temperature", "Seasonality")),
         dir_star = ifelse(direct_p < 0.05, "direct *", "direct n.s.")) %>%
  pivot_longer(c(`Direct effect`, `Via taxon richness`, `Via prey niche area`), names_to = "path", values_to = "eff") %>%
  mutate(path = factor(path, levels = c("Direct effect", "Via taxon richness", "Via prey niche area")))
tot <- f3 %>% distinct(driver, outcome_f, total, dir_star)
fig3 <- ggplot(f3, aes(eff, outcome_f)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_col(aes(fill = path), width = 0.6) +
  geom_point(data = tot, aes(total, outcome_f), shape = 18, size = 3.2, colour = "black", inherit.aes = FALSE) +
  geom_text(data = tot, aes(x = 0.012, y = outcome_f, label = dir_star), hjust = 0, size = 2.6, colour = "grey25", inherit.aes = FALSE) +
  scale_fill_manual(values = c(`Direct effect` = COL_DIRECT, `Via taxon richness` = COL_RICH, `Via prey niche area` = COL_PREY), name = NULL) +
  facet_wrap(~driver, nrow = 1) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.42)), breaks = c(-0.3, -0.2, -0.1, 0)) +
  labs(x = "Standardised effect of moving poleward (path-traced, piecewise SEM)", y = NULL,
       title = "Poleward effects split into direct and diversity-mediated paths",
       subtitle = "Sampling effort controlled. Bars = decomposition of total effect (diamond); text = direct path significance")
save_fig(fig3, "Fig3_direct_vs_mediated", 180, 80)

# ==========================================================================
# Fig 4  SEM path coefficients (Latitude as driver, sampling effort controlled)
# ==========================================================================
cl <- read_csv("sem_nsample_coefficients.csv", show_col_types = FALSE) %>% filter(driver == "Latitude", model == "+ sampling effort")
RESP_LAB <- c(Richness_z = "Taxon richness", PreyArea_z = "Prey niche area", Connectance_z = "Connectance",
              zNODF_z = "Nestedness (z-NODF)", Strength_z = "Mean link strength")
PRED_LAB <- c(Exog_z = "Latitude", Richness_z = "Taxon richness", PreyArea_z = "Prey niche area", N_z = "Sampling effort (covariate)")
f4 <- cl %>% filter(Response %in% names(RESP_LAB), Predictor %in% names(PRED_LAB)) %>%
  mutate(resp = factor(RESP_LAB[Response], levels = RESP_LAB), pred = factor(PRED_LAB[Predictor], levels = rev(PRED_LAB)),
         grp = factor(ifelse(Response %in% c("Richness_z", "PreyArea_z"), "Mediators", "Network structure"), levels = c("Mediators", "Network structure")),
         st = stars(P.Value), lab = sprintf("%.2f%s", Std.Estimate, st))
fig4 <- ggplot(f4, aes(resp, pred, fill = Std.Estimate)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = lab, fontface = ifelse(P.Value < 0.05, "bold", "plain")), size = 3.1) +
  scale_fill_gradient2(low = DIV_NEG, mid = DIV_MID, high = DIV_POS, midpoint = 0, limits = c(-1, 1), oob = squish, name = "Std. path\ncoefficient") +
  facet_grid(~grp, scales = "free_x", space = "free_x") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "Predictor", title = "Path coefficients of the piecewise SEM (Latitude as driver)",
       subtitle = "Rows predict columns; blank = path not in the model; bold* = p < 0.05; n = 474 food webs") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right", strip.text = element_text(face = "bold"))
fig4 <- fig4 + scale_x_discrete(position = "top", labels = function(v) gsub(" \\(", "\n(", gsub("Mean link", "Mean link\n", v)))
save_fig(fig4, "Fig4_sem_paths", 180, 85)


# ==========================================================================
# Supplementary
# ==========================================================================
# S1 threshold sensitivity
load("IsoFresh_Metrics_v2_thresholdsweep.RData")
s1 <- threshold_sweep_df %>% mutate(sig = p.value < 0.05, sgn = ifelse(term == "Mean_Annual_Temp", -1, 1), est = estimate * sgn)
figS1 <- ggplot(s1, aes(threshold, est)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_line(colour = COL_LENTIC) + geom_point(aes(shape = sig), size = 2.6, colour = COL_LENTIC, fill = "white") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  labs(x = "Global similarity threshold used to define links", y = "Standardised effect of cooler temperature on Connectance",
       title = "The Connectance-climate relationship is robust to the link threshold")
save_fig(figS1, "FigS1_threshold_sensitivity", 120, 80)

# S2 bootstrap mediation
bm <- read_csv("mediation_bootstrap_nsample.csv", show_col_types = FALSE) %>%
  mutate(driver = factor(DRIVER_LAB[predictor], levels = DRIVER_LAB), sgn = flip(as.character(driver)),
         est = indirect_effect * sgn, lo = pmin(ci_lwr * sgn, ci_upr * sgn), hi = pmax(ci_lwr * sgn, ci_upr * sgn), sig = as.logical(sig),
         metric_f = factor(METRIC_LAB[metric], levels = rev(METRIC_LAB)),
         med = factor(ifelse(mediator == "Log_Richness", "Taxon richness", "Prey niche area"), levels = c("Taxon richness", "Prey niche area")))
figS2 <- ggplot(bm, aes(est, metric_f, colour = med, shape = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.2, fill = "white", position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(`Taxon richness` = COL_RICH, `Prey niche area` = COL_DIRECT), name = "Mediator") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "95% CI excludes 0", `FALSE` = "includes 0"), name = NULL) +
  facet_wrap(~driver) +
  labs(x = "Indirect effect of moving poleward (standardised a x b; cluster-bootstrap 95% CI)", y = NULL,
       title = "Bootstrap mediation by taxon richness and prey niche area", subtitle = "Sampling effort (log individuals) included in both paths")
save_fig(figS2, "FigS2_bootstrap_mediation", 180, 80)

# S3 all predictors x all metrics heatmap (S10, pooled ecosystems)
PRED_GRP <- c(Abs_Latitude = "Climate", Mean_Annual_Temp = "Climate", Temp_Seasonality = "Climate",
              Log_Richness = "Community", Log_N_Prey = "Community", Prey_Niche_Area = "Community", N_Modules = "Community",
              BodySize_Ratio = "Fish traits", BodySize_MaxGap = "Fish traits", Mean_Omnivory_Index = "Fish traits",
              Omnivory_FoodSeTroph = "Fish traits", Omnivory_DietSeTroph = "Fish traits", Log_Lake_Area = "Habitat")
PRED_NAME <- c(Abs_Latitude = "Latitude", Mean_Annual_Temp = "Temperature", Temp_Seasonality = "Seasonality",
               Log_Richness = "Taxon richness", Log_N_Prey = "No. prey taxa", Prey_Niche_Area = "Prey niche area", N_Modules = "No. modules",
               BodySize_Ratio = "Size ratio", BodySize_MaxGap = "Size max gap", Mean_Omnivory_Index = "Omnivory v1",
               Omnivory_FoodSeTroph = "Omnivory (FoodSeTroph)", Omnivory_DietSeTroph = "Omnivory (DietSeTroph)", Log_Lake_Area = "Lake area")
s3 <- battery %>% filter(population == "S10", ecosystem == "All", predictor %in% names(PRED_GRP), metric %in% names(METRIC_LAB)) %>%
  mutate(pred = factor(PRED_NAME[predictor], levels = rev(unique(PRED_NAME))), grp = factor(PRED_GRP[predictor], levels = c("Climate", "Community", "Fish traits", "Habitat")),
         metric_f = factor(METRIC_LAB[metric], levels = METRIC_LAB), st = stars(q.value))
figS3 <- ggplot(s3, aes(metric_f, pred, fill = estimate)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f%s", estimate, st)), size = 2.5) +
  scale_fill_gradient2(low = DIV_NEG, mid = DIV_MID, high = DIV_POS, limits = c(-0.8, 0.8), oob = squish, name = "Std. effect") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") + scale_x_discrete(position = "top", guide = guide_axis(angle = 30)) +
  labs(x = NULL, y = NULL, title = "All candidate drivers vs. all network metrics", subtitle = "S >= 10, pooled ecosystems; * q<0.05, ** q<0.01, *** q<0.001 (BH-FDR)") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "right")
save_fig(figS3, "FigS3_all_predictors_heatmap", 180, 150)

# S4 interaction effects (all FDR n.s.)
si <- read_csv("seasonality_interaction.csv", show_col_types = FALSE) %>% filter(population == "S10") %>%
  transmute(test = "Temperature x seasonality", metric, ecosystem, est = interaction_est, q = interaction_q)
bi <- read_csv("bodysize_temperature_interaction.csv", show_col_types = FALSE) %>% filter(population == "S10") %>%
  transmute(test = paste0("Temperature x body-size (", ifelse(size_metric == "SDsize", "SD", "ratio"), ")"), metric = "Connectance", ecosystem, est = interaction_est, q = interaction_q)
s4 <- bind_rows(si, bi) %>% mutate(ecosystem = factor(ecosystem, levels = c("All", "Lentic", "Lotic")),
                                   row = ifelse(test == "Temperature x seasonality", METRIC_LAB[metric], test),
                                   row = factor(row, levels = rev(unique(row))), sig = q < 0.05)
figS4 <- ggplot(s4, aes(est, row, colour = ecosystem, shape = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(size = 2.2, fill = "white", position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(All = COL_ALL, Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = "Ecosystem") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), labels = c(`TRUE` = "q < 0.05", `FALSE` = "n.s."), name = "BH-FDR") +
  labs(x = "Standardised interaction coefficient", y = NULL, title = "Moderation of climate effects: no interaction survives FDR correction",
       subtitle = "Seasonality x mean temperature (all six metrics) and body-size x temperature (Connectance); S >= 10")
save_fig(figS4, "FigS4_interactions", 180, 100)
message("done -> ", OUT)
