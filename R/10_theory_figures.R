# Figures for the theory-driven tests (reads theory_results/*.csv and the
# per-web analysis table written by theory_tests.R / theory_sem.R).
suppressMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork); library(lme4); library(lmerTest)})
OUT <- "manuscript_figs"; dir.create(OUT, showWarnings = FALSE)
R <- function(f) read_csv(file.path("theory_results", f), show_col_types = FALSE)
d <- readRDS("theory_results/theory_analysis_df.rds")
COL_INC <- "#1c5cab"; COL_DEC <- "#e34948"; COL_NONE <- "grey60"; COL_LENTIC <- "#2a78d6"; COL_LOTIC <- "#eb6834"
theme_set(theme_classic(base_size = 9) + theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
  plot.title = element_text(face = "bold", size = 9.5), plot.tag = element_text(face = "bold", size = 11), legend.position = "bottom",
  axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3)))
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white"); message("saved ", name)
}
re <- "(1 | Region/Source_DOI)"
lmm_line <- function(dat, y, x, covs = "Log_N_Ind", n = 60) {
  dd <- dat %>% filter(!is.na(.data[[y]]), !is.na(.data[[x]]), if_all(all_of(covs), ~ !is.na(.x)), !is.na(Region), !is.na(Source_DOI))
  m <- suppressMessages(suppressWarnings(lmer(as.formula(paste(y, "~", x, "+", paste(covs, collapse = "+"), "+", re)), dd)))
  nd <- tibble(!!x := seq(min(dd[[x]]), max(dd[[x]]), length.out = n)); for (cv in covs) nd[[cv]] <- mean(dd[[cv]])
  X <- model.matrix(as.formula(paste("~", x, "+", paste(covs, collapse = "+"))), nd); V <- as.matrix(vcov(m))
  nd$fit <- as.numeric(X %*% fixef(m)); se <- sqrt(rowSums((X %*% V) * X)); nd$lwr <- nd$fit - 1.96 * se; nd$upr <- nd$fit + 1.96 * se
  list(band = nd, coef = summary(m)$coefficients[x, ], n = nrow(dd))
}

# ==========================================================================
# Fig 5  Test of Petchey, Brose & Rall (2010)
# ==========================================================================
ps <- R("T1_petchey_predicted_sign_from_activation_energies.csv") %>%
  mutate(group = factor(recode(group, "Freshwater ectotherm vertebrates (fish)" = "Freshwater fish", "Freshwater invertebrates" = "Freshwater invertebrates", "All taxa pooled" = "All taxa pooled", "Marine ectotherm vertebrates" = "Marine fish"), levels = rev(c("Freshwater fish", "Freshwater invertebrates", "All taxa pooled", "Marine fish"))), lab = sprintf("%.2f / %.2f", E_A, E_H),
         pred = factor(predicted_effect_of_warming_on_connectance, levels = c("increase", "~none", "decrease")))
p5a <- ggplot(ps, aes(E_sum, group, fill = pred)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") + geom_col(width = 0.6) +
  geom_text(aes(x = ifelse(E_sum > 0, E_sum + 0.03, 0.03), label = lab), hjust = 0, size = 2.4, colour = "grey25") +
  scale_fill_manual(values = c(increase = COL_INC, `~none` = COL_NONE, decrease = COL_DEC), name = "Warming should") +
  scale_x_continuous(limits = c(-0.25, 1.3), breaks = c(0, 0.5, 1)) +
  labs(x = expression(E[A] + E[H] ~ "(eV);  labels:" ~ E[A] / E[H]), y = NULL, title = "Petchey et al. 2010 prediction\nwith Rall et al. 2012 activation energies") +
  theme(legend.position = "bottom")
t1 <- R("T1_temperature_connectance_raw_and_arrhenius.csv") %>% filter(grepl("^Connectance per deg C", test)) %>%
  mutate(lab = ifelse(grepl("total", test), "Total", "Direct\n(richness, prey area fixed)"), lab = factor(lab, levels = rev(unique(lab))))
p5b <- ggplot(t1, aes(est, lab)) +
  annotate("rect", xmin = 0.0005, xmax = 0.01, ymin = -Inf, ymax = Inf, fill = COL_INC, alpha = 0.12) +
  annotate("text", x = 0.0052, y = 2.55, label = "Range predicted by Petchey et al. 2010", size = 2.3, colour = COL_INC) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, orientation = "y", linewidth = 0.6) + geom_point(size = 2.4) +
  scale_x_continuous(limits = c(-0.001, 0.0105), breaks = c(0, 0.005, 0.01)) + scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(x = "Connectance change per +1 deg C", y = NULL, title = "Observed slope\n(within the predicted range)")
t2 <- R("T2_temperature_slope_by_size_tertile.csv") %>% mutate(tertile = factor(tertile, levels = c("Similar sizes", "Intermediate", "Dissimilar sizes")))
t2i <- R("T2_temperature_x_bodysize.csv") %>% filter(term == "Tz:Bz", size_metric == "SD_logMass", ecosystem == "All")
p5c <- ggplot(t2, aes(tertile, est)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0, linewidth = 0.6) + geom_point(size = 2.4, colour = COL_INC) +
  annotate("text", x = 3.45, y = max(t2$upr) + 0.08, hjust = 1, vjust = 1, size = 2.4,
           label = sprintf("interaction %.2f (p = %.2f)", t2i$est, t2i$p)) +
  scale_x_discrete(labels = c("Similar", "Intermediate", "Dissimilar")) +
  labs(x = "Fish body-size dispersion (tertile)", y = "Std. temperature effect\non connectance",
       title = "Similar body sizes:\nmore temperature-sensitive?")
lnS <- lmm_line(d %>% mutate(lnS = log(Richness)), "lnS", "InvKT")
mte_int <- mean(log(d$Richness), na.rm = TRUE) + 0.65 * mean(d$InvKT, na.rm = TRUE)
t1e <- R("T1_temperature_connectance_raw_and_arrhenius.csv") %>% filter(grepl("Richness", test))
p5d <- ggplot(d, aes(InvKT, log(Richness))) +
  geom_point(aes(colour = Ecosystem_Grp), alpha = 0.3, size = 0.8, stroke = 0) +
  geom_abline(intercept = mte_int, slope = -0.65, linetype = "dashed", colour = "grey30") +
  geom_ribbon(data = lnS$band, aes(x = InvKT, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "grey30", alpha = 0.2) +
  geom_line(data = lnS$band, aes(InvKT, fit), inherit.aes = FALSE, linewidth = 0.8) +
  annotate("text", x = 4.2, y = 4.5, hjust = 0, vjust = 1, size = 2.3, label = sprintf("Observed E = %.2f eV\n(95%% CI %.2f-%.2f)\nMTE: E = 0.65 eV (dashed)", t1e$est, t1e$lwr, t1e$upr)) +
  scale_x_reverse() + scale_colour_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = NULL, guide = "none") +
  coord_cartesian(ylim = c(2.2, 4.6)) +
  labs(x = "1/kT (centred at 15 deg C; warmer to right)", y = "ln(number of taxa)",
       title = "Richness vs temperature:\nshallower than metabolic theory")
fig5 <- (p5a | p5b) / (p5c | p5d) + plot_annotation(tag_levels = "a")
save_fig(fig5, "Fig5_petchey_test", 180, 150)

# ==========================================================================
# Fig 6  Isotopic niche geometry: expansion vs packing, link-species scaling
# ==========================================================================
t4 <- R("T4_niche_packing_expansion.csv")
sl <- t4 %>% filter(model == "scaling with richness (log-log)")
p6a <- ggplot(d, aes(Log_S, Log_TA)) + geom_point(alpha = 0.3, size = 0.8, stroke = 0, colour = "grey30") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = FALSE) +
  annotate("text", x = 1.02, y = max(d$Log_TA, na.rm = TRUE), hjust = 0, vjust = 1, size = 2.5,
           label = sprintf("slope = %.2f\n(1 = pure expansion)", sl$est[sl$response == "Log_TA"])) +
  scale_x_continuous(labels = function(v) round(10^v)) +
  labs(x = "Number of taxa (log scale)", y = "log10 isotopic total area (per mil^2)", title = "Niche expansion")
p6b <- ggplot(d, aes(Log_S, Log_MNND)) + geom_point(alpha = 0.3, size = 0.8, stroke = 0, colour = "grey30") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "black", linewidth = 0.7, se = FALSE) +
  annotate("text", x = 1.02, y = min(d$Log_MNND, na.rm = TRUE), hjust = 0, vjust = 0, size = 2.5,
           label = sprintf("slope = %.2f\n(0 = pure expansion; -0.5 = pure packing)", sl$est[sl$response == "Log_MNND"])) +
  scale_x_continuous(labels = function(v) round(10^v)) +
  labs(x = "Number of taxa (log scale)", y = "log10 mean nearest-neighbour distance", title = "Niche packing")
cl <- t4 %>% filter(model == "climate effect | richness, effort") %>%
  mutate(driver = recode(term, Abs_Latitude = "Latitude", Mean_Annual_Temp = "Temperature", Temp_Seasonality = "Seasonality"),
         response = recode(response, Packing_Ratio = "Packing ratio (obs/null)", Log_TA = "Total area", Log_MNND = "Nearest-neighbour distance"),
         driver = factor(driver, levels = c("Latitude", "Temperature", "Seasonality")))
p6c <- ggplot(cl, aes(est, response, colour = driver)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = lwr, xmax = upr), width = 0, orientation = "y", position = position_dodge(0.6)) +
  geom_point(size = 2, shape = 21, fill = "white", position = position_dodge(0.6)) +
  scale_colour_manual(values = c(Latitude = "grey20", Temperature = COL_LOTIC, Seasonality = COL_LENTIC), name = NULL) +
  labs(x = "Std. climate effect, richness held constant", y = NULL, title = "Climate does not alter isotopic geometry\nbeyond richness (all n.s.)")
t3 <- R("T3_link_species_scaling_x_latitude.csv")
dl <- d %>% filter(is.finite(Log_L)) %>% mutate(LatT = cut(Abs_Latitude, quantile(Abs_Latitude, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE,
                                                        labels = c("Low lat.", "Mid lat.", "High lat.")))
p6d <- ggplot(dl, aes(Log_S, Log_L, colour = LatT)) + geom_point(alpha = 0.25, size = 0.7, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  scale_colour_manual(values = c("#eb6834", "grey50", "#1c5cab"), name = NULL) +
  annotate("text", x = 1, y = max(dl$Log_L), hjust = 0, vjust = 1, size = 2.5,
           label = sprintf("L ~ S^%.2f (2 = constant C)\nS x latitude, lentic: %.3f (p = %.2f)\nS x latitude, lotic: %.3f (p = %.2f)",
                           t3$est[t3$term == "Log_S_c" & t3$ecosystem == "All"],
                           t3$est[t3$term == "Log_S_c:Lat_z" & t3$ecosystem == "Lentic"], t3$p[t3$term == "Log_S_c:Lat_z" & t3$ecosystem == "Lentic"],
                           t3$est[t3$term == "Log_S_c:Lat_z" & t3$ecosystem == "Lotic"], t3$p[t3$term == "Log_S_c:Lat_z" & t3$ecosystem == "Lotic"])) +
  scale_x_continuous(labels = function(v) round(10^v)) + scale_y_continuous(labels = function(v) round(10^v)) +
  labs(x = "Number of taxa (log scale)", y = "Number of links (log scale)", title = "Link-species scaling vs latitude")
fig6 <- (p6a | p6b) / (p6c | p6d) + plot_annotation(tag_levels = "a")
save_fig(fig6, "Fig6_niche_geometry", 180, 150)

# ==========================================================================
# Fig 7  Generalism pathway and integrated decomposition
# ==========================================================================
th <- lmm_line(d, "Fish_Trophic_Height", "Mean_Annual_Temp", covs = c("Log_S", "Log_N_Ind"))
p7a <- ggplot(d, aes(Mean_Annual_Temp, Fish_Trophic_Height)) +
  geom_point(aes(colour = Ecosystem_Grp), alpha = 0.3, size = 0.8, stroke = 0) +
  geom_ribbon(data = th$band, aes(x = Mean_Annual_Temp, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "grey30", alpha = 0.2) +
  geom_line(data = th$band, aes(Mean_Annual_Temp, fit), inherit.aes = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = NULL) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + coord_cartesian(ylim = c(-1, 2.8)) +
  annotate("text", x = -10, y = 2.75, hjust = 0, vjust = 1, size = 2.5, label = sprintf("slope = %.3f per deg C, p = %.3f\n(richness and effort held constant)", th$coef[1], th$coef[5])) +
  labs(x = "Mean annual temperature (deg C)", y = "Fish trophic height above other taxa\n(d15N difference / 3.4)", title = "Warm-water fish sit lower in the web\n(more omnivory)")
cn <- lmm_line(d, "Connectance", "Fish_Trophic_Height", covs = c("Log_S", "Log_N_Ind"))
p7b <- ggplot(d, aes(Fish_Trophic_Height, Connectance)) +
  geom_point(aes(colour = Ecosystem_Grp), alpha = 0.3, size = 0.8, stroke = 0) +
  geom_ribbon(data = cn$band, aes(x = Fish_Trophic_Height, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = "grey30", alpha = 0.2) +
  geom_line(data = cn$band, aes(Fish_Trophic_Height, fit), inherit.aes = FALSE, linewidth = 0.8) +
  scale_colour_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = NULL) + guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) + coord_cartesian(xlim = c(-1, 2.8)) +
  annotate("text", x = 2.75, y = 0.93, hjust = 1, vjust = 1, size = 2.5, label = sprintf("slope p = %.3f\n(richness, effort held constant)", cn$coef[5])) +
  labs(x = "Fish trophic height above other taxa", y = "Connectance", title = "Lower fish trophic height,\nhigher isotopic overlap")
t6 <- R("T6_theory_sem_connectance.csv") %>% mutate(sgn = ifelse(driver == "Temperature", -1, 1)) %>%
  pivot_longer(c(direct, via_richness, via_prey_area, via_fish_trophic_height), names_to = "path", values_to = "v") %>%
  mutate(v = v * sgn, total = total * sgn, driver = factor(driver, levels = c("Latitude", "Temperature", "Seasonality")),
         path = factor(path, levels = c("direct", "via_richness", "via_prey_area", "via_fish_trophic_height"),
                       labels = c("Direct (foraging / unmeasured)", "Via taxon richness", "Via prey niche area", "Via fish trophic height")))
p7c <- ggplot(t6, aes(v, driver, fill = path)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey40") +
  geom_col(width = 0.6) + geom_point(data = t6 %>% distinct(driver, total), aes(total, driver), inherit.aes = FALSE, shape = 18, size = 3.2) +
  scale_fill_manual(values = c("#eb6834", "#2a78d6", "#9cc4f1", "#7a5195"), name = NULL) + guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Std. effect of moving poleward on connectance (diamond = total)", y = NULL,
       title = "Decomposition of the poleward decline in connectance (piecewise SEM, n = 473)")
fig7 <- (p7a | p7b) / p7c + plot_layout(heights = c(1, 0.8)) + plot_annotation(tag_levels = "a")
save_fig(fig7, "Fig7_generalism_pathway", 180, 150)
