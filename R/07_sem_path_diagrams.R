# Path diagrams of the piecewise SEM (one per climate driver: latitude, mean
# annual temperature, temperature seasonality), sampling effort as a covariate in
# every equation. Reads sem_nsample_coefficients.csv (sem_nsample.R).
# Sampling effort = log10 of the summed IsoFresh N_Sample (samples analysed per taxon).
suppressMessages({library(dplyr); library(readr); library(ggplot2)})
coef_all <- read_csv("sem_nsample_coefficients.csv", show_col_types = FALSE) %>% filter(model == "+ sampling effort")
DRIVERS <- list(
  Latitude    = list(label = "Latitude\n(absolute)",                         title = "latitude",    file = "Fig8_sem_diagram"),
  Temperature = list(label = "Mean annual\ntemperature",                      title = "temperature", file = "FigS10_sem_diagram_temperature"),
  Seasonality = list(label = "Temperature\nseasonality\n(SD of monthly means)", title = "temperature seasonality", file = "FigS11_sem_diagram_seasonality"))
for (drv in names(DRIVERS)) {
info <- DRIVERS[[drv]]
cf <- coef_all %>% filter(driver == drv)

nodes <- tribble(~id, ~label, ~x, ~y, ~fill,
  "Exog_z",        info$label,                            0.0,  0.0, "#f1e4d3",
  "N_z",           "Sampling effort\n(log no. of samples)", 0.0,  2.7, "#eeeeee",
  "Richness_z",    "Taxon richness",                      2.9,  2.7, "#dbe8f8",
  "PreyArea_z",    "Prey niche area",                     2.9, -2.7, "#dbe8f8",
  "Connectance_z", "Connectance",                         5.9,  1.9, "#ffffff",
  "zNODF_z",       "Nestedness\n(z-NODF)",                5.8,  0.0, "#ffffff",
  "Strength_z",    "Mean link\nstrength",                 5.9, -1.9, "#ffffff")
HW <- 0.78; HH <- 0.4
# show every structural path; for the effort covariate only its path to richness (others in caption)
edges <- cf %>% filter(Predictor != "N_z" | Response == "Richness_z") %>%
  left_join(nodes %>% select(id, x1 = x, y1 = y), by = c("Predictor" = "id")) %>%
  left_join(nodes %>% select(id, x2 = x, y2 = y), by = c("Response" = "id")) %>%
  mutate(dx = x2 - x1, dy = y2 - y1,
         t0 = pmin(HW / pmax(abs(dx), 1e-9), HH / pmax(abs(dy), 1e-9)),
         xs = x1 + t0 * dx, ys = y1 + t0 * dy, xe = x2 - t0 * dx, ye = y2 - t0 * dy,
         sig = P.Value < 0.05, sign = ifelse(Std.Estimate > 0, "positive", "negative"),
         lw = 0.3 + 3 * abs(Std.Estimate),
         lab = sprintf("%.2f%s", Std.Estimate, ifelse(P.Value < 0.001, "***", ifelse(P.Value < 0.01, "**", ifelse(P.Value < 0.05, "*", "")))),
         pos = 0.5)
# nudge label positions on crowded edges
edges <- edges %>% mutate(pos = case_when(Predictor == "Exog_z" & Response == "Connectance_z" ~ 0.3,
                                          Predictor == "Richness_z" & Response == "PreyArea_z" ~ 0.25,
                                          Predictor == "Exog_z" & Response == "Strength_z" ~ 0.3,
                                          Predictor == "Exog_z" & Response == "zNODF_z" ~ 0.68,
                                          Predictor == "Richness_z" & Response == "Strength_z" ~ 0.25,
                                          Predictor == "PreyArea_z" & Response == "Connectance_z" ~ 0.25,
                                          Predictor == "Richness_z" & Response == "zNODF_z" ~ 0.3,
                                          Predictor == "PreyArea_z" & Response == "zNODF_z" ~ 0.3,
                                          TRUE ~ pos),
                          lx = xs + pos * (xe - xs), ly = ys + pos * (ye - ys))

p <- ggplot() +
  geom_segment(data = edges, aes(x = xs, y = ys, xend = xe, yend = ye, colour = sign, linetype = sig, linewidth = lw),
               arrow = arrow(length = unit(2.2, "mm"), type = "closed"), lineend = "butt", linejoin = "mitre") +
  geom_tile(data = nodes, aes(x, y, fill = fill), width = 2 * HW, height = 2 * HH, colour = "grey30", linewidth = 0.4) +
  geom_text(data = nodes, aes(x, y, label = label), size = 2.9, lineheight = 0.9) +
  geom_label(data = edges, aes(lx, ly, label = lab, colour = sign), size = 2.5, linewidth = 0, label.padding = unit(0.6, "mm"),
             fill = "white", fontface = ifelse(edges$sig, "bold", "plain"), show.legend = FALSE) +
  scale_fill_identity() + scale_linewidth_identity() +
  scale_colour_manual(values = c(positive = "#1c5cab", negative = "#d03b3a"), name = NULL, labels = c(negative = "Negative path", positive = "Positive path")) +
  scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "22"), labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."), name = NULL) +
  coord_equal(xlim = c(-0.95, 6.85), ylim = c(-3.3, 3.3), expand = FALSE) +
  labs(title = paste0("Piecewise SEM: ", info$title, ", taxon richness and network structure"),
       subtitle = "Standardised path coefficients; arrow width proportional to |coefficient|; n = 474 food webs") +
  theme_void(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"), plot.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)), linetype = guide_legend(override.aes = list(linewidth = 0.7)))
ggsave(file.path("manuscript_figs", paste0(info$file, ".png")), p, width = 170, height = 155, units = "mm", dpi = 300, bg = "white")
ggsave(file.path("manuscript_figs", paste0(info$file, ".pdf")), p, width = 170, height = 155, units = "mm")
message("saved ", info$file)
}
