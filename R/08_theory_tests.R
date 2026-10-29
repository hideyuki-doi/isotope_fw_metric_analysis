# ==========================================================================
# Theory-driven tests of the IsoFresh climate - network-structure results
#
# T1  Petchey, Brose & Rall (2010) temperature-ADBM + metabolic theory:
#       raw connectance slope per deg C (total and "fixed community" direct),
#       Boltzmann-Arrhenius apparent activation energies of connectance and
#       richness (MTE richness prediction E ~ 0.65 eV, Allen et al. 2002),
#       sign predicted from measured activation energies (Rall et al. 2012).
# T2  Petchey et al. (2010) Fig. 6: sensitivity of connectance to temperature
#       is larger when body sizes are more similar -> T x size-dispersion.
# T3  Latitude-niche-breadth hypothesis as tested by Cirtwill et al. (2015):
#       does the link-species scaling exponent change with latitude?
# T4  Niche packing vs niche expansion (MacArthur; Pigot et al. 2016):
#       isotopic total area (expansion) and nearest-neighbour spacing (packing)
#       against richness; packing ratio vs a uniform-in-hull null; climate.
# T5  Generalism / latitudinal omnivory gradient (Gonzalez-Bergonzoni et al.
#       2012) and energy-channel coupling (Rooney et al. 2006): fish trophic
#       height above other taxa and fish d13C centrality.
#
# All mixed models: (1 | Region/Source_DOI), log10 individuals measured as a
# covariate (sampling effort), S >= 10 webs. Outputs -> theory_results/
# ==========================================================================
suppressMessages({library(dplyr); library(tidyr); library(readr); library(lme4); library(lmerTest); library(parallel); library(sp)})
RAW  <- Sys.getenv("ISOFRESH_RAW", "data")
DATA <- Sys.getenv("ISOFRESH_DATA", ".")
N_NULL <- as.integer(Sys.getenv("ISOFRESH_N_NULL", "199")); N_CORES <- as.integer(Sys.getenv("ISOFRESH_CORES", "8")); SEED <- 42
OUT <- "theory_results"; dir.create(OUT, showWarnings = FALSE)
KB <- 8.617333e-5   # Boltzmann constant, eV/K

load(file.path(DATA, "IsoFresh_Metrics_v2.RData"))
rd <- function(f) read_delim(file.path(RAW, f), delim = ";", show_col_types = FALSE)
fish_data <- rd("IsoFresh_Isotope_Fish.csv"); other_data <- rd("IsoFresh_Isotope_Other.csv")

# taxon-level isotope means, built exactly as in the main pipeline
mk <- function(x, type) x %>% mutate(GIG = coalesce(Scientific_name, Common_name, paste0("Unidentified_", row_number()))) %>%
  select(FoodWeb_ID, GIG, d13C, d15N) %>% mutate(Type = type)
mean_isotop <- bind_rows(mk(fish_data, "Fish"), mk(other_data, "Other")) %>% filter(!is.na(d13C), !is.na(d15N)) %>%
  group_by(FoodWeb_ID, GIG, Type) %>% summarize(d15N = mean(d15N), d13C = mean(d13C), .groups = "drop")
effort <- bind_rows(fish_data %>% select(FoodWeb_ID, N_Sample, d13C, d15N), other_data %>% select(FoodWeb_ID, N_Sample, d13C, d15N)) %>%
  filter(!is.na(d13C), !is.na(d15N)) %>% group_by(FoodWeb_ID) %>%
  summarize(N_Individuals = ifelse(all(is.na(N_Sample)), NA_real_, sum(N_Sample, na.rm = TRUE)), .groups = "drop")

# ---- T4 / T5 per-web isotope-space metrics ---------------------------------
hull_area <- function(x, y) { h <- chull(x, y); xh <- x[h]; yh <- y[h]; 0.5 * abs(sum(xh * c(yh[-1], yh[1]) - c(xh[-1], xh[1]) * yh)) }
mnnd <- function(x, y) { d <- as.matrix(dist(cbind(x, y))); diag(d) <- Inf; mean(apply(d, 1, min)) }
# null: hull vertices fixed, interior taxa placed uniformly at random inside the hull
packing_null <- function(x, y, n_null, seed) {
  h <- chull(x, y); n_int <- length(x) - length(h)
  if (n_int < 2) return(c(NA_real_, NA_real_))
  hx <- x[h]; hy <- y[h]; set.seed(seed)
  sims <- replicate(n_null, {
    px <- numeric(0); py <- numeric(0)
    while (length(px) < n_int) {
      cx <- runif(4 * n_int, min(hx), max(hx)); cy <- runif(4 * n_int, min(hy), max(hy))
      ok <- point.in.polygon(cx, cy, hx, hy) == 1
      px <- c(px, cx[ok]); py <- c(py, cy[ok])
    }
    mnnd(c(hx, px[1:n_int]), c(hy, py[1:n_int]))
  })
  c(mean(sims), sd(sims))
}
webs <- split(mean_isotop, mean_isotop$FoodWeb_ID)
webs <- webs[names(webs) %in% final_df$FoodWeb_ID[final_df$Include_S10]]
cat(sprintf("Computing isotope-space metrics for %d webs (%d null reps, %d cores)...\n", length(webs), N_NULL, N_CORES))
iso_metrics <- bind_rows(mclapply(names(webs), function(id) {
  w <- webs[[id]]; S <- nrow(w)
  TA <- hull_area(w$d13C, w$d15N); obs <- mnnd(w$d13C, w$d15N)
  nl <- packing_null(w$d13C, w$d15N, N_NULL, SEED + match(id, names(webs)))
  fsh <- w %>% filter(Type == "Fish"); oth <- w %>% filter(Type == "Other")
  th <- cen <- NA_real_
  if (nrow(fsh) >= 2 && nrow(oth) >= 3) {
    th <- (mean(fsh$d15N) - mean(oth$d15N)) / 3.4
    cr <- diff(range(oth$d13C)); mid <- mean(range(oth$d13C))
    if (cr > 0) cen <- mean(abs(fsh$d13C - mid)) / (cr / 2)
  }
  tibble(FoodWeb_ID = id, TA = TA, MNND = obs, MNND_null = nl[1], MNND_null_sd = nl[2],
         Packing_Ratio = obs / nl[1], Packing_z = (obs - nl[1]) / nl[2],
         CR = diff(range(w$d13C)), NR = diff(range(w$d15N)), N_Fish = nrow(fsh), N_Other = nrow(oth),
         Fish_Trophic_Height = th, Fish_C_Centrality = cen)
}, mc.cores = N_CORES))
write_csv(iso_metrics, file.path(OUT, "iso_space_metrics.csv"))

d <- final_df %>% filter(Include_S10) %>% left_join(effort, by = "FoodWeb_ID") %>% left_join(iso_metrics, by = "FoodWeb_ID") %>%
  mutate(Log_S = log10(Richness), Log_N_Ind = log10(N_Individuals),
         L = round(Connectance * S * (S - 1) / 2), Log_L = log10(L), Mean_Degree = Connectance * (S - 1),
         InvKT = 1 / (KB * (Mean_Annual_Temp + 273.15)) - 1 / (KB * 288.15),   # centred at 15 deg C
         Log_TA = log10(TA), Log_MNND = log10(MNND), SD_logMass = 3 * SD_log10_Length)
saveRDS(d, file.path(OUT, "theory_analysis_df.rds"))
cat(sprintf("Webs: %d; with effort: %d; with packing ratio: %d; with fish metrics: %d\n", nrow(d), sum(!is.na(d$Log_N_Ind)),
            sum(!is.na(d$Packing_Ratio)), sum(!is.na(d$Fish_Trophic_Height))))

z <- function(x) as.numeric(scale(x))
re <- "(1 | Region/Source_DOI)"
fitlmm <- function(formula, dat) suppressMessages(suppressWarnings(lmer(as.formula(formula), data = dat)))
coef_row <- function(m, term) { s <- summary(m)$coefficients; ci <- s[term, "Estimate"] + c(-1.96, 1.96) * s[term, "Std. Error"]
  tibble(term = term, est = s[term, "Estimate"], se = s[term, "Std. Error"], lwr = ci[1], upr = ci[2], p = s[term, "Pr(>|t|)"], n = nobs(m)) }
cc <- function(dat, vars) dat %>% filter(if_all(all_of(c(vars, "Region", "Source_DOI")), ~ !is.na(.x)))

# ==========================================================================
# T1  Temperature - connectance in raw units and Arrhenius form
# ==========================================================================
t1 <- list()
dd <- cc(d, c("Connectance", "Mean_Annual_Temp", "Log_N_Ind"))
m <- fitlmm(paste("Connectance ~ Mean_Annual_Temp + Log_N_Ind +", re), dd)
t1[[1]] <- coef_row(m, "Mean_Annual_Temp") %>% mutate(test = "Connectance per deg C (total)")
dd2 <- cc(d, c("Connectance", "Mean_Annual_Temp", "Log_N_Ind", "Log_S", "Prey_Niche_Area"))
m <- fitlmm(paste("Connectance ~ Mean_Annual_Temp + Log_S + Prey_Niche_Area + Log_N_Ind +", re), dd2)
t1[[2]] <- coef_row(m, "Mean_Annual_Temp") %>% mutate(test = "Connectance per deg C (direct; richness & prey area fixed)")
m <- fitlmm(paste("Mean_Degree ~ Mean_Annual_Temp + Log_N_Ind +", re), dd)
t1[[3]] <- coef_row(m, "Mean_Annual_Temp") %>% mutate(test = "Mean degree (links per taxon) per deg C (total)")
m <- fitlmm(paste("Mean_Degree ~ Mean_Annual_Temp + Log_S + Prey_Niche_Area + Log_N_Ind +", re), dd2)
t1[[4]] <- coef_row(m, "Mean_Annual_Temp") %>% mutate(test = "Mean degree per deg C (direct)")
# Arrhenius: ln(Y) = c - E * (1/kT)  => E = -slope
for (v in c("Connectance", "Richness")) {
  ddv <- cc(d, c(v, "InvKT", "Log_N_Ind")); ddv$lnY <- log(ddv[[v]])
  m <- fitlmm(paste("lnY ~ InvKT + Log_N_Ind +", re), ddv)
  r <- coef_row(m, "InvKT"); t1[[length(t1) + 1]] <- r %>% mutate(est = -est, lwr0 = -upr, upr = -lwr, lwr = lwr0, test = paste0("Apparent activation energy of ", v, " (eV)")) %>% select(-lwr0)
}
ddv <- cc(d, c("Connectance", "InvKT", "Log_N_Ind", "Log_S", "Prey_Niche_Area")); ddv$lnY <- log(ddv$Connectance)
m <- fitlmm(paste("lnY ~ InvKT + Log_S + Prey_Niche_Area + Log_N_Ind +", re), ddv)
r <- coef_row(m, "InvKT"); t1[[length(t1) + 1]] <- r %>% mutate(est = -est, lwr0 = -upr, upr = -lwr, lwr = lwr0, test = "Apparent activation energy of Connectance, direct (eV)") %>% select(-lwr0)
t1 <- bind_rows(t1) %>% select(test, everything())
write_csv(t1, file.path(OUT, "T1_temperature_connectance_raw_and_arrhenius.csv"))

# Petchey et al. (2010): ln(lambda*H) shifts by -(E_A + E_H) * d(1/kT). Sign of E_A + E_H decides the sign of the
# temperature-connectance relationship (|E_A| < |E_H| -> warming increases connectance). Values: Rall et al. 2012 Table 1.
tq <- quantile(d$Mean_Annual_Temp, c(0.05, 0.95), na.rm = TRUE)
dinv <- 1 / (KB * (tq[1] + 273.15)) - 1 / (KB * (tq[2] + 273.15))
pred_sign <- tribble(~group, ~E_A, ~E_H,
  "Freshwater ectotherm vertebrates (fish)", 0.31, -0.46,
  "Freshwater invertebrates", 0.58, -0.56,
  "All taxa pooled", 0.42, -0.30,
  "Marine ectotherm vertebrates", 0.93, -0.09) %>%
  mutate(E_sum = E_A + E_H, dlog10_lambdaH = E_sum * dinv / log(10),
         predicted_effect_of_warming_on_connectance = case_when(E_sum < -0.02 ~ "increase", E_sum > 0.02 ~ "decrease", TRUE ~ "~none"),
         temp_range = sprintf("%.1f to %.1f deg C (5-95%% of webs)", tq[1], tq[2]))
write_csv(pred_sign, file.path(OUT, "T1_petchey_predicted_sign_from_activation_energies.csv"))

# ==========================================================================
# T2  Temperature x body-size dispersion (Petchey et al. 2010 Fig. 6)
# ==========================================================================
t2 <- list()
for (sz in c("SD_logMass", "BodySize_Ratio")) for (eco in c("All", "Lentic", "Lotic")) {
  de <- if (eco == "All") d else d %>% filter(Ecosystem_Grp == eco)
  dd <- cc(de, c("Connectance", "Mean_Annual_Temp", sz, "Log_N_Ind", "Log_S"))
  if (nrow(dd) < 40) next
  dd <- dd %>% mutate(Cz = z(Connectance), Tz = z(Mean_Annual_Temp), Bz = z(if (sz == "BodySize_Ratio") log10(.data[[sz]]) else .data[[sz]]), Sz = z(Log_S), Nz = z(Log_N_Ind))
  m <- fitlmm(paste("Cz ~ Tz * Bz + Sz + Nz +", re), dd)
  t2[[length(t2) + 1]] <- bind_rows(coef_row(m, "Tz"), coef_row(m, "Tz:Bz")) %>% mutate(size_metric = sz, ecosystem = eco)
}
t2 <- bind_rows(t2) %>% group_by(term) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
write_csv(t2, file.path(OUT, "T2_temperature_x_bodysize.csv"))
# simple slopes by size-dispersion tertile (pooled ecosystems)
dd <- cc(d, c("Connectance", "Mean_Annual_Temp", "SD_logMass", "Log_N_Ind", "Log_S")) %>%
  mutate(tert = cut(SD_logMass, quantile(SD_logMass, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = c("Similar sizes", "Intermediate", "Dissimilar sizes")))
t2b <- bind_rows(lapply(levels(dd$tert), function(tl) { x <- dd %>% filter(tert == tl) %>% mutate(Cz = z(Connectance), Tz = z(Mean_Annual_Temp), Sz = z(Log_S), Nz = z(Log_N_Ind))
  m <- fitlmm(paste("Cz ~ Tz + Sz + Nz +", re), x); coef_row(m, "Tz") %>% mutate(tertile = tl, SD_range = sprintf("%.2f-%.2f", min(x$SD_logMass), max(x$SD_logMass))) }))
write_csv(t2b, file.path(OUT, "T2_temperature_slope_by_size_tertile.csv"))

# ==========================================================================
# T3  Link-species scaling and latitude (Cirtwill et al. 2015)
# ==========================================================================
t3 <- list()
for (eco in c("All", "Lentic", "Lotic")) {
  de <- if (eco == "All") d else d %>% filter(Ecosystem_Grp == eco)
  dd <- cc(de, c("Log_L", "Log_S", "Abs_Latitude", "Log_N_Ind")) %>% filter(is.finite(Log_L)) %>%
    mutate(Lat_z = z(Abs_Latitude), Nz = z(Log_N_Ind), Log_S_c = Log_S - mean(Log_S))
  m <- fitlmm(paste("Log_L ~ Log_S_c * Lat_z + Nz +", re), dd)
  t3[[length(t3) + 1]] <- bind_rows(coef_row(m, "Log_S_c"), coef_row(m, "Lat_z"), coef_row(m, "Log_S_c:Lat_z")) %>% mutate(ecosystem = eco)
}
t3 <- bind_rows(t3); write_csv(t3, file.path(OUT, "T3_link_species_scaling_x_latitude.csv"))

# ==========================================================================
# T4  Niche packing vs expansion
# ==========================================================================
t4 <- list()
dd <- cc(d, c("Log_TA", "Log_MNND", "Log_S", "Log_N_Ind")) %>% filter(is.finite(Log_TA), is.finite(Log_MNND))
for (y in c("Log_TA", "Log_MNND")) { m <- fitlmm(paste(y, "~ Log_S + Log_N_Ind +", re), dd); t4[[length(t4) + 1]] <- coef_row(m, "Log_S") %>% mutate(response = y, model = "scaling with richness (log-log)") }
# climate effects on packing and expansion, richness held constant
for (y in c("Packing_Ratio", "Log_TA", "Log_MNND")) for (x in c("Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality")) {
  dd <- cc(d, c(y, x, "Log_S", "Log_N_Ind")) %>% filter(is.finite(.data[[y]])) %>% mutate(Yz = z(.data[[y]]), Xz = z(.data[[x]]), Sz = z(Log_S), Nz = z(Log_N_Ind))
  m <- fitlmm(paste("Yz ~ Xz + Sz + Nz +", re), dd)
  t4[[length(t4) + 1]] <- coef_row(m, "Xz") %>% mutate(term = x, response = y, model = "climate effect | richness, effort")
}
# does packing explain connectance and absorb the direct climate effect?
dd <- cc(d, c("Connectance", "Packing_Ratio", "Log_S", "Prey_Niche_Area", "Log_N_Ind", "Mean_Annual_Temp", "Abs_Latitude", "Temp_Seasonality")) %>%
  mutate(Cz = z(Connectance), Pz = z(Packing_Ratio), Sz = z(Log_S), Az = z(Prey_Niche_Area), Nz = z(Log_N_Ind))
for (x in c("Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality")) {
  dd$Xz <- z(dd[[x]])
  m0 <- fitlmm(paste("Cz ~ Xz + Sz + Az + Nz +", re), dd); m1 <- fitlmm(paste("Cz ~ Xz + Pz + Sz + Az + Nz +", re), dd)
  t4[[length(t4) + 1]] <- coef_row(m0, "Xz") %>% mutate(term = x, response = "Connectance", model = "direct climate effect (no packing)")
  t4[[length(t4) + 1]] <- coef_row(m1, "Xz") %>% mutate(term = x, response = "Connectance", model = "direct climate effect | packing")
  if (x == "Abs_Latitude") t4[[length(t4) + 1]] <- coef_row(m1, "Pz") %>% mutate(term = "Packing_Ratio", response = "Connectance", model = "packing effect on connectance")
}
t4 <- bind_rows(t4) %>% group_by(model) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
write_csv(t4, file.path(OUT, "T4_niche_packing_expansion.csv"))

# ==========================================================================
# T5  Generalism: fish trophic height and energy-channel coupling
# ==========================================================================
t5 <- list()
for (y in c("Fish_Trophic_Height", "Fish_C_Centrality")) for (x in c("Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality")) {
  dd <- cc(d, c(y, x, "Log_S", "Log_N_Ind")) %>% mutate(Yz = z(.data[[y]]), Xz = z(.data[[x]]), Sz = z(Log_S), Nz = z(Log_N_Ind))
  m <- fitlmm(paste("Yz ~ Xz + Sz + Nz +", re), dd)
  t5[[length(t5) + 1]] <- coef_row(m, "Xz") %>% mutate(term = x, response = y, model = "climate effect | richness, effort")
}
dd <- cc(d, c("Connectance", "Fish_Trophic_Height", "Fish_C_Centrality", "Log_S", "Prey_Niche_Area", "Log_N_Ind", "Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality")) %>%
  mutate(Cz = z(Connectance), THz = z(Fish_Trophic_Height), CEz = z(Fish_C_Centrality), Sz = z(Log_S), Az = z(Prey_Niche_Area), Nz = z(Log_N_Ind))
m <- fitlmm(paste("Cz ~ THz + CEz + Sz + Az + Nz +", re), dd)
t5[[length(t5) + 1]] <- bind_rows(coef_row(m, "THz"), coef_row(m, "CEz")) %>% mutate(response = "Connectance", model = "fish generalism -> connectance")
for (x in c("Abs_Latitude", "Mean_Annual_Temp", "Temp_Seasonality")) {
  dd$Xz <- z(dd[[x]])
  m0 <- fitlmm(paste("Cz ~ Xz + Sz + Az + Nz +", re), dd); m1 <- fitlmm(paste("Cz ~ Xz + THz + CEz + Sz + Az + Nz +", re), dd)
  t5[[length(t5) + 1]] <- coef_row(m0, "Xz") %>% mutate(term = x, response = "Connectance", model = "direct climate effect (fish subset, no generalism)")
  t5[[length(t5) + 1]] <- coef_row(m1, "Xz") %>% mutate(term = x, response = "Connectance", model = "direct climate effect | fish generalism")
}
t5 <- bind_rows(t5) %>% group_by(model) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup()
write_csv(t5, file.path(OUT, "T5_fish_generalism.csv"))

for (f in c("T1_temperature_connectance_raw_and_arrhenius", "T1_petchey_predicted_sign_from_activation_energies", "T2_temperature_x_bodysize",
            "T2_temperature_slope_by_size_tertile", "T3_link_species_scaling_x_latitude", "T4_niche_packing_expansion", "T5_fish_generalism")) {
  cat("\n=====", f, "=====\n"); x <- read_csv(file.path(OUT, paste0(f, ".csv")), show_col_types = FALSE)
  print(as.data.frame(x %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))), row.names = FALSE)
}
cat("\nTHEORY TESTS COMPLETE\n")
