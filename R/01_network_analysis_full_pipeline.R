# ==========================================================================
# IsoFresh network analysis -- FULL PIPELINE (single file, all phases)
#
# Runs everything in one `Rscript` call: network/null-model computation,
# omnivory v2, threshold sensitivity, the statistical battery + Fig0-5,
# and the bootstrap mediation + piecewiseSEM (run for BOTH Temperature and
# Latitude as the exogenous driver, per the 2026-09-01 request to foreground
# Latitude as the primary, sampling-time-independent integrative predictor,
# since WorldClim's Mean_Annual_Temp is a 30-year climatological normal with
# no record of each pond's actual sampling month/season -- IsoFresh has no
# such field, confirmed by inspection, so this is a real, uncorrectable
# limitation of using it as a precise per-sample "cause").
#
# CHECKPOINTING: phases 1-3 are the only slow-ish/expensive ones (Phase 1 is
# the big one, likely 1.5-3h; Phase 5's bootstrap is the other big one,
# likely 1.5-3h). Each checks whether its output file already exists and
# skips recomputation unless you force it:
#   ISOFRESH_FORCE_RECOMPUTE=1   forces phases 1-3 to redo even if present
#   ISOFRESH_FORCE_REBOOT=1      forces the Phase 5 bootstrap to redo
#   ISOFRESH_TEST_N=<n>          Phase 1 smoke-test on n ponds only
#   ISOFRESH_N_SIMUL=<n>         Phase 1 null-model replicates (default 199)
#   ISOFRESH_N_BOOT=<n>          Phase 5 bootstrap replicates (default 999)
# So a first full run does everything (long); a later re-run (e.g. after
# tweaking a figure) reuses the saved RData/CSVs from 1-3 and Phase 5's
# bootstrap, and just redoes the fast parts (Phase 4 figures, Phase 5 SEM).
# ==========================================================================

suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(igraph); library(rnetcarto)
  library(vegan); library(parallel); library(geodata); library(terra); library(rfishbase)
  library(ggplot2); library(lme4); library(lmerTest); library(cowplot); library(corrplot)
  library(maps); library(viridis); library(piecewiseSEM); library(car)
})
# NOTE: car::recode() masks dplyr::recode() -- every recode() call below is
# namespaced dplyr::recode() to avoid silently picking up the wrong one.

FORCE_RECOMPUTE <- Sys.getenv("ISOFRESH_FORCE_RECOMPUTE", "0") == "1"
FORCE_REBOOT    <- Sys.getenv("ISOFRESH_FORCE_REBOOT", "0") == "1"
TEST_N   <- as.integer(Sys.getenv("ISOFRESH_TEST_N", "0"))
N_SIMUL  <- as.integer(Sys.getenv("ISOFRESH_N_SIMUL", "199"))
N_BOOT   <- as.integer(Sys.getenv("ISOFRESH_N_BOOT", "999"))
N_CORES  <- max(1, parallel::detectCores() - 1)
SEED     <- 123
zscore   <- function(x) as.numeric(scale(x))
dir.create("figures_v2", showWarnings = FALSE)

COL_LENTIC <- "#2a78d6"; COL_LOTIC <- "#eb6834"
DIV_LOW <- "#1c5cab"; DIV_MID <- "#f0efec"; DIV_HIGH <- "#e34948"
theme_set(theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.2),
                                            strip.background = element_rect(fill = "grey92", color = NA),
                                            plot.title = element_text(size = 10, hjust = 0.5), legend.position = "bottom"))

# ==========================================================================
# PHASE 1: Network metrics (expensive: netcarto + null models)
# ==========================================================================
cat("\n================================================================\nPHASE 1: Network metrics\n================================================================\n")
if (!file.exists("IsoFresh_Metrics_v2.RData") || FORCE_RECOMPUTE) {
local({

fish_data  <- read_delim("IsoFresh_Isotope_Fish.csv",  delim = ";", show_col_types = FALSE)
other_data <- read_delim("IsoFresh_Isotope_Other.csv", delim = ";", show_col_types = FALSE)
eco_data   <- read_delim("IsoFresh_Ecosystem.csv",     delim = ";", show_col_types = FALSE)
fw_data    <- read_delim("IsoFresh_FoodWeb.csv",       delim = ";", show_col_types = FALSE)

# 38% of "Other" and 62/8036 fish records have NA Scientific_name; using it
# alone as node ID both pools distinct unidentified taxa into one phantom
# node (dplyr groups all NAs together) and crashes rnetcarto's C code on a
# literal NA dimname. Common_name is populated for ~100% of these rows.
fish_sub  <- fish_data  %>% mutate(GIG = coalesce(Scientific_name, Common_name, paste0("Unidentified_", row_number()))) %>%
  select(FoodWeb_ID, GIG, d13C, d15N) %>% mutate(Type = "Fish")
other_sub <- other_data %>% mutate(GIG = coalesce(Scientific_name, Common_name, paste0("Unidentified_", row_number()))) %>%
  select(FoodWeb_ID, GIG, d13C, d15N) %>% mutate(Type = "Other")
isotop_all <- bind_rows(fish_sub, other_sub) %>% filter(!is.na(d13C) & !is.na(d15N))
mean_isotop <- isotop_all %>% group_by(FoodWeb_ID, GIG, Type) %>%
  summarize(d15N = mean(d15N), d13C = mean(d13C), .groups = "drop")

# Resource ("prey") space per pond
chull_area <- function(x, y) {
  if (length(x) < 3) return(NA_real_)
  h <- chull(x, y); xh <- x[h]; yh <- y[h]
  0.5 * abs(sum(xh * c(yh[-1], yh[1]) - c(xh[-1], xh[1]) * yh))
}
prey_space_df <- mean_isotop %>% filter(Type == "Other") %>% group_by(FoodWeb_ID) %>%
  summarize(N_Prey = n(), Prey_Niche_Area = chull_area(d13C, d15N), .groups = "drop")

# Site-level covariates: temperature (WorldClim), lake area (HydroLAKES), Source_DOI
clim_dir <- "worldclim_data"; dir.create(clim_dir, showWarnings = FALSE)
bio <- worldclim_global(var = "bio", res = 10, path = clim_dir)
pts <- eco_data %>% select(FoodWeb_ID, Longitude, Latitude)
temp_vals <- terra::extract(bio[[c(1, 4, 12)]], pts[, c("Longitude", "Latitude")])
names(temp_vals) <- c("ID", "Mean_Annual_Temp", "Temp_Seasonality", "Annual_Precip")
eco_data <- bind_cols(eco_data, temp_vals %>% select(-ID))
eco_data$Temp_Seasonality <- eco_data$Temp_Seasonality / 100

lakes_dir <- "hydrolakes_data"
lakes_zip <- file.path(lakes_dir, "HydroLAKES_polys_v10_shp.zip")
lakes_shp <- file.path(lakes_dir, "HydroLAKES_polys_v10_shp", "HydroLAKES_polys_v10.shp")
dir.create(lakes_dir, showWarnings = FALSE)
if (!file.exists(lakes_shp) && !file.exists(lakes_zip)) {
  message("Downloading HydroLAKES v10 (~820MB, one-time)...")
  tryCatch(download.file("https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip",
                          lakes_zip, mode = "wb", quiet = TRUE),
           error = function(e) warning("HydroLAKES download failed: ", conditionMessage(e)))
}
if (!file.exists(lakes_shp) && file.exists(lakes_zip)) unzip(lakes_zip, exdir = lakes_dir)
if (file.exists(lakes_shp)) {
  lakes <- vect(lakes_shp)
  eco_pts <- vect(eco_data, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
  lake_hit <- terra::extract(lakes[, "Lake_area"], eco_pts)$Lake_area
  miss <- which(is.na(lake_hit) & eco_data$Ecosystem %in% c("Lentic", "Both"))
  if (length(miss) > 0) {
    nn <- nearest(eco_pts[miss], lakes)
    within_2km <- nn$distance <= 2000
    lake_area_vals <- values(lakes)$Lake_area
    lake_hit[miss[within_2km]] <- lake_area_vals[nn$to_id[within_2km]]
  }
  lake_hit[eco_data$Ecosystem == "Lotic"] <- NA
  eco_data$Lake_Area_km2 <- lake_hit
  cat(sprintf("Lake area: %d/%d Lentic/Both sites matched to HydroLAKES\n",
              sum(!is.na(eco_data$Lake_Area_km2[eco_data$Ecosystem != "Lotic"])), sum(eco_data$Ecosystem != "Lotic")))
} else {
  eco_data$Lake_Area_km2 <- NA_real_
  warning("HydroLAKES shapefile not found -- Lake_Area_km2 will be all NA")
}

study_lookup <- fw_data %>% select(FoodWeb_ID, Source_DOI) %>% distinct(FoodWeb_ID, .keep_all = TRUE)
eco_data <- eco_data %>% left_join(study_lookup, by = "FoodWeb_ID")
mean_isotop <- mean_isotop %>%
  left_join(eco_data %>% select(FoodWeb_ID, Latitude, Longitude, Ecosystem, Region,
                                 Mean_Annual_Temp, Temp_Seasonality, Annual_Precip, Lake_Area_km2, Source_DOI),
            by = "FoodWeb_ID")
cat(sprintf("Temperature extraction: %d/%d sites\n", sum(!is.na(eco_data$Mean_Annual_Temp)), nrow(eco_data)))

# Fish body-size structure + Omnivory Index proxy v1
codes <- fish_data$FishBase_code %>% unique() %>% na.omit()
fb <- tryCatch(rfishbase::species(SpecCode = codes) %>% select(SpecCode, Length),
               error = function(e) { warning("rfishbase lookup failed: ", conditionMessage(e)); NULL })
diet_tbl <- tryCatch(rfishbase::diet(SpecCode = codes) %>% filter(!is.na(seTroph)) %>%
                        group_by(SpecCode) %>% summarize(Omnivory_Index = mean(seTroph), .groups = "drop"),
                      error = function(e) { warning("diet() lookup failed: ", conditionMessage(e)); NULL })
if (!is.null(fb)) {
  fish_traits <- fish_data %>% select(FoodWeb_ID, FishBase_code) %>%
    left_join(fb, by = c("FishBase_code" = "SpecCode")) %>%
    left_join(diet_tbl, by = c("FishBase_code" = "SpecCode"))
  bodysize_df <- fish_traits %>% filter(!is.na(Length) & Length > 0) %>% group_by(FoodWeb_ID) %>%
    summarize(n_fish_with_length = n(), Mean_log10_Length = mean(log10(Length)), SD_log10_Length = sd(log10(Length)),
              BodySize_Ratio = max(Length) / min(Length),
              BodySize_MaxGap = { lg <- sort(log10(Length)); if (length(lg) >= 2 && diff(range(lg)) > 0) max(diff(lg)) / diff(range(lg)) else NA_real_ },
              .groups = "drop")
  omniv_df <- fish_traits %>% filter(!is.na(Omnivory_Index)) %>% group_by(FoodWeb_ID) %>%
    summarize(Mean_Omnivory_Index = mean(Omnivory_Index), .groups = "drop")
  bodysize_df <- bodysize_df %>% full_join(omniv_df, by = "FoodWeb_ID")
  cat(sprintf("Body-size: %d/%d ponds; Omnivory v1: %d ponds\n",
              sum(!is.na(bodysize_df$Mean_log10_Length)), length(unique(mean_isotop$FoodWeb_ID)),
              sum(!is.na(bodysize_df$Mean_Omnivory_Index))))
} else {
  bodysize_df <- tibble(FoodWeb_ID = character(), n_fish_with_length = integer(), Mean_log10_Length = double(),
                         SD_log10_Length = double(), BodySize_Ratio = double(), BodySize_MaxGap = double(),
                         Mean_Omnivory_Index = double())
}

pond_list <- split(mean_isotop, mean_isotop$FoodWeb_ID)
pond_list <- pond_list[sapply(pond_list, nrow) >= 3]
cat(sprintf("%d/%d ponds retained (S>=3)\n", length(pond_list), length(unique(mean_isotop$FoodWeb_ID))))
if (TEST_N > 0) {
  set.seed(SEED)
  pond_list <- pond_list[sample(names(pond_list), min(TEST_N, length(pond_list)))]
  cat(sprintf("[TEST MODE] restricted to %d ponds\n", length(pond_list)))
}
sim_list <- lapply(pond_list, function(sub_data) {
  dd <- as.matrix(dist(sub_data[, c("d13C", "d15N")]))
  rownames(dd) <- colnames(dd) <- sub_data$GIG
  max_d <- max(dd, na.rm = TRUE); if (max_d == 0) max_d <- 1
  sim <- 1 - dd / max_d; diag(sim) <- 0; sim
})
pooled_sim <- unlist(lapply(sim_list, function(m) m[upper.tri(m)]))
GLOBAL_THRESHOLD <- mean(pooled_sim[pooled_sim > 0])
cat(sprintf("Global absolute similarity threshold: %.4f\n", GLOBAL_THRESHOLD))
saveRDS(pooled_sim, "pooled_similarity_values.rds")

perm_sym <- function(m) {
  n <- nrow(m); ut <- upper.tri(m); vals <- m[ut]
  m2 <- matrix(0, n, n); m2[ut] <- sample(vals); m2 <- m2 + t(m2); dimnames(m2) <- dimnames(m); m2
}
iterfac_obs <- function(S) max(10, min(100, round(2000 / S)))
ITERFAC_NULL <- 5
n_simul_for <- function(S, n_simul) if (S <= 30) n_simul else if (S <= 50) max(99, round(n_simul / 2)) else max(49, round(n_simul / 4))
strip_names <- function(m) { dimnames(m) <- NULL; m }

compute_pond_metrics <- function(id, sim, threshold = GLOBAL_THRESHOLD, n_simul = N_SIMUL, seed = SEED) {
  set.seed(seed); S <- nrow(sim); n_simul <- n_simul_for(S, n_simul); real_names <- rownames(sim)
  simf <- sim; simf[simf < threshold] <- 0
  g <- graph_from_adjacency_matrix(simf, mode = "undirected", weighted = TRUE)
  L <- length(E(g)); conn_val <- L / (S * (S - 1) / 2)
  ifac <- iterfac_obs(S)
  mod_result <- tryCatch(netcarto(strip_names(simf), seed = seed, iterfac = ifac, bipartite = FALSE), error = function(e) NULL)
  q_val <- if (is.null(mod_result)) NA else mod_result[[2]]
  node_df <- if (!is.null(mod_result)) { dd <- mod_result[[1]]; dd$name <- real_names[as.integer(as.character(dd$name))]; dd$FoodWeb_ID <- id; dd } else NULL
  str_val <- mean(strength(g, weights = E(g)$weight), na.rm = TRUE)
  eig_val <- mean(eigen_centrality(g, weights = E(g)$weight)$vector, na.rm = TRUE)
  dist_filtered <- (1 - sim); dist_filtered[simf == 0] <- 0
  g_dist <- graph_from_adjacency_matrix(dist_filtered, mode = "undirected", weighted = TRUE)
  clo_val <- mean(closeness(g_dist, weights = E(g_dist)$weight, normalized = TRUE), na.rm = TRUE)
  bet_val <- mean(betweenness(g_dist, weights = E(g_dist)$weight, normalized = TRUE), na.rm = TRUE)
  nodf_obs <- tryCatch(nestednodf(strip_names(simf), order = TRUE, weighted = TRUE)$statistic["NODF"], error = function(e) NA)
  nodf_null <- replicate(n_simul, { mp <- perm_sym(sim); mp[mp < threshold] <- 0
    tryCatch(nestednodf(strip_names(mp), order = TRUE, weighted = TRUE)$statistic["NODF"], error = function(e) NA) })
  z_nodf <- (nodf_obs - mean(nodf_null, na.rm = TRUE)) / sd(nodf_null, na.rm = TRUE)
  q_null <- replicate(n_simul, { mp <- perm_sym(sim); mp[mp < threshold] <- 0
    tryCatch(suppressWarnings(netcarto(strip_names(mp), bipartite = FALSE, iterfac = ITERFAC_NULL)[[2]]), error = function(e) NA) })
  z_q <- (q_val - mean(q_null, na.rm = TRUE)) / sd(q_null, na.rm = TRUE)
  list(metrics = data.frame(FoodWeb_ID = id, S = S, Connectance = conn_val, Modularity = q_val, z_Q = z_q, z_NODF = z_nodf,
                             Mean_Strength = str_val, Mean_Eigen = eig_val, Mean_Closeness = clo_val, Mean_Betweenness = bet_val),
       node = node_df)
}

cat(sprintf("Running %d ponds on %d cores, n_simul=%d (capped for large S) ...\n", length(pond_list), N_CORES, N_SIMUL))
t0 <- Sys.time()
ordered_ids <- names(sim_list)[order(-sapply(sim_list, nrow))]
results <- mclapply(ordered_ids, function(id) tryCatch(compute_pond_metrics(id, sim_list[[id]]), error = function(e) NULL),
                     mc.cores = N_CORES, mc.preschedule = FALSE)
names(results) <- ordered_ids
cat(sprintf("Done in %.1f minutes\n", as.numeric(Sys.time() - t0, units = "mins")))
failed <- names(results)[sapply(results, is.null)]
if (length(failed) > 0) cat(sprintf("%d ponds failed\n", length(failed)))

metrics_results <- bind_rows(lapply(results, function(r) if (!is.null(r)) r$metrics))
node_metrics_results <- bind_rows(lapply(results, function(r) if (!is.null(r)) r$node))
n_modules_df <- node_metrics_results %>% group_by(FoodWeb_ID) %>% summarize(N_Modules = n_distinct(module), .groups = "drop")
richness_df <- mean_isotop %>% group_by(FoodWeb_ID) %>% summarize(Richness = n(), .groups = "drop")

final_df <- metrics_results %>%
  left_join(eco_data %>% select(FoodWeb_ID, Latitude, Longitude, Ecosystem, Region, Source_DOI,
                                 Mean_Annual_Temp, Temp_Seasonality, Annual_Precip, Lake_Area_km2), by = "FoodWeb_ID") %>%
  left_join(richness_df, by = "FoodWeb_ID") %>% left_join(prey_space_df, by = "FoodWeb_ID") %>%
  left_join(bodysize_df, by = "FoodWeb_ID") %>% left_join(n_modules_df, by = "FoodWeb_ID") %>%
  filter(!is.na(Latitude)) %>%
  mutate(Abs_Latitude = abs(Latitude), Ecosystem_Grp = ifelse(Ecosystem == "Both", "Lentic", Ecosystem), Include_S10 = S >= 10)

save(final_df, node_metrics_results, GLOBAL_THRESHOLD, pooled_sim, file = "IsoFresh_Metrics_v2.RData")
cat(sprintf("Saved IsoFresh_Metrics_v2.RData: N=%d (%d pass S>=10)\n", nrow(final_df), sum(final_df$Include_S10)))
})
} else {
  cat("SKIPPED -- IsoFresh_Metrics_v2.RData exists (ISOFRESH_FORCE_RECOMPUTE=1 to redo)\n")
}
load("IsoFresh_Metrics_v2.RData")

# ==========================================================================
# PHASE 2: Omnivory Index v2 (FishBase ecology(), better coverage)
# ==========================================================================
cat("\n================================================================\nPHASE 2: Omnivory Index v2\n================================================================\n")
if (!file.exists("IsoFresh_Omnivory_v2.RData") || FORCE_RECOMPUTE) {
local({
  fish_data <- read_delim("IsoFresh_Isotope_Fish.csv", delim = ";", show_col_types = FALSE)
  codes <- fish_data$FishBase_code %>% unique() %>% na.omit()
  eco_tbl <- tryCatch(rfishbase::ecology() %>% filter(SpecCode %in% codes) %>%
                         select(SpecCode, FoodTroph, FoodSeTroph, DietTroph, DietSeTroph) %>%
                         group_by(SpecCode) %>%
                         summarize(across(c(FoodTroph, FoodSeTroph, DietTroph, DietSeTroph), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
                         mutate(across(everything(), ~ ifelse(is.nan(.), NA, .))),
                       error = function(e) { warning("ecology() lookup failed: ", conditionMessage(e)); NULL })
  if (is.null(eco_tbl)) {
    omnivory_v2_df <- tibble(FoodWeb_ID = character(), Omnivory_FoodSeTroph = double(), Omnivory_DietSeTroph = double())
  } else {
    fish_omniv <- fish_data %>% select(FoodWeb_ID, FishBase_code) %>% left_join(eco_tbl, by = c("FishBase_code" = "SpecCode"))
    omnivory_v2_df <- fish_omniv %>% group_by(FoodWeb_ID) %>%
      summarize(n_fish_FoodSeTroph = sum(!is.na(FoodSeTroph)), Omnivory_FoodSeTroph = mean(FoodSeTroph, na.rm = TRUE),
                n_fish_DietSeTroph = sum(!is.na(DietSeTroph)), Omnivory_DietSeTroph = mean(DietSeTroph, na.rm = TRUE), .groups = "drop") %>%
      mutate(across(c(Omnivory_FoodSeTroph, Omnivory_DietSeTroph), ~ ifelse(is.nan(.), NA, .)))
    cat(sprintf("Omnivory v2: %d/%d ponds FoodSeTroph; %d DietSeTroph\n",
                sum(omnivory_v2_df$n_fish_FoodSeTroph > 0), nrow(omnivory_v2_df), sum(omnivory_v2_df$n_fish_DietSeTroph > 0)))
  }
  save(omnivory_v2_df, file = "IsoFresh_Omnivory_v2.RData")
})
} else {
  cat("SKIPPED -- IsoFresh_Omnivory_v2.RData exists\n")
}
load("IsoFresh_Omnivory_v2.RData")
final_df <- final_df %>% left_join(omnivory_v2_df, by = "FoodWeb_ID")
have_omniv_v2 <- TRUE

# ==========================================================================
# PHASE 3: Threshold-sensitivity sweep
# ==========================================================================
cat("\n================================================================\nPHASE 3: Threshold sensitivity\n================================================================\n")
if (!file.exists("IsoFresh_Metrics_v2_thresholdsweep.RData") || FORCE_RECOMPUTE) {
local({
  fish_data  <- read_delim("IsoFresh_Isotope_Fish.csv",  delim = ";", show_col_types = FALSE)
  other_data <- read_delim("IsoFresh_Isotope_Other.csv", delim = ";", show_col_types = FALSE)
  eco_data   <- read_delim("IsoFresh_Ecosystem.csv",     delim = ";", show_col_types = FALSE)
  fw_data    <- read_delim("IsoFresh_FoodWeb.csv",       delim = ";", show_col_types = FALSE)
  fish_sub  <- fish_data  %>% mutate(GIG = coalesce(Scientific_name, Common_name, paste0("Unidentified_", row_number()))) %>%
    select(FoodWeb_ID, GIG, d13C, d15N) %>% mutate(Type = "Fish")
  other_sub <- other_data %>% mutate(GIG = coalesce(Scientific_name, Common_name, paste0("Unidentified_", row_number()))) %>%
    select(FoodWeb_ID, GIG, d13C, d15N) %>% mutate(Type = "Other")
  mean_isotop <- bind_rows(fish_sub, other_sub) %>% filter(!is.na(d13C) & !is.na(d15N)) %>%
    group_by(FoodWeb_ID, GIG, Type) %>% summarize(d15N = mean(d15N), d13C = mean(d13C), .groups = "drop")
  pond_list <- split(mean_isotop, mean_isotop$FoodWeb_ID); pond_list <- pond_list[sapply(pond_list, nrow) >= 3]
  sim_list <- lapply(pond_list, function(sub_data) {
    dd <- as.matrix(dist(sub_data[, c("d13C", "d15N")])); max_d <- max(dd, na.rm = TRUE); if (max_d == 0) max_d <- 1
    sim <- 1 - dd / max_d; diag(sim) <- 0; sim
  })
  sweep_pooled_sim <- unlist(lapply(sim_list, function(m) m[upper.tri(m)]))
  have_temp <- FALSE
  if (dir.exists("worldclim_data")) {
    bio <- tryCatch(worldclim_global(var = "bio", res = 10, path = "worldclim_data"), error = function(e) NULL)
    if (!is.null(bio)) { eco_data$Mean_Annual_Temp <- terra::extract(bio[[1]], eco_data[, c("Longitude", "Latitude")])[, 2]; have_temp <- TRUE }
  }
  study_lookup <- fw_data %>% select(FoodWeb_ID, Source_DOI) %>% distinct(FoodWeb_ID, .keep_all = TRUE)
  eco_data <- eco_data %>% left_join(study_lookup, by = "FoodWeb_ID")
  thresholds <- c(quantile(sweep_pooled_sim[sweep_pooled_sim > 0], c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)), mean = mean(sweep_pooled_sim[sweep_pooled_sim > 0]))
  richness <- sapply(sim_list, nrow)
  sweep_rows <- lapply(names(thresholds), function(th_name) {
    th <- thresholds[[th_name]]
    conn <- sapply(sim_list, function(m) { S <- nrow(m); L <- sum(m >= th) / 2; L / (S * (S - 1) / 2) })
    dd <- tibble(FoodWeb_ID = names(sim_list), Connectance = conn, S = richness) %>%
      left_join(eco_data %>% select(FoodWeb_ID, Latitude, Region, Source_DOI, any_of("Mean_Annual_Temp")), by = "FoodWeb_ID") %>%
      filter(S >= 10, !is.na(Region), !is.na(Source_DOI))
    dd$Abs_Latitude <- abs(dd$Latitude)
    pred <- if (have_temp && "Mean_Annual_Temp" %in% names(dd) && !all(is.na(dd$Mean_Annual_Temp))) "Mean_Annual_Temp" else "Abs_Latitude"
    dd <- dd %>% filter(!is.na(.data[[pred]]))
    dd$pred_z <- as.numeric(scale(dd[[pred]])); dd$conn_z <- as.numeric(scale(dd$Connectance))
    m <- tryCatch(lmer(conn_z ~ pred_z + (1 | Region/Source_DOI), data = dd), error = function(e) NULL)
    if (is.null(m)) return(tibble(threshold = th, threshold_label = th_name, term = pred, estimate = NA, p.value = NA, n = nrow(dd)))
    s <- summary(m)$coefficients
    tibble(threshold = th, threshold_label = th_name, term = pred, estimate = s["pred_z", "Estimate"], p.value = s["pred_z", "Pr(>|t|)"], n = nrow(dd))
  })
  threshold_sweep_df <- bind_rows(sweep_rows)
  save(threshold_sweep_df, file = "IsoFresh_Metrics_v2_thresholdsweep.RData")
  cat("Saved IsoFresh_Metrics_v2_thresholdsweep.RData\n")
})
} else {
  cat("SKIPPED -- threshold sweep RData exists\n")
}

# ==========================================================================
# PHASE 4: Statistical battery + Fig0-5
# ==========================================================================
cat("\n================================================================\nPHASE 4: Statistical battery + figures\n================================================================\n")

if (file.exists("IsoFresh_Metrics_v2_thresholdsweep.RData")) {
  load("IsoFresh_Metrics_v2_thresholdsweep.RData")
  p_sweep <- ggplot(threshold_sweep_df, aes(x = threshold, y = estimate)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
    geom_line(color = COL_LENTIC) + geom_point(aes(color = p.value < 0.05), size = 2) +
    scale_color_manual(values = c(`TRUE` = COL_LOTIC, `FALSE` = "grey60"), name = "p<0.05") +
    labs(x = "Global similarity threshold", y = "Standardized slope: Connectance ~ predictor",
         title = "Robustness of the Connectance-climate effect to threshold choice")
  ggsave("figures_v2/threshold_sensitivity.png", p_sweep, width = 6, height = 4, dpi = 300)
}

fit_lmm <- function(df, metric, pred, extra_terms = NULL) {
  dd <- df %>% filter(!is.na(.data[[metric]]), !is.na(.data[[pred]]), !is.na(Region), !is.na(Source_DOI))
  if (nrow(dd) < 15 || n_distinct(dd$Region) < 2) return(NULL)
  rhs <- paste(c(pred, extra_terms, "(1 | Region/Source_DOI)"), collapse = " + ")
  tryCatch(lmer(as.formula(paste(metric, "~", rhs)), data = dd, REML = TRUE), error = function(e) NULL)
}
tidy_lmm <- function(model, pred) {
  if (is.null(model)) return(tibble(estimate = NA, std.error = NA, p.value = NA, n = NA))
  s <- summary(model)$coefficients
  if (!pred %in% rownames(s)) return(tibble(estimate = NA, std.error = NA, p.value = NA, n = NA))
  tibble(estimate = s[pred, "Estimate"], std.error = s[pred, "Std. Error"], p.value = s[pred, "Pr(>|t|)"], n = nobs(model))
}

METRICS <- c("Connectance", "Modularity", "z_Q", "z_NODF", "Mean_Strength", "Mean_Closeness")
final_df <- final_df %>% mutate(Log_Lake_Area = log10(Lake_Area_km2 + 1), Log_Richness = log10(Richness), Log_N_Prey = log10(N_Prey + 1))
PREDICTOR_META <- tibble::tribble(
  ~predictor,          ~label,
  "Abs_Latitude",      "Absolute latitude",
  "Mean_Annual_Temp",  "Mean annual temperature (deg C)",
  "Temp_Seasonality",  "Temperature seasonality (SD of monthly temp, deg C)",
  "Log_Lake_Area",     "Lake area (log10 km2 + 1) -- ecosystem size",
  "Log_N_Prey",        "Prey/baseline taxon count (log10 + 1) -- prey number",
  "Prey_Niche_Area",   "Prey isotopic niche area -- prey space",
  "Mean_Omnivory_Index","Omnivory Index proxy (FishBase diet seTroph)",
  "BodySize_Ratio",    "Body-size ratio (max/min fish length)",
  "BodySize_MaxGap",   "Body-size max gap (log-size continuity)",
  "Log_Richness",      "Richness (log10 species count)",
  "N_Modules",         "Guild count (netcarto modules)"
)
if (have_omniv_v2) {
  PREDICTOR_META <- bind_rows(PREDICTOR_META, tibble::tribble(
    ~predictor,               ~label,
    "Omnivory_FoodSeTroph",  "Omnivory Index v2 (FishBase FoodSeTroph, ~64% coverage)",
    "Omnivory_DietSeTroph",  "Omnivory Index v3 (FishBase DietSeTroph, diet studies only)"
  ))
}
PREDICTORS <- PREDICTOR_META$predictor
POPULATIONS <- list(full = final_df, S10 = filter(final_df, Include_S10))
ECOSYSTEMS <- c("All", "Lentic", "Lotic")

run_battery <- function(pop_name, pop_df) {
  rows <- list()
  for (eco in ECOSYSTEMS) {
    d0 <- if (eco == "All") pop_df else filter(pop_df, Ecosystem_Grp == eco)
    for (metric in METRICS) for (pred in PREDICTORS) {
      dd <- d0
      if (!metric %in% names(dd) || !pred %in% names(dd)) next
      dd[[paste0(metric, "_z")]] <- zscore(dd[[metric]]); dd[[paste0(pred, "_z")]] <- zscore(dd[[pred]])
      m <- fit_lmm(dd, paste0(metric, "_z"), paste0(pred, "_z"))
      rows[[length(rows) + 1]] <- tibble(population = pop_name, ecosystem = eco, metric = metric, predictor = pred,
                                          !!!tidy_lmm(m, paste0(pred, "_z")))
    }
  }
  bind_rows(rows)
}
battery <- bind_rows(lapply(names(POPULATIONS), function(nm) run_battery(nm, POPULATIONS[[nm]]))) %>%
  filter(!is.na(estimate)) %>% mutate(q.value = p.adjust(p.value, method = "BH"))
write.csv(battery, "figures_v2/model_battery_results.csv", row.names = FALSE)
cat(sprintf("Model battery: %d models; %d sig raw p<0.05; %d sig BH-FDR\n", nrow(battery), sum(battery$p.value < 0.05), sum(battery$q.value < 0.05)))

THERMAL_PREDICTORS <- c("Abs_Latitude", "Mean_Annual_Temp")

mediation_results <- list()
for (nm in names(POPULATIONS)) for (metric in METRICS) for (pred in THERMAL_PREDICTORS) {
  dd <- POPULATIONS[[nm]] %>% filter(!is.na(.data[[metric]]), !is.na(.data[[pred]]), !is.na(Log_Richness), !is.na(Prey_Niche_Area))
  if (nrow(dd) < 20) next
  dd$m_z <- zscore(dd[[metric]]); dd$p_z <- zscore(dd[[pred]]); dd$S_z <- zscore(dd$Log_Richness); dd$PA_z <- zscore(dd$Prey_Niche_Area)
  b0 <- tidy_lmm(fit_lmm(dd, "m_z", "p_z"), "p_z")
  b2 <- tidy_lmm(fit_lmm(dd, "m_z", "p_z", extra_terms = "S_z + PA_z"), "p_z")
  mediation_results[[length(mediation_results) + 1]] <- tibble(population = nm, metric = metric, predictor = pred, n = nrow(dd),
    beta_alone = b0$estimate, p_alone = b0$p.value, beta_with_Richness_PreyArea = b2$estimate, p_with_Richness_PreyArea = b2$p.value,
    pct_attenuation = 100 * (1 - abs(b2$estimate) / abs(b0$estimate)))
}
mediation_results <- bind_rows(mediation_results)
write.csv(mediation_results, "figures_v2/mediation_analysis.csv", row.names = FALSE)

prey_area_gradient <- bind_rows(lapply(names(POPULATIONS), function(nm) bind_rows(lapply(THERMAL_PREDICTORS, function(pred) {
  dd <- POPULATIONS[[nm]] %>% filter(!is.na(Prey_Niche_Area), !is.na(.data[[pred]]))
  dd$pa_z <- zscore(dd$Prey_Niche_Area); dd$p_z <- zscore(dd[[pred]])
  tibble(population = nm, predictor = pred, !!!tidy_lmm(fit_lmm(dd, "pa_z", "p_z"), "p_z"))
}))))
write.csv(prey_area_gradient, "figures_v2/prey_niche_area_gradient.csv", row.names = FALSE)

scaling_results <- bind_rows(lapply(names(POPULATIONS), function(nm) {
  dd <- POPULATIONS[[nm]] %>% filter(!is.na(Connectance), !is.na(S), S >= 3)
  dd$L <- dd$Connectance * dd$S * (dd$S - 1) / 2; dd <- dd %>% filter(L > 0)
  dd$logL <- log10(dd$L); dd$logS <- log10(dd$S)
  base_slope <- tidy_lmm(fit_lmm(dd, "logL", "logS"), "logS")$estimate
  bind_rows(lapply(THERMAL_PREDICTORS, function(pred) {
    d2 <- dd %>% filter(!is.na(.data[[pred]])); d2$pred_z <- zscore(d2[[pred]])
    m <- fit_lmm(d2, "logL", "logS", extra_terms = "pred_z + logS:pred_z")
    if (is.null(m)) return(tibble(population = nm, predictor = pred, base_scaling_exponent = base_slope, interaction_est = NA, interaction_p = NA, n = nrow(d2)))
    s <- summary(m)$coefficients; int_term <- grep(":", rownames(s), value = TRUE)
    tibble(population = nm, predictor = pred, base_scaling_exponent = base_slope,
           interaction_est = if (length(int_term)) s[int_term, "Estimate"] else NA,
           interaction_p = if (length(int_term)) s[int_term, "Pr(>|t|)"] else NA, n = nobs(m))
  }))
}))
write.csv(scaling_results, "figures_v2/link_species_scaling.csv", row.names = FALSE)

bodysize_predictors <- c(SDsize = "SD_log10_Length", Ratio = "BodySize_Ratio")
bodysize_mod <- list()
for (nm in names(POPULATIONS)) for (eco in c("All", "Lentic", "Lotic")) {
  d0 <- POPULATIONS[[nm]]; if (eco != "All") d0 <- filter(d0, Ecosystem_Grp == eco)
  for (sz_label in names(bodysize_predictors)) {
    sz_col <- bodysize_predictors[[sz_label]]
    dd <- d0 %>% filter(!is.na(Connectance), !is.na(Mean_Annual_Temp), !is.na(.data[[sz_col]]), n_fish_with_length >= 3)
    if (nrow(dd) < 20) next
    dd$Temp_z <- zscore(dd$Mean_Annual_Temp); dd$Size_z <- zscore(dd[[sz_col]])
    m <- fit_lmm(dd, "Connectance", "Temp_z", extra_terms = "Size_z + Temp_z:Size_z")
    if (is.null(m)) next
    s <- summary(m)$coefficients; int_term <- grep(":", rownames(s), value = TRUE)
    bodysize_mod[[length(bodysize_mod) + 1]] <- tibble(population = nm, ecosystem = eco, size_metric = sz_label, n = nobs(m),
      interaction_est = if (length(int_term)) s[int_term, "Estimate"] else NA, interaction_p = if (length(int_term)) s[int_term, "Pr(>|t|)"] else NA)
  }
}
# bind_rows(list()) on an empty accumulator (e.g. too few ponds pass the
# nrow>=20 guard in every population x ecosystem x size-metric slice, as
# happens on small test/subset runs) yields a 0-column tibble; mutate() then
# fails with "object not found" since the column it adjusts never existed.
bodysize_mod <- bind_rows(bodysize_mod)
if (nrow(bodysize_mod) > 0) bodysize_mod <- bodysize_mod %>% mutate(interaction_q = p.adjust(interaction_p, method = "BH"))
write.csv(bodysize_mod, "figures_v2/bodysize_temperature_interaction.csv", row.names = FALSE)

# ---- Temperature seasonality: main effect (already covered via battery/Fig0
# since Temp_Seasonality is now in PREDICTOR_META) PLUS the more novel test --
# does climate VARIABILITY moderate the mean-temperature~metric relationship?
# (i.e. is the thermal sensitivity of network structure itself a function of
# how seasonal the site's climate is -- relevant to climate-change framing,
# since variability/seasonality shifts are as much a part of warming as the
# mean). Tested across all 6 metrics x 3 ecosystem splits x 2 populations,
# matching the body-size-interaction design above so it gets the same
# FDR treatment.
seasonality_interaction <- list()
for (nm in names(POPULATIONS)) for (eco in c("All", "Lentic", "Lotic")) {
  d0 <- POPULATIONS[[nm]]; if (eco != "All") d0 <- filter(d0, Ecosystem_Grp == eco)
  for (metric in METRICS) {
    dd <- d0 %>% filter(!is.na(.data[[metric]]), !is.na(Mean_Annual_Temp), !is.na(Temp_Seasonality))
    if (nrow(dd) < 20) next
    dd$m_z <- zscore(dd[[metric]]); dd$Temp_z <- zscore(dd$Mean_Annual_Temp); dd$Seas_z <- zscore(dd$Temp_Seasonality)
    m <- fit_lmm(dd, "m_z", "Temp_z", extra_terms = "Seas_z + Temp_z:Seas_z")
    if (is.null(m)) next
    s <- summary(m)$coefficients; int_term <- grep(":", rownames(s), value = TRUE)
    seasonality_interaction[[length(seasonality_interaction) + 1]] <- tibble(
      population = nm, ecosystem = eco, metric = metric, n = nobs(m),
      main_seasonality_est = if ("Seas_z" %in% rownames(s)) s["Seas_z", "Estimate"] else NA,
      main_seasonality_p = if ("Seas_z" %in% rownames(s)) s["Seas_z", "Pr(>|t|)"] else NA,
      interaction_est = if (length(int_term)) s[int_term, "Estimate"] else NA,
      interaction_p = if (length(int_term)) s[int_term, "Pr(>|t|)"] else NA)
  }
}
seasonality_interaction <- bind_rows(seasonality_interaction)
if (nrow(seasonality_interaction) > 0) {
  seasonality_interaction <- seasonality_interaction %>%
    mutate(main_seasonality_q = p.adjust(main_seasonality_p, method = "BH"), interaction_q = p.adjust(interaction_p, method = "BH"))
}
write.csv(seasonality_interaction, "figures_v2/seasonality_interaction.csv", row.names = FALSE)
cat("\n--- Temperature seasonality: main effect + interaction with mean temperature (S10, All) ---\n")
if (nrow(seasonality_interaction) > 0) {
  print(as.data.frame(seasonality_interaction %>% filter(population == "S10", ecosystem == "All") %>%
                         select(metric, main_seasonality_est, main_seasonality_q, interaction_est, interaction_q)))
} else {
  cat("(no seasonality_interaction rows -- too few ponds passed the nrow>=20 guard in this run)\n")
}

safe_step <- function(label, expr) tryCatch(expr, error = function(e) message(sprintf("[SKIPPED] %s: %s", label, conditionMessage(e))))

one_relationship_panel <- function(df, metric, pred, label, batt_row) {
  dd <- df %>% filter(!is.na(.data[[metric]]), !is.na(.data[[pred]]), !is.na(Region), !is.na(Source_DOI))
  if (nrow(dd) < 15) return(NULL)
  m <- tryCatch(lmer(as.formula(paste(metric, "~", pred, "+ (1 | Region/Source_DOI)")), data = dd, REML = TRUE), error = function(e) NULL)
  is_sig <- !is.null(batt_row) && !is.na(batt_row$q.value) && batt_row$q.value < 0.05
  line_col <- if (is_sig) "#e34948" else "grey55"
  p <- ggplot(dd, aes(x = .data[[pred]], y = .data[[metric]]))
  if (!is.null(m)) {
    xr <- range(dd[[pred]], na.rm = TRUE); newd <- tibble(x = seq(xr[1], xr[2], length.out = 60)); names(newd) <- pred
    fx <- fixef(m); X <- cbind(1, newd[[pred]]); V <- as.matrix(vcov(m))
    newd$fit <- as.numeric(X %*% fx); newd$se <- sqrt(diag(X %*% V %*% t(X)))
    newd$lwr <- newd$fit - 1.96 * newd$se; newd$upr <- newd$fit + 1.96 * newd$se
    p <- p + geom_ribbon(data = newd, aes(x = .data[[pred]], ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = line_col, alpha = 0.15) +
      geom_line(data = newd, aes(x = .data[[pred]], y = fit), inherit.aes = FALSE, color = line_col, linewidth = 1)
  }
  qtxt <- if (is.null(batt_row) || is.na(batt_row$q.value)) "n/a" else sprintf("%.3f", batt_row$q.value)
  p + geom_point(aes(color = Ecosystem_Grp), alpha = 0.4, size = 1.2) +
    scale_color_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), guide = "none") +
    labs(x = label, y = NULL, title = sprintf("%s (q=%s, n=%d)", label, qtxt, nrow(dd))) +
    theme(plot.title = element_text(size = 8, color = if (is_sig) "#e34948" else "grey30"), axis.title.x = element_text(size = 8))
}
make_predictor_grid <- function(metric) {
  dd <- filter(final_df, Include_S10)
  panels <- lapply(seq_len(nrow(PREDICTOR_META)), function(i) {
    pred <- PREDICTOR_META$predictor[i]; label <- PREDICTOR_META$label[i]
    br <- battery %>% filter(population == "S10", ecosystem == "All", metric == !!metric, predictor == !!pred)
    one_relationship_panel(dd, metric, pred, label, if (nrow(br)) br[1, ] else NULL)
  })
  plot_grid(plotlist = Filter(Negate(is.null), panels), ncol = 2)
}
for (m in c("Connectance", "Modularity", "z_NODF")) {
  safe_step(paste0("Fig0_", m), {
    g <- make_predictor_grid(m)
    title <- ggdraw() + draw_label(sprintf("All requested factors vs. %s (S>=10 ponds; red = BH-FDR q<0.05)", m), fontface = "bold", size = 13) +
      theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(sprintf("figures_v2/Fig0_all_predictors_vs_%s.png", m),
           plot_grid(title, g, ncol = 1, rel_heights = c(0.06, 1)) + theme(plot.background = element_rect(fill = "white", color = NA)),
           width = 9, height = 14, dpi = 300, bg = "white")
  })
}

safe_step("Fig1_forest_plot", {
  fp_df <- battery %>% filter(population == "S10", ecosystem != "All") %>%
    mutate(metric = factor(metric, levels = METRICS),
           predictor_label = dplyr::recode(predictor, Abs_Latitude = "Absolute latitude", Mean_Annual_Temp = "Mean annual temperature"),
           sig = q.value < 0.05)
  if (nrow(fp_df) == 0) stop("no S>=10 ecosystem-split models available")
  p1 <- ggplot(fp_df, aes(x = estimate, y = metric, color = ecosystem)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), orientation = "y", width = 0.15, position = position_dodge(width = 0.5)) +
    geom_point(aes(shape = sig), size = 2.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c(Lentic = COL_LENTIC, Lotic = COL_LOTIC), name = "Ecosystem") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH-FDR q<0.05") +
    facet_wrap(~predictor_label) +
    labs(x = "Standardized effect (95% CI)", y = NULL, title = "Standardized effects of latitude and temperature on network metrics",
         subtitle = "S >= 10 ponds; (1|Region/Source_DOI); filled = survives FDR correction")
  ggsave("figures_v2/Fig1_forest_plot.png", p1, width = 8, height = 5, dpi = 300)
})

safe_step("Fig2_correlation_matrix", {
  cor_vars <- c("Connectance", "Modularity", "z_Q", "z_NODF", "Mean_Strength", "Mean_Closeness",
                "Abs_Latitude", "Mean_Annual_Temp", "Richness", "N_Prey", "Mean_Omnivory_Index", "BodySize_Ratio", "BodySize_MaxGap", "Lake_Area_km2")
  cor_data <- final_df %>% filter(Include_S10) %>% select(any_of(cor_vars)) %>% drop_na()
  if (nrow(cor_data) <= 20) stop(sprintf("only %d complete rows (need >20)", nrow(cor_data)))
  cor_mat <- cor(cor_data, method = "spearman"); p_mat <- cor.mtest(cor_data)$p
  png("figures_v2/Fig2_correlation_matrix.png", width = 1000, height = 1000, res = 120)
  corrplot(cor_mat, method = "color", type = "upper", order = "hclust", col = colorRampPalette(c(DIV_LOW, DIV_MID, DIV_HIGH))(200),
           p.mat = p_mat, sig.level = 0.05, insig = "blank", addCoef.col = "black", tl.col = "black", tl.srt = 45, diag = FALSE,
           title = "Spearman correlations (S>=10 ponds; blanks = p>0.05, uncorrected)", mar = c(0, 0, 2, 0))
  dev.off()
})

safe_step("Fig3_global_maps", {
  world_map <- map_data("world")
  map_list <- lapply(c("Connectance", "Modularity", "z_Q", "Mean_Annual_Temp"), function(metric) {
    dd <- final_df %>% filter(Include_S10, !is.na(.data[[metric]]), !is.na(Latitude), !is.na(Longitude))
    ggplot() + geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey88", color = "white", linewidth = 0.1) +
      geom_point(data = dd, aes(x = Longitude, y = Latitude, color = .data[[metric]]), alpha = 0.8, size = 1.6) +
      scale_color_viridis_c(option = "viridis") + coord_fixed(1.3) + theme_void(base_size = 10) +
      labs(title = metric, color = NULL) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11), legend.position = "bottom", legend.key.width = unit(1, "cm"))
  })
  ggsave("figures_v2/Fig3_global_maps.png", plot_grid(plotlist = map_list, ncol = 2), width = 12, height = 8, dpi = 300)
})

TERTILE_COLORS <- c(Narrow = "#86b6ef", Medium = "#2a78d6", Wide = "#104281")
make_bodysize_interaction_fig <- function(size_col, size_label, out_name) {
  dd <- final_df %>% filter(Include_S10, !is.na(Connectance), !is.na(Mean_Annual_Temp), n_fish_with_length >= 3, !is.na(.data[[size_col]]), !is.na(Ecosystem_Grp))
  if (nrow(dd) < 15) stop(sprintf("only %d ponds (need >=15)", nrow(dd)))
  brks <- quantile(dd[[size_col]], c(0, 1/3, 2/3, 1), na.rm = TRUE)
  dd$SizeGroup <- cut(dd[[size_col]], breaks = unique(brks), labels = NULL, include.lowest = TRUE)
  levels(dd$SizeGroup) <- c("Narrow", "Medium", "Wide")[seq_along(levels(dd$SizeGroup))]
  p <- ggplot(dd, aes(x = Mean_Annual_Temp, y = Connectance, color = SizeGroup, fill = SizeGroup)) +
    geom_point(alpha = 0.4, size = 1.3) + geom_smooth(method = "lm", linewidth = 1) +
    scale_color_manual(values = TERTILE_COLORS, name = paste(size_label, "tertile")) +
    scale_fill_manual(values = TERTILE_COLORS, name = paste(size_label, "tertile")) +
    facet_wrap(~Ecosystem_Grp) +
    labs(x = "Mean annual temperature (deg C)", y = "Connectance", title = "Does body-size structure moderate the temperature-connectance relationship?",
         subtitle = sprintf("(Petchey et al. 2010 Fig. 6 analog, by %s)", size_label))
  ggsave(out_name, p, width = 9, height = 5, dpi = 300)
}
safe_step("Fig4a", make_bodysize_interaction_fig("SD_log10_Length", "size SD", "figures_v2/Fig4a_bodysize_SD_temperature_interaction.png"))
safe_step("Fig4b", make_bodysize_interaction_fig("BodySize_Ratio", "size ratio", "figures_v2/Fig4b_bodysize_ratio_temperature_interaction.png"))

safe_step("Fig5_modularity_vs_zQ", {
  fp <- battery %>% filter(population == "S10", ecosystem == "All", metric %in% c("Modularity", "z_Q"), predictor %in% c("Abs_Latitude", "Mean_Annual_Temp")) %>%
    mutate(predictor_label = dplyr::recode(predictor, Abs_Latitude = "Absolute latitude", Mean_Annual_Temp = "Mean annual temperature"),
           metric_label = dplyr::recode(metric, Modularity = "Raw Modularity (Q)", z_Q = "Null-corrected (z_Q)"),
           metric_label = factor(metric_label, levels = c("Raw Modularity (Q)", "Null-corrected (z_Q)")), sig = q.value < 0.05)
  if (nrow(fp) == 0) stop("no Modularity/z_Q rows")
  p5 <- ggplot(fp, aes(x = estimate, y = metric_label)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), orientation = "y", width = 0.1) +
    geom_point(aes(color = sig), size = 3) + scale_color_manual(values = c(`TRUE` = "#e34948", `FALSE` = "grey55"), guide = "none") +
    facet_wrap(~predictor_label) +
    labs(x = "Standardized effect (95% CI)", y = NULL, title = "Does S-correction change the modularity-climate relationship?",
         subtitle = "Same model spec, S>=10 ponds; red = BH-FDR q<0.05")
  ggsave("figures_v2/Fig5_modularity_vs_zQ_nullcorrection.png", p5, width = 9, height = 3.8, dpi = 300)
})

# ---- Fig 8: does climate variability (seasonality) moderate the mean-
#      temperature ~ Connectance relationship? Tertile split (ordinal
#      sequential ramp, same convention as the body-size tertile figure),
#      faceted by ecosystem. ----
safe_step("Fig8_seasonality_interaction", {
  dd <- final_df %>% filter(Include_S10, !is.na(Connectance), !is.na(Mean_Annual_Temp), !is.na(Temp_Seasonality), !is.na(Ecosystem_Grp))
  if (nrow(dd) < 15) stop(sprintf("only %d ponds (need >=15)", nrow(dd)))
  brks <- quantile(dd$Temp_Seasonality, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  dd$SeasGroup <- cut(dd$Temp_Seasonality, breaks = unique(brks), labels = NULL, include.lowest = TRUE)
  levels(dd$SeasGroup) <- c("Low seasonality", "Medium seasonality", "High seasonality")[seq_along(levels(dd$SeasGroup))]
  p8 <- ggplot(dd, aes(x = Mean_Annual_Temp, y = Connectance, color = SeasGroup, fill = SeasGroup)) +
    geom_point(alpha = 0.4, size = 1.3) + geom_smooth(method = "lm", linewidth = 1) +
    scale_color_manual(values = setNames(TERTILE_COLORS, c("Low seasonality", "Medium seasonality", "High seasonality")), name = "Temp. seasonality tertile") +
    scale_fill_manual(values = setNames(TERTILE_COLORS, c("Low seasonality", "Medium seasonality", "High seasonality")), name = "Temp. seasonality tertile") +
    facet_wrap(~Ecosystem_Grp) +
    labs(x = "Mean annual temperature (deg C)", y = "Connectance",
         title = "Does climate variability moderate the temperature-connectance relationship?",
         subtitle = "Tertile split on BIO4 (temperature seasonality, SD of monthly means)")
  ggsave("figures_v2/Fig8_seasonality_interaction.png", p8, width = 9, height = 5, dpi = 300)
})

# ==========================================================================
# PHASE 5: Bootstrap mediation + piecewiseSEM (Temperature AND Latitude) + VIF
# ==========================================================================
cat("\n================================================================\nPHASE 5: Mediation bootstrap + SEM + VIF\n================================================================\n")
d <- final_df %>% filter(Include_S10) %>% mutate(Log_Richness = log10(Richness))

fit_path <- function(data, response, predictor, extra = NULL) {
  rhs <- paste(c(predictor, extra, "(1 | Region/Source_DOI)"), collapse = " + ")
  tryCatch(lmer(as.formula(paste(response, "~", rhs)), data = data, REML = TRUE), error = function(e) NULL)
}
get_coef <- function(model, term) {
  if (is.null(model)) return(NA_real_)
  s <- summary(model)$coefficients
  if (!term %in% rownames(s)) return(NA_real_)
  s[term, "Estimate"]
}
MEDIATORS <- c("Log_Richness", "Prey_Niche_Area")
boot_indirect <- function(metric, predictor, mediator, base_df) {
  dd <- base_df %>% filter(!is.na(.data[[metric]]), !is.na(.data[[predictor]]), !is.na(.data[[mediator]]))
  if (nrow(dd) < 30 || n_distinct(dd$Source_DOI) < 10) return(NULL)
  dd$m_z <- zscore(dd[[metric]]); dd$p_z <- zscore(dd[[predictor]]); dd$med_z <- zscore(dd[[mediator]])
  point_a <- get_coef(fit_path(dd, "med_z", "p_z"), "p_z")
  point_b <- get_coef(fit_path(dd, "m_z", "p_z", extra = "med_z"), "med_z")
  study_ids <- unique(dd$Source_DOI); set.seed(SEED)
  boot_vals <- replicate(N_BOOT, {
    resampled_studies <- sample(study_ids, length(study_ids), replace = TRUE)
    bd <- bind_rows(lapply(seq_along(resampled_studies), function(i) {
      rows <- dd %>% filter(Source_DOI == resampled_studies[i]); rows$Source_DOI <- paste0(rows$Source_DOI, "__rep", i); rows
    }))
    a <- get_coef(fit_path(bd, "med_z", "p_z"), "p_z"); b <- get_coef(fit_path(bd, "m_z", "p_z", extra = "med_z"), "med_z")
    a * b
  })
  boot_vals <- boot_vals[!is.na(boot_vals)]
  if (length(boot_vals) < N_BOOT * 0.5) return(NULL)
  ci <- quantile(boot_vals, c(0.025, 0.975))
  tibble(metric = metric, predictor = predictor, mediator = mediator, n = nrow(dd), path_a = point_a, path_b = point_b,
         indirect_effect = point_a * point_b, ci_lwr = ci[1], ci_upr = ci[2], n_boot_ok = length(boot_vals),
         sig = (ci[1] > 0 & ci[2] > 0) | (ci[1] < 0 & ci[2] < 0))
}

boot_file <- "figures_v2/mediation_bootstrap.csv"
if (!file.exists(boot_file) || FORCE_REBOOT) {
  cat(sprintf("Running cluster-bootstrap mediation (%d reps, %d tests) -- slow part...\n", N_BOOT, length(METRICS) * length(THERMAL_PREDICTORS) * length(MEDIATORS)))
  t0 <- Sys.time()
  mediation_boot <- list()
  for (metric in METRICS) for (pred in THERMAL_PREDICTORS) for (med in MEDIATORS) {
    res <- tryCatch(boot_indirect(metric, pred, med, d), error = function(e) NULL)
    if (!is.null(res)) mediation_boot[[length(mediation_boot) + 1]] <- res
    cat(sprintf("  done: %s ~ %s via %s (%.1f min elapsed)\n", metric, pred, med, as.numeric(Sys.time() - t0, units = "mins")))
  }
  mediation_boot <- bind_rows(mediation_boot)
  write.csv(mediation_boot, boot_file, row.names = FALSE)
  cat(sprintf("Bootstrap mediation done in %.1f minutes\n", as.numeric(Sys.time() - t0, units = "mins")))
} else {
  cat("SKIPPED -- mediation_bootstrap.csv exists (ISOFRESH_FORCE_REBOOT=1 to redo)\n")
  mediation_boot <- read_csv(boot_file, show_col_types = FALSE)
}

safe_step("Fig6_mediation_bootstrap", {
  mb <- mediation_boot %>% mutate(predictor_label = dplyr::recode(predictor, Abs_Latitude = "Absolute latitude", Mean_Annual_Temp = "Mean annual temperature"))
  p6 <- ggplot(mb, aes(x = indirect_effect, y = metric, color = mediator)) + geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = ci_lwr, xmax = ci_upr), orientation = "y", width = 0.15, position = position_dodge(width = 0.5)) +
    geom_point(aes(shape = sig), size = 2.5, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c(Log_Richness = "#2a78d6", Prey_Niche_Area = "#eb6834"), name = "Mediator") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "95% CI excludes 0") + facet_wrap(~predictor_label) +
    labs(x = "Indirect effect (standardized a x b, 95% cluster-bootstrap CI)", y = NULL,
         title = "Bootstrap mediation: indirect effects via Richness and Prey Niche Area",
         subtitle = sprintf("S>=10 ponds; bootstrap resampled at the Source_DOI level (n=%d reps)", N_BOOT))
  ggsave("figures_v2/Fig6_mediation_bootstrap.png", p6, width = 9, height = 5, dpi = 300)
})

# ---- piecewiseSEM, parameterized so it can run for EITHER Temperature or
# Latitude as the exogenous driver -- run for both and compare, per the
# 2026-09-01 request to foreground Latitude as the primary, integrative,
# sampling-time-independent predictor (WorldClim's Mean_Annual_Temp is a
# 30-year climatological normal; IsoFresh records no sampling month/season,
# confirmed by inspection, so it cannot be matched to actual conditions at
# collection time the way a fixed spatial covariate like latitude can be).
build_sem <- function(exog_col, exog_label, tag) {
  sem_df <- d %>% filter(!is.na(Connectance), !is.na(z_NODF), !is.na(Mean_Strength), !is.na(.data[[exog_col]]),
                          !is.na(Log_Richness), !is.na(Prey_Niche_Area), !is.na(Region), !is.na(Source_DOI)) %>%
    mutate(Exog_z = zscore(.data[[exog_col]]), Richness_z = zscore(Log_Richness), PreyArea_z = zscore(Prey_Niche_Area),
           Connectance_z = zscore(Connectance), zNODF_z = zscore(z_NODF), Strength_z = zscore(Mean_Strength))
  cat(sprintf("\nSEM (%s) sample size: n=%d\n", exog_label, nrow(sem_df)))

  m_rich <- lmer(Richness_z ~ Exog_z + (1 | Region/Source_DOI), data = sem_df)
  m_prey <- lmer(PreyArea_z ~ Exog_z + Richness_z + (1 | Region/Source_DOI), data = sem_df)
  m_conn <- lmer(Connectance_z ~ Exog_z + Richness_z + PreyArea_z + (1 | Region/Source_DOI), data = sem_df)
  m_nodf <- lmer(zNODF_z ~ Exog_z + Richness_z + PreyArea_z + (1 | Region/Source_DOI), data = sem_df)
  m_str  <- lmer(Strength_z ~ Exog_z + Richness_z + PreyArea_z + (1 | Region/Source_DOI), data = sem_df)
  sem_mod <- psem(m_rich, m_prey, m_conn, m_nodf, m_str, data = as.data.frame(sem_df))
  sem_summary <- summary(sem_mod)
  sink(sprintf("figures_v2/piecewiseSEM_summary_%s.txt", tag)); print(sem_summary); sink()

  coef_tbl <- sem_summary$coefficients
  names(coef_tbl)[names(coef_tbl) == ""] <- "Sig"
  write.csv(coef_tbl, sprintf("figures_v2/piecewiseSEM_coefficients_%s.csv", tag), row.names = FALSE)

  get_path <- function(response, predictor) {
    row <- coef_tbl %>% filter(Response == response, Predictor == predictor)
    if (nrow(row) == 0) return(NA_real_)
    row$Std.Estimate[1]
  }
  trace_effects <- function(outcome) {
    a_exog <- get_path("Richness_z", "Exog_z"); a_prey_direct <- get_path("PreyArea_z", "Exog_z"); a_prey_via_rich <- get_path("PreyArea_z", "Richness_z")
    b_rich <- get_path(outcome, "Richness_z"); b_prey <- get_path(outcome, "PreyArea_z"); direct <- get_path(outcome, "Exog_z")
    ind_via_rich <- a_exog * b_rich; ind_via_prey_direct <- a_prey_direct * b_prey; ind_via_rich_then_prey <- a_exog * a_prey_via_rich * b_prey
    total_indirect <- sum(ind_via_rich, ind_via_prey_direct, ind_via_rich_then_prey, na.rm = TRUE)
    tibble(driver = exog_label, outcome = outcome, direct_effect = direct, indirect_via_Richness = ind_via_rich,
           indirect_via_PreyArea_direct = ind_via_prey_direct, indirect_via_Richness_then_PreyArea = ind_via_rich_then_prey,
           total_indirect = total_indirect, total_effect = direct + total_indirect)
  }
  sem_effects <- bind_rows(lapply(c("Connectance_z", "zNODF_z", "Strength_z"), trace_effects))
  write.csv(sem_effects, sprintf("figures_v2/piecewiseSEM_traced_effects_%s.csv", tag), row.names = FALSE)

  safe_step(paste0("Fig7_", tag), {
    nodes <- tibble(name = c(exog_label, "Richness", "Prey_Niche_Area", "Connectance", "z_NODF", "Mean_Strength"),
                     x = c(0, 2.2, 2.2, 4.6, 4.6, 4.6), y = c(0, 1.9, -1.9, 2.6, 0, -2.6),
                     half_w = c(max(0.55, nchar(exog_label) * 0.09), 0.62, 0.85, 0.62, 0.55, 0.68), half_h = 0.28)
    edges <- coef_tbl %>% transmute(from = Predictor, to = Response, estimate = Std.Estimate, p = P.Value) %>%
      mutate(from = dplyr::recode(from, Exog_z = exog_label, Richness_z = "Richness", PreyArea_z = "Prey_Niche_Area"),
             to = dplyr::recode(to, Connectance_z = "Connectance", zNODF_z = "z_NODF", Strength_z = "Mean_Strength", PreyArea_z = "Prey_Niche_Area", Richness_z = "Richness")) %>%
      left_join(nodes, by = c("from" = "name")) %>% rename(x0 = x, y0 = y, hw0 = half_w, hh0 = half_h) %>%
      left_join(nodes, by = c("to" = "name")) %>% rename(x1 = x, y1 = y, hw1 = half_w, hh1 = half_h) %>%
      mutate(sig = p < 0.05, sign = ifelse(estimate > 0, "Positive", "Negative"), sig_label = ifelse(sig, "p < 0.05", "n.s."))
    point_at_box_edge <- function(cx, cy, ox, oy, hw, hh) {
      dx <- ox - cx; dy <- oy - cy
      t_x <- if (dx != 0) hw / abs(dx) else Inf; t_y <- if (dy != 0) hh / abs(dy) else Inf
      t <- min(t_x, t_y, 1); c(cx + dx * t, cy + dy * t)
    }
    edges <- edges %>% rowwise() %>% mutate(
      start = list(point_at_box_edge(x0, y0, x1, y1, hw0, hh0)), end = list(point_at_box_edge(x1, y1, x0, y0, hw1, hh1)),
      xs = start[1], ys = start[2], xe = end[1], ye = end[2],
      curve = dplyr::case_when(from == exog_label ~ 0.05, from == "Richness" ~ 0.18, from == "Prey_Niche_Area" ~ -0.18, TRUE ~ 0.1)
    ) %>% ungroup()
    p_sem <- ggplot()
    for (cv in sort(unique(edges$curve))) {
      p_sem <- p_sem + geom_curve(data = filter(edges, curve == cv),
                                   aes(x = xs, y = ys, xend = xe, yend = ye, color = sign, linewidth = abs(estimate), linetype = sig_label),
                                   curvature = cv, arrow = arrow(length = unit(0.22, "cm"), type = "closed"), alpha = 0.85)
    }
    edges <- edges %>% mutate(mx = (xs + xe) / 2, my = (ys + ye) / 2, dx = xe - xs, dy = ye - ys, lx = mx - dy * curve * 0.5, ly = my + dx * curve * 0.5)
    p_sem <- p_sem +
      geom_label(data = edges, aes(x = lx, y = ly, label = sprintf("%.2f%s", estimate, ifelse(sig, "*", ""))), size = 3.1, label.size = 0, fill = "white", alpha = 0.9) +
      geom_label(data = nodes, aes(x = x, y = y, label = gsub("_", " ", name)), size = 4.2, fontface = "bold", fill = "grey95", label.size = 0.3, label.padding = unit(0.4, "lines")) +
      scale_color_manual(values = c(Positive = "#2a78d6", Negative = "#e34948"), name = "Path sign") +
      scale_linewidth(range = c(0.3, 2.2), guide = "none") +
      scale_linetype_manual(values = c("p < 0.05" = "solid", "n.s." = "22"), name = "Significance") +
      scale_x_continuous(expand = expansion(mult = 0.12)) + scale_y_continuous(expand = expansion(mult = 0.12)) +
      coord_cartesian(clip = "off") + theme_void(base_size = 12) +
      labs(title = sprintf("piecewiseSEM path model: %s -> Richness/Prey niche -> network structure", exog_label),
           subtitle = sprintf("S>=10 ponds, n=%d; standardized path coefficients (* = p<0.05); (1|Region/Source_DOI) in every component model", nrow(sem_df))) +
      theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10), legend.position = "bottom", plot.margin = margin(15, 25, 15, 15))
    ggsave(sprintf("figures_v2/Fig7_piecewiseSEM_path_diagram_%s.png", tag), p_sem, width = 12, height = 8, dpi = 300, bg = "white")
  })
  list(coefficients = coef_tbl, effects = sem_effects, n = nrow(sem_df))
}

# Each build_sem() call is wrapped so one driver's model-fitting failure
# (e.g. too few observations relative to the number of Region/Source_DOI
# groups, which only arises on tiny test/subset runs -- the real ~510-pond
# dataset has far more rows than groups) doesn't halt the whole script and
# lose the other drivers' results, or anything already computed upstream
# (the bootstrap mediation CSV is already saved by this point regardless).
safe_build_sem <- function(exog_col, exog_label, tag) {
  tryCatch(build_sem(exog_col, exog_label, tag), error = function(e) {
    message(sprintf("[SKIPPED SEM] %s: %s", exog_label, conditionMessage(e))); NULL })
}
sem_temp <- safe_build_sem("Mean_Annual_Temp", "Temperature", "temperature")
sem_lat  <- safe_build_sem("Abs_Latitude", "Latitude", "latitude")
sem_seas <- safe_build_sem("Temp_Seasonality", "Temp_Seasonality", "seasonality")

sem_compare <- bind_rows(sem_temp$effects, sem_lat$effects, sem_seas$effects)
write.csv(sem_compare, "figures_v2/piecewiseSEM_driver_comparison.csv", row.names = FALSE)
cat("\n--- SEM comparison: Temperature- vs Latitude- vs Seasonality-driven (standardized total effects) ---\n")
if (nrow(sem_compare) > 0) {
  print(as.data.frame(sem_compare %>% select(driver, outcome, direct_effect, total_indirect, total_effect)))
} else {
  cat("(no SEM driver succeeded -- see [SKIPPED SEM] messages above)\n")
}

# ---- VIF check on the z_Q suppression-effect model, both drivers ----
safe_step("VIF_temperature", {
  vif_df <- d %>% filter(!is.na(z_Q), !is.na(Mean_Annual_Temp), !is.na(Log_Richness), !is.na(Prey_Niche_Area))
  vif_df$Temp_z <- zscore(vif_df$Mean_Annual_Temp); vif_df$S_z <- zscore(vif_df$Log_Richness); vif_df$PA_z <- zscore(vif_df$Prey_Niche_Area)
  vif_vals <- car::vif(lmer(zscore(z_Q) ~ Temp_z + S_z + PA_z + (1 | Region/Source_DOI), data = vif_df))
  cat("\n--- VIF check (z_Q ~ Temp + Richness + PreyArea) ---\n"); print(vif_vals)
  write.csv(data.frame(term = names(vif_vals), VIF = vif_vals), "figures_v2/vif_check_zQ_model_temperature.csv", row.names = FALSE)
})

safe_step("VIF_latitude", {
  vif_df2 <- d %>% filter(!is.na(z_Q), !is.na(Abs_Latitude), !is.na(Log_Richness), !is.na(Prey_Niche_Area))
  vif_df2$Lat_z <- zscore(vif_df2$Abs_Latitude); vif_df2$S_z <- zscore(vif_df2$Log_Richness); vif_df2$PA_z <- zscore(vif_df2$Prey_Niche_Area)
  vif_vals2 <- car::vif(lmer(zscore(z_Q) ~ Lat_z + S_z + PA_z + (1 | Region/Source_DOI), data = vif_df2))
  cat("\n--- VIF check (z_Q ~ Latitude + Richness + PreyArea) ---\n"); print(vif_vals2)
  write.csv(data.frame(term = names(vif_vals2), VIF = vif_vals2), "figures_v2/vif_check_zQ_model_latitude.csv", row.names = FALSE)
})

cat("\n================================================================\nALL PHASES COMPLETE -- see figures_v2/\n================================================================\n")
