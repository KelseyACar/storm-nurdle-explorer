# SCRIPT 8: Site Characteristics — Manufacturer Proximity & MP Redistribution
  # Purpose: Do site characteristics predict post-storm MP magnitude/direction?

# Load Libraries
# Read Data
# Config
# Helper FNs
  # Section 1: Build site-level predictors
  # Section 2: Build storm sequence & time-since-previous storm analysis df
  # Section 3: Build main analysis df
  # Section 4: Site predictory models
    # A) Magnitude: CPLM
    # B) Direction: CLMM (ordinal)
  # Section 5: Repeat storm exposure plot
  # Section 6: Sampling resolution diagnositc
  # Section 7: scatterplots of pcs_chang vs site predictors


# NOTE: Sites 6 & 7 are river-connected, not coastal
# NOTE: sediment type / exposure level not available in current data; need an external shoreline classification source (e.g. NOAA ESI) added later.


# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr', 'tidyr', 'ggplot2', 'here', 'cplm', 'ordinal')
source(here("functions.R"))


# Read Data ---------------------------------------------------------------
sites <- read.csv(here("output","summary_tables", "distinct_sites.csv")) 

manufacturers <- read.csv(here("output", "manufacturer_df.csv"))

mp_change <- read.csv(here("output", "summary_tables", "MP_delta_bySite.csv"))

phase_abund <- read.csv(here("output", "summary_tables", "site_abund_byStormPhase.csv"))

master_df <- read.csv(here("output", "master_df.csv"))

regional_abund_bySiteStorm <- read.csv(here("output", "summary_tables", "regional_abund_bySiteStorm.csv"))

# Config ------------------------------------------------------------------
strand <- "total"
stat <- "median"
alpha <- 0.05

pcs_delta_col   <- paste0("pcs_change_", stat, "_", strand)
change_type_col <- paste0("change_type_", stat, "_", strand)
mp_change[["pcs_change"]]  <- mp_change[[pcs_delta_col]]
mp_change[["change_type"]] <- mp_change[[change_type_col]]



# Helper Functions --------------------------------------------------------


# Helper fn to get Haversine distance (km) between two lat/long points
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  lat1 <- lat1 * pi / 180; lat2 <- lat2 * pi / 180
  dlat <- lat2 - lat1; dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}


# nearest manufacturer distance per site
site_mfg_dist <- sites %>%
  rowwise() %>%
  mutate(
    dist_to_manufacturer_km = min(haversine_km(latitude, longitude, manufacturers$lat, manufacturers$long))
  ) %>%
  ungroup() %>%
  select(site_id, site_name, dist_to_manufacturer_km)

write.csv(site_mfg_dist, here("output", "site_level", "stats", "site_manufacturer_distance.csv"), row.names = FALSE)



# Section 1: Build site-level predictors -------------------------------------

# 1a) manufacturer distance
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  lat1 <- lat1 * pi / 180; lat2 <- lat2 * pi / 180
  dlat <- lat2 - lat1; dlon <- (lon2 - lon1) * pi / 180
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}

site_mfg_dist <- sites %>%
  rowwise() %>%
  mutate(dist_to_manufacturer_km = min(haversine_km(latitude, longitude, manufacturers$lat, manufacturers$long))) %>%
  ungroup() %>%
  select(site_id, dist_to_manufacturer_km)

# 1b) baseline / control-period abundance per site
baseline_abund <- phase_abund %>%
  filter(fine_storm_phase == "control") %>%
  select(site_id, baseline_abundance = paste0(stat, "_abundance_", strand), baseline_n_surveys = n_surveys)

# 1c) repeat storm exposure -- how many storms were sampled at this site
sampling_effort <- regional_abund_bySiteStorm %>%
  group_by(site_id) %>%
  summarise(total_surveys = sum(n_surveys, na.rm = TRUE), .groups = "drop")

site_predictors <- sites %>%
  select(site_id) %>%
  left_join(site_mfg_dist, by = "site_id") %>%
  left_join(baseline_abund, by = "site_id") %>%
  left_join(storm_exposure, by = "site_id") %>%
  left_join(sampling_effort, by = "site_id") %>%
  mutate(
    dist_to_manufacturer_km_z = as.numeric(scale(dist_to_manufacturer_km)),
    baseline_abundance_z = as.numeric(scale(baseline_abundance)),
    n_storms_site_z = as.numeric(scale(n_storms_site))
  )

write.csv(site_predictors, here("output", "site_level", "stats", "site_characteristics_predictors.csv"), row.names = FALSE)



# Section 2: Build storm sequence & time-since-previous storm -------------------------------------------------------

storm_dates <- master_df %>%
  filter(survey_type == "storm_influenced") %>%
  distinct(site_id, name, date_storm_nearest_site) %>%
  mutate(date_storm_nearest_site = as.Date(date_storm_nearest_site)) %>%
  filter(site_id != 7)

storm_sequence <- mp_change %>%
  select(site_id, name, pcs_change, change_type) %>%
  left_join(storm_dates, by = c("site_id", "name")) %>%
  left_join(storm_exposure, by = "site_id") %>%
  filter(repeat_exposure) %>%
  group_by(site_id) %>%
  arrange(date_storm_nearest_site, .by_group = TRUE) %>%
  mutate(
    storm_sequence = row_number(),
    # NA for each site's first storm -- there's no "previous" storm to measure a gap from
    days_since_prev_storm = as.numeric(date_storm_nearest_site - lag(date_storm_nearest_site))
  ) %>%
  ungroup()

# gap table to merge into site_base -> only storms with a real previous-storm gap carry a value;
# a site's first storm (and any site with only 1 storm total) will be NA here and get excluded
# automatically by the predictor loop's own NA filter in Section 4, not by dropping rows now
storm_gap <- storm_sequence %>%
  select(site_id, name, days_since_prev_storm)



# Section 3: Build main analysis df --------------------------------------
site_base <- mp_change %>%
  select(site_id, name, pcs_change, change_type) %>%
  left_join(site_predictors, by = "site_id") %>%
  left_join(storm_gap, by = c("site_id", "name")) %>%
  mutate(
    abs_pcs_change = abs(pcs_change),
    change_type_ord = factor(change_type, levels = c("depletion", "no_change", "accumulation"), ordered = TRUE),
    days_since_prev_storm_z = as.numeric(scale(days_since_prev_storm))
  )


# Section 4: Predictor models --------------------------------------------
# one predictor per model, same framework as script 7's storm-metric models

site_predictors_z <- c("dist_to_manufacturer_km_z", "baseline_abundance_z", "n_storms_site_z", "days_since_prev_storm_z")
site_predictors_z <- c("dist_to_manufacturer_km_z", "baseline_abundance_z", "n_storms_site_z", "days_since_prev_storm_z")
magnitude_results <- data.frame()
direction_results <- data.frame()

for (z_col in site_predictors_z) {
  p <- sub("_z$", "", z_col)
  df <- site_base %>% filter(!is.na(.data[[z_col]]), !is.na(abs_pcs_change))
  
  # magnitude -- Tweedie GLMM (cplm) in place of LMER: abs_pcs_change is non-negative
  # with exact zeros (no_change interactions), violating LMER's normality assumption
  mag_model <- tryCatch(
    cpglmm(reformulate(c(z_col, "(1 | site_id)", "(1 | name)"), response = "abs_pcs_change"), data = df),
    error = function(e) NULL
  )
  if (!is.null(mag_model)) {
    cf <- summary(mag_model)@coefs
    if (z_col %in% rownames(cf)) {
      magnitude_results <- rbind(magnitude_results, data.frame(
        predictor = p,
        n_obs = nrow(df),
        estimate_z = cf[z_col, "Estimate"],
        std_error_z = cf[z_col, "Std. Error"],
        statistic = cf[z_col, "t value"],
        p_value = 2 * pnorm(abs(cf[z_col, "t value"]), lower.tail = FALSE)  # normal approximation; cplm does not report df-based p-values
      ))
    }
  }
  
  # direction -- unchanged (CLMM remains the correct model for the ordinal outcome)
  dir_model <- tryCatch(
    ordinal::clmm(reformulate(c(z_col, "(1 | site_id)", "(1 | name)"), response = "change_type_ord"), data = df),
    error = function(e) NULL
  )
  if (!is.null(dir_model)) {
    cf <- coef(summary(dir_model))
    if (z_col %in% rownames(cf)) {
      direction_results <- rbind(direction_results, data.frame(
        predictor = p,
        n_obs = nrow(df),
        estimate_z = cf[z_col, "Estimate"],
        std_error_z = cf[z_col, "Std. Error"],
        statistic = cf[z_col, "z value"],
        p_value = cf[z_col, "Pr(>|z|)"]
      ))
    }
  }
}

magnitude_results$p_adj <- p.adjust(magnitude_results$p_value, method = "BH")
magnitude_results$sig_adj <- magnitude_results$p_adj < alpha
direction_results$p_adj <- p.adjust(direction_results$p_value, method = "BH")
direction_results$sig_adj <- direction_results$p_adj < alpha

write.csv(magnitude_results, here("output", "site_level", "stats", "site_characteristics", "site_characteristics_magnitude_models.csv"), row.names = FALSE)
write.csv(direction_results, here("output", "site_level", "stats", "site_characteristics", "site_characteristics_direction_models.csv"), row.names = FALSE)

# Section 5: Repeat storm exposure visual for model above ------
# does response change across a site's 1st, 2nd, 3rd sampled storm?
# both plots below are descriptive/exploratory; the actual test of whether days_since_prev_storm predicts magnitude/direction lives in Section 4's loop above


# 5a) per-site trajectory across its own storms, x-axis = actual date (not ordinal rank)
# ordinal rank (1st/2nd/3rd storm) can't be read meaningfully since it ignores how much time
# actually passed between storms at a site; real dates fix that. free x/y scales since each
# site has its own date range and its own magnitude range

# first, drop control columns & sites with only 2 storms
storm_sequence <- storm_sequence %>% 
  filter(!is.na(name)) %>%
  filter(n_storms_site > 3)

p_sequence <- ggplot(storm_sequence, aes(x = date_storm_nearest_site, y = pcs_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  facet_wrap(~ site_id, scales = "free") +
  labs(title = "MP Change Across Sequential Storms — Repeat-Exposure Sites",
       x = "Date of Storm (nearest approach to site)",
       y = "Pcs Change (post - pre)") +
  theme_apa(x.font.size = 9, y.font.size = 9, facet.title.size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "none", plot.title = element_text(size = 16))
p_sequence

ggsave(here("output", "site_level", "graphs", "site_characteristics", "repeat_exposure_sequence.pdf"),
       plot = p_sequence, width = 12, height = 10, device = cairo_pdf)


# 5b) does the size of the gap since the previous storm relate to MP change? 
# one shared plot across all sites/storm-pairs (not faceted), since "days since previous storm" is a
# comparable axis across sites in a way that ordinal rank never was
gap_plot_data <- storm_gap %>%
  filter(!is.na(days_since_prev_storm)) %>%
  left_join(mp_change %>% select(site_id, name, pcs_change, change_type), by = c("site_id", "name"))

p_gap <- ggplot(gap_plot_data, aes(x = days_since_prev_storm, y = pcs_change)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = change_type), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = colors) +
  labs(title = "MP Change vs. Time Since Previous Storm",
       x = "Days Since Previous Storm at This Site",
       y = "Pcs Change (post - pre)", color = "Response") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "none", plot.title = element_text(size = 16))
p_gap


ggsave(here("output", "site_level", "graphs", "site_characteristics", "days_since_prev_storm.pdf"),
       plot = p_gap, width = 12, height = 10, device = cairo_pdf)


# Section 6: Sampling resolution diagnostic ------------------------------
# not a formal model, a check to support the text discussion of sampling-resolution limits
# (17 of 32 sites have only 1 storm sampled)

resolution_check <- site_base %>%
  filter(!is.na(total_surveys), !is.na(abs_pcs_change)) %>%
  select(site_id, name, total_surveys, abs_pcs_change, pcs_change, change_type)


# plot
p_resolution <- ggplot(resolution_check, aes(x = total_surveys, y = abs_pcs_change)) +
  geom_point(aes(color = change_type), size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(values = colors) +
  labs(title = "MP Change Magnitude vs. Site Sampling Effort",
       x = "Total Surveys at Site", y = "|Pcs Change| (post - pre)") +
  theme_apa(x.font.size = 13, y.font.size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "none", plot.title = element_text(size = 16))
p_resolution

ggsave(here("output", "site_level", "graphs", "site_characteristics", "sampling_resolution_check.pdf"),
       plot = p_resolution, width = 12, height = 10, device = cairo_pdf)


resolution_cor <- suppressWarnings(
  cor.test(resolution_check$total_surveys, resolution_check$abs_pcs_change, method = "spearman")
)

write.csv(data.frame(rho = resolution_cor$estimate, p_value = resolution_cor$p.value),
          here("output", "site_level", "stats", "sampling_resolution_correlation.csv"), row.names = FALSE)


# Section 7: Scatter plots — pcs_change vs site-level predictors ----------------

plot_predictor_scatter <- function(df, x_var, x_label, filename) {
  p <- ggplot(df, aes(x = .data[[x_var]], y = pcs_change, color = change_type)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
    scale_color_manual(values = colors) +
    labs(title = paste0("MP Delta vs. ", x_label),
         x = x_label, y = "Pcs Change (post - pre)", color = "Response") +
    theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13)
  
  ggsave(here("output", "site_level", "graphs", "site_characteristics", filename),
         plot = p, width = 12, height = 10, device = cairo_pdf)
  p
}

p_mfg <- plot_predictor_scatter(site_base, "dist_to_manufacturer_km", "Distance to Nearest Manufacturer (km)", "scatter_manufacturer_dist.pdf")
p_base <- plot_predictor_scatter(site_base, "baseline_abundance", "Baseline (Control) Abundance", "scatter_baseline_abundance.pdf")
p_nstorm <- plot_predictor_scatter(site_base, "n_storms_site", "Number of Storms Sampled at Site", "scatter_n_storms_site.pdf")
