# SCRIPT 6: Strandline-Specific Storm Response Analysis
# Purpose: Test whether MPs at a single strandline change pre v post storm and at different response windows
# Purpose: Test whether MPs at a single strandline change pre v post storm and at different response windows
# Output shape: long format -- one row per (group x metric), metric in {mean, median, max}.
# LMER is run for mean and median (not max, matching original scope), as its own row (test = "lmer").

# Load Libraries
# Read Data
# Config
# Pre-Process
# Section 1: Storm Phase - GoM level (mean, median, max; primary tests + LMER)
# Section 2: Storm Phase - Storm-level (mean, median, max; primary tests + LMER)
# Section 3: Storm Phase - Site-level (mean, median, max; primary tests + LMER)
# Section 4: Temporal Window Response (mean, median, max; primary tests, no LMER)
# Section 5: Uni-directional Post-storm (storm & site level) -- unaffected by metric looping
# Section 6: Robustness Checks (window sensitivity, outlier trim) -- run per metric

# Load Libraries ----------------------------------------------------------
pacman::p_load('dplyr', 'tidyr', 'here', 'lme4', 'lmerTest', 'emmeans', 'pbkrtest', 'jtools')
source(here("functions.R"))


# Read Data ---------------------------------------------------------------

# pre/post means & change
mp_change <- read.csv(here("output", "summary_tables", "MP_delta_bySite.csv"))

# abundance by phase (pre & post storm) — produced by Script 2, control surveys dropped
abund_storm_phase <- read.csv(here('output','summary_tables','siteStorm_abund_byStormPhase.csv')) %>%
  filter(fine_storm_phase != "control")

# Site/storm/fine_response means (pre & post responses) — produced by Script 2
abund_storm_response <- read.csv(here('output','summary_tables', 'site_temp_window_abund.csv')) %>%
  filter(fine_storm_phase != "control")

# master_df used only in Section 5 robustness checks (needs fine_days, complete_storm_phase, linked_id)
master_df <- read.csv(here("output", "master_df.csv")) %>%
  mutate(date_nurdle = as.Date(date_nurdle),
         date_storm = as.Date(date_storm))


# Config ------------------------------------------------------------------
strand <- "old" # options: old, new, total
metrics <- c("mean", "median") # max is always run in addition to this, required different pre-processing, see helper fns
alpha = 0.05

# draw out relevant columns (for each strand, in case analysis of others is desired)
strand_cols <- list(
  total = list(median = "median_abundance_total", mean = "mean_abundance_total", max = "max_abundance_total", raw = "weighted_total_amount"),
  old = list(median = "median_abundance_old", mean = "mean_abundance_old", max = "max_abundance_old", raw = "weighted_old_strandline"),
  new = list(median = "median_abundance_new", mean = "mean_abundance_new", max = "max_abundance_new", raw = "weighted_new_strandline")
)


# configure column names
col <- strand_cols[[strand]]


# Helper: run paired Wilcoxon/t-test across mean, median, AND max in one call ---
# long_data: long-format df (one row per site/storm/phase) with a fine_storm_phase column
run_all_metrics_paired <- function(long_data) {
  all_metrics <- c(metrics, "max")
  bind_rows(lapply(all_metrics, function(m) {
    wide <- pivot_pre_post(long_data, col[[m]])
    stats <- run_paired_tests(wide$v_post, wide$v_pre)
    stats$metric <- m
    stats
  }))
}

# Helper: run LMER across mean and median only (max LMER was never part of the original storm/site-level scope)
run_all_metrics_lmer <- function(long_data, random_formula) {
  bind_rows(lapply(metrics, function(m) {
    d <- long_data %>%
      filter(fine_storm_phase %in% c("pre", "post")) %>%
      mutate(mp_abundance = .data[[col[[m]]]],
             fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")))
    res <- run_lmer(d, random_formula)
    data.frame(metric = m, lmer_estimate = res$estimate, lmer_p = res$p_value, lmer_sig = res$sig)
  }))
}



# Pre-Process -------------------------------------------------------------
# mp_change: build abundance_pre/post + change_type for EACH metric (mean, median), long-stacked
paired_data <- bind_rows(lapply(metrics, function(m) {
  mp_change %>%
    filter(!is.na(name)) %>%
    mutate(
      metric = m,
      abundance_pre = .data[[paste0(col[[m]], "_pre")]],
      abundance_post = .data[[paste0(col[[m]], "_post")]],
      pcs_change = .data[[paste0("pcs_change_", m, "_", strand)]],
      change_type = .data[[paste0("change_type_", m, "_", strand)]]
    )
}))


# Section 1) Storm Phase: GoM level ---------------------------------------------

# A) Normality check on raw deltas per metric & max
GoM_normality <- bind_rows(lapply(metrics, function(m) {
  d <- paired_data %>% filter(metric == m)
  sw <- shapiro.test(d$pcs_change)
  data.frame(metric = m, shapiro_W = sw$statistic, shapiro_p = sw$p.value)
}))
  
max_wide_gom <- pivot_pre_post(abund_storm_phase, col$max)

sw_max <- shapiro.test(max_wide_gom$v_post - max_wide_gom$v_pre)

GoM_normality <- bind_rows(GoM_normality, data.frame(metric = "max", shapiro_W = sw_max$statistic, shapiro_p = sw_max$p.value)) %>%
  mutate(normal = shapiro_p > alpha, test = "shapiro", level = "GoM", strand = strand)

# paired tests, all metrics + max, in one call
GoM_paired <- run_all_metrics_paired(abund_storm_phase) %>%
  mutate(test = "paired_wilcox_ttest", level = "GoM", strand = strand)

# n_accumulation/depletion/no_change per metric
GoM_change_counts <- paired_data %>%
  group_by(metric) %>%
  summarise(n_accumulation = sum(change_type == "accumulation"),
            n_depletion = sum(change_type == "depletion"),
            n_no_change = sum(change_type == "no_change"), .groups = "drop")

GoM_paired <- GoM_paired %>% left_join(GoM_change_counts, by = "metric")

# LMER, mean + median (site & storm random effects)
GoM_lmer <- run_all_metrics_lmer(abund_storm_phase, "(1 | site_id) + (1 | name)") %>%
  mutate(test = "lmer", level = "GoM", strand = strand)

GoM_results <- bind_rows(GoM_normality, GoM_paired, GoM_lmer)

write.csv(GoM_results, here("output", "GoM_level", "stats", "pre_post", strand, paste0("GoM_pre_post_results_", strand, ".csv")), row.names = FALSE)



# Section 2) Storm Phase - Storm-level -------------------------------
# Pre vs. post across all sites linked to each storm
storms <- unique(paired_data$name)

per_storm_results <- bind_rows(lapply(storms, function(storm) {
  
  d_phase <- abund_storm_phase[abund_storm_phase$name == storm, ]
  if (nrow(d_phase) == 0) return(NULL)
  
  paired_res <- run_all_metrics_paired(d_phase) %>% mutate(name = storm)
  lmer_res <- run_all_metrics_lmer(d_phase, "(1 | site_id)") %>% mutate(name = storm)
  
  d_change <- paired_data %>% filter(name == storm)
  change_counts <- d_change %>%
    group_by(metric) %>%
    summarise(n_accumulation = sum(change_type == "accumulation", na.rm = TRUE),
              n_depletion = sum(change_type == "depletion", na.rm = TRUE),
              n_no_change = sum(change_type == "no_change", na.rm = TRUE), .groups = "drop")
  
  bind_rows(
    paired_res %>% mutate(test = "paired_wilcox_ttest") %>% left_join(change_counts, by = "metric"),
    lmer_res %>% mutate(test = "lmer")
  )
}))

# FDR adj -- grouped by metric, since each metric is its own family of comparisons across storms
per_storm_results <- per_storm_results %>%
  group_by(metric, test) %>%
  mutate(
    wilcox_p_adj = if ("wilcox_p" %in% names(.)) p.adjust(wilcox_p, method = "BH") else NA_real_,
    wilcox_sig_adj = wilcox_p_adj < alpha,
    ttest_p_adj = if ("ttest_p" %in% names(.)) p.adjust(ttest_p, method = "BH") else NA_real_,
    lmer_p_adj = if ("lmer_p" %in% names(.)) p.adjust(lmer_p, method = "BH") else NA_real_,
    lmer_sig_adj = lmer_p_adj < alpha
  ) %>%
  ungroup()

write.csv(per_storm_results, here("output", "storm_level", "stats", "pre_post", strand, paste0("per_storm_paired_stats_", strand, ".csv")), row.names = FALSE)



# Section 3: Storm Phase - Site-level --------------------------------
sites <- unique(paired_data$site_id)

per_site_results <- bind_rows(lapply(sites, function(sid) {
  cat("Processing site:", sid, "\n")
  
  d_phase <- abund_storm_phase[abund_storm_phase$site_id == sid, ]
  if (nrow(d_phase) == 0) return(NULL)
  
  paired_res <- run_all_metrics_paired(d_phase) %>% mutate(site_id = gsub("_cluster", "", sid))
  lmer_res <- run_all_metrics_lmer(d_phase, "(1 | name)") %>% mutate(site_id = gsub("_cluster", "", sid))
  
  d_change <- paired_data %>% filter(site_id == sid)
  change_counts <- d_change %>%
    group_by(metric) %>%
    summarise(n_accumulation = sum(change_type == "accumulation", na.rm = TRUE),
              n_depletion = sum(change_type == "depletion", na.rm = TRUE),
              n_no_change = sum(change_type == "no_change", na.rm = TRUE),
              n_storms = n(), .groups = "drop")
  
  bind_rows(
    paired_res %>% mutate(test = "paired_wilcox_ttest") %>% left_join(change_counts, by = "metric"),
    lmer_res %>% mutate(test = "lmer")
  )
}))

per_site_results <- per_site_results %>%
  group_by(metric, test) %>%
  mutate(
    wilcox_p_adj = if ("wilcox_p" %in% names(.)) p.adjust(wilcox_p, method = "BH") else NA_real_,
    wilcox_sig_adj = wilcox_p_adj < alpha,
    ttest_p_adj = if ("ttest_p" %in% names(.)) p.adjust(ttest_p, method = "BH") else NA_real_,
    lmer_p_adj = if ("lmer_p" %in% names(.)) p.adjust(lmer_p, method = "BH") else NA_real_,
    lmer_sig_adj = lmer_p_adj < alpha
  ) %>%
  ungroup()

write.csv(per_site_results, here("output", "site_level", "stats", "pre_post", strand, paste0("per_site_paired_stats_", strand, ".csv")), row.names = FALSE)


# Section 4: Temporal Window Response ------------------------------------------------
# only site x storm pairs where both pre & post surveys fall in the same fine_response window are included. n in run_paired_tests reflects this.

response_levels <- c("acute", "subacute", "extended")

# A) GoM level per window
GoM_resp <- bind_rows(lapply(response_levels, function(level) {
  level_data <- abund_storm_response[abund_storm_response$fine_response == level, ]
  run_all_metrics_paired(level_data) %>% mutate(response_level = level)
}))

write.csv(GoM_resp, here("output", "GoM_level", "stats", "resp_window", strand, paste0("GoM_resp_window_stats_", strand, ".csv")), row.names = FALSE)


# B) Storm-level per window
storm_resp <- bind_rows(lapply(response_levels, function(level) {
  level_data <- abund_storm_response[abund_storm_response$fine_response == level, ]
  bind_rows(lapply(storms, function(storm) {
    d <- level_data[level_data$name == storm, ]
    if (nrow(d) == 0) return(NULL)
    run_all_metrics_paired(d) %>% mutate(response_level = level, name = storm)
  }))
}))

storm_resp <- storm_resp %>%
  group_by(response_level, metric) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p, method = "BH"),
    wilcox_sig_adj = wilcox_p_adj < alpha
  ) %>%
  ungroup()

write.csv(storm_resp, here("output", "storm_level", "stats", "resp_window", strand, paste0("per_storm_resp_window_stats_", strand, ".csv")), row.names = FALSE)


# C) Site-level per window
site_resp <- bind_rows(lapply(response_levels, function(level) {
  level_data <- abund_storm_response[abund_storm_response$fine_response == level, ]
  bind_rows(lapply(sites, function(sid) {
    d <- level_data[level_data$site_id == sid, ]
    if (nrow(d) == 0) return(NULL)
    run_all_metrics_paired(d) %>% mutate(response_level = level, site_id = gsub("_cluster", "", sid))
  }))
}))

site_resp <- site_resp %>%
  group_by(response_level, metric) %>%
  mutate(
    wilcox_p_adj = p.adjust(wilcox_p, method = "BH"),
    wilcox_sig_adj = wilcox_p_adj < alpha
  ) %>%
  ungroup()

write.csv(site_resp, here("output", "site_level", "stats", "resp_window", strand, paste0("per_site_resp_window_stats_", strand, ".csv")), row.names = FALSE)


# Section 5: Uni-directional Post-storm (storm & site level) ---------------------------------------------
# unaffected by metric looping -- change_type/direction is categorical, uses paired_data filtered to one metric
# NOTE: pick which metric's change_type classification drives this test

primary_metric_for_direction <- "median"

unidata_storm <- paired_data %>%
  filter(metric == primary_metric_for_direction) %>%
  distinct(site_id, name, .keep_all = TRUE)

storm_unidirectional <- run_unidirectional_tests(unidata_storm, group_col = "name", group_label = "storm")

write.csv(storm_unidirectional, here("output", "storm_level", "stats", "rate_and_dir", strand, paste0("unidirectional_storm_resp_", strand, ".csv")), row.names = FALSE)

unidata_site <- paired_data %>%
  filter(metric == primary_metric_for_direction) %>%
  distinct(site_id, name, .keep_all = TRUE)

site_unidirectional <- run_unidirectional_tests(unidata_site, group_col = "site_id", group_label = "site") %>%
  filter(!is.na(group))

write.csv(site_unidirectional, here("output", "site_level", "stats", "rate_and_dir", strand, paste0("site_unidirectional_", strand, ".csv")), row.names = FALSE)


# Section 6: Robustness Checks --------------------------------------------
# A) Window sensitivity (fine_days cutoffs)
day_cutoffs <- c(4, 10, 28)

window_sensitivity <- bind_rows(lapply(metrics, function(m) {
  bind_rows(lapply(day_cutoffs, function(cut) {
    w_data <- master_df %>%
      filter(survey_type == "storm_influenced",
             fine_storm_phase %in% c("pre", "post"),
             abs(fine_days) <= cut) %>%
      distinct(linked_id, date_nurdle, fine_storm_phase, .keep_all = TRUE) %>%
      mutate(mp_abundance = .data[[col$raw]]) %>%
      group_by(linked_id, fine_storm_phase) %>%
      summarise(stat_val = if (m == "median") median(mp_abundance, na.rm = TRUE) else mean(mp_abundance, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = fine_storm_phase, values_from = stat_val) %>%
      filter(!is.na(pre), !is.na(post))

    stats <- run_paired_tests(w_data$post, w_data$pre)
    cbind(data.frame(metric = m, cutoff_days = cut), stats)
  }))
}))

write.csv(window_sensitivity, here("output", "GoM_level", "stats", "pre_post", strand, paste0("GoM_window_sensitivity_", strand, ".csv")), row.names = FALSE)


# Outlier trim, per metric
outlier_trim <- bind_rows(lapply(metrics, function(m) {
  d <- paired_data %>% filter(metric == m)
  pre_q95 <- quantile(d$abundance_pre, 0.95, na.rm = TRUE)
  post_q95 <- quantile(d$abundance_post, 0.95, na.rm = TRUE)
  trimmed <- d %>% filter(abundance_pre <= pre_q95, abundance_post <= post_q95)

  trim_stats <- run_paired_tests(trimmed$abundance_post, trimmed$abundance_pre)
  original_stats <- GoM_paired %>% filter(metric == m)

  data.frame(
    metric = m,
    n_original = nrow(d),
    n_remaining = nrow(trimmed),
    original_p = original_stats$wilcox_p,
    trimmed_p = trim_stats$wilcox_p,
    result_stable = (original_stats$wilcox_p < alpha) == (trim_stats$wilcox_p < alpha)
  )
}))

write.csv(outlier_trim, here("output", "GoM_level", "stats", "pre_post", strand, paste0("GoM_outlier_trim_", strand, ".csv")), row.names = FALSE)


# Test sensitivity to extreme values (robustness of primary tests) by repeating primary paired test after outlier trip (top 5% of pre and post-storm abundance values)
  # results will either be stable or not stable compared to pre-trim
outlier_trim_summary <- outlier_trim %>%
  mutate(
    n_trimmed = n_original - n_remaining,
    original_p_fmt = format.pval(original_p, digits = 3),
    trimmed_p_fmt = format.pval(trimmed_p, digits = 3),
    sentence = sprintf(
      "To assess sensitivity to extreme values, we repeated the primary paired test after excluding the top 5%% of pre- and post-storm abundance values. Results were %s: the %s-based test remained %s after trimming (original p = %s, trimmed p = %s).",
      ifelse(result_stable, "stable", "not stable"),
      metric,
      ifelse(trimmed_p < 0.05, "significant", "non-significant"),
      original_p_fmt, trimmed_p_fmt
    )
  )

write.csv(outlier_trim_summary, here("output", "GoM_level", "stats", "pre_post", strand, paste0("GoM_outlier_trim_summary_", strand, ".csv")), row.names = FALSE)
