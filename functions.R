 # Poster Theme ------------------------------------------------------------
theme_poster <- function() {
  theme_minimal(base_size = 16) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "gray40", linewidth = 0.4),
      axis.text = element_text(size = 14, color = "gray20"),
      axis.title = element_text(size = 16, face = "bold", color = "gray10"),
      legend.background = element_blank(),
      legend.key = element_rect(fill = "white", color = NA),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = "right",
      plot.title = element_text(size = 18, face = "bold", hjust = 0, color = "gray10"),
      plot.subtitle = element_text(size = 14, hjust = 0, color = "gray30"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Color Palette -----------------------------------------------------------
colors <- c(
  "accumulation" = "#D97B7B",
  "depletion" = "#5B7C99",
  "no_change" = "#A8A8A8",
  "control" = "black",
  "pre" = "#ffcccc",
  "storm_influenced" = "#ff0000",
  "storm" = "#ff0000",
  "post" = "#996666",
  "acute" = "#FFCC00",
  "subacute" = "#FF9900",
  "extended" = "#993300",
  "old" = "#33cc66",
  "new" = "#3fffff"
)




# Summary Table FNs --------------------------------------------------------------

# 1) Abundance summaries across strandlines for any table structure. 
  # To use: call df %>% distinct(..., .keep_all = TRUE) %>% group_by(summary variable) %>% summarize_abundance()
summarise_abundance <- function(df) { 
  df %>%
    summarise(
      # cumulative
    sum_abundance_total = round(sum(weighted_total_amount, na.rm = TRUE), 0),
    max_abundance_total = round(max(weighted_total_amount, na.rm = TRUE), 0),
    min_abundance_total = round(min(weighted_total_amount, na.rm = TRUE), 0),
    mean_abundance_total = round(mean(weighted_total_amount, na.rm = TRUE), 0),
    se_total = round(sd(weighted_total_amount, na.rm = TRUE) / sqrt(n()), 0), # mean +/- se
    median_abundance_total = round(median(weighted_total_amount, na.rm = TRUE), 0),
    q1_abundance_total = round(quantile(weighted_total_amount, 0.25, na.rm = TRUE), 0), # median +/- IQR
    q3_abundance_total = round(quantile(weighted_total_amount, 0.75, na.rm = TRUE), 0),
      # new strandline
    sum_abundance_new = round(sum(weighted_new_strandline, na.rm = TRUE), 0),
    max_abundance_new = round(max(weighted_new_strandline, na.rm = TRUE), 0),
    min_abundance_new = round(min(weighted_new_strandline, na.rm = TRUE), 0),
    mean_abundance_new = round(mean(weighted_new_strandline, na.rm = TRUE), 0),
    se_new = round(sd(weighted_new_strandline, na.rm = TRUE) / sqrt(n()), 0),
    median_abundance_new = round(median(weighted_new_strandline, na.rm = TRUE), 0),
    q1_abundance_new = round(quantile(weighted_new_strandline, 0.25, na.rm = TRUE), 0), 
    q3_abundance_new = round(quantile(weighted_new_strandline, 0.75, na.rm = TRUE), 0),
    # old strandline
    sum_abundance_old = round(sum(weighted_old_strandline, na.rm = TRUE), 0),
    max_abundance_old = round(max(weighted_old_strandline, na.rm = TRUE), 0),
    min_abundance_old = round(min(weighted_old_strandline, na.rm = TRUE), 0),
    mean_abundance_old = round(mean(weighted_old_strandline, na.rm = TRUE), 0),
    se_old = round(sd(weighted_old_strandline, na.rm = TRUE) / sqrt(n()), 0),
    median_abundance_old = round(median(weighted_old_strandline, na.rm = TRUE), 0),
    q1_abundance_old = round(quantile(weighted_old_strandline, 0.25, na.rm = TRUE), 0), 
    q3_abundance_old = round(quantile(weighted_old_strandline, 0.75, na.rm = TRUE), 0),
    n_surveys = n_distinct(id),
    .groups = "drop")
}





# Strandline Functions ----------------------------------------------------

# Fn to label directional change for all strandline options
label_direction <- function(df) {
  df %>% mutate(
    direction_old = case_when(delta_old > 0 ~ "accumulation", delta_old < 0 ~ "depletion", TRUE ~ "no_change"),
    direction_new = case_when(delta_new > 0 ~ "accumulation", delta_new < 0 ~ "depletion", TRUE ~ "no_change"),
    direction_total = case_when(delta_total > 0 ~ "accumulation", delta_total < 0 ~ "depletion", TRUE ~ "no_change")
  )
}


# Fn to rename strand-specific columns to generic names for use in analysis
prep_strand <- function(df, strand) {
  df[["rate_of_change"]] <- df[[paste0("rate_", strand)]]
  df[["direction"]] <- df[[paste0("direction_", strand)]]
  df
}



# Fn to pivot abundance dfs to wide pre-post for a given value column
pivot_pre_post <- function(df, value_col, phase_col = "fine_storm_phase") {
  if (nrow(df) == 0) return(data.frame(site_id = character(), name = character(),
                                       v_pre = numeric(), v_post = numeric()))
  
  wide <- df %>%
    filter(.data[[phase_col]] %in% c("pre", "post")) %>%
    select(site_id, name, phase = all_of(phase_col), val = all_of(value_col)) %>%
    pivot_wider(names_from = phase, values_from = val, names_prefix = "v_")
  
  if (!"v_pre" %in% names(wide)) wide$v_pre  <- NA_real_ # handle sites or storms without these data
  if (!"v_post" %in% names(wide)) wide$v_post <- NA_real_
  
  wide %>% filter(!is.na(v_pre), !is.na(v_post))
}



# reshape strandline data columns to long format
# value var is strand value column  
# se_var can be NA if the stat has no SE (like max values)

pivot_strand_long <- function(data, value_var, se_var = NA, value_name = "value"){
  
  out <- data %>%
    pivot_longer(cols = all_of(value_var), names_to = "strandline_type", values_to = "value") %>%
    mutate(strandline_type = case_when(
      strandline_type == value_var[1] ~ "new",
      strandline_type == value_var[2] ~ "old"
    ))
  
  if (!is.null(se_var[1])) {
    out <- out %>%
      mutate(se = case_when(
        strandline_type == "new" ~ .data[[se_var[1]]],
        strandline_type == "old" ~ .data[[se_var[2]]]
      ))
  }
  out 
}



# Statistical Functions ---------------------------------------------------


# 1) Fn to run Wilcoxon (rate) + Fisher (contingency < 5) /chi-square (contingency > 5) (direction), return one-row data frame
# based on normality, Welch's t-test if normal, Wilcoxon otherwise (rate of change)

run_rate_dir_tests <- function(df, group_label) {
  n_storm <- sum(df$survey_type == "storm_influenced", na.rm = TRUE)
  n_ctrl <- sum(df$survey_type == "control",  na.rm = TRUE)
  
  if (n_storm < 3 | n_ctrl < 3 | length(unique(df$survey_type[!is.na(df$rate_of_change)])) < 2) { # not powered lower than 3, so returns NAs
    return(data.frame(
      group = group_label,
      n_storm = n_storm,
      n_control = n_ctrl,
      median_storm = NA_real_,
      median_ctrl = NA_real_,
      pct_accum_storm = NA_real_,
      pct_deplete_storm = NA_real_,
      pct_unchanged_storm = NA_real_,
      pct_accum_ctrl = NA_real_,
      pct_deplete_ctrl = NA_real_,
      pct_unchanged_ctrl = NA_real_,
      shapiro_p_storm = NA_real_,
      shapiro_p_ctrl = NA_real_,
      rate_W_or_t = NA_real_,
      rate_p = NA_real_,
      rate_sig = NA,
      rate_test = NA_character_,
      rate_df = NA_real_,
      levene_p = NA_character_,
      levene_sig = NA,
      levene_F = NA_real_,
      levene_df_group = NA_real_,
      levene_df_residual = NA_real_,
      direction_p = NA_real_,
      direction_sig = NA,
      direction_test = NA_character_,
      direction_stat = NA_real_,
      direction_df = NA_real_
    ))
  }
  
  storm_rows <- df[df$survey_type == "storm_influenced", ]
  ctrl_rows  <- df[df$survey_type == "control", ]
  
  # normality check on rate of change for each group
    # shapiro.test req n>3 (already a req of function run)
    # lower power at small n, but test results are stored for transparency
  sw_storm <- tryCatch(
    shapiro.test(storm_rows$rate_of_change),
    error = function(e) list(p.value = 0)) # Shapiro fails when variance = 0 (i.e. values same), this defaults to Wilcoxon when shapiro fails (identical values are definitively non-normal)
  sw_ctrl <- tryCatch(
    shapiro.test(ctrl_rows$rate_of_change),
    error = function(e) list(p.value = 0))
  both_normal <- sw_storm$p.value > 0.05 & sw_ctrl$p.value > 0.05
  
  
  # rate test
  if (both_normal) {
    rate_result <- t.test(rate_of_change ~ survey_type, data = df, var.equal = FALSE)
    rate_stat <- rate_result$statistic
    rate_type <- "Welch t-test"
  } else {
    rate_result <- wilcox.test(rate_of_change ~ survey_type, data = df, exact = FALSE)
    rate_stat <- rate_result$statistic
    rate_type <- "Wilcoxon rank sum"
  }
  
  # variance test: Levene's test for equality of variance b/w storm and control rate distributions
  # tests whether storm conditions increase or decrease variance relative to control
  levene_result <- tryCatch({
    if (var(storm_rows$rate_of_change, na.rm = TRUE) == 0 | # tests whether zero variance first, if either group var = 0, NA is returned
        var(ctrl_rows$rate_of_change, na.rm = TRUE) == 0) {
      data.frame("Pr(>F)" = NA_real_, check.names = FALSE)
    } else {
      car::leveneTest(rate_of_change ~ as.factor(survey_type), data = df)
    }
  }, error = function(e) data.frame("Pr(>F)" = NA_real_, check.names = FALSE))
  
  levene_p <- levene_result[["Pr(>F)"]][1]
  levene_F <- if (!is.null(levene_result[["F value"]])) levene_result[["F value"]][1] else NA_real_
  levene_df_group <- if (!is.null(levene_result$Df)) levene_result$Df[1] else NA_real_
  levene_df_residual <- if (!is.null(levene_result$Df)) levene_result$Df[2] else NA_real_
  
  # direction test
  contingency <- table(df$survey_type, df$direction)
  
  # dir test
    # tryCatch handles cases that produce errors due to 1 row returns (i.e. unidirectional responses)
  dir_result <- tryCatch({
    if (any(contingency < 5)) {
      list(
        test = fisher.test(contingency, simulate.p.value = TRUE, B = 5000),
        type = "Fisher")
    } else {
      list(
        test = chisq.test(contingency),
        type = "Chi-square"
      )
    }
  }, error = function(e) {
    list(test = list(p.value = NA), type = "insufficient_variation")
  })
  
  direction_stat <- if (dir_result$type == "Chi-square") unname(dir_result$test$statistic) else NA_real_
  direction_df <- if (dir_result$type == "Chi-square") unname(dir_result$test$parameter) else NA_real_
  
  # define storm and ctrl rows
  n_storm_total <- nrow(storm_rows)
  n_ctrl_total  <- nrow(ctrl_rows)
  
  # append df with results
  data.frame(
    group = group_label,
    n_storm = n_storm,
    n_control = n_ctrl,
    median_storm = median(storm_rows$rate_of_change, na.rm = TRUE),
    median_ctrl = median(ctrl_rows$rate_of_change,  na.rm = TRUE),
    pct_accum_storm = mean(as.character(storm_rows$direction) == "accumulation", na.rm = TRUE) * 100,
    pct_deplete_storm = mean(as.character(storm_rows$direction) == "depletion",    na.rm = TRUE) * 100,
    pct_unchanged_storm = mean(as.character(storm_rows$direction) == "no_change",  na.rm = TRUE) * 100,
    pct_accum_ctrl = mean(as.character(ctrl_rows$direction)  == "accumulation", na.rm = TRUE) * 100,
    pct_deplete_ctrl = mean(as.character(ctrl_rows$direction)  == "depletion",    na.rm = TRUE) * 100,
    pct_unchanged_ctrl = mean(as.character(ctrl_rows$direction) == "no_change",   na.rm = TRUE) * 100,
    rate_W_or_t = rate_stat,
    rate_p = rate_result$p.value,
    rate_sig = rate_result$p.value < 0.05,
    rate_test = rate_type,
    rate_df = ifelse(rate_type == "Welch t-test", rate_result$parameter, NA_real_), # no df for Wilcoxon, but Welch's will print
    levene_p = levene_p,
    levene_sig = !is.na(levene_p) & levene_p < 0.05,
    levene_F = levene_F,
    levene_df_group = levene_df_group,
    levene_df_residual = levene_df_residual,
    shapiro_p_storm = sw_storm$p.value,
    shapiro_p_ctrl = sw_ctrl$p.value,
    direction_p = dir_result$test$p.value,
    direction_sig = dir_result$test$p.value < 0.05,
    direction_test = dir_result$type,
    direction_stat = direction_stat,
    direction_df = direction_df
  )
}


# 3) Fn to test whether indiv sites or storms display unidirection changes following storms
run_unidirectional_tests <- function(df, group_col, group_label) {
  
  groups <- unique(df[[group_col]]) # can use with site_id or storm name
  
  results <- do.call(rbind, lapply(groups, function(g) {
    sub <- df[df[[group_col]] == g, ]
    
    n_accum <- sum(sub$change_type == "accumulation", na.rm = TRUE)
    n_deplete <- sum(sub$change_type == "depletion", na.rm = TRUE)
    n_no_change <- sum(sub$change_type == "no_change", na.rm = TRUE)
    n_total <- nrow(sub)
    n_directional <- n_accum + n_deplete # excludes no_change from binomial test, per methods
    
    if (n_directional < 3) return(NULL) # not powered if < 3 storms linked to the site (or vice versa), returns NAs
    
    # unidirectional change test
    bt <- binom.test(n_accum, n_directional, p = 0.5)
    
    data.frame(
      group = g,
      group_type = group_label,
      n = n_total,
      n_accumulation = n_accum,
      n_depletion = n_deplete,
      pct_accumulation = n_accum / n_directional * 100,
      pct_no_change = n_no_change / n_total * 100,
      dominant_direction_of_change = case_when(
        n_accum == n_directional ~ "always_accumulation",
        n_accum == 0 ~ "always_depletion",
        n_accum == n_directional / 2 ~ NA_character_, # exact tie b/w accum. and dep. 
        n_accum > n_directional / 2 ~ "mostly_accumulation",
        TRUE ~ "mostly_depletion"
      ),
      binom_p = bt$p.value,
      binom_sig = bt$p.value < 0.05
    )
  }))
  
  results[!sapply(results, is.null), ]
}



# Fn to run normality + Wilcoxon + t-test + effect sizes on a paired pre-post vector
# requires n >= 3 pairs to run tests, otherwise fails silently and na produced
run_paired_tests <- function(post_vals, pre_vals) {
  deltas <- post_vals - pre_vals
  n <- sum(!is.na(deltas))
  
  if (n < 3) {  # not powered at < 3 
    return(data.frame(
      strand = strand,
      n_pairs = n,
      mean_pre = NA_real_,
      mean_post = NA_real_,
      median_pre = NA_real_,
      median_post = NA_real_,
      median_delta = NA_real_,
      shapiro_p = NA_real_,
      normal = NA,
      wilcox_W = NA_real_,
      wilcox_p = NA_real_,
      wilcox_sig = NA,
      ttest_t = NA_real_,
      ttest_p = NA_real_,
      ttest_df = NA_real_,
      ttest_sig = NA,
      r_rb = NA_real_, # effect size for Wilcoxon Ranked tests
      cohens_d = NA_real_ # effect size for t-test
    ))
  }
  
  shap <- tryCatch(
    shapiro.test(deltas),
    error = function(e) list(p.value = NA_real_))
  w <- wilcox.test(post_vals, pre_vals, paired = TRUE, exact = FALSE)
  tt <- t.test(post_vals, pre_vals, paired = TRUE)
  ttest_df <- unname(tt$parameter)
  W_stat <- w$statistic
  r_rb <- 1 - (2 * W_stat) / (n * (n + 1) / 2 + W_stat)  # if +: post tends to be > pre, if -: post tends to be < pre, 0 = no systemic diff
  cohens_d <- mean(deltas, na.rm = TRUE) / sd(deltas, na.rm = TRUE)
  
  data.frame(
    strand = strand,
    n_pairs = n,
    mean_pre = mean(pre_vals, na.rm = TRUE),
    mean_post = mean(post_vals, na.rm = TRUE),
    median_pre = median(pre_vals, na.rm = TRUE),
    median_post = median(post_vals, na.rm = TRUE),
    median_delta = median(deltas, na.rm = TRUE),
    shapiro_p = shap$p.value,
    normal = !is.na(shap$p.value) & shap$p.value > alpha,
    wilcox_W = W_stat,
    wilcox_p = w$p.value,
    wilcox_sig = w$p.value < alpha,
    ttest_t = tt$statistic,
    ttest_p = tt$p.value,
    ttest_df = ttest_df,
    ttest_sig = tt$p.value < alpha,
    r_rb = round(r_rb, 3),
    cohens_d = round(cohens_d, 3)
  )
}


# Fn to run LMER and return contrast (pre vs. post)
# choosing LMER bc data have a nested, non-indep. structure, i.e. sies can be linked to multiple storms and vice versa
# fixed effect = storm phase
# random intercepts: site + storm (GoM), site (storm-level), storm (site-level)
run_lmer <- function(df, random_formula) {
  if (nrow(df) < 6) {  # requires 3 pairs, phase_long has 2 rows (pre & post) per pair
    return(data.frame(estimate = NA_real_, p_value = NA_real_, sig = NA)) #set same n <3 power as other tests
  }
  
  tryCatch({
  lmer_data <- df %>% mutate(mp_log = log1p(mp_abundance)) # log transform bc LMER assumes normally dist. data
  fit <- lmer(as.formula(paste("mp_log ~ fine_storm_phase +", random_formula)), 
              data = lmer_data)
  # return NA if model fit is singular, do not go on to emmeans
  if (isSingular(fit)) return(data.frame(estimate = NA_real_, p_value = NA_real_, sig = NA))
  # get comparisons
  emm <- emmeans(fit, ~ fine_storm_phase, type = "response")
  contrast <- as.data.frame(summary(pairs(emm, adjust = "holm"), infer = c(TRUE, TRUE)))
  data.frame(estimate = contrast$estimate, p_value = contrast$p.value,
             sig = contrast$p.value < alpha)
  }, error = function(e) {
    data.frame(estimate = NA_real_, p_value = NA_real_, sig = NA)
  })
}





# Storm quadrant FN
# test storm quadrant data for sig patterns (Fisher's exact if <5 obs, otherwise Chi-square)
# Monte Carlo simulation
# randomly generate tables under null hypoth. of independence to see what proportion (the p value) of them produces a test stat as extreme or more exereme than results.
quad_test <- if (any(quad_contingency < 5)) { 
  list(test = "Fisher", result = fisher.test(quad_contingency, simulate.p.value = TRUE, B = 5000)) 
} else {
  list(test = "Chi-square", result = chisq.test(quad_contingency)) 
}



# Plot Functions ----------------------------------------------------------

# 0) choose plot helper (boxplot versus violin for spread based on n - useful for rate of change plots)
choose_plot_type <- function(data, group_var, threshold = 10) {
n_per_group <- table(data[[group_var]])
if (min(n_per_group) >= threshold) "violin" else "boxplot"
}


# 1) Trend plots, can be: line,, box, or violin
build_trend_plot <- function(data, x_var, var, se_var = NA, q1_var = NA, q3_var = NA,
                             title, y_label, plot_type = "line",
                             group_var = "fine_storm_phase", facet_var = NULL, 
                             colors = NULL, log_scale = TRUE, base_size = 14) {
  
  has_group <- !is.null(group_var) && group_var %in% names(data)
  has_se <- plot_type == "line" && !is.na(se_var)
  has_iqr <- plot_type == "line" && !is.na(q1_var) && !is.na(q3_var)
  
  base_aes <- if (has_group) {
    aes(x = .data[[x_var]], y = .data[[var]], color = .data[[group_var]],
        fill = .data[[group_var]], group = .data[[group_var]])
  } else {
    aes(x = .data[[x_var]], y = .data[[var]])
  }
  
  p <- ggplot(data, base_aes)
  
  p <- switch(plot_type,
              "line" = {
                if (has_se) {
                  p <- p + geom_ribbon(aes(ymin = .data[[var]] - .data[[se_var]],
                                           ymax = .data[[var]] + .data[[se_var]]),
                                       alpha = 0.2, color = NA)
                } else if (has_iqr) {
                  p <- p + geom_ribbon(aes(ymin = .data[[q1_var]], ymax = .data[[q3_var]]),
                                       alpha = 0.2, color = NA)
                }
                p + geom_line(linewidth = 1) + geom_point(size = 2)
              },
              "violin"  = p + geom_violin(alpha = 0.5) + 
                geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), alpha = 0.2) +
                stat_summary(fun = median, geom = "point", size = 4, color = "black", position = position_dodge(width = 0.9)),
              "boxplot" = p + geom_boxplot(alpha = 0.7, outlier.alpha = 0.5, width = 0.5),
              stop("plot_type must be 'line', 'violin', or 'boxplot'")
  )
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)))
  }
    
  p <- p +
    labs(title = title, x = NULL, y = y_label, fill = NULL, color = NULL) +
    theme_apa() +
    theme(
      text = element_text(family = "sans", size = base_size),
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size - 1, angle = 45, hjust = 1),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2)
    )
  
  if (log_scale) {
    p <- p + scale_y_continuous(trans = "log1p", labels = scales::comma)
  } else {
    p <- p + scale_y_continuous(labels = scales::comma)
  }
  
  if (has_group && !is.null(colors)) {
    p <- p + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
  }
  
  p
}


# 2) Bar/col plot: comparisons across a categorical x 
build_bar_plot <- function(data, x_var, y_var, fill_var = NULL, se_var = NA, q1_var = NA, q3_var = NA,
                           title, y_label, facet_var = NULL, facet_formula = NULL,
                           facet_scales = "fixed", colors = NULL, base_size = 14, hline_zero = FALSE) {
 
  dodge_width <- 0.7
  
  mapping <- aes(x = .data[[x_var]], y = .data[[y_var]])
  if (!is.null(fill_var)) mapping <- modifyList(mapping, aes(fill = .data[[fill_var]]))
  
  p <- ggplot(data, mapping) +
    geom_col(position = position_dodge2(width = 0.7, preserve = "single"), width = 0.7, color = "black")
  
  # mean_se
  if (!is.na(se_var)) {
    p <- p + geom_errorbar(aes(ymin = .data[[y_var]] - .data[[se_var]], ymax = .data[[y_var]] + .data[[se_var]]),
                           position = position_dodge(width = dodge_width), width = 0.2, alpha = 0.4)
  } 
  
  # median_IQR
  else if (!is.na(q1_var) && !is.na(q3_var)) {
    p <- p + geom_errorbar(aes(ymin = .data[[q1_var]], ymax = .data[[q3_var]]),
                           position = position_dodge(width = dodge_width), width = 0.2, alpha = 0.4)
  }
  
  if (hline_zero) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6)
  }
  
  # facet formula
  if (!is.null(facet_formula)) {
    p <- p + facet_grid(as.formula(facet_formula), scales = facet_scales)
  } else if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = facet_scales)
  }
  
  is_change_var <- grepl("pcs_change|_change$", y_var)
  
  p <- p +
    labs(title = title, x = NULL, y = y_label, fill = NULL) +
    theme_apa() +
    theme(
      text = element_text(family = "sans", size = base_size),
      axis.title = element_text(size = base_size + 2),
      axis.text = element_text(size = base_size - 1, angle = 45, hjust = 1),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      plot.title = element_text(size = base_size + 2)
    )
  
  p <- if (is_change_var) {
    p + scale_y_continuous(labels = scales::comma)
  } else {
    p + scale_y_continuous(trans = "log1p", labels = scales::comma)
  }
  
  if (!is.null(fill_var) && !is.null(colors)) {
    p <- p + scale_fill_manual(values = colors)
  }
  
  p
}



# per site timeseries plots of control & storm-influenced data 
build_site_timeseries <- function(site, master_df, colors, raw_var, window_days = 28) {
  
  site_label <- master_df %>% filter(site_id == site) %>% pull(site_name) %>% first()
  
  site_control <- master_df %>% filter(site_id == site, survey_type == "control") %>% 
    distinct(date_nurdle, .keep_all = TRUE) %>% mutate(date_nurdle = as.Date(date_nurdle))
  
  site_storm <- master_df %>% filter(site_id == site, survey_type == "storm_influenced") %>% 
    distinct(date_nurdle, name, fine_storm_phase, .keep_all = TRUE) %>% 
    mutate(date_nurdle = as.Date(date_nurdle))
  
  if (nrow(site_control) == 0 && nrow(site_storm) == 0) return(NULL)
  
  storm_markers <- site_storm %>% group_by(name) %>% summarise(event_date = as.Date(first(date_storm_nearest_site)), .groups = "drop")
  
  site_all <- bind_rows(site_control, site_storm) %>% arrange(date_nurdle) %>% distinct(date_nurdle, .keep_all = TRUE)
  
  p <- ggplot() +
    geom_vline(data = storm_markers, aes(xintercept = event_date), linetype = "dashed", color = colors["storm"], linewidth = 0.6, alpha = 0.4) +
    geom_rect(data = storm_markers, aes(xmin = event_date - days(window_days), xmax = event_date + days(window_days), ymin = 0, ymax = Inf), fill = colors["storm"], alpha = 0.05) +
    geom_text(data = storm_markers, aes(x = event_date, y = Inf, label = name), angle = 90, vjust = -0.3, hjust = 1.1, size = 4, color = colors["storm"]) +
    geom_line(data = site_all, aes(x = date_nurdle, y = .data[[raw_var]]), color = "gray50", linewidth = 0.5, alpha = 0.7) +
    geom_point(data = site_control, aes(x = date_nurdle, y = .data[[raw_var]]), color = colors["control"], shape = 16, size = 2.5, alpha = 0.8) +
    geom_point(data = site_storm, aes(x = date_nurdle, y = .data[[raw_var]], color = survey_type), shape = 17, size = 2.5, alpha = 0.9) +
    scale_color_manual(values = colors) + scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
    labs(title = paste0(site, " : ", site_label), subtitle = "Circle = control, Triangle = storm-influenced, Vertical line = storm event", x = NULL, y = "Abundance (pcs/10min/person)") +
    scale_y_continuous(trans = "log1p", labels = scales::comma) + 
    theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13))
  list(plot = p, site_label = site_label)
}



# 4) FN to add storm tracks to ggplot (weighted by max storm category)
add_storm_tracks <- function(data) {
  list(
    geom_sf(data = data, 
            aes(linewidth = factor(max_category)), 
            color = "red", 
            alpha = 0.6),
    scale_linewidth_manual(
      values = c("-3" = 0.1, "0" = 0.3, "1" = 0.5, "2" = 0.8, 
                 "3" = 1.0, "4" = 1.5),
      name = "Storm Category"
    )
  )
}



# 5) FN to plot strandlines per individual storms
plot_storm_strandlines <- function(storm_name, data) {
  
  storm_data <- data %>% filter(name == storm_name)
  
  if(nrow(storm_data) < 3) {
    cat("  - Skipping", storm_name, "- insufficient data\n")
    return(NULL)
  }
  
  p <- ggplot(storm_data, aes(x = site_id, y = diff, fill = change_type)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
    facet_wrap(~ strand, ncol = 1) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = paste("Storm", storm_name, ": Strandline Changes"),
      x = "Cluster",
      y = "Change in Abundance (post - pre)",
      fill = "Change Type"
    ) +
    theme_apa() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}





# Tide functions (scrap) ----------------------------------------------------------

# Function to characterize site based on geography
# characterize_site_spatial <- function(cluster_lon, cluster_lat) {
#   
#   # Known Gulf harbors/ports (major ones) # Used Claude AI for this info -> confirm accuaracy! 1.3.26
#   gulf_harbors <- tibble(
#     name = c("Galveston", "Houston Ship Channel", "Corpus Christi", 
#              "Port Aransas", "South Padre Island", "Brownsville",
#              "Tampa", "St. Petersburg", "Clearwater",
#              "Pensacola", "Mobile", "New Orleans", "Gulfport", "Biloxi",
#              "Panama City", "Apalachicola", "Cedar Key"),
#     lon = c(-94.79, -94.88, -97.39, -97.07, -97.16, -97.15,
#             -82.46, -82.63, -82.80,
#             -87.21, -88.04, -90.07, -89.09, -88.89,
#             -85.66, -84.98, -83.03),
#     lat = c(29.31, 29.73, 27.80, 27.84, 26.07, 25.90,
#             27.95, 27.77, 27.97,
#             30.40, 30.69, 29.95, 30.37, 30.40,
#             30.16, 29.73, 29.14)
#   )
#   
#   # Major river mouths
#   gulf_rivers <- tibble(
#     name = c("Mississippi", "Atchafalaya", "Brazos", "Trinity", "Sabine",
#              "Rio Grande", "Apalachicola", "Mobile", "Pearl"),
#     lon = c(-89.15, -91.34, -95.40, -94.72, -93.87,
#             -97.15, -84.98, -88.04, -89.63),
#     lat = c(29.26, 29.54, 28.93, 29.73, 29.73,
#             25.96, 29.73, 30.69, 30.20)
#   )
#   
#   # Calculate distances (km) -> note, my other measures are in meters (storms) but i can convert on the backside or change here if desired. 1.3.26
#   dist_to_harbor <- min(distHaversine(
#     c(cluster_lon, cluster_lat),
#     gulf_harbors[, c("lon", "lat")]
#   )) / 1000
#   
#   dist_to_river <- min(distHaversine(
#     c(cluster_lon, cluster_lat),
#     gulf_rivers[, c("lon", "lat")]
#   )) / 1000
#   
#   # Classify
#   site_type <- case_when(
#     dist_to_harbor < 5 ~ "harbor",
#     dist_to_river < 3 ~ "river_mouth",
#     dist_to_harbor < 15 ~ "near_harbor",
#     TRUE ~ "open_coast"
#   )
#   
#   exposure <- case_when(
#     site_type %in% c("harbor", "near_harbor") ~ "protected",
#     site_type == "river_mouth" ~ "intermediate",
#     TRUE ~ "exposed"
#   )
#   
#   tibble(
#     site_type = site_type,
#     exposure = exposure,
#     dist_to_harbor_km = dist_to_harbor,
#     dist_to_river_km = dist_to_river
#   )
# }

# Function to fetch tide predictions for a station and date range
# fetch_tide_data <- function(station_id, begin_date, end_date) {
#   
#   # Format dates for API
#   begin_str <- format(begin_date, "%Y%m%d")
#   end_str <- format(end_date, "%Y%m%d")
#   
#   # Build API URL for tide predictions (not observations)
#   url <- paste0(
#     "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?",
#     "product=predictions&", # using predicted rather than observed tides to minimize NAs when measurements were not taken during storms, etc. 
#     "application=NOS.COOPS.TAC.WL&",
#     "begin_date=", begin_str, "&",
#     "end_date=", end_str, "&",
#     "datum=MLLW&", # Mean Lower Low Water "baseline" to measure tide height from -> standard for navigation and NOAA charts (water is never lower than charts show). Other options include MSL = Mean Sea Level or NAVD88 - North American Vertical Datum 1988 (land elevation based -> maybe this one?)
#     "station=", station_id, "&",
#     "time_zone=gmt&",
#     "units=metric&",
#     "interval=hilo&",  # High/Low tides only
#     "format=json"
#   )
#   
#   # fetch data
#   resp <- GET(url)
#   content_data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
#   
#   # Check if we got predictions
#   if(is.null(content_data$predictions)) {
#     cat("! No data for station", station_id, "\n")
#     return(NULL)
#   }
#   
#   # Convert to tibble
#   tide_data <- as_tibble(content_data$predictions) %>%
#     mutate(
#       station_id = station_id,
#       datetime = as.POSIXct(t, format = "%Y-%m-%d %H:%M", tz = "UTC"),
#       date = as.Date(datetime),
#       tide_height_m = as.numeric(v),
#       tide_type = type  # H = high, L = low
#     ) %>%
#     select(station_id, datetime, date, tide_height_m, tide_type)
#   
#   cat("Station", station_id, ":", nrow(tide_data), "observations\n")
#   
#   return(tide_data)
# }


# Function to find nearest 3 stations
# find_nearest_stations <- function(cluster_lon, cluster_lat, stations, n = 3) {
#   
#   # create coordinate matrix
#   station_coords <- cbind(stations$lon, stations$lat)
#   
#   # calculate distances
#   dists <- distHaversine(
#     c(cluster_lon, cluster_lat),
#     station_coords
#   )
#   
#   # get n nearest indices
#   ord <- order(dists)[1:n]
#   
#   # create result tibble explicitly
#   result <- tibble(
#     station_id = as.character(stations$id[ord]),
#     station_name = as.character(stations$name[ord]),
#     dist_km = dists[ord] / 1000,
#     rank = 1:n
#   )
#   
#   return(result)
# }

