# SCRIPT 5: MP Redistribution Visualizations
# Purpose: Visualize MP redistribution patterns from Script 4 analysis

# Load Libraries
# Read Data
# Config: strand and stat
# Section 1: Rate & Direction of Change - Regional (Levene's result)
# Section 2: Rate of Change per Storm: storm pairs vs control at linked sites
# Section 3: Storm-Specific Post-storm Directional Heterogeneity
# Section 4: Rate of Change per Site: storm pairs vs control at each site
# Section 5: Site-specific Post Storm directional heterogeneity
# Section 6: Storm Variance / Levene's Results - only sign. findings

# Prep Strand

# Load Libraries ----------------------------------------------------------
pacman::p_load("dplyr", "ggplot2", "ggrepel", "patchwork", "here")
source(here("functions.R"))


# Config ------------------------------------------------------------------
# set strand: must match strand used in Script 4: ("old" | "new" | "total")
strand <- "total" 
stat <- "median" # -> sections 2 & 5 only, all other variance is raw values


#  Read Data -------------------------------------------------------------
master_df <- read.csv(here("output", "master_df.csv")) %>%
  mutate(date_nurdle = as.Date(date_nurdle))

survey_pairs <- read.csv(here("output", "summary_tables", paste0("survey_pairs_", strand, ".csv"))) %>%
  mutate(date_nurdle = as.Date(date_nurdle))

mp_change <- read.csv(here("output", "summary_tables", "MP_delta_bySite.csv"))

per_storm_rate_dir <- read.csv(here("output", "storm_level", "stats", "rate_and_dir", strand, paste0("per_storm_rate_dir_", strand, ".csv")))


# Prep Strand -------------------------------------------------------------

# draw out relevant columns (for each strand, in case analysis of others is desired)
survey_pairs <- prep_strand(survey_pairs, strand)

change <- mp_change %>%
  mutate(direction = .data[[paste0("change_type_", stat, "_", strand)]])

# set max days b/w surveys to match storm-influenced limit (+/- 28d = 56)
max_days = 56

# subset data from survey pairs
plot_pairs <- survey_pairs %>%
  filter(survey_type == "storm_influenced" | (survey_type == "control" & days_between <= max_days))

storm_names <- unique(survey_pairs$name[survey_pairs$survey_type == "storm_influenced"])
site_ids <- unique(change$site_id)

# Section 1: Rate & Direction of Change - Regional  -----------------------------------------

# 1) Direction of change — storm vs. control (all pairs, w/ max days cap)
p_direction <- ggplot(plot_pairs, aes(x = survey_type, fill = direction)) +
  geom_bar(position = "dodge", color = "white", linewidth = 0.3, width = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = colors) +
  labs(title = paste0("Direction of Change: Control vs. Storm (", strand, " strandline)"),
      x = "Period Type", y = "Number of Survey Pairs", fill = "Direction") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(size = 13), plot.title = element_text(size = 16))

p_direction

ggsave(here("output", "GoM_level", "graphs", "rate_and_dir", strand, paste0("GoM_dir_change_", strand, ".pdf")),
       plot = p_direction, width = 8, height = 6, dpi = 300)



# 2) Rate of change boxplot — storm vs. control (all pairs, w/ max days cap)
gom_levene_p <- read.csv(here("output", "GoM_level", "stats", "rate_and_dir", strand, paste0("GoM_rate_dir_", strand, ".csv")))$levene_p


p_rate <- build_trend_plot(plot_pairs, x_var = "survey_type", var = "rate_of_change", group_var = "survey_type",
                           plot_type = "boxplot", colors = colors, log_scale = FALSE,
                           title = paste0("Rate of Change: Control vs. Storm (", strand, " strandline)"),
                           y_label = "Rate of Change (nurdles/day)")
p_rate


# compute SD per group directly from the data driving p_rate
sd_summary <- plot_pairs %>%
  group_by(survey_type) %>%
  summarise(sd_rate = round(sd(rate_of_change, na.rm = TRUE), 1), .groups = "drop") %>%
  mutate(x_pos = as.numeric(factor(survey_type, levels = c("control", "storm_influenced"))),
         label = paste0("SD = ", sd_rate))

p_rate <- p_rate +
  geom_text(data = sd_summary, aes(x = x_pos, y = min(plot_pairs$rate_of_change, na.rm = TRUE) * 1.05, label = label),
            inherit.aes = FALSE, size = 3.5, fontface = "italic")

p_rate
ggsave(here("output", "GoM_level", "graphs", "rate_and_dir", strand, paste0("GoM_rate_change_", strand, ".pdf")),
       plot = p_rate, width = 8, height = 6, dpi = 300)



# Section 2: Rate of Change - PerStorm --------
for (storm in storm_names) {
  
  storm_sub <- survey_pairs %>% filter(survey_type == "storm_influenced", name == storm)
  linked_sites <- unique(storm_sub$site_id)
  ctrl_sub <- survey_pairs %>% filter(survey_type == "control", site_id %in% linked_sites, days_between <= max_days)
  d <- bind_rows(storm_sub, ctrl_sub)
  if (nrow(d) == 0) next
  
  p <- build_trend_plot(d, x_var = "survey_type", var = "rate_of_change", group_var = "survey_type",
                        plot_type = "boxplot", colors = colors,
                        title = paste0("Rate of Change, Storm: ", storm, " (", strand, " strandline)"),
                        y_label = "Rate of Change (nurdles/day)")
  
  ggsave(here("output", "storm_level", "graphs", "rate_and_dir", strand, paste0("rate_change_", storm, "_", strand, ".pdf")),
         plot = p, width = 8, height = 6, dpi = 300)
}




# Section 3: Storm-Specific Post-storm Directional Heterogeneity (w N linked Sites)  --------------
for (storm in storm_names) {
  
  storm_sub <- survey_pairs %>% filter(survey_type == "storm_influenced", name == storm)
  linked_sites <- unique(storm_sub$site_id)
  ctrl_sub <- survey_pairs %>% filter(survey_type == "control", site_id %in% linked_sites, days_between <= max_days)
  d <- bind_rows(storm_sub, ctrl_sub)
  if (nrow(d) == 0) next
  
  p <- build_trend_plot(d, x_var = "survey_type", var = "rate_of_change", group_var = "survey_type",
                        plot_type = "boxplot", colors = colors, log_scale = FALSE,
                        title = paste0("Rate of Change, Storm: ", storm, " (", strand, " strandline)"),
                        y_label = "Rate of Change (nurdles/day)")
  
  # SD labels per group, same pattern as the GoM-level plot
  sd_summary <- d %>%
    group_by(survey_type) %>%
    summarise(sd_rate = round(sd(rate_of_change, na.rm = TRUE), 1), .groups = "drop") %>%
    mutate(x_pos = as.numeric(factor(survey_type, levels = c("control", "storm_influenced"))),
           label = paste0("SD = ", sd_rate))
  
  p <- p +
    geom_text(data = sd_summary, aes(x = x_pos, y = min(d$rate_of_change, na.rm = TRUE) * 1.05, label = label),
              inherit.aes = FALSE, size = 3, fontface = "italic")
  
  ggsave(here("output", "storm_level", "graphs", "rate_and_dir", strand, paste0("rate_change_", storm, "_", strand, ".pdf")),
         plot = p, width = 8, height = 6, dpi = 300)
}




# Section 4:  Rate of Change by Site (storm pairs vs. control at t --------
for (site in unique(survey_pairs$site_id)) {
  
  d <- survey_pairs %>% filter(site_id == site)
  if (nrow(d) == 0 || n_distinct(d$survey_type) < 2) next
  site_label <- unique(master_df$site_name[master_df$site_id == site])[1]
  site_number <- gsub("_cluster", "", site)
  
  p <- build_trend_plot(d, x_var = "survey_type", var = "rate_of_change", group_var = "survey_type",
                        plot_type = "boxplot", colors = colors, log_scale = FALSE,
                        title = paste0("Rate of Change, Site ", site_number, ": ", site_label, " (", strand, " strandline)"),
                        y_label = "Rate of Change (nurdles/day)")
  
  # SD labels per group, same pattern as storm/GoM level
  sd_summary <- d %>%
    group_by(survey_type) %>%
    summarise(sd_rate = round(sd(rate_of_change, na.rm = TRUE), 1), .groups = "drop") %>%
    mutate(x_pos = as.numeric(factor(survey_type, levels = c("control", "storm_influenced"))),
           label = paste0("SD = ", sd_rate))
  
  p <- p +
    geom_text(data = sd_summary, aes(x = x_pos, y = min(d$rate_of_change, na.rm = TRUE) * 1.05, label = label),
              inherit.aes = FALSE, size = 3, fontface = "italic")
  
  ggsave(here("output", "site_level", "graphs", "rate_and_dir", strand, paste0("rate_change_site_", site_number, "_", strand, ".pdf")),
         plot = p, width = 8, height = 6, dpi = 300)
}


# Section 5: Site-Specific Post-storm Directional Heterogeneity ----------------------
site_ids <- unique(change$site_id)

for (site in site_ids) {
  
  d <- change %>% filter(site_id == site)
  if (nrow(d) == 0) next
  
  d$direction <- as.character(d$direction)
  site_label <- unique(d$site_name)[1]
  site_number <- gsub("_cluster", "", site)
  
  p_sites <- ggplot(d, aes(x = direction, fill = direction)) +
    geom_bar(color = "white", linewidth = 0.3, width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_discrete(limits = c("accumulation", "depletion", "no_change")) +
    scale_fill_manual(values = colors) +
    labs(title = paste0("Site ", site_number, ": ", site_label, " (", stat, ", ", strand, " strandline)"),
         x = "Post-storm Directional Change", y = "Number of Storms", fill = "Direction") +
    theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
          legend.position = "bottom", plot.title = element_text(size = 16))
  
  ggsave(here("output", "site_level", "graphs", "directional", strand, paste0(site_number, "_response_", strand, ".pdf")),
         plot = p_sites, width = 8, height = 5, dpi = 300)
}




# Section 6: Storm Variance / Levene's Results ----------------------------------

# Storm-level: which storms show significant variance inflation (BH-adjusted)
sig_storms <- per_storm_rate_dir %>% filter(levene_sig_adj == TRUE) %>% pull(group)

storm_variance_facets <- do.call(rbind, lapply(storm_names, function(storm) {
  storm_sub <- survey_pairs %>% filter(survey_type == "storm_influenced", name == storm)
  linked_sites <- unique(storm_sub$site_id)
  ctrl_sub <- survey_pairs %>% filter(survey_type == "control", site_id %in% linked_sites, days_between <= max_days) %>%
    mutate(name = storm)
  bind_rows(storm_sub, ctrl_sub) %>%
    mutate(levene_flag = if_else(storm %in% sig_storms, "Significant (p < .05)", "Not Significant"))
}))

p_variance_facets <- ggplot(storm_variance_facets, aes(x = survey_type, y = rate_of_change, fill = survey_type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  facet_wrap(~ name + levene_flag, labeller = label_wrap_gen(width = 20)) +
  scale_fill_manual(values = colors) +
  labs(title = paste0("Rate of Change Spread by Storm, Flagged by Levene's Significance (", strand, " strandline)"),
       x = NULL, y = "Rate of Change (nurdles/day)") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "none", plot.title = element_text(size = 16))

ggsave(here("output", "storm_level", "graphs", "rate_and_dir", strand, paste0("variance_byStorm_levene_flagged_", strand, ".pdf")),
       plot = p_variance_facets, width = 16, height = 12, dpi = 300)

# highlight-only version: just the significant storms, easier to read for a defense slide
p_variance_sig <- storm_variance_facets %>% filter(name %in% sig_storms) %>%
  ggplot(aes(x = survey_type, y = rate_of_change, fill = survey_type)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
  facet_wrap(~ name) +
  scale_fill_manual(values = colors) +
  labs(title = paste0("Storms with Significant Variance Inflation (Levene's p < .05, BH-adjusted, ", strand, " strandline)"),
       x = NULL, y = "Rate of Change (nurdles/day)") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.position = "none", plot.title = element_text(size = 16))

# significance stars, pulled directly from adjusted p-values
levene_labels <- per_storm_rate_dir %>%
  filter(group %in% sig_storms) %>%
  rename(name = group) %>%
  mutate(label = case_when(
    levene_p_adj < 0.001 ~ "***",
    levene_p_adj < 0.01  ~ "**",
    levene_p_adj < 0.05  ~ "*",
    TRUE ~ ""
  ))

p_variance_sig <- p_variance_sig +
  geom_text(data = levene_labels, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 5, fontface = "bold")

ggsave(here("output", "storm_level", "graphs", "rate_and_dir", strand, paste0("variance_significantStorms_", strand, ".pdf")),
       plot = p_variance_sig, width = 10, height = 8, dpi = 300)

