# SCRIPT 3: Linked DF Visualizations - control & storm-influenced
  # Purpose: EDA plots for figures
  # calls on build_trend-plot and build_bar_plot from Functions.r (line plots, bar plots, violin plots, box plots)

# Read Data: full df & summary dfs created in script 2  
# Config: strandline (total, new, old) & stat to plot (mean_se, median_IQR, sum, max) 
# Section 1: Regional Plots
  # cntrl vs storm surveys: overall MP, yearly MP (w/ and w/out sampling effort), yearly MP by storm phase, variance (raw abundance w/ median)
  # per all storms in region: MP
  # yearly temporal window: pre-post facet
# Section 2: Per storm
  # control v storm: raw variance w/ median
  # pre-post abund 
  # abund per linked site
  # change type by # sites
  # temp window abund: pre-post facet
# Section 3: By Site Plots
  # control v storm: raw variance w/ median
  # pre-post abund 
  # abund per linked storms
  # change type by # storms
  # temp window abund: pre-post facet
# Section 4: Per-site full time series plots (control + storm data)



# Load Libraries ----------------------------------------------------------
install.packages("rnaturalearthdata")
install.packages("remotes")
install.packages("jtools") # APA format package
install.packages("scales") # allows for breaking ggplot2 internal breaks and scales
remotes::install_github("ropensci/rnaturalearthhires")
pacman::p_load('lubridate', 'ggplot2',"ggrepel", "ggpattern", 'dplyr','sf','rnaturalearth','here','cowplot', 'scales', 'gridExtra','tidygeocoder', 'jtools') 

# Load functions
source(here("functions.R"))  


# Read Data -------------------------------------------------------------------
master_df <- read.csv(here('output', 'master_df.csv')) %>% 
  filter(fine_storm_phase != "during") %>%
  mutate(site_id = gsub("_cluster", "", site_id),
         site_id = factor(site_id, levels = unique(site_id[order(as.numeric(site_id))])),
         name = factor(name, levels = sort(unique(name))),
         survey_type = factor(survey_type, levels = c("control", "storm_influenced")),
         fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")),
         fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")))


regional_yearly_abund <- read.csv(here("output", "summary_tables", "regional_yearly_abund.csv")) 


abund_byCondition <- read.csv(here('output', 'summary_tables', 'regional_abund_bySurveyCondtion.csv')) %>%
  mutate(survey_type =factor(survey_type, levels = c("control", "storm_influenced")))


abund_byStormPhase <- read.csv(here('output', 'summary_tables', 'regional_abund_byStormPhase.csv')) %>%
  mutate(fine_storm_phase =factor(fine_storm_phase, levels = c("pre", "post")))


yearly_abund_byCondition <- read.csv(here('output', 'summary_tables', 'regional_yearly_abund_bySurveyCondition.csv'))  %>%
  mutate(survey_type = factor(survey_type, levels = c("control", "storm_influenced")))


yearly_abund_byStormPhase <- read.csv(here('output', 'summary_tables', 'regional_yearly_abund_byStormPhase.csv')) %>%
  filter(fine_storm_phase != 'during') %>%
  mutate(fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")))


abund_byStorm <- read.csv(here('output', 'summary_tables', 'regional_abund_byStorm.csv')) %>%
  mutate(name = factor(name, levels = sort(unique(name))),
         fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")))

site_abund_byStormPhase <- read.csv(here('output', 'summary_tables', 'site_abund_byStormPhase.csv')) %>%
  filter(fine_storm_phase != 'during') %>%
  mutate(fine_storm_phase = factor(fine_storm_phase, levels = c("control", "pre", "post")),
         site_id = gsub("_cluster", "", site_id),
         site_id = factor(site_id, levels = unique(site_id[order(as.numeric(site_id))])))


abund_bySiteStorm <- read.csv(here('output', 'summary_tables', 'regional_abund_bySiteStorm.csv')) %>% 
  mutate(site_id = gsub("_cluster", "", site_id),
         site_id = factor(site_id, levels = unique(site_id[order(as.numeric(site_id))])),
         name = factor(name, levels = sort(unique(name))))

mp_change <- read.csv(here('output', 'summary_tables', 'MP_delta_bySite.csv')) %>% 
  mutate(site_id = gsub("_cluster", "", site_id),
         site_id = factor(site_id, levels = unique(site_id[order(as.numeric(site_id))])),
         name = factor(name, levels = sort(unique(name))))


yearly_temp_resp <- read.csv(here('output', 'summary_tables', 'regional_yearly_temp_window_abund.csv')) %>% 
  filter(fine_storm_phase != "control" & fine_storm_phase != "during") %>%
  mutate(fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")),
         fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")))


storm_temp_resp <- read.csv(here('output', 'summary_tables', 'storm_temp_window_abund.csv')) %>% 
  filter(fine_storm_phase != "control" & fine_storm_phase != "during") %>%
  mutate(fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")),
         fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")),
         name = factor(name, levels = sort(unique(name))))


site_temp_resp <- read.csv(here('output', 'summary_tables', 'site_temp_window_abund.csv')) %>% 
  filter(fine_storm_phase != "control" & fine_storm_phase != "during") %>%
  mutate(site_id = gsub("_cluster", "", site_id), 
         site_id = factor(site_id, levels = unique(site_id[order(as.numeric(site_id))])),
         fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")),
         fine_storm_phase = factor(fine_storm_phase, levels = c("pre", "post")))


# Config --------------------------------------------------------------------
strand <- "total"   # options: total, new, old
stat <- "median"  # options: mean, median, max


# Derive column names from config (same pattern as scripts 4/6)
var <- paste0(stat, "_abundance_", strand)
se_var <- if (stat == "mean") paste0("se_", strand) else NA
q1_var <- if (stat == "median") paste0("q1_abundance_", strand) else NA
q3_var <- if (stat == "median") paste0("q3_abundance_", strand) else NA
y_label <- paste0(tools::toTitleCase(stat), " abundance (pcs/10min/person)")
change_stat <- if (stat == "max") "median" else stat   # no max-based change type exists; falls back to median
change_var <- paste0("change_type_", change_stat, "_", strand)
raw_var <- if (strand == "total") "weighted_total_amount" else paste0("weighted_", strand, "_strandline")
site_ids <- unique(master_df$site_id)
storm_names <- unique(abund_byStorm$name)


pcs_change_col <- paste0("pcs_change_", stat, "_", strand)
change_type_col <- paste0("change_type_", stat, "_", strand)

mp_change <- mp_change %>% 
  mutate(pcs_change_plot = .data[[pcs_change_col]],
         change_type = .data[[change_type_col]])

# Regional plots ----------------------------------------------------------
out_dir <- here('output', 'GoM_level', 'graphs', 'EDA', var, strand)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# standardized sum abundance over time (all cumulative pellet counts in region over years)
regional_yearly_abund <- regional_yearly_abund %>%
  filter(year != "2025") # partial year in dataset

p <- ggplot(regional_yearly_abund, aes(x = factor(year), y = sum_abundance_total, group = 1)) +
  geom_vline(xintercept = which(regional_yearly_abund$year %in% c(2019, 2020)),
             linetype = "dashed", color = "grey70") +
  annotate("text", x = which(regional_yearly_abund$year == 2019), y = 100000,
           label = "Formosa\nruling", size = 3.5, color = "grey40", hjust = -0.1) +
  annotate("text", x = which(regional_yearly_abund$year == 2020), y = 100000,
           label = "COVID-19\nonset", size = 3.5, color = "grey40", hjust = -0.1) +
  geom_line(color = "#4C4C4C", linewidth = 1) +
  geom_point(aes(color = year == 2019), size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "#D97B7B", "FALSE" = "#4C4C4C")) +
  geom_text(data = regional_yearly_abund %>% filter(sum_abundance_total > 0),
            aes(label = n_surveys, y = sum_abundance_total),
            vjust = -0.8, size = 4) +
  scale_x_discrete(name = "Year") +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                     labels = scales::comma) +
  labs(title = "Regional Yearly Total Pellet Abundance with Sampling Effort",
       y = "Sum Standardized Pellet Abundance") +
  theme_apa(x.font.size = 13, y.font.size = 13) +
  theme(plot.title = element_text(size = 16))
p
ggsave(here("output", "GoM_level", "graphs", "EDA", "regional_yearly_abund.pdf"),
       plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)




# pre-post median abundances by region (all storms pooled)
p <- build_bar_plot(abund_byStormPhase, x_var = "fine_storm_phase", y_var = var, fill_var = "fine_storm_phase", se_var = se_var, 
                    title = "Regional Abundance by Storm Phase", y_label = y_label, colors = colors)

p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                            breaks = c(0, 10, 100, 1000, 10000, 30000, 50000, 100000),
                            labels = scales::comma)
p
ggsave(file.path(out_dir, paste0("regional_", var, "_byStormPhase.pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)



# cntrl vs storm: overall MP
p <- build_bar_plot(abund_byCondition, x_var = "survey_type", y_var = var, fill_var = "survey_type", se_var = se_var, 
                    title = "Regional Abundance by Survey Condition", y_label = y_label, colors = colors)

p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000, 100000),
                     labels = scales::comma)
p
ggsave(file.path(out_dir, paste0("regional_", var, "_byCondition.pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)


# cntrl vs storm: yearly MP w/ sampling effort bar plot
p <- build_bar_plot(yearly_abund_byCondition, x_var = "year", y_var = var, fill_var = "survey_type", se_var = se_var,
                    title = "Regional Yearly Abundance with Sampling Effort", y_label = y_label, colors = colors)

# add sampling n
p <- p + geom_text(data = yearly_abund_byCondition %>% filter(!is.na(.data[[var]]), .data[[var]] > 0),
                   aes(label = n_surveys, y = .data[[var]]),
                   position = position_dodge2(width = 0.7, preserve = "single"),
                   vjust = -0.6, size = 4)

p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                     labels = scales::comma)
p
ggsave(file.path(out_dir, paste0("regional_yearly_", var, "_with_n.pdf")), plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)


# cntrl vs storm:  yearly MP no sampling effort
p <- build_trend_plot(yearly_abund_byCondition, x_var = "year", var = var, se_var = se_var, 
                      q1_var = q1_var, q3_var = q3_var, title = "Yearly Abundance by Survey Condition", y_label = y_label, 
                      plot_type = "line", group_var = "survey_type", colors = colors)
p
ggsave(file.path(out_dir, paste0("regional_yearly_", var, ".pdf")), plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)


# cntrl vs storm surveys: yearly MP by storm phase,
p <- build_trend_plot(yearly_abund_byStormPhase, x_var = "year", var = var, se_var = se_var, 
                      q1_var = q1_var, q3_var = q3_var, title = "Yearly Abundance by Storm Phase", y_label = y_label, 
                      plot_type = "line", group_var = "fine_storm_phase", colors = colors)

p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                            breaks = c(0, 10, 100, 1000, 10000, 30000, 50000, 100000),
                            labels = scales::comma)
p
ggsave(file.path(out_dir, paste0("regional_yearly_", var, ".pdf")), plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)


# cntrl vs storm surveys: overall variance violin (point is median)
p <- build_trend_plot(master_df %>% distinct(id, .keep_all = TRUE), x_var = "survey_type", var = raw_var, colors = colors,
                      group_var = "survey_type", plot_type = "violin", 
                      title = "Abundance Spread: Control vs. Storm-Influenced", y_label = y_label)
p <- p + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000, 100000),
                     labels = scales::comma)
p
ggsave(file.path(out_dir, paste0('regional_variance', var, '_violin_perCondition.pdf')), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)


# cntrl vs storm phase: variance violin per storm phase (point is median)
p <- build_trend_plot(master_df %>% distinct(id, .keep_all = TRUE), x_var = "survey_type", var = raw_var, colors = colors,
                      plot_type = "violin", title = "Abundance Spread", y_label = y_label)
# scale y-axis
p <- p + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                     labels = scales::comma)

p
ggsave(file.path(out_dir, paste0('regional_variance', var, '_violin_perStormPhase.pdf')), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)



# all storms: MP in region
p <- build_bar_plot(abund_byStorm %>% filter(fine_storm_phase %in% c("pre", "post")), x_var = "name", y_var = var, fill_var = "fine_storm_phase", 
                    se_var = se_var, q1_var = q1_var, q3_var = q3_var, title = "Abundance by Storm: Pre vs. Post", 
                    y_label = y_label, colors = colors)
p
ggsave(file.path(out_dir, paste0("regional_", var, "_byStorm.pdf")), plot = p, width = 10, height = 6, dpi = 300, device = cairo_pdf)


# all storms: change type proportions
change_by_storm <- mp_change %>% filter(!is.na(name)) %>% count(name, .data[[change_var]]) %>% rename(change_type = 2)

p <- ggplot(change_by_storm, aes(x = name, y = n, fill = change_type)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = colors) +
  labs(title = "Change Type Proportion by Storm", x = NULL, y = "Proportion of Site-Storm Interactions", fill = "Change Type") +
  theme_apa(x.font.size = 13, y.font.size = 13, legend.font.size = 13, facet.title.size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(file.path(out_dir, paste0("regional_", change_var, "_proportions.pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)


# all storms overall change magnitude (pcs_change, 0 baseline)
regional_pcs_by_storm <- mp_change %>%
  filter(!is.na(name)) %>%
  group_by(name) %>%
  summarise(pcs_change_plot = if (stat == "median") median(pcs_change_plot, na.rm = TRUE) else mean(pcs_change_plot, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(name = factor(name, levels = sort(unique(name), decreasing = TRUE)),
         change_type = if_else(pcs_change_plot >= 0, "accumulation", "depletion"))

p <- build_bar_plot(regional_pcs_by_storm, x_var = "name", y_var = "pcs_change_plot", fill_var = "change_type",
                    title = "Overall Directional Change per Storm",
                    y_label = "Pcs Change (post - pre)", colors = colors, hline_zero = TRUE) 
p
# log scale 
p <- p +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10), 
                     labels = scales::comma,
                     breaks = c(-2000, -1000, -500, -100, 0, 100, 500, 1000, 2000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11))
p

ggsave(file.path(out_dir, paste0("change_magnitude_byStorms_", change_var, stat, ".pdf")),
       plot = p, width = 8, height = max(6, length(unique(regional_pcs_by_storm$name)) * 0.3), device = cairo_pdf)



# all storms count of directional change type 
storm_direction_counts <- mp_change %>%
  filter(!is.na(name), !is.na(change_type)) %>%
  count(name, change_type, name = "n_sites") %>%
  mutate(name = factor(name, levels = sort(unique(name), decreasing = TRUE)),
         change_type = factor(change_type, levels = c("depletion", "no_change", "accumulation")))

p <- ggplot(storm_direction_counts, aes(x = name, y = n_sites, fill = change_type)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = colors) +
  coord_flip() +
  labs(title = "Directional Change Instances per Storm",
       x = NULL, y = "Number of Sites", fill = "Change Type") +
  theme_apa(x.font.size = 11, y.font.size = 11, legend.font.size = 11)

p
ggsave(file.path(out_dir, paste0('perStorm_number_sites_change_type', var, '.pdf')), plot = p, width = 12, height = 8, dpi = 300, device = cairo_pdf)


# response window: pre vs post dodged, post solid / pre faded, faceted by year (no control data)
p <- build_bar_plot(yearly_temp_resp,x_var = "fine_response", y_var = var, fill_var = "fine_response",
                    se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                    facet_formula = "fine_storm_phase ~ year",
                    title = "Regional Temporal Window Abundance: Pre vs Post",
                    y_label = y_label, colors = colors)
p <- p + 
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                     labels = scales::comma)
p
ggsave(file.path(out_dir, paste0('regional_response_window_', var, '.pdf')), plot = p, width = 12, height = 8, dpi = 300, device = cairo_pdf)




# Per Storm (no control data)---------------------------------------------------------------
out_dir <- here('output', 'storm_level', 'graphs', 'EDA', var, strand)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


for (storm in storm_names) {
  
  # control vs storm surveys per linked site per storm
  linked_sites <- unique(master_df$site_id[master_df$name == storm])
  
  d_raw <- master_df %>% filter(site_id %in% linked_sites) %>% distinct(id, .keep_all = TRUE)
  p <- build_trend_plot(d_raw, x_var = "survey_type", var = raw_var, group_var = "survey_type", plot_type = "violin", 
                        title = paste("Abundance Spread, Storm:", storm), y_label = "Abundance (pcs/10min/person)", colors = colors)
  p <- p + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                       labels = scales::comma)
  ggsave(file.path(out_dir, paste0("variance_", storm, ".pdf")), plot = p, width = 6, height = 6, dpi = 300, device = cairo_pdf)
  
  # abundance per storm (pre vs post, this storm only, no site breakdown)
  if (sum(abund_byStorm$name == storm) > 0) {
    p <- build_bar_plot(abund_byStorm %>% filter(name == storm), x_var = "fine_storm_phase", y_var = var, se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                        title = paste("Overall Abundance, Storm:", storm), y_label = y_label, fill_var = "fine_storm_phase", colors = colors)
    p <- p + 
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir,  paste0("overall_", var, "_", storm, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # abundance per site linked to this storm
  if (sum(abund_bySiteStorm$name == storm) > 0) {
    p <- build_bar_plot(abund_bySiteStorm %>% filter(name == storm), x_var = "site_id", y_var = var, se_var = se_var, q1_var = q1_var, q3_var = q3_var, 
                        title = paste("Abundance by Linked Site, Storm:", storm), y_label = y_label)
    p <- p + geom_col(fill = "#3B7C8C", position = position_dodge2(width = 0.7, preserve = "single"), width = 0.7)
    p <- p + 
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir,  paste0("linkedSites_", var, "_", storm, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # change type per storm (# sites per change type)
  d_change <- mp_change %>% filter(name == storm) %>% count(.data[[change_var]]) %>% rename(change_type = 1)
  if (nrow(d_change) > 0) {
    p <- build_bar_plot(d_change, x_var = "change_type", y_var = "n", fill_var = "change_type", 
                        title = paste("Change Type, Storm:", storm), y_label = "Number of Sites", colors = colors)
    p <- p + 
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0(change_var, "_", storm, ".pdf")), plot = p, width = 6, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # change type for all sites linked to each storm
  d <- mp_change %>% filter(name == storm, !is.na(pcs_change_plot))
  if (nrow(d) == 0) next
  d <- d %>% mutate(site_id = factor(site_id, levels = sort(unique(as.numeric(site_id)))))
  
  p <- build_bar_plot(d, x_var = "site_id", y_var = "pcs_change_plot", fill_var = "change_type",
                      title = paste0("Directional Change by Site, Storm: ", storm),
                      y_label = "Pcs Change (post - pre)", colors = colors, hline_zero = TRUE)
  p <- p + 
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                       labels = scales::comma)
  ggsave(file.path(out_dir, paste0("pcs_change_by_site_", storm, "_", var, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  
  
  # temporal window abundance, per storm, pre vs post faceted
  if (sum(storm_temp_resp$name == storm) > 0) {
    p <- build_bar_plot(storm_temp_resp %>% filter(name == storm), x_var = "fine_response", y_var = var, 
                        fill_var = "fine_response",
                        se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                        facet_var = "fine_storm_phase",
                        title = paste("Temporal Window Abundance: Pre vs Post, Storm:", storm),
                        y_label = y_label, colors = colors)
    p <- p + 
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0("response_window_", var, "_", storm, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
}




# Per Site ---------------------------------------------------------------
out_dir <- here('output', 'site_level', 'graphs', 'EDA', var, strand)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

phase_levels <- c("control", "pre", "post")   # single source of truth for ordering, reused below

for (site in site_ids) {
  
  # control vs storm surveys at this site
  d_raw <- master_df %>%
    filter(site_id == site, survey_type != "before") %>%
    distinct(id, .keep_all = TRUE)
  if (nrow(d_raw) > 0) {
    p <- build_trend_plot(d_raw, x_var = "survey_type", var = raw_var, group_var = "survey_type", plot_type = "violin",
                          title = paste("Abundance Spread, Site:", site), y_label = "Abundance (pcs/10min/person)", colors = colors)
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0(var, "variance_site_", site, ".pdf")), plot = p, width = 6, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # abundance at this site (control, pre and post)
  d_phase <- site_abund_byStormPhase %>%
    filter(site_id == site, fine_storm_phase != "before") %>%
    mutate(fine_storm_phase = factor(fine_storm_phase, levels = phase_levels))
  if (nrow(d_phase) > 0) {
    p <- build_bar_plot(d_phase, x_var = "fine_storm_phase", y_var = var, fill_var = "fine_storm_phase", se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                        title = paste("Overall Abundance, Site:", site), y_label = y_label, colors = colors)
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0("overall_", var, "_site_", site, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # abundance per storm linked to this site
  if (sum(abund_bySiteStorm$site_id == site) > 0) {
    p <- build_bar_plot(abund_bySiteStorm %>% filter(site_id == site), x_var = "name", y_var = var, se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                        title = paste("Abundance by Linked Storm, Site:", site), y_label = y_label)
    p <- p + geom_col(fill = "#3B7C8C", position = position_dodge2(width = 0.7, preserve = "single"), width = 0.7)
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir,  paste0("linkedStorms_", var, "_site_", site, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  # change type per site (# storms per change type)
  d_change <- mp_change %>% filter(site_id == site) %>% count(.data[[change_var]]) %>% rename(change_type = 1)
  if (nrow(d_change) > 0) {
    p <- build_bar_plot(d_change, x_var = "change_type", y_var = "n", fill_var = "change_type",
                        title = paste("Change Type, Site:", site), y_label = "Number of Storms", colors = colors)
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0(change_var, "_site_", site, ".pdf")), plot = p, width = 6, height = 6, dpi = 300, device = cairo_pdf)
  }
  
  
  # change type all storms linked to each site
  d <- mp_change %>% filter(site_id == site, !is.na(pcs_change_plot), !is.na(name))
  if (nrow(d) == 0) next
  d <- d %>% mutate(name = factor(name, levels = sort(unique(name))))
  
  p <- build_bar_plot(d, x_var = "name", y_var = "pcs_change_plot", fill_var = "change_type",
                      title = paste0("Directional Change by Storm, Site: ", site, " (", stat, ", ", strand, ")"),
                      y_label = "Pcs Change (post - pre)", colors = colors, hline_zero = TRUE) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13))
  p <- p +
    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                       breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                       labels = scales::comma)
  ggsave(file.path(out_dir, paste0("pcs_change_by_storm_site_", site, "_", stat, ".pdf")),
         plot = p, width = 8, height = 6)
  
  
  # temporal window abundance, this site, pre vs post faceted
  d_temp_resp <- site_temp_resp %>%
    filter(site_id == site, fine_storm_phase != "before") %>%
    mutate(fine_storm_phase = factor(fine_storm_phase, levels = phase_levels))
  if (nrow(d_temp_resp) > 0) {
    p <- build_bar_plot(d_temp_resp, x_var = "fine_response", y_var = var, fill_var = "fine_response",
                        se_var = se_var, q1_var = q1_var, q3_var = q3_var,
                        facet_var = "fine_storm_phase",
                        title = paste("Temporal Window Abundance: Pre vs Post, Site:", site),
                        y_label = y_label, colors = colors)
    p <- p +
      scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 30000, 50000),
                         labels = scales::comma)
    ggsave(file.path(out_dir, paste0("response_window_site_", var, "_", site, ".pdf")), plot = p, width = 8, height = 6, dpi = 300, device = cairo_pdf)
  }
}



# Section 4: Per-Site Timeseries ------------------------------------------

for (site in site_ids) {
  
  result <- build_site_timeseries(site, master_df, colors, raw_var)
  if (is.null(result)) next
  ggsave(file.path(out_dir, paste0(site, "_", strand, "_", result$site_label, ".pdf")), plot = result$plot, width = 14, height = 6, dpi = 300, device = cairo_pdf)
  
}



# Section 5: Spatial Visualizations ---------------------------------------

# build regional map
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- ne_states(country = c("United States of America", "Mexico"), returnclass = "sf")

# fixed map extent -- wide enough to show full Florida + more of Mexico/each state
lon_limits <- c(-98, -80)
lat_limits <- c(17, 31)

# ---- storm tracks, segmented so intensity can vary along a single storm's path ----
sf_df <- master_df %>% filter(!is.na(longitude_storm), !is.na(latitude_storm), !is.na(cluster_lon), !is.na(cluster_lat))

make_segment <- function(lon1, lat1, lon2, lat2) {
  st_linestring(matrix(c(lon1, lon2, lat1, lat2), ncol = 2))
}

storm_tracks_prepped <- sf_df %>%
  filter(!is.na(usa_sshs), !is.na(longitude_storm)) %>%
  arrange(name, date_storm) %>%
  group_by(name) %>%
  mutate(next_lon = lead(longitude_storm), next_lat = lead(latitude_storm)) %>%
  filter(!is.na(next_lon)) %>%  # drops each storm's final point, which has no "next" to segment to
  ungroup()

segment_geoms <- mapply(make_segment,
                        storm_tracks_prepped$longitude_storm,
                        storm_tracks_prepped$latitude_storm,
                        storm_tracks_prepped$next_lon,
                        storm_tracks_prepped$next_lat,
                        SIMPLIFY = FALSE)

storm_segments_sf <- storm_tracks_prepped %>%
  mutate(geometry = st_sfc(segment_geoms, crs = 4326)) %>%
  st_as_sf() %>%
  select(name, date_storm, usa_sshs, geometry)

# check what range of categories actually appears at the segment level before trusting the scale below
unique(storm_segments_sf$usa_sshs)

# 4) FN to add storm tracks to ggplot (weighted by category at each segment, not storm max)
add_storm_tracks <- function(data) {
  list(
    geom_sf(data = data,
            aes(linewidth = factor(usa_sshs)),
            color = "red",
            alpha = 0.6),
    scale_linewidth_manual(
      values = c("-5" = 0.05, "-4" = 0.1, "-3" = 0.1, "-2" = 0.15, "-1" = 0.2, "0" = 0.3,
                 "1" = 0.5, "2" = 0.8, "3" = 1.0, "4" = 1.5, "5" = 2.0),
      labels = c("-5" = "Unknown", "-4" = "Post-tropical", "-3" = "Disturbance",         # all possibilities from IBTrACS
                 "-2" = "Subtropical", "-1" = "Tropical Depression", "0" = "Tropical Storm",
                 "1" = "Hurricane 1", "2" = "Hurricane 2", "3" = "Hurricane 3",
                 "4" = "Hurricane 4", "5" = "Hurricane 5"),
      name = "Storm Category"
    )
  )
}

build_base_map <- function(world, states, storm_tracks_data) {
  ggplot() + geom_sf(data = world, fill = "gray95", color = "black", linewidth = 0.2) +
    geom_sf(data = states, fill = NA, color = "gray60", linewidth = 0.15) + add_storm_tracks(storm_tracks_data)
}

# MP manufacturers: read & spatial vector
manu_df <- read.csv(here('output', 'manufacturer_df.csv')) %>% filter(lat >= lat_limits[1] & lat <= lat_limits[2], long >= lon_limits[1] & long <= lon_limits[2])
manu_sf <- st_as_sf(manu_df, coords = c("long", "lat"), crs = 4326)

# pellet spatial vector
nurdle_sf <- mp_change %>%
  filter(!is.na(name)) %>%
  distinct(site_id, name, cluster_lon, cluster_lat, .keep_all = TRUE) %>%
  st_as_sf(coords = c("cluster_lon", "cluster_lat"), crs = 4326)

nurd <- nurdle_sf %>% mutate(change_type = .data[[change_var]])

cluster_lookup <- master_df %>%
  distinct(site_id, site_name, cluster_lon, cluster_lat) %>%
  arrange(site_id) %>%
  rename(lon = cluster_lon, lat = cluster_lat) %>%
  mutate(cluster_label = paste0(site_id, ": ", site_name))


#### A) Regional Map Plot no Manufacturers ####
p_all <- build_base_map(world, states, storm_segments_sf) +
  geom_point(data = cluster_lookup, aes(x = lon, y = lat, shape = "Survey Site"), fill = 'lightblue', size = 4, alpha = 0.4,
             position = position_jitter(width = 0.2, height = 0.1)) +
  scale_shape_manual(values = c("Survey Site" = 21), name = NULL) +
  coord_sf(xlim = lon_limits, ylim = lat_limits, expand = FALSE, clip = "on") +
  labs(title = "A) Regional Synthesized Data", x = "Longitude", y = "Latitude") +
  theme_apa(legend.use.title = TRUE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.margin = margin(0.15, 0.3, 0.4, 0.15, "cm"))

p_all

ggsave(here('output', 'GoM_level', 'graphs', 'maps', 'synthesized_data.pdf'), plot = p_all, width = 6, height = 7, dpi = 300, device = cairo_pdf)

# save legend separately
storm_legend <- get_legend(p_all)
ggsave(here('output', 'GoM_level', 'graphs', 'maps', 'storm_legend.pdf'), plot = storm_legend, width = 5, height = 6, dpi = 300, device = cairo_pdf)

cluster_legend <- cluster_lookup %>%
  mutate(ID = site_id) %>%
  select(ID, Location = site_name) %>%
  arrange(as.numeric(ID))

table_grob <- tableGrob(cluster_legend, rows = NULL, 
                        theme = ttheme_minimal(core = list(fg_params = list(cex = 1.0, fontface = 'plain'), bg_params = list(fill = 'white')),
                                               colhead = list(fg_params = list(cex = 1.4, fontface = "bold"), bg_params = list(fill = 'grey90'))))

ggsave(here('output', 'GoM_level', 'graphs', 'maps', 'site_legend.pdf'), plot = table_grob, width = 5, height = 6, dpi = 300, device = cairo_pdf)

# save w/out legend to decide how to piece together later
p_all_noLegend <- p_all + theme(legend.position = "none")
p_all_noLegend
ggsave(here('output', 'GoM_level', 'graphs', 'maps', 'synthesized_data_noLegend.pdf'), plot = p_all_noLegend, width = 6, height = 7, dpi = 300, device = cairo_pdf)



#### B) Regional Map w/ Manufacturers ####
p_all_manu <- p_all + geom_sf(data = manu_sf, shape = 22, fill = "black", size = 2, alpha = 0.6) + geom_sf(data = states, fill = NA, color = "gray60", linewidth = 0.15) +
  scale_shape_manual(values = c("Known Nurdle Facility" = 22, "Survey Site" = 21), name = NULL) +
  coord_sf(xlim = lon_limits, ylim = lat_limits, expand = FALSE, clip = "on") +
  labs(title = "Regional Synthesized Data", x = "Longitude", y = "Latitude") +
  theme_apa(legend.use.title = TRUE) +
  theme(legend.position = 'right', legend.box = 'vertical', axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
        legend.key.size = unit(0.5, "cm"), plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"))

ggsave(here('output', 'GoM_level', 'graphs', 'maps', 'synthesized_data_plus_manufacturers.pdf'), plot = p_all_manu, width = 12, height = 8, dpi = 300, device = cairo_pdf)



#### C) Maps per storm: drops sites not linked to particular storm, sites colored by post-storm change direction ####
# alias the change type column for nurd so it matches script structure
nurd <- nurd %>% mutate(change_type = .data[[change_type_col]])

# get storm tracks w/ ROCI info for mapping a buffer of spatial extent
storm_tracks_full <- read.csv(here("output", "GOM_storms.csv")) %>%
  filter(name %in% storm_names) %>%
  select(name, date, longitude, latitude, usa_roci) %>%
  filter(!is.na(longitude), !is.na(latitude), is.finite(usa_roci), usa_roci > 0) %>%
  mutate(roci_m = usa_roci * 1852)   # usa_roci is nautical miles; confirmed against master_df's roci_m

# buffer each track point by its ROCI, then union per storm into one continuous swath
storm_buffers <- storm_tracks_full %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(3857) %>%   # project to meters before buffering -- st_buffer's dist is in the
  # current CRS's units, and 4326 (lat/long degrees) would buffer wrong
  st_buffer(dist = .$roci_m)

roci_buffers <- storm_buffers %>%
  group_by(name) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%   # group_by(name) first, or
  # st_union() drops "name" entirely
  st_transform(4326)   # back to lat/long to match the other map layers

# C) Maps per storm: drops sites not linked to particular storm
for (storm in storm_names) {
  
  m <- storm_segments_sf %>% filter(name == storm)
  n <- nurd %>% filter(name == storm)
  r <- roci_buffers %>% filter(name == storm)
  
  if (nrow(m) == 0) next
  
  map1 <- build_base_map(world, states, m) +
    geom_sf(data = r, fill = "#a6b8c4", alpha = 0.25, color = NA) +
    geom_sf(data = n, aes(fill = change_type), shape = 21, size = 4, alpha = 0.8) +
    geom_text_repel(data = n %>% mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]), aes(x = lon, y = lat, label = site_id), size = 4, force = 2) +
    scale_fill_manual(values = colors, name = "MP Directional Change") +
    coord_sf(xlim = lon_limits, ylim = lat_limits, expand = FALSE) +
    labs(title = paste("Storm:", storm), x = "Longitude", y = "Latitude") +
    theme_apa(legend.use.title = TRUE, x.font.size = 13, y.font.size = 13, legend.font.size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13), plot.title = element_text(size = 14))
  
  ggsave(here("output", "storm_level", "graphs", "maps", paste0("map_of_", storm, ".pdf")), plot = map1, width = 12, height = 8, dpi = 300, device = cairo_pdf)
  
}
