# This script makes the output directory for my thesis project

# Run once to create the full output directory structure.
# source at the top of all scripts (except functions) to ensure dirs exist before writing

library(here)

dirs <- c(
  
  # Shared intermediate tables
  "output/summary_tables",
  
  # GoM level ---------------------------------------------------------------
  "output/GoM_level/stats/paired_surveys",
  
  "output/GoM_level/stats/rate_and_dir/old",
  "output/GoM_level/stats/rate_and_dir/new",
  "output/GoM_level/stats/rate_and_dir/total",
  
  "output/GoM_level/stats/pre_post/old",
  "output/GoM_level/stats/pre_post/new",
  "output/GoM_level/stats/pre_post/total",
  
  "output/GoM_level/stats/resp_window/old",
  "output/GoM_level/stats/resp_window/new",
  "output/GoM_level/stats/resp_window/total",
  
  "output/GoM_level/stats/storm_metrics",
  "output/GoM_level/stats/storm_metrics/old",
  "output/GoM_level/stats/storm_metrics/new",
  "output/GoM_level/stats/storm_metrics/total",
  
  "output/GoM_level/graphs/rate_and_dir/old",
  "output/GoM_level/graphs/rate_and_dir/new",
  "output/GoM_level/graphs/rate_and_dir/total",
  
  "output/GoM_level/graphs/pre_post/old",
  "output/GoM_level/graphs/pre_post/new",
  "output/GoM_level/graphs/pre_post/total",
  
  "output/GoM_level/graphs/resp_window/old",
  "output/GoM_level/graphs/resp_window/new",
  "output/GoM_level/graphs/resp_window/total",
  
  "output/GoM_level/graphs/storm_metrics/old",
  "output/GoM_level/graphs/storm_metrics/new",
  "output/GoM_level/graphs/storm_metrics/total",
  
  "output/GoM_level/graphs/EDA",
  "output/GoM_level/graphs/maps",
  "output/GoM_level/graphs/maps/by_storm",
  
  # Storm level -------------------------------------------------------------
  "output/storm_level/stats/rate_and_dir/old",
  "output/storm_level/stats/rate_and_dir/new",
  "output/storm_level/stats/rate_and_dir/total",
  
  "output/storm_level/stats/pre_post/old",
  "output/storm_level/stats/pre_post/new",
  "output/storm_level/stats/pre_post/total",
  
  "output/storm_level/stats/resp_window/old",
  "output/storm_level/stats/resp_window/new",
  "output/storm_level/stats/resp_window/total",
  
  "output/storm_level/stats/storm_metrics",
  "output/storm_level/stats/storm_metrics/old",
  "output/storm_level/stats/storm_metrics/new",
  "output/storm_level/stats/storm_metrics/total",
  
  "output/storm_level/graphs/rate_and_dir/old",
  "output/storm_level/graphs/rate_and_dir/new",
  "output/storm_level/graphs/rate_and_dir/total",
  
  "output/storm_level/graphs/pre_post/old",
  "output/storm_level/graphs/pre_post/new",
  "output/storm_level/graphs/pre_post/total",
  
  "output/storm_level/graphs/resp_window/old",
  "output/storm_level/graphs/resp_window/new",
  "output/storm_level/graphs/resp_window/total",
  
  "output/storm_level/graphs/storm_metrics/old",
  "output/storm_level/graphs/storm_metrics/new",
  "output/storm_level/graphs/storm_metrics/total",
  
  "output/storm_level/graphs/EDA",
  
  # Site level --------------------------------------------------------------
  "output/site_level/stats/rate_and_dir/old",
  "output/site_level/stats/rate_and_dir/new",
  "output/site_level/stats/rate_and_dir/total",
  
  "output/site_level/stats/pre_post/old",
  "output/site_level/stats/pre_post/new",
  "output/site_level/stats/pre_post/total",
  
  "output/site_level/stats/resp_window/old",
  "output/site_level/stats/resp_window/new",
  "output/site_level/stats/resp_window/total",
  
  "output/site_level/stats/storm_metrics",
  "output/site_level/stats/storm_metrics/old",
  "output/site_level/stats/storm_metrics/new",
  "output/site_level/stats/storm_metrics/total",
  
  "output/site_level/graphs/rate_and_dir/old",
  "output/site_level/graphs/rate_and_dir/new",
  "output/site_level/graphs/rate_and_dir/total",
  
  "output/site_level/graphs/pre_post/old",
  "output/site_level/graphs/pre_post/new",
  "output/site_level/graphs/pre_post/total",
  
  "output/site_level/graphs/resp_window/old",
  "output/site_level/graphs/resp_window/new",
  "output/site_level/graphs/resp_window/total",
  
  "output/site_level/graphs/directional/old",
  "output/site_level/graphs/directional/new",
  "output/site_level/graphs/directional/total",
  
  "output/site_level/graphs/storm_metrics/old",
  "output/site_level/graphs/storm_metrics/new",
  "output/site_level/graphs/storm_metrics/total",
  
  "output/site_level/graphs/EDA"
)

for (d in dirs) {
  dir.create(here(d), recursive = TRUE, showWarnings = FALSE)
}

message("Output directories created.")

