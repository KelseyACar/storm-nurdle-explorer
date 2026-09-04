# This script builds the shiny app Server using the SuperZip example as a framework
  # built from the structure of: https://github.com/rstudio/shiny-examples/tree/main/063-superzip-example
  # Run: shiny::runApp()
  # to print linenumber of errors: source("server.R")

# Understanding 
  # This is where the logic lives!
  # server reads input$ values, does the work, and writes to output$ placeholders created in ui.R
  # it runs continuously, reacting to user input changes
  # renderLeaflet() - runs once to build base map canvas
  # leafletProxy() - modifies existing map w/out full redraw (layer updates go here)
  # reactive() - reusable computation that recalculates when inputs change, called with () line a fn
  # observe() - runs automatically when inputs change, no return value
  # reactiveVal() - stores a single value that can be updated and read reactively
  # observeEvent() - like observe but only triggers on one specific input

# ui.R declares a placeholder by:
  # plotOutput("my_plot") (example)
# server.R fills the placeholder by:
  # output$my_plot <- renderPlot({...})



# to call specific line number for an error run: source("server.R")
# to parse syntax errors parse("serve.R")

library(leaflet)
library(RColorBrewer)
library(scales)
library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)


# survey data is ordered by weighted_total_amount so higher count sites render on top (from SuperZip's centile ordering)
surveydata <- master_df[order(master_df$weighted_total_amount), ] 

# storms in chronological order so site-level change plots match timeseries plots
storm_date_order <- master_df %>%
  filter(!is.na(name)) %>%
  distinct(name, full_storm_start) %>%
  arrange(full_storm_start)


# server function starts

function(input, output, session) {

  # ── Reactive helpers ──────────────────────────────────────────────────────
  
  # translates input$strandline selection to actual column name in master_df
  strandCol <- reactive({
    input$strandline %||% "weighted_total_amount" # %||% is the null coalescing operator from rlang -> it means use left side but if NULL use right side
  })
  
  # TRUE when a specific storm is selected (not "All storms") - drive site color changes (state 1 vs state 2/3)
  # drives the State 1 vs State 2/3 switch throughout the server
  storm_active <- reactive({
    val <- input$active_storm
    !is.null(val) && !is.na(val) && val != "All Paired Storms"
  })
  
  # extracts just the storm name string (drops year) for filtering
  sel_storm_name <- reactive({
    strsplit(input$active_storm, " ")[[1]][1]
  })
  
  # stores the layerId of the currently clicked site marker
  # NULL until a site is clicked; cleared when storm selection changes
  selected_site <- reactiveVal(NULL)
  
  # resolves change type col, no stat consideration needed bc median is default on map, but mean and median selectable in site explorer tab
  change_type_col <- reactive({
    switch(input$strandline,
           "weighted_total_amount"   = "change_type_median_total",
           "weighted_new_strandline" = "change_type_median_new",
           "weighted_old_strandline" = "change_type_median_old"
    )
  })
  
  # resolves pcs_change col from mp_change based on strandline + stat
  pcs_change_col <- function(strand_str, stat_str) {
    paste0("pcs_change_", stat_str, "_", strand_str)
  }
  
  # resolves abundance col from summary dfs based on strand + stat
  abund_col <- function(strand_str, stat_str) {
    paste0(stat_str, "_abundance_", strand_str)
  }
  
  # resolves SE/IQR stat error cols
  se_col <- function(strand_str) paste0("se_", strand_str)
  q1_col <- function(strand_str) paste0("q1_abundance_", strand_str)
  q3_col <- function(strand_str) paste0("q3_abundance_", strand_str)
  
  
  # ── Base map ──────────────────────────────────────────────────────────────

  # Create base map once w/ GoM bounding box (from global.R) - only static elements go here (tiles, view, scale)
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      setView(lng = GoM_lon, lat = GoM_lat, zoom = GoM_zoom) %>%
      addScaleBar(position = "bottomright")
  })
  
  
  
  # ── Survey site markers ───────────────────────────────────────────────────
  # This observer is responsible for maintaining the circles and legend according to the variables the user has chosen to map to color and size.
 
  # state 1: all sites shown in a yellowish color
  # state 2/3 (storm selected): only sites that interact w/ selected storm are show, color of site reactive to post-storm change type
   observe({
     if (!isTRUE(storm_active())) {
       # state 1
       site_points <- surveydata %>%
         distinct(site_id, site_name, cluster_lon, cluster_lat) %>%
         mutate(fill = "#ffc500")

       print(sapply(site_points, class))   # ← add here
     } else {
       # state 2/3: change type color interactive with strandline selection
       ct_col <- change_type_col()
        site_points <- mp_change %>%
         filter(name == sel_storm_name()) %>%
         distinct(site_id, cluster_lon, cluster_lat, change_type = .data[[ct_col]]) %>%
         filter(!is.na(change_type)) %>%
         left_join(surveydata %>% distinct(site_id, site_name), by = "site_id") %>%
          mutate(fill = ifelse(!is.na(change_type) & change_type %in% names(colors),
                               colors[change_type],
                               "#A8A8A8"))
     }

    # marker shape, sizes, etc
    leafletProxy("map") %>%
      clearGroup("sites") %>%
      addCircleMarkers(
        lng = site_points$cluster_lon,
        lat = site_points$cluster_lat,
        group = "sites",
        radius = 6,
        color = "#333333",
        weight = 0.8,
        fillColor = site_points$fill,
        fillOpacity = 0.8,
        stroke = TRUE,
        layerId = site_points$site_id,
        label = paste(site_points$site_id, site_points$site_name, sep = " | ")
      )
  })
  
  
  
  # ── Map legend ────────────────────────────────────────────────────────────
    # separate observe so legend updates independently from markers
    # State 1: no legend (neutral site color, so needs no key)
    # State 2/3: accumulation / depletion / no_change legend - site color changes
  observe({
    proxy <- leafletProxy("map")
    if (!isTRUE(storm_active())) {
      proxy %>% removeControl("colorLegend")
    } else {
      proxy %>% addLegend(
        "bottomleft",
        colors = unname(colors[c("accumulation", "depletion", "no_change")]),
        labels = c("Accumulation", "Depletion", "No change"),
        title = "Post-storm MP Redistribution",
        layerId = "colorLegend"
      )
    }
  })
  
  
  
  # ── Nurdle facilities markers ───────────────────────────────────────────────────
  # this is a separate layer toggled by input$show_manufacturers checkbox
  observe ({
    proxy <- leafletProxy("map")
    if (input$show_manufacturers) {
      proxy %>%
        clearGroup("manufacturers") %>%
        addCircleMarkers(
          data = manu_df,
          lng = ~long,
          lat = ~lat,
          group = "manufacturers",
          radius = 4,
          color = "black",
          weight = 1,
          fillColor = "black",
          fillOpacity = 0.85,
          stroke = TRUE,
          popup = ~paste0("<b>", plant_name, "</b><br>", city, ", ", state),
          label = ~plant_name
        )
    } else {
      proxy %>% clearGroup("manufacturers")
    }
  })
  
  
  
  # ── Storm track layer ─────────────────────────────────────────────────────
    # "All Paired Storms": all track lines, no per-point markers to reduce clutter
    # specific storm: track line & per-day point markers w/ data popup
  observe({
    proxy <- leafletProxy("map")
    
    if (input$show_tracks)  { 
      if (!isTRUE(storm_active())) {
       
         # All Paired Storms: draw all tracks
        proxy %>%
          clearGroup("tracks") %>%
          addPolylines(
            data = storm_tracks_cat, 
            color = "#ff0000", 
            weight = ~track_weight,
            group = "tracks", 
            label = ~paste(name, "| Category: ", max_category)
            )
        
      } else {
        # draw selected storm only
      sel_name <-sel_storm_name()
      track <- storm_tracks_cat %>% filter(name == sel_name)
      points <- storm_points %>% filter(name == sel_name)
      
      proxy %>% 
        clearGroup("tracks") %>%
        addPolylines(
          data = track, 
          color = "#ff0000",
          weight = ~track_weight,
          group = "tracks",
          label = ~paste(name, "|Category: ", max_category)
          ) %>%
        addCircleMarkers(
          data = points,
          lng = ~lon,
          lat = ~lat,
          group = "tracks",
          radius = 4,
          color = "darkred",
          fillColor = ~colorNumeric("YlOrRd", domain = NULL)(as.numeric(as.Date(date_storm))), # track pts colored by storm date so we can visualize start, end and direction
          fillOpacity = 0.9,
          stroke = TRUE,
          weight = 1,
          popup = ~paste0(
          "<b>", name, "-", date_storm, "</b><br>", # <b> = bolds text, </b> = stops bold, <br> = line break
          "Status: ", usa_status, "<br>",
          "Saffir-Simpson Category: ", usa_sshs, "<br>",
          "Max ROCI: ", roci, "km<br>",
          "Wind: ", usa_wind, "kts<br>",
          "Pressure: ", usa_pres, "mb<br>",
          "Speed: ", storm_speed, "kts")
        )
      } 
      
      } else {
      proxy %>% clearGroup("tracks")
  }
    })

    
  # ── Site popup ─────────────────────────────────────────────────────
  # reactive to both selected_site() and input$strandline
  # separate observe so popup re-renders when strandline changes without requiring re-click
  observe({
    req(selected_site())
    
    sid <- selected_site()
    site <- site_popup_lookup[site_popup_lookup$site_id == sid, ]
    if (nrow(site) == 0) return()
    
    # pick correct mean columns based on strandline selection
    ctrl_col <- switch(input$strandline,
                       "weighted_total_amount" = "ctrl_mean_total",
                       "weighted_new_strandline" = "ctrl_mean_new",
                       "weighted_old_strandline" = "ctrl_mean_old"
    )
    pre_col  <- "storm_mean_pre" # total mean strandline only in site pop on map 
    post_col <- "storm_mean_post"
    
    fmt <- function(col) {
      v <- site[[col]]
      if (!is.null(v) && length(v) > 0 && !is.na(v)) sprintf("%.1f pcs", v) else "No data"
      }
    
    content <- as.character(tagList(
      tags$h4(site$site_name),
      tags$strong(paste("Site ID:", site$site_id)), tags$br(),
      tags$hr(),
      paste("Control mean:", fmt(ctrl_col)), tags$br(),
      paste("Pre-storm mean:", fmt(pre_col)), tags$br(),
      paste("Post-storm mean:", fmt(post_col)), tags$br(),
      tags$em("See Storm/Site Explorer tabs for full plots, median values, and strandline resolution.")
    ))
    
    leafletProxy("map") %>%
      clearPopups() %>%
      addPopups(lng = site$cluster_lon, lat = site$cluster_lat, content, layerId = sid)
  })
  

  # ── Click observer ────────────────────────────────────────────────────────
  # clears old popup, stores clicked site_id in selected_site(), fires popup
  # isolate() reads event$id without creating a reactive dependency on it
  observe({
    event <- input$map_shape_click
    if (is.null(event)) return()
    isolate({ selected_site(event$id) })
  })
  
  # clear selected site when active_storm changes so site plot resets
  observeEvent(input$active_storm, { selected_site(NULL) })

  
  
  # ── Tab 2: Regional Explorer ──────────────────────────────────────────────
    # plots from EDA scripts in project
  
  # 1) Yearly abundance by survey type (control, pre, post)
  output$plot_reg_yearly <- renderPlot({
    s    <- input$reg_strand
    stat <- input$reg_stat
    var  <- abund_col(s, stat)
    se   <- se_col(s)
    
    if (!var %in% names(combined_yearly)) return(NULL)
    
    ggplot(regional_byStormPhase,
           aes(x = year, y = .data[[var]], color = fine_storm_phase,
               fill = fine_storm_phase, group = fine_storm_phase)) +
      geom_ribbon(aes(ymin = .data[[var]] - .data[[se]],
                      ymax = .data[[var]] + .data[[se]]),
                  alpha = 0.2, color = NA, na.rm = TRUE) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 50000),
                         labels = comma) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values  = colors) +
      labs(title = "Yearly Abundance by Survey Type",
           x = NULL, y = paste(stat, "abundance (pcs/10min/person)"),
           color = "fine_storm_phase", fill = "fine_storm_phase") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # 2) Abundance by survey condition (control vs storm-influenced pooled across all storms)
  output$plot_regional_condition <- renderPlot({
    s    <- input$reg_strand
    stat <- input$reg_stat
    var  <- abund_col(s, stat)
    se   <- se_col(s)
    
    if (!var %in% names(regional_byStormPhase)) return(NULL)
    
    ggplot(regional_byCondition,
           aes(x = survey_type, y = .data[[var]], fill = survey_type)) +
      geom_col(width = 0.6) +
      geom_errorbar(aes(ymin = .data[[var]] - .data[[se]],
                        ymax = .data[[var]] + .data[[se]]),
                    width = 0.2) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 50000),
                         labels = comma) +
      scale_fill_manual(values = colors) +
      scale_x_discrete(labels = c("control" = "baseline conditions", "storm_influenced" = "storm-influenced conditions")) +
      labs(title = "Regional Abundance by Survey Conditions",
           x = NULL, y = paste(stat, "abundance (pcs/10min/person)")) +
      theme_bw(base_size = 18) +
      theme(legend.position = "none")
  })
  
  # 3) Abundance per storm pre vs post
  output$plot_reg_by_storm <- renderPlot({
    s    <- input$reg_strand
    stat <- input$reg_stat
    var  <- abund_col(s, stat)
    se   <- se_col(s)
    
    if (!var %in% names(storm_abund_byPhase)) return(NULL)
    
    ggplot(storm_abund_byPhase %>% filter(fine_storm_phase %in% c("pre", "post")),
           aes(x = name, y = .data[[var]], fill = fine_storm_phase)) +
      geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
      geom_errorbar(aes(ymin = .data[[var]] - .data[[se]],
                        ymax = .data[[var]] + .data[[se]]),
                    position = position_dodge(width = 0.7), width = 0.3, alpha = 0.4) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 50000),
                         labels = comma) +
      scale_fill_manual(values = colors) +
      labs(title = "Abundance by Storm: Pre vs. Post",
           x = NULL, y = paste(stat, "abundance (pcs/10min/person)"), fill = "Phase") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # 4) Directional change count per storm (# sites accumulating/depleting)
  output$plot_reg_direction_count <- renderPlot({
    s    <- input$reg_strand
    stat <- input$reg_stat
    ct   <- paste0("change_type_", stat, "_", s)
    
    if (!ct %in% names(mp_change)) return(NULL)
    
    storm_direction_counts <- mp_change %>%
      filter(!is.na(name), !is.na(.data[[ct]])) %>%
      count(name, change_type = .data[[ct]], name = "n_sites") %>%
      mutate(name = factor(name, levels = sort(unique(as.character(name)), decreasing = TRUE)),
             change_type = factor(change_type, levels = c("accumulation", "no_change", "depletion")))
    
    ggplot(storm_direction_counts, aes(x = name, y = n_sites, fill = change_type)) +
      geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
      scale_fill_manual(values = colors) +
      coord_flip() +
      labs(title = "Directional Change Instances per Storm",
           x = NULL, y = "Number of Sites", fill = "Change Type") +
      theme_bw(base_size = 18)
  })
  
  # 5) Rate of change boxplot
  output$plot_reg_rate <- renderPlot({
    req(nrow(survey_pairs_total) > 0)
    s        <- input$reg_strand
    rate_col <- paste0("rate_", s)
    if (!rate_col %in% names(survey_pairs_total)) return(NULL)
    
    ggplot(survey_pairs_total, aes(x = survey_type, y = .data[[rate_col]], fill = survey_type)) +
      geom_boxplot(outlier.size = 1.5, linewidth = 0.6, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
      scale_fill_manual(values = colors) +
      scale_x_discrete(labels = c("control" = "Control", "storm_influenced" = "Storm-influenced")) +
      labs(title = paste("Rate of Change: Control vs. Storm-influenced —", s, "strand"),
           x = NULL, y = "rate of change (nurdles/day)") +
      theme_bw(base_size = 18) +
      theme(legend.position = "none")
  })
  
  # 6) Direction of change bar chart
  output$plot_reg_direction <- renderPlot({
    req(nrow(survey_pairs_total) > 0)
    s       <- input$reg_strand
    dir_col <- paste0("direction_", s)
    if (!dir_col %in% names(survey_pairs_total)) return(NULL)
    
    ggplot(survey_pairs_total, aes(x = survey_type, fill = .data[[dir_col]])) +
      geom_bar(position = "fill", color = "white", linewidth = 0.3, width = 0.7) +
      scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
      scale_fill_manual(values = colors) +
      scale_x_discrete(labels = c("control" = "Control", "storm_influenced" = "Storm-influenced")) +
      labs(title = paste("Proportion of Directional Change under Control vs. Storm-influenced Survey Conditions —", s, "strand"),
           x = NULL, y = "number of survey pairs", fill = "Direction") +
      theme_bw(base_size = 18)
  })
  
  
  # ── Tab 3: Storm Explorer ─────────────────────────────────────────────────
  
  # Main plot: prompt OR storm x site abundance OR response window
  output$plot_storm_exp_main <- renderPlot({
    storm_sel <- input$storm_exp_storm
    stat      <- input$storm_exp_stat %||% "median"
    s         <- "total"   # storm explorer always shows total in main plot; strandline in bottom plot
    
    if (storm_sel == "All Paired Storms") {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "Select a storm to explore redistribution at each paired site.",
                        size = 5, color = "gray50", hjust = 0.5) +
               theme_void())
    }
    
    sel_name <- strsplit(storm_sel, " ")[[1]][1]
    var <- abund_col(s, stat)
    se  <- se_col(s)
    
    if (input$storm_exp_view == "phase") {
      # abundance per linked site, pre vs post
      plot_df <- abund_bySiteStorm %>%
        filter(name == sel_name)
      
      if (!var %in% names(plot_df) || nrow(plot_df) == 0) return(NULL)
      
      ggplot(plot_df, aes(x = site_id, y = .data[[var]], fill = fine_storm_phase)) +
        geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
        geom_errorbar(aes(ymin = .data[[var]] - .data[[se]],
                          ymax = .data[[var]] + .data[[se]]),
                      position = position_dodge(width = 0.7), width = 0.3, alpha = 0.4) +
        scale_y_continuous(trans = pseudo_log_trans(base = 10),
                           breaks = c(0, 10, 100, 1000, 10000, 50000), labels = comma) +
        scale_fill_manual(values = colors) +
        labs(title = paste("MP Abundance — All Paired Sites:", sel_name),
             x = "Paired Sites", y = paste(stat, "MP abundance (pcs/10min/person)"), fill = "Phase") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      # response window by site
      plot_df <- storm_temp_window %>%
        filter(name == sel_name) %>%
        mutate(fine_storm_phase = factor(as.character(fine_storm_phase), levels = c("pre", "post")),
               fine_response = factor(fine_response, levels = c("acute", "subacute", "extended")))
      
      if (nrow(plot_df) == 0) return(NULL)
      
      ggplot(plot_df, aes(x = fine_response, y = .data[[var]], fill = fine_response)) +
        facet_wrap(~ fine_storm_phase, labeller = nurdle_labeller) +
        geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
        geom_errorbar(aes(ymin = .data[[var]] - .data[[se]],
                          ymax = .data[[var]] + .data[[se]]),
                      position = position_dodge(width = 0.7), width = 0.3, alpha = 0.4) +
        scale_y_continuous(trans = pseudo_log_trans(base = 10),
                           breaks = c(0, 10, 100, 1000, 10000, 50000), labels = comma) +
        scale_fill_manual(values = colors) +
        labs(title = paste("MP Temporal Response:", sel_name),
             x = "Response Window", y = "mean MP abundance (pcs/10min/person)", fill = "Response Window") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
  
  # Pieces changed per storm
  output$plot_storm_pcs <- renderPlot({
    storm_sel <- input$storm_exp_storm
    req(storm_sel != "All Paired Storms")
    sel_name <- strsplit(storm_sel, " ")[[1]][1]
    
    stat <- input$storm_exp_stat %||% "median"
    s    <- "total"
    pcs  <- pcs_change_col(s, stat)
    ct   <- paste0("change_type_", stat, "_", s)
    
    if (!pcs %in% names(mp_change) || !ct %in% names(mp_change)) return(NULL)
    
    plot_df <- mp_change %>%
      filter(name == sel_name, !is.na(.data[[pcs]])) %>%
      mutate(change_type = .data[[ct]])
    
    if (nrow(plot_df) == 0) return(NULL)
    
    ggplot(plot_df, aes(x = site_id, y = .data[[pcs]], fill = change_type)) +
      geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(-5000, -1000, -100, 0, 100, 1000, 5000), labels = comma) +
      scale_fill_manual(values = colors) +
      labs(title = paste("Directional MP Redistribution:", sel_name),
           x = NULL, y = "MP change (pieces)", fill = "Change Type") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Strandline breakdown
  output$plot_storm_exp_strand <- renderPlot({
    storm_sel <- input$storm_exp_storm
    req(storm_sel != "All Paired Storms")
    req(input$storm_exp_view)
    sel_name <- strsplit(storm_sel, " ")[[1]][1]
    stat     <- input$storm_exp_stat %||% "median"
    
    if (input$storm_exp_view == "phase") {
      
      plot_df <- abund_bySiteStorm %>%
        filter(name == sel_name) %>%
        pivot_longer(
          cols = c(abund_col("new", stat), abund_col("old", stat)),
          names_to = "strandline_type", values_to = "strand_abund"
        ) %>%
        mutate(
          se = if_else(grepl("new", strandline_type), .data[[se_col("new")]], .data[[se_col("old")]]),
          strandline_type = if_else(grepl("new", strandline_type), "new", "old")
        )
      
      if (nrow(plot_df) == 0) return(NULL)
      
      ggplot(plot_df, aes(x = site_id, y = strand_abund, fill = fine_storm_phase)) +
        facet_wrap(~ strandline_type, labeller = nurdle_labeller) +
        geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
        geom_errorbar(aes(ymin = strand_abund - se, ymax = strand_abund + se),
                      position = position_dodge(width = 0.7), width = 0.3, alpha = 0.4) +
        scale_y_continuous(trans = pseudo_log_trans(base = 10),
                           breaks = c(0, 10, 100, 1000, 10000, 50000), labels = comma) +
        scale_fill_manual(values = colors) +
        labs(title = paste("Abundance per Strandline:", sel_name),
             x = "Paired Sites", y = paste(stat, "abundance (pcs/10min/person)"), fill = "Phase") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    } else {
      # response window strandline breakdown
      plot_df <- storm_temp_window %>%
        filter(name == sel_name) %>%
        mutate(fine_storm_phase = factor(as.character(fine_storm_phase), levels = c("pre", "post")),
               fine_response    = factor(fine_response, levels = c("acute", "subacute", "extended"))) %>%
        pivot_longer(
          cols      = c(abund_col("new", "mean"), abund_col("old", "mean")),
          names_to  = "strandline_type", values_to = "strand_abund"
        ) %>%
        mutate(
          se              = if_else(grepl("new", strandline_type), .data[[se_col("new")]], .data[[se_col("old")]]),
          strandline_type = if_else(grepl("new", strandline_type), "new", "old")
        )
      
      if (nrow(plot_df) == 0) return(NULL)
      
      ggplot(plot_df, aes(x = fine_response, y = strand_abund, fill = fine_response)) +
        facet_grid(strandline_type ~ fine_storm_phase, labeller = nurdle_labeller) +
        geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
        geom_errorbar(aes(ymin = strand_abund - se, ymax = strand_abund + se),
                      position = position_dodge(width = 0.7), width = 0.3, alpha = 0.4) +
        scale_fill_manual(values = colors) +
        labs(title = paste("Strandline Response Window:", sel_name),
             x = "Response Window", y = "mean abundance (pcs/10min/person)", fill = "Response Window") +
        theme_bw(base_size = 18) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
  })
  
  # Storm metrics card
  output$storm_metrics_card <- renderUI({
    req(input$storm_exp_storm != "All Paired Storms")
    sel_name <- strsplit(input$storm_exp_storm, " ")[[1]][1]
    m <- GoM_level_storm_info %>% filter(name == sel_name)
    if (nrow(m) == 0) return(NULL)
    
    tags$div(style = "background:#f8f9fa; border-left: 4px solid #cc6600;
                      padding: 12px 16px; border-radius: 4px; margin-bottom: 10px;",
             tags$h5(style = "margin-top:0;",
                     paste(m$name, "\u2014", m$primary_usa_status, "| Category", m$primary_usa_sshs)),
             fluidRow(
               column(3, tags$b("Duration:"),               paste(m$duration_days, "days")),
               column(3, tags$b("Active:"),                 paste(m$full_storm_start, "\u2192", m$full_storm_end)),
               column(3, tags$b("Proximal survey sites:"),  m$n_survey_sites),
               column(3, tags$b("Max category:"),           m$max_sshs)
             ),
             tags$br(),
             fluidRow(
               column(3, tags$b("Max wind:"),    paste(m$max_wind_kts, "kts")),
               column(3, tags$b("Mean wind:"),   paste(round(m$mean_wind_kts, 1), "kts")),
               column(3, tags$b("Min pressure:"), paste(m$min_usa_pres, "mb")),
               column(3, tags$b("Mean pressure:"), paste(round(m$mean_usa_pres, 1), "mb"))
             ),
             tags$br(),
             fluidRow(
               column(3, tags$b("Max speed:"),  paste(m$max_speed_kts, "kts")),
               column(3, tags$b("Mean speed:"), paste(round(m$mean_speed_kts, 1), "kts"))
             )
    )
  })
  
  
  # ── Tab 4: Site Explorer ──────────────────────────────────────────────────
  
  output$site_exp_label <- renderUI({
    req(input$site_exp_site)
    sname <- master_df %>% filter(site_id == input$site_exp_site) %>%
      pull(site_name) %>% first()
    tags$p(tags$strong(sname), style = "font-size: 13px; margin-top: 4px;")
  })
  
  # Full time series
  output$plot_site_timeseries <- renderPlot({
    
    if (input$site_exp_site == "") {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "Select a site to explore redistribution across paired storms.",
                        size = 5, color = "gray50", hjust = 0.5) +
               theme_void())
    }
    
    site       <- input$site_exp_site
    s_col      <- input$site_exp_strand
    site_label <- master_df %>% filter(site_id == site) %>% pull(site_name) %>% first()
    
    site_control <- master_df %>%
      filter(site_id == site, survey_type == "control") %>%
      distinct(date_nurdle, .keep_all = TRUE) %>%
      mutate(date_nurdle = as.Date(date_nurdle))
    
    site_storm <- master_df %>%
      filter(site_id == site, survey_type == "storm_influenced") %>%
      distinct(date_nurdle, name, fine_storm_phase, .keep_all = TRUE) %>%
      mutate(date_nurdle = as.Date(date_nurdle))
    
    if (nrow(site_control) == 0 && nrow(site_storm) == 0) return(NULL)
    
    storm_markers <- site_storm %>%
      group_by(name) %>%
      summarise(event_date = as.Date(first(date_storm_nearest_site)), .groups = "drop")
    
    site_all <- bind_rows(site_control, site_storm) %>%
      arrange(date_nurdle) %>%
      distinct(date_nurdle, .keep_all = TRUE)
    
    ggplot() +
      geom_rect(data = storm_markers,
                aes(xmin = event_date - days(28), xmax = event_date + days(28),
                    ymin = -Inf, ymax = Inf),
                fill = colors["post"], alpha = 0.05) +
      geom_vline(data = storm_markers,
                 aes(xintercept = event_date),
                 linetype = "dashed", color = colors["post"], linewidth = 0.6, alpha = 0.5) +
      geom_text(data = storm_markers,
                aes(x = event_date, y = Inf, label = name),
                angle = 90, vjust = -0.3, hjust = 1.1, size = 2.8, color = colors["post"]) +
      geom_line(data = site_all,
                aes(x = date_nurdle, y = .data[[s_col]]),
                color = "gray50", linewidth = 0.5, alpha = 0.7) +
      geom_point(data = site_control,
                 aes(x = date_nurdle, y = .data[[s_col]]),
                 color = colors["control"], shape = 16, size = 2.5, alpha = 0.8) +
      geom_point(data = site_storm,
                 aes(x = date_nurdle, y = .data[[s_col]], color = fine_storm_phase),
                 shape = 17, size = 2.5, alpha = 0.9) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(0, 10, 100, 1000, 10000, 50000), labels = comma) +
      scale_color_manual(values = colors) +
      scale_x_date(date_labels = "%b %Y", date_breaks = "6 months") +
      labs(title    = paste0(site, ": ", site_label),
           subtitle = "Circle = control  |  Triangle = storm-influenced  |  Dashed line = storm event",
           x = NULL, y = "MP abundance (pcs/10min/person)", color = "Phase") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Pieces changed per site across all paired storms
  output$plot_site_pcs_change <- renderPlot({
    req(input$site_exp_site)
    
    site <- input$site_exp_site
    stat <- input$site_exp_stat %||% "median"
    s    <- switch(input$site_exp_strand,
                   "weighted_total_amount"   = "total",
                   "weighted_new_strandline" = "new",
                   "weighted_old_strandline" = "old")
    
    pcs <- pcs_change_col(s, stat)
    ct  <- paste0("change_type_", stat, "_", s)
    
    if (!pcs %in% names(mp_change)) return(NULL)
    
    plot_df <- mp_change %>%
      filter(site_id == site, !is.na(.data[[pcs]])) %>%
      mutate(change_type = .data[[ct]]) %>%
      left_join(storm_date_order %>% select(name, full_storm_start), by = "name") %>%
      mutate(name = factor(name,
                           levels = storm_date_order$name[storm_date_order$name %in% name]))
    
    if (nrow(plot_df) == 0) {
      return(ggplot() +
               annotate("text", x = 0.5, y = 0.5,
                        label = "No storm interaction data for this site",
                        size = 5, color = "gray50") +
               theme_void())
    }
    
    ggplot(plot_df, aes(x = name, y = .data[[pcs]], fill = change_type)) +
      geom_col(position = position_dodge(width = 0.7, preserve = "single"), width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.6) +
      scale_y_continuous(trans = pseudo_log_trans(base = 10),
                         breaks = c(-5000, -1000, -100, 0, 100, 1000, 5000), labels = comma) +
      scale_fill_manual(values = colors) +
      labs(title = paste("Directional MP Redistribution Across Linked Storms:", site, "-", stat),
           x = NULL, y = "MP change (pieces)", fill = "Change Type") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Storm closest-approach metrics table
  output$closest_metrics_table <- DT::renderDataTable({
    req(input$site_exp_site)
    
    closest_metrics %>%
      filter(site_id == input$site_exp_site) %>%
      left_join(storm_date_order %>% select(name, full_storm_start), by = "name") %>%
      mutate(name = factor(name,
                           levels = storm_date_order$name[storm_date_order$name %in% name])) %>%
      arrange(name) %>%
      select(
        Storm = name,
        Status = usa_status,
        `Cat.` = usa_sshs,
        `Wind (kts)` = usa_wind,
        `Pres. (mb)` = usa_pres,
        `Min dist (m)` = min_dist_to_track_m,
        `RMW (m)` = rmw_m,
        `ROCI (m)` = roci_m,
        `Speed (kts)` = storm_speed,
        `Dir. (°)` = storm_dir,
        `Land (km)` = dist2land
      )
  }, options = list(pageLength = 10, dom = "t", scrollX = TRUE),
  rownames = FALSE, width = "100%")
  
  
  # ── Tab 5: Raw Survey Data Explorer ──────────────────────────────────────────
  
  strandCol_data <- reactive({ input$strandline %||% "weighted_total_amount" })
  
  observe({
    states <- if (length(input$country) == 0) character(0) else {
      filter(master_df, country %in% input$country) %>%
        `$`("state") %>% unique() %>% sort()
    }
    stillSelected <- isolate(input$states[input$states %in% states])
    updateSelectizeInput(session, "states", choices = states,
                         selected = stillSelected, server = TRUE)
  })
  
  observe({
    storms <- if (length(input$states) == 0) character(0) else {
      master_df %>%
        filter(state %in% input$states, !is.na(name)) %>%
        mutate(storm_label = paste(name, storm_year)) %>%
        `$`("storm_label") %>% unique() %>% sort()
    }
    stillSelected <- isolate(input$storms_table[input$storms_table %in% storms])
    updateSelectizeInput(session, "storms_table", choices = storms,
                         selected = stillSelected, server = TRUE)
  })
  
  observe({
    conditions <- if (length(input$states) == 0) character(0) else {
      df <- master_df %>% filter(state %in% input$states)
      if (length(input$storms_table) > 0)
        df <- df %>% filter(paste(name, storm_year) %in% input$storms_table)
      df %>% `$`("fine_storm_phase") %>% as.character() %>% unique() %>% na.omit() %>% sort()
    }
    stillSelected <- isolate(input$conditions_table[input$conditions_table %in% conditions])
    updateSelectizeInput(session, "conditions_table", choices = conditions,
                         selected = stillSelected, server = TRUE)
  })
  
  output$survey_table <- DT::renderDataTable({
    df <- master_df %>%
      distinct(site_id, date_nurdle, fine_storm_phase, name, .keep_all = TRUE)
    
    if (length(input$country) > 0)
      df <- df %>% filter(country %in% input$country)
    if (length(input$states) > 0)
      df <- df %>% filter(state %in% input$states)
    if (length(input$storms_table) > 0)
      df <- df %>% filter(paste(name, storm_year) %in% input$storms_table)
    if (length(input$conditions_table) > 0)
      df <- df %>% filter(as.character(fine_storm_phase) %in% input$conditions_table)
    
    df <- df %>%
      filter(.data[[strandCol_data()]] >= input$minCount,
             .data[[strandCol_data()]] <= input$maxCount) %>%
      select(site_id, site_name, state, country, date_nurdle,
             fine_storm_phase, name, storm_year,
             weighted_total_amount, weighted_new_strandline, weighted_old_strandline)
    
    DT::datatable(df, options = list(pageLength = 25, dom = "tp"), escape = FALSE, rownames = FALSE)
  })
  
  
} # end server function
