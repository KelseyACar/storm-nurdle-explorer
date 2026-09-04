# This script builds the shiny app UI using the SuperZip example as a framework
  # built from the structure of: https://github.com/rstudio/shiny-examples/tree/main/063-superzip-example
  
# Understanding
  # This file is a description of the layout: what to draw and where.
    # no logic
    # every Input widget here creates an input$name variable that the server can read
    # every Output widget here is a placeholder that the server fills

# ui.R declares a placeholder by:
  # plotOutput("my_plot") (example)
# server.R fills the placeholder by:
  # output$my_plot <- renderPlot({...})


# app layout:
  # full-screen map (leaflet)
  # floating control panel
  # tabs inside the panel to display different layers of data

# every single thing is commented while I learn to use this, can make a clean copy for github / defense (4.12.26)

# to call specific line number for an error run: source("ui.R")
# to parse syntax errors parse("ui.R")

library(leaflet)


navbarPage("Impacts of marine storms on primary microplastics: A Gulf of Mexico study", id="nav", # 1st arg is app title  # title still undecided 8.18.26
           
           
# Page layout -------------------------------------------------------------
# ── Tab 1: Interactive Map ──────────────────────────────────────────
tabPanel("Interactive map", 
         div(class="outer", # full-screen
             
             tags$head(
               includeCSS("styles.css"), # from SuperZip repo - autozoom to location
               includeScript("gomap.js") # from SuperZip repo
             ),
             
             # map canvas
             leafletOutput("map", width="100%", height="100%"),
             
             # floating control panel (structure & positioning args from SuperZip repo) 
             absolutePanel(id = "controls", class = "panel panel-default", fixed = TRUE,
                           draggable = TRUE, top = 60, left = "auto", right = 20, bottom = "auto",
                           width = 330, height = "auto",
                           
                h4("Nurdle Survey Explorer"), # panel title
                           
                # ── Layer toggles ────────────────────────────────────────────────────
                checkboxInput("show_manufacturers", "Show known nurdle facilities", value = FALSE),
                checkboxInput("show_tracks", "Show storm tracks", value = FALSE),
                           
                # single storm selector drives track layer & site color (accum/dep) & panel plots (always visible)
                selectInput("active_storm", "Select storm:", choices = storm_choices),
                           
                # strandline selector drives popup values & panel plots
                selectInput("strandline", "Strandline", 
                            choices = c(
                            "Total MP (both strandlines)" = "weighted_total_amount",
                            "New strandline" = "weighted_new_strandline",
                            "Old strandline" = "weighted_old_strandline")),
             
                tags$p(tags$em("Hover over storm to see name and max storm category. Use explorer tabs for full plots."),
                       style = "color: gray; font-size: 12px; margin-top: 6px;")
              ), # end absolutePanel
         
              tags$div(id = "cite",
                       # give text a white box to live in
                       style = " background: rgba(255,255,255,0.85); 
                                padding: 6px 10px; border-radius: 4px; 
                                font-size: 11px; color: #333; 
                                max-width: 900px; z-index: 1000;",
              "Citizen-science coastal nurdle surveys (nurdlepatrol.org/en/map),
              Storm track data (NOAA's International Best Track Archive for Climate Stewardship (IBTrACS) data, accessed on [Feb 22, 2025]")
         
  ) # end div.outer
), # end Tab 1


# Tab 2: About Section ----------------------------------------------------
  tabPanel("About",
    fluidRow(
      column(8, offset = 2,
             tags$br(),
             h2("Gulf of Mexico Explorer"),
             tags$p(style = "font-size: 16px;",
                "Nurdles (pre-production plastic pellets) are a ubiquitous form of microplastic (MP) pollution with deleterious environmental impacts, yet are not globally regulated.",
                tags$br(),
                "Under typical conditions, MPs are transported along coastlines via wind, waves, ocean currents, and riverine input.",
                tags$br(),
                "Tropical storms may mobilize and redistribute beached nurdles beyond these standard pathways, yet spatiotemporal patterns of redistribution remain poorly understood.",
                tags$br(),
                "The Gulf of Mexico is a hub for plastics manufacturing and is dominated by the Atlantic Hurricane season.",
                tags$br(),
                "This app presents findings from a master's thesis investigating whether tropical marine storms redistribute beached nurdles along Gulf Coast shorelines."
                ),
             
             h4("Research Questions"),
             tags$ul(
               tags$li("Do tropical storms redistribute previously beached nurdles?"),
               tags$li("Are regional patterns of redistribution consistent across individual storms and sites?"),
               tags$li("Do specific storm characteristics or site attributes drive directional and magnitude changes in post-storm nurdle abundance?")
             ),
             
             h4("Key Findings"),
             tags$ul(
               tags$li("Together, regional, storm, and site-scale analyses provide consistent evidence that marine storms redistribute beached MPs in the Gulf"),
               tags$li("However, redistribution is not characterized by a uniform increase or decrease in abundance."),
               tags$li("Instead, storms produce highly heterogeneous responses, with localized accumulation and depletion occurring simultaneously across spatial and temporal scales, reflecting the interplay of coastal characteristics and storm metrics"), 
               tags$li("The stochastic nature of both storm-driven outcomes and citizen sampling highlight the importance of regulatory framework to classify nurdles as pollutants and mitigate pellet loss at the source of pollution."),  
               ),
             
             h4("Data Sources"),
             tags$ul(
               tags$li(tags$b("Nurdle Patrol"), " — citizen science nurdle survey data (",
                       tags$a("nurdlepatrol.org", href = "https://www.nurdlepatrol.org", target = "_blank"), ")"),
               tags$li(tags$b("IBTrACS"), " — NOAA's NCEI International Best Track Archive for Climate Stewardship storm track data (",
                       tags$a("v4.01, Knapp et al., 2010", href = "https://www.ncei.noaa.gov/products/international-best-track-archive", target = "_blank"), ")"),
               tags$li(tags$b("Plastics Manufacturers"), " — plastic manufacturers known to handle nurdles in the Gulf (",
                       tags$a("Beyond Plastics, 2024", href = "https://www.beyondplastics.org", target = "_blank"), ")")
             ),
             
             h4("App Structure"),
             tags$ul(
               tags$li(tags$b("Interactive Map"), " - Explore survey sites and storm tracks spatially"),
               tags$li(tags$b("Regional Explorer"), " - Gulf-wide nurdle abundance trends and regional redistribution patterns"),
               tags$li(tags$b("Storm Explorer"), " - Site-level data for each paired storm"),
               tags$li(tags$b("Site Explorer"), " - Full survey history and storm response for each site in dataset"),
               tags$li(tags$b("Raw Data"), " - Browse and filter raw survey records for sites paired to at least one tropical storm in the region")
             ),
             
             h4("Key Terms"),
             tags$ul(
               tags$li(tags$b("Survey condition"), " - control (non-storm, baseline conditions) vs storm-influenced (storm within ± 28-days of survey)"),
               tags$li(tags$b("Storm phase"), " - pre-storm vs post-storm survey relative to storm's closest approach to survey site"),
               tags$li(tags$b("Response window"), " - storm-influenced surveys graded by onset and dissipation of storm: acute (± 0-4 days), subacute (± 5-10 days), extended (± 11-28 days)"),
               tags$li(tags$b("Strandline"), " - where on beach survey ocurred: new = most recent tide line, where sand is wet; old = established wrack line; total = reported pellet count, regardless of strandline")
             ),
             
             h4("How to Get Involved"),
             tags$ul(
               tags$li("Participate in Nurdle Patrol data collection along river beds and coastlines — visit ",
                       tags$a("nurdlepatrol.org", href = "https://www.nurdlepatrol.org", target = "_blank"),
                       " to learn the simple sampling protocol and submit data reports."),
               tags$li("Lobby your local politicians to classify nurdles as pollutants and reduce your own plastic use by choosing reusable and sustainable products where possible.")
             ),
             
             tags$br(),
             tags$p(style = "color: gray; font-size: 12px;",
                    "Thesis project by Kelsey Carter- Arizona State University, Biological Data Science.",
                    tags$br(),
                    "Survey data © Nurdle Patrol contributors.",
                    tags$br(),
                    "App layout adapted from the SuperZip Shiny Example (RStudio/Posit, 2023).",
                    tags$br(),
                    "Developed iteratively with assistance from Claude (Anthropic, 2025).",
                    tags$br(),
                    tags$a("Github", href = "https://github.com/KelseyACar", target = "_blank"),
                    " | Publication pending"
                     )
      ) # closes columns
    ) # closes fluid row
  ), # end tab 2



           
  # ── Tab 3: Regional Explorer ──────────────────────────────────────────────
  tabPanel("Regional Explorer",
    sidebarLayout(
      sidebarPanel(width = 3,
        radioButtons("reg_stat", "Statistic:",
          choices = c("Mean" = "mean", "Median" = "median"),
          selected = "median", inline = TRUE
        ),
          
        hr(),
         
        radioButtons("reg_strand", "Strandline:",
          choices  = c("Total" = "total", "New" = "new", "Old" = "old"),
          selected = "total", inline = TRUE
        ),
                                   
        hr(),
                                   
        helpText("Scroll down to see all regional plots.")
                                   
      ), # end sidebarPanel
                      
      mainPanel(width = 9,
                h4("Yearly Abundance"),
                plotOutput("plot_reg_yearly",       height = 340),
                hr(),
                
                h4("Abundance by Storm Phase"),
                plotOutput("plot_regional_condition",     height = 300),
                hr(),
                
                h4("Abundance per Storm: Pre vs. Post"),
                plotOutput("plot_reg_by_storm",     height = 320),
                hr(),
                
                h4("Directional Change Instances per Storm"),
                plotOutput("plot_reg_direction_count", height = 400),
                hr(),
                
                h4("Rate of Change: Control vs. Storm-influenced"),
                plotOutput("plot_reg_rate",         height = 300),
                hr(),
                
                h4("Direction of Change: Control vs. Storm-influenced"),
                plotOutput("plot_reg_direction",    height = 300)
                
      ) # end mainPanel
    ) # end sidebarLayout
  ), # end Tab 3                     
           
           
           
  # ── Tab 4: Storm Explorer ───────────────────────────────────
  tabPanel("Storm Explorer", 
    sidebarLayout(
      sidebarPanel(width = 3,
                                   
        #storm selector            
        selectInput("storm_exp_storm", "Select storm:", choices = storm_choices),
             
                      
        # view & metrics (phase/response window) only relevant when a storm is selected
        conditionalPanel("input.storm_exp_storm != 'All Paired Storms'",
          hr(),
          radioButtons("storm_exp_view", "View:",       # radioButtons allows an item selection from a list
            choices = c("Abundance (pre/post)" = "phase",
                        "Response window" = "response"),
            selected = "phase")
        ), # close conditionalPanel
                      
        # mean or median stat selector
        conditionalPanel("input.storm_exp_view == 'phase'",
          radioButtons("storm_exp_stat", "Statistic",
            choices  = c("Mean" = "mean", "Median" = "median"),
            selected = "median", inline = TRUE)
          )

        ), # end sidebarPanel
                    
        mainPanel(width = 9,
          plotOutput("plot_storm_exp_main", height = 380), # MP abundance
          conditionalPanel("input.storm_exp_storm != 'All Paired Storms'",
            hr(),
            plotOutput("plot_storm_pcs", height = 300), # dir change
            hr(),
            plotOutput("plot_storm_exp_strand", height = 340), # strand line
            hr(),
            uiOutput("storm_metrics_card")) # storm metrics
        ) # end mainPanel
      ) # end sidebarLayout
    ), # end Tab 4
                      
                      
# ── Tab 5: Site Explorer ───────────────────────────────────
tabPanel("Site Explorer",
  sidebarLayout(
    sidebarPanel(width = 3,
                 
      selectInput("site_exp_site", "Select site:",
                  choices = c("Select a site" = "", site_choice_vec)
      ),

      radioButtons("site_exp_stat", "Statistic:",
        choices  = c("Mean" = "mean", "Median" = "median"),
        selected = "median", inline = TRUE
      ),

      selectInput("site_exp_strand", "Strandline:",
        choices = c(
          "Total MP (both strandlines)" = "weighted_total_amount",
          "New strandline"              = "weighted_new_strandline",
          "Old strandline"              = "weighted_old_strandline")
      ),

      uiOutput("site_exp_label")

    ), # end sidebarPanel

    mainPanel(width = 9,
      plotOutput("plot_site_timeseries", height = 380),
      hr(),
      plotOutput("plot_site_pcs_change", height = 320),
      hr(),
      h4("Storm Closest-Approach Metrics"),
      helpText("Meteorological conditions on the day each paired storm was nearest this site."),
      DT::dataTableOutput("closest_metrics_table")
    ) # end mainPanel

  ) # end sidebarLayout
), # end Tab 5
           
          
  # ── Tab 6: Nurdle Survey Data Explorer ───────────────────────────────────
  tabPanel("Nurdle Survey Data Explorer",
           
    fluidRow(
      # countries
      column(3,
        selectInput("country", "Country",
        c("All Countries" = "", "United States", "Mexico"),
        multiple = TRUE)
      ),
                      
      #states - cascade off country
      column(3,
        conditionalPanel("input.country",
          selectInput("states", "States",
            c("All States" = "", structure(
              c("TX", "LA", "MS", "AL", "FL", "Veracruz"),
              names = c("Texas", "Louisiana", "Mississippi", "Alabama", "Florida", "Veracruz")
            )),
            multiple = TRUE)
      )
     ),
                      
     # storms - cascade off states
     column(3,
      conditionalPanel("input.states",
        selectInput("storms_table", "Storm",
          c("All Storms" = "", storm_choices),
          multiple = TRUE)
                                              )
      ),
                      
      # conditions - cascade off storms
      column(3,
        conditionalPanel("input.states",
          selectInput("conditions_table", "Phase",
            c("All Conditons" = "", "control", "pre", "post"),
            multiple=TRUE)
        )
      )
    ), # close fluidRow 1
                    
    fluidRow(
      column(2, numericInput("minCount", "Min nurdles", min=0, max=40000, value=0)),
      column(2, numericInput("maxCount", "Max nurdles", min=0, max=40000, value=10000))
    ), # end fluidRow 2
                    
    hr(),
    DT::dataTableOutput("survey_table")
                    
  ) # end tabPanel 6
           
) # end navbarPage
