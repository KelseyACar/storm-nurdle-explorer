# Script 0: NOAA International Best Track record Archive of Climate Stewardship -------------------------------------------------------------
# Purpose: reads in the latest version of NOAA NCEI's IBTrACS data, then pre-processes, reduces, and cleans it for synthesis with Nurdle Patrol data
# Load Libraries
# Load Data
  # full IBTrACS csv from 1980 on
# Pre-Process 
  # Temporal Filter: reduce df to 2 years before nurdle data starts (2018)
  # Spatial Filter: reduce df to storms in the sub-basin Gulf of Mexico 'GM'
  # Clean df: drop cols w/ no observations, duplicated cols, clean col names, set data types and reorder df

# Load Libraries ---------------------------------------------------------------
install.packages('pacman')
pacman::p_load('lubridate', 'dplyr', 'here')
source(here("functions.R"))

# Call Data ----------------------------------------------------
# read CSV downloaded from https://www.ncei.noaa.gov/products/international-best-track-archive
storms <- read.csv(here('IBTrACS','ibtracs.since1980.list.v04r01.csv'))
#View(storms)


# Pre-Process ------------------------------------------------------------------

# 1) Temporal Filter: 2016, 2 years before nurdle df
storms <- storms %>% filter(ISO_TIME >= '2016-11-01 00:00:00')

# 2) Spatial Filter: Gulf of Mexico subbasin
storms <- storms %>% filter(SUBBASIN == 'GM')

# 3) df cleaning

#convert to df
storms <- as.data.frame(storms)
class(storms) # check conversion 

# drop all blank and all NA columns
result <- sapply(storms, function(x) all(is.na(x) | x == " ")) # ID blank and NA cols
storms <- storms[,!result] # drops results (blank cols)
# head(storms)

# drop duplicate cols w/ diff reporting structure
storms <- storms[, !names(storms) %in% c('USA_LAT', 'USA_LON','WMO_WIND', "WMO_PRES")]

# rename date, lat and lon columns to match nurdle df
storms <- storms %>% rename(date = ISO_TIME, Latitude = LAT, Longitude = LON)

# convert date to POSIXct
storms$date <- as.POSIXct(storms$date, format = "%Y-%m-%d %H:%M:%S")

# reorder cols with leading col names
leading_cols <- c("SID", "NAME", "date", "Longitude", "Latitude", "NATURE", 
                  "USA_STATUS", "USA_WIND", "USA_PRES", "STORM_SPEED", 
                  "STORM_DIR", "DIST2LAND", "USA_SEAHGT")

# remaining columns
remaining_cols <- setdiff(names(storms), leading_cols)

# reorder df
storms <- storms[, c(leading_cols, remaining_cols)]
#head(storms)


# clean column names with janitor
storms <- storms %>% clean_names()


# Detect the full storm_start and storm_end dates so full extend of storm is captured in linked df even if MP surveys don't occur across full storm (to look at full effect of storm later)
storms <- storms %>%
  group_by(sid) %>% # use sid because some storm names repeat
  mutate(
    full_storm_start = min(date),
    full_storm_end = max(date)
  )
#View(storms)

# write to csv
write.csv(storms, file = here('output', 'GOM_storms.csv'), row.names = FALSE)


