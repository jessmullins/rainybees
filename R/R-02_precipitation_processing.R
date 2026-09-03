# =============================================================================
# R-02_precipitation_processing.R
#
# Purpose : Build annual winter precipitation totals and lagged predictors
#           from the daily Montgomery Field record, and map study plots and
#           weather stations.
#
#           Winter year Y is defined as 1 November of year Y-1 through
#           31 March of year Y. November and December are therefore assigned
#           to the FOLLOWING calendar year; January, February and March keep
#           their own. WPy0 ("same") is the winter total for the survey year;
#           WPy0-1 ("previous1") is the preceding winter, and so on to
#           previous6.
#
# Inputs  : data/raw/precipitation_manual.csv  — daily PRCP, Montgomery Field
#                                                (NOAA CDO download, inches)
#           data/raw/site_stations.csv         — plot-to-station assignments
#
# Outputs : data/cleaned/precipyrs.csv         — annual winter totals + lags
#           figures/stations_map.png           — Fig. S2
#
# API key : Section 7 requires a free Stadia Maps key: https://stadiamaps.com/
#           Do not commit your key to GitHub.
#
# Author  : Jessica Mullins
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(tibble)
library(ggplot2)
library(ggmap)
library(ggrepel)

raw_dir     <- "data/raw"
cleaned_dir <- "data/cleaned"
figures_dir <- "figures"

MONTGOMERY <- "SAN DIEGO MONTGOMERY FIELD, CA US"

# NOAA CDO serves daily precipitation in inches; all analyses use mm.
IN_TO_MM <- 25.4

# =============================================================================
# SECTION 1 — Daily record to winter-year totals
# =============================================================================

daily <- read.csv(file.path(raw_dir, "precipitation_manual.csv"),
                  stringsAsFactors = FALSE) %>%
  mutate(DATE = as.Date(DATE),
         PRCP = as.numeric(PRCP) * IN_TO_MM)

winter_totals <- daily %>%
  filter(month(DATE) %in% c(11, 12, 1, 2, 3)) %>%
  mutate(year = if_else(month(DATE) %in% c(11, 12),
                        year(DATE) + 1L,
                        year(DATE))) %>%
  group_by(year) %>%
  summarise(same = sum(PRCP, na.rm = TRUE), .groups = "drop") %>%
  arrange(year)

# =============================================================================
# SECTION 2 — Lagged winter precipitation (WPy0-1 through WPy0-6)
# =============================================================================

precipyrs <- winter_totals %>%
  mutate(
    previous1 = lag(same, 1),
    previous2 = lag(same, 2),
    previous3 = lag(same, 3),
    previous4 = lag(same, 4),
    previous5 = lag(same, 5),
    previous6 = lag(same, 6)
  )

write.csv(precipyrs,
          file.path(cleaned_dir, "precipyrs.csv"),
          row.names = FALSE)

message("Wrote ", file.path(cleaned_dir, "precipyrs.csv"),
        " (", nrow(precipyrs), " winter years: ",
        min(precipyrs$year), "-", max(precipyrs$year), ")")

# =============================================================================
# SECTION 7 — Map of study plots and weather stations (Fig. S2)
#
# Requires a free Stadia Maps API key: https://stadiamaps.com/
# Do not commit your key to GitHub.
# =============================================================================

register_stadiamaps(key = "YOUR TOKEN HERE")

site_stations <- read.csv(file.path(raw_dir, "site_stations.csv"))

# --- Combine legacy-alias site labels for the map -------------------------
# SWS1/SWS1A share identical coordinates (same physical site, never
# collapsed in R-01), and ECR1->ECR4 / ECR2->ECR5 are legacy site-code
# renames (R-01 Section 3). Whichever code(s) actually appear in
# site_stations.csv, relabel to a combined "old/new" form and de-duplicate
# so only ONE point/label renders per physical location, regardless of
# which specific code(s) this file happens to use.
combine_alias_label <- function(label) {
  dplyr::case_when(
    label %in% c("SWS1", "SWS1A")     ~ "SWS1/SWS1A",
    label %in% c("ECR1", "ECR4")      ~ "ECR1/4",
    label %in% c("ECR2", "ECR5")      ~ "ECR2/5",
    TRUE                              ~ label
  )
}

plot_points <- site_stations %>%
  distinct(fieldNumber, Type, Latitude, Longitude) %>%
  mutate(PointType = Type,
         label      = combine_alias_label(fieldNumber)) %>%
  distinct(label, Latitude, Longitude, PointType, .keep_all = TRUE)  # collapse duplicate-coordinate aliases to one point

station_coords <- tibble::tibble(
  NAME      = c(
    MONTGOMERY,
    "CHULA VISTA, CA US",
    "SAN DIEGO INTERNATIONAL AIRPORT, CA US",
    "SAN DIEGO SEAWORLD, CA US",
    "EL CAJON, CA US",
    "SAN DIEGO MIRAMAR NAS, CA US"
  ),
  Latitude  = c(32.81444, 32.6400, 32.73361, 32.7631, 32.7384, 32.8711),
  Longitude = c(-117.13639, -117.0800, -117.18306, -117.2299, -116.91, -117.1719)
)

# mutate() rather than transmute() so a label column survives for the
# Montgomery Field point; the other stations stay unlabelled.
weather_points <- station_coords %>%
  mutate(
    PointType = if_else(NAME == MONTGOMERY, "Montgomery Field", "Weather Station"),
    label     = if_else(NAME == MONTGOMERY, "Montgomery Field", NA_character_)
  ) %>%
  select(Latitude, Longitude, PointType, label)

map_points <- bind_rows(
  plot_points %>% select(Latitude, Longitude, PointType, label),
  weather_points
)

buffer <- 0.02
bbox <- c(
  left   = min(map_points$Longitude, na.rm = TRUE) - buffer,
  bottom = min(map_points$Latitude,  na.rm = TRUE) - buffer,
  right  = max(map_points$Longitude, na.rm = TRUE) + buffer,
  top    = max(map_points$Latitude,  na.rm = TRUE) + buffer
)

sd_map <- get_stadiamap(bbox = bbox, zoom = 11,
                        maptype = "stamen_terrain_background")

mappy <- ggmap(sd_map, darken = c(0.4, "white")) +
  geom_point(
    data = map_points,
    aes(x = Longitude, y = Latitude, shape = PointType, color = PointType),
    size = 3, stroke = 1
  ) +
  geom_text_repel(
    data = filter(map_points, PointType %in% c("Fragment", "Reserve", "Montgomery Field")),
    aes(x = Longitude, y = Latitude, label = label),
    size = 3, box.padding = 0.3, point.padding = 0.2,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_shape_manual(values = c("Fragment"         = 16, "Reserve"          = 17,
                                "Weather Station"  = 17, "Montgomery Field" = 17)) +
  scale_color_manual(values = c("Fragment"         = "gray52",
                                "Reserve"          = "springgreen4",
                                "Weather Station"  = "blue",
                                "Montgomery Field" = "orange")) +
  labs(title    = "Sampling plots and weather stations",
       subtitle = "Coastal San Diego County",
       x = "Longitude", y = "Latitude",
       color = "Type", shape = "Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(figures_dir, "stations_map.png"),
       mappy, width = 6, height = 6, dpi = 300)