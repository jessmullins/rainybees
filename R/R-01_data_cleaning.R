## =============================================================================
## R-01_data_cleaning.R
##
## Purpose : Clean and prepare raw GBIF bee occurrence data for all downstream
##           analyses. Joins site metadata, harmonizes species names, joins
##           species traits, calculates sampling effort, builds the community
##           abundance matrix, and runs iNEXT coverage-based rarefaction.
##
##           The iNEXT coverage target is saved as target_coverage.rds and
##           loaded by downstream analysis scripts to ensure identical
##           standardization across all analyses.
##
## Inputs  : data/raw/occurrences.csv    — GBIF occurrence records
##           data/raw/siteInfo.csv       — site-level metadata
##           data/raw/traits.csv         — species traits (ITD, diet, sociality,
##                                         range size)
##
## Outputs : data/cleaned/df_analysis.csv       — filtered, trait-joined records
##           data/cleaned/comm_abun.csv          — site_year × species matrix
##           data/cleaned/comm_rare.csv          — coverage-rarefied matrix
##           data/cleaned/meta.csv               — one row per site_year + effort
##           data/cleaned/rarefied_n.csv         — per-assemblage rarefied n
##           data/cleaned/target_coverage.rds    — shared coverage target scalar
##
## Author  : Jessica Mullins
## Date    : April 2026
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(iNEXT)
})

## =============================================================================
## SECTION 0 — File paths + settings
## =============================================================================

raw_dir     <- "data/raw"
cleaned_dir <- "data/cleaned"
figures_dir <- "figures"
results_dir <- "results"
dir.create(cleaned_dir, recursive = TRUE, showWarnings = FALSE)

gbif_file   <- file.path(raw_dir, "occurrences.csv")
sites_file  <- file.path(raw_dir, "siteInfo.csv")
traits_file <- file.path(raw_dir, "traits.csv")

set.seed(816)   # for reproducible rarefaction

## =============================================================================
## SECTION 1 — Canonical species name function
##
## Converts raw species strings into a clean "Genus_species" key used as
## column names in the community matrix.
##
## Rules applied in order:
##   1. Trim whitespace
##   2. Strip "_SIM" suffix (tree simulation artefact)
##   3. Drop trailing periods
##   4. Replace non-alphanumeric/non-underscore chars with "_"
##   5. Split on "_"; remove empty tokens
##   6. Drop open-nomenclature qualifiers: "cf", "aff", "nr"
##   7. Return "Genus_species" — single tokens returned as-is (genus only)
## =============================================================================

canon_name <- function(x) {
  x <- trimws(x)
  x <- gsub("_SIM$",         "", x, ignore.case = TRUE)
  x <- gsub("\\.+$",         "", x)
  x <- gsub("[^A-Za-z0-9_]", "_", x)

  parts <- strsplit(x, "_")

  out <- vapply(parts, function(tokens) {
    tokens <- tokens[tokens != ""]
    if (length(tokens) == 0) return(NA_character_)
    keep   <- !tolower(tokens) %in% c("cf", "aff", "nr")
    tokens <- tokens[keep]
    if (length(tokens) == 0) return(NA_character_)
    if (length(tokens) >= 2) {
      paste(tokens[1], tolower(tokens[2]), sep = "_")
    } else {
      tokens[1]
    }
  }, character(1))

  out
}

## =============================================================================
## SECTION 2 — Load raw data
## =============================================================================

df_gbif  <- read.csv(gbif_file)
df_sites <- read.csv(sites_file)

## =============================================================================
## SECTION 3 — Join site metadata + harmonize legacy site aliases
##
## Three sites were renamed between survey years. Collapse each pair to the
## canonical current site code so all years for the same physical location
## share a single fieldNumber.
## =============================================================================

df_gbif <- df_gbif %>%
  left_join(df_sites %>% select(fieldNumber, Type), by = "fieldNumber")

df_gbif <- df_gbif %>%
  mutate(fieldNumber = case_when(
    fieldNumber == "ECR1" ~ "ECR4",   # legacy alias -> current code
    fieldNumber == "ECR2" ~ "ECR5",
    fieldNumber == "MTS1" ~ "MTS1A",
    TRUE ~ fieldNumber
  ))

## =============================================================================
## SECTION 4 — Filter lists
## =============================================================================

# Sites excluded: irregular sampling, urban matrix outliers, short time series
l_noplots <- c("SWS10", "CrCyn", "MTLA", "SWI1", "MTS1A-Urb",
               "UCSDN", "TRR1",  "TRS1", "SCR")

# Survey years included in the analysis
l_years   <- c("2011", "2012", "2015", "2016", "2022", "2023", "2024")

# Pan-trap records only — exclude malaise, net, etc.
l_sample  <- c("Fluorescent blue bowl", "White bowl trap",
               "Fluorescent yellow bowl", "Bowl trap")

# Habitat types retained
l_type    <- c("Reserve", "Fragment")

## =============================================================================
## SECTION 5 — Filter occurrences
##
## Removes:
##   - non-pan-trap records
##   - out-of-scope survey years
##   - excluded / irregular plots
##   - honeybee (Apis mellifera)
##   - records with missing site or habitat-type info
## =============================================================================

df_gbif2 <- df_gbif %>%
  mutate(
    year      = factor(year, levels = l_years),
    Count     = as.numeric(individualCount),
    site_year = paste(fieldNumber, year, sep = "_")
  ) %>%
  filter(
    samplingProtocol %in% l_sample,
    year             %in% l_years,
    !fieldNumber     %in% l_noplots,
    scientificName   != "Apis mellifera",
    !is.na(fieldNumber),
    Type             %in% l_type
  )

## =============================================================================
## SECTION 6 — Harmonize species names
##
## Map synonyms and subspecific names in GBIF to currently accepted names
## used in the trait file, so joins don't lose data due to name mismatches.
## =============================================================================

df_gbif2 <- df_gbif2 %>%
  mutate(scientificName = case_when(
    scientificName == "Agapostemon subtilior"         ~ "Agapostemon texanus",
    scientificName == "Dianthidium pudicum consimile" ~ "Dianthidium pudicum",
    scientificName == "Neopasites cressoni"           ~ "Biastes cressoni",
    scientificName == "Neopasites mojavensis"         ~ "Biastes mojavensis",
    scientificName == "Eucera frater albopilosa"      ~ "Eucera frater",
    scientificName == "Hoplitis albifrons maura"      ~ "Hoplitis albifrons",
    scientificName == "Hoplitis fulgida platyura"     ~ "Hoplitis fulgida",
    scientificName == "Hylaeus episcopalis metzi"     ~ "Hylaeus episcopalis",
    scientificName == "Hylaeus mesillae cressoni"     ~ "Hylaeus mesillae",
    scientificName == "Melissodes communis alopex"    ~ "Melissodes communis",
    scientificName == "Osmia montana quadriceps"      ~ "Osmia montana",
    scientificName == "Perdita claypolei australior"  ~ "Perdita claypolei",
    scientificName == "Perdita interrupta interrupta" ~ "Perdita interrupta",
    TRUE ~ scientificName
  ))

## =============================================================================
## SECTION 7 — Join species traits
##
## Joins trait columns from traits.csv by scientificName:
##   itd_species_mean — intertegular distance (species mean; body size proxy)
##   itd              — intertegular distance (individual measurement)
##   diet           — dietary specialisation (Oligolectic / Polylectic)
##   excavator        — nest excavation behaviour
##   nestloc          — nest location (Ground, Cavity, etc.)
##   soc              — sociality (Solitary, Social, etc.)
##   range_km2        — total range area (km², including Canada/Mexico from Chesshire et al. 2023 Ecography)
##
## Species missing traits are retained — NAs are handled downstream.
## =============================================================================

df_traits <- read.csv(traits_file)

df_gbif2 <- df_gbif2 %>%
  left_join(
    df_traits %>%
      select(scientificName, itd, diet, excavator, nestloc, soc, range_km2),
    by = "scientificName"
  )

missing_range <- df_gbif2 %>%
  filter(is.na(range_km2), grepl(" ", scientificName)) %>%
  distinct(scientificName) %>%
  arrange(scientificName)

if (nrow(missing_range) > 0) {
  message("\n--- Species missing trait data ---")
  print(missing_range)
} else {
  message("All species matched to trait data successfully.")
}

## =============================================================================
## SECTION 8 — Sampling effort: trap-days per site_year
##
## trap_days = sampling_rounds × traps
## Trap counts differ by era:
##   2011–2012: 30 traps per round
##   2015+    : 15 traps per round
## =============================================================================

df_sampling_rounds <- tibble(
  site_year = c(
    "CFS1_2012",  "ECR4_2012",   "ECR5_2012",  "ECR4_2015",  "ECR4_2016",
    "ECR4_2022",  "ECR4_2023",   "ECR4_2024",  "ECR5_2015",  "ECR5_2016",
    "ECR5_2022",  "ECR5_2023",   "ECR5_2024",  "MTE2_2011",  "MTE2_2012",
    "MTE3_2015",  "MTE3_2016",   "MTE3_2022",  "MTE3_2023",  "MTE3_2024",
    "MTE4_2015",  "MTE4_2016",   "MTE4_2022",  "MTE4_2023",  "MTE4_2024",
    "MTI2_2011",  "MTI2_2012",   "MTI2_2022",  "MTI2_2023",  "MTI2_2024",
    "MTLB1_2011", "MTS1A_2011",  "MTS1A_2012", "MTS1A_2015", "MTS1A_2016",
    "MTS1A_2022", "MTS1A_2023",  "MTS1A_2024", "MTS2_2012",  "MTS2_2015",
    "MTS2_2016",  "MTS2_2022",   "MTS2_2023",  "MTS2_2024",  "MTS3_2012",
    "MTS6_2012",  "MTS6_2015",   "MTS6_2016",  "MTS6_2022",  "MTS6_2023",
    "MTS6_2024",  "MTS7_2012",   "MTS7_2022",  "MTS7_2023",  "MTS7_2024",
    "SWEA_2011",  "SWEA_2012",   "SWEA_2015",  "SWEA_2016",  "SWEA_2022",
    "SWEA_2023",  "SWEA_2024",   "SWI2_2011",  "SWI2_2012",  "SWI4_2015",
    "SWI4_2016",  "SWI4_2022",   "SWI4_2023",  "SWI4_2024",  "SWS1_2011",
    "SWS1_2012",  "SWS1_2015",   "SWS1_2016",  "SWS1A_2022", "SWS1A_2023",
    "SWS1A_2024", "SWS3_2011",   "SWS3_2012",  "SWS3_2015",  "SWS3_2016",
    "SWS3_2022",  "SWS3_2023",   "SWS3_2024",  "TRL1_2012"
  ),
  sampling_rounds = c(
     5,  4,  4,  7,  6,  8,  6,  6,  7,  6,  8,  6,  6,  9,  6,
     7,  6,  8,  6,  6,  7,  6,  8,  6,  6,  9,  6,  8,  6,  6,
     9,  9,  5,  7,  6,  8,  5,  6,  5,  7,  6,  8,  6,  6,  5,
     5,  7,  6,  8,  6,  6,  5,  8,  6,  6, 11,  6,  7,  7,  8,
     6,  6, 10,  6,  7,  7,  8,  6,  6, 12,  5,  7,  6,  8,  5,
     6, 12,  5,  7,  6,  8,  6,  6,  5   # TRL1_2012: 5 rounds
  )
)

df_traps <- df_sampling_rounds %>%
  mutate(
    year  = as.integer(str_extract(site_year, "\\d{4}")),
    traps = if_else(year %in% c(2011L, 2012L), 30L, 15L)
  ) %>%
  select(site_year, traps)

df_gbif2 <- df_gbif2 %>%
  left_join(df_sampling_rounds, by = "site_year") %>%
  left_join(df_traps,           by = "site_year") %>%
  mutate(trap_days = sampling_rounds * traps)

## =============================================================================
## SECTION 9 — Trim to analysis columns + add canonical name
## =============================================================================

df_analysis <- df_gbif2 %>%
  select(
    catalogNumber, family, genus, specificEpithet, scientificName,
    year, eventDate, fieldNumber, sex, individualCount, samplingProtocol,
    decimalLatitude, decimalLongitude, Type, site_year, itd, diet,
    excavator, nestloc, soc, range_km2,
    sampling_rounds, traps, trap_days
  ) %>%
  mutate(
    Count  = as.numeric(individualCount),
    sp_can = canon_name(gsub(" ", "_", scientificName))
  )

## =============================================================================
## SECTION 10 — Build community abundance matrix (site_year × species)
## =============================================================================

comm_abun <- df_analysis %>%
  filter(!is.na(sp_can), !is.na(site_year), !is.na(Count)) %>%
  group_by(site_year, sp_can) %>%
  summarise(abundance = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from  = sp_can,
    values_from = abundance,
    values_fill = 0
  ) %>%
  as.data.frame() %>%
  tibble::column_to_rownames("site_year")

cat("Community matrix dimensions:",
    nrow(comm_abun), "assemblages x", ncol(comm_abun), "species\n")

## =============================================================================
## SECTION 11 — iNEXT: coverage-based rarefaction
##
## Step 1: Build abundance list from community matrix
## Step 2: Run iNEXT to get observed coverage per assemblage
## Step 3: Set target_coverage = minimum observed coverage across all sites
## Step 4: For each assemblage, find the sample size (n) that achieves
##         target_coverage via interpolation from iNEXT size-based curves
## Step 5: Rarefy each assemblage to its coverage-matched n (999 iterations)
##
## target_coverage.rds is loaded by downstream scripts (R-05 through R-08)
## to ensure identical standardization across all analyses.
##
## NOTE: iNEXT can take several minutes on large datasets.
## =============================================================================

cat("\nBuilding iNEXT abundance list...\n")

comm_mat <- as.matrix(comm_abun)

abundance_list <- lapply(1:nrow(comm_mat), function(i) {
  x <- comm_mat[i, ]
  as.numeric(x[x > 0])
})
names(abundance_list) <- rownames(comm_mat)

cat("Running iNEXT on", length(abundance_list), "assemblages...\n")
cat("(This may take several minutes)\n")

inext_out <- iNEXT(abundance_list, q = 0, datatype = "abundance")

# Extract observed coverage per assemblage
observed_coverage <- inext_out$DataInfo %>%
  select(Assemblage, SC) %>%
  rename(site_year = Assemblage, coverage = SC)

target_coverage <- min(observed_coverage$coverage)

cat("\n--- Coverage summary across assemblages ---\n")
print(summary(observed_coverage$coverage))
cat("Target coverage (minimum observed):", round(target_coverage, 4), "\n")
cat("Driven by site:",
    observed_coverage$site_year[which.min(observed_coverage$coverage)], "\n\n")

saveRDS(target_coverage, file.path(cleaned_dir, "target_coverage.rds"))

## =============================================================================
## SECTION 12 — Per-assemblage rarefied sample sizes
##
## For each assemblage find the sample size that achieves target_coverage
## via interpolation from iNEXT size-based curves.
## =============================================================================

rarefied_n <- map(names(abundance_list), function(sy) {
  est <- inext_out$iNextEst$size_based %>%
    filter(Assemblage == sy, Method == "Rarefaction")
  closest <- est %>%
    slice_min(abs(SC - target_coverage), n = 1)
  tibble(
    site_year         = sy,
    rarefied_n        = as.integer(round(closest$m)),
    achieved_coverage = closest$SC
  )
}) %>%
  list_rbind()

cat("Rarefied sample size summary:\n")
print(summary(rarefied_n$rarefied_n))

write.csv(rarefied_n,
          file.path(cleaned_dir, "rarefied_n.csv"),
          row.names = FALSE)

## =============================================================================
## SECTION 13 — Coverage-rarefied community matrix
##
## Rarefies each assemblage to its coverage-matched n using 999 iterations
## and averages the results.
## =============================================================================

rarefy_community <- function(abund_vec, n_target, n_iter = 999) {
  species <- names(abund_vec)
  pool    <- rep(species, abund_vec)
  results <- replicate(n_iter, {
    samp <- sample(pool, size = min(n_target, length(pool)), replace = FALSE)
    table(factor(samp, levels = species))
  }, simplify = FALSE)
  avg_abund <- Reduce("+", results) / n_iter
  avg_abund[avg_abund > 0]
}

cat("\nRarefying community matrix (999 iterations per assemblage)...\n")

rarefied_communities <- map2(
  abundance_list,
  rarefied_n$rarefied_n[match(names(abundance_list), rarefied_n$site_year)],
  ~ rarefy_community(.x, .y, 999)
)

all_species <- unique(unlist(lapply(rarefied_communities, names)))

comm_rare <- do.call(rbind, lapply(
  names(rarefied_communities), function(sy) {
    x   <- rarefied_communities[[sy]]
    vec <- setNames(rep(0, length(all_species)), all_species)
    vec[names(x)] <- as.numeric(x)
    vec
  }
))
rownames(comm_rare) <- names(rarefied_communities)

cat("Rarefied matrix dimensions:",
    nrow(comm_rare), "assemblages x", ncol(comm_rare), "species\n")

## =============================================================================
## SECTION 14 — Build site metadata table
## =============================================================================

meta <- df_analysis %>%
  group_by(site_year) %>%
  summarise(
    fieldNumber = first(fieldNumber),
    year        = first(year),
    Type        = first(Type),
    trap_days   = suppressWarnings(first(trap_days)),
    N_raw       = sum(Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  left_join(observed_coverage, by = "site_year") %>%
  left_join(rarefied_n,        by = "site_year")

## =============================================================================
## SECTION 15 — Save all outputs
## =============================================================================

write.csv(df_analysis,
          file.path(cleaned_dir, "df_analysis.csv"),   row.names = FALSE)
write.csv(comm_abun,
          file.path(cleaned_dir, "comm_abun.csv"),     row.names = TRUE)
write.csv(comm_rare,
          file.path(cleaned_dir, "comm_rare.csv"),     row.names = TRUE)
write.csv(meta,
          file.path(cleaned_dir, "meta.csv"),          row.names = FALSE)

# target_coverage.rds already saved in Section 11

message("\n--- Cleaning complete ---")
message("Assemblages (site_year):     ", nrow(comm_abun))
message("Species (full matrix):       ", ncol(comm_abun))
message("Species (rarefied matrix):   ", ncol(comm_rare))
message("Target coverage:             ", round(target_coverage, 4))
message("Total records (df_analysis): ", nrow(df_analysis))
message("Files written to:            ", cleaned_dir)
