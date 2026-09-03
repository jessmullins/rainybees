#### Spatial autocorrelation test: does site location predict unexplained
#### variance in bee richness models, after accounting for Type x precipitation?
####
#### Reviewer concern: reserve sites are spatially clustered relative to
#### fragment sites; are site-type estimates confounded with spatial location?
####
#### Approach:
#### 1) Fit the richness model as before (Type * WPy0-1 + Type * WPy0),
####    now WITH year as a random effect to match the primary model used
####    throughout the rest of the paper (WPy0/WPy0-1 are year-level
####    variables, so (1|year) corrects for pseudoreplication across sites
####    sampled in the same year -- see main Methods).
#### 2) Extract site-level random-effect intercepts (BLUPs) -- one value per
####    Plot, which is the right unit for a spatial test (site-years at the
####    same site are not independent spatial replicates)
#### 3) Test Moran's I on those BLUPs using great-circle distance between sites
#### 4) Visualize with a spline correlogram (distance vs. residual correlation)
#### 5) Refit as a GLS with an explicit spatial correlation structure and
####    compare AIC / Type effect to the original lmer, as a robustness check
####    (GLS uses one row per SITE, averaged across years, so no separate
####    year-random-effect correction is needed for this specific check)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(car)
  library(ape)        # Moran.I
  library(ncf)         # spline.correlog
  library(nlme)        # gls with spatial correlation structures
  library(geosphere)   # great-circle (Haversine) distances
})

figures_dir <- "figures"
results_dir <- "results"
cleaned_dir <- "data/cleaned"
raw_dir     <- "data/raw"

#### 0) Load site coordinates -----------------------------------------------
# siteInfo.csv already uses canonical fieldNumber codes for most sites (post
# R-01 Section 3 alias collapsing: ECR1->ECR4, ECR2->ECR5, MTS1->MTS1A), and
# critically already contains BOTH "SWS1" and "SWS1A" as separate rows with
# identical coordinates (32.7496, -117.0316) -- these were never collapsed
# in R-01, so df_analysis$fieldNumber (-> Plot, -> the model's random-effect
# levels) still treats them as two site codes. We match that here: join
# directly on Plot == fieldNumber with no collapsing, so each retains its
# own BLUP. This does introduce one zero-distance pair (SWS1 <-> SWS1A) into
# the distance matrix, which the inverse-distance weighting below already
# guards against (w[!is.finite(w)] <- 0).
site_meta <- read.csv(file.path(raw_dir, "siteInfo.csv")) %>%
  rename(Plot = fieldNumber) %>%
  distinct(Plot, Latitude, Longitude, .keep_all = TRUE)

#### 1) Refit richness model (same specification as main analysis, +year) ---

df_analysis <- read.csv(file.path(cleaned_dir, "df_analysis.csv"))
target_cov  <- readRDS(file.path(cleaned_dir, "target_coverage.rds"))
precipitation_yearly <- read.csv(file.path(cleaned_dir, "precipyrs.csv")) %>% drop_na()
precipitation_yearly$year <- as.factor(precipitation_yearly$year)

mean_same  <- mean(precipitation_yearly$same,      na.rm = TRUE)
mean_prev  <- mean(precipitation_yearly$previous1, na.rm = TRUE)

allbees <- df_analysis %>%
  rename(Plot = fieldNumber) %>%
  mutate(year = as.factor(year), Count = as.numeric(Count),
         plot_year = paste(Plot, year, sep = "_"))

comm <- allbees %>%
  group_by(plot_year, Plot, year, Type, scientificName) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  pivot_wider(names_from = scientificName, values_from = n, values_fill = 0)

meta <- comm %>% dplyr::select(plot_year, Plot, year, Type)
species_cols <- setdiff(names(comm), c("plot_year","Plot","year","Type"))

abund_list <- apply(comm[, species_cols, drop = FALSE], 1, function(z) {
  z <- as.numeric(z); z[z > 0]
})
names(abund_list) <- comm$plot_year

est <- iNEXT::estimateD(abund_list, q = c(0,1,2), datatype = "abundance",
                        base = "coverage", level = target_cov) %>% as_tibble()

long_cov <- est %>%
  transmute(plot_year = Assemblage, Order.q, SC, qD) %>%
  left_join(meta, by = "plot_year") %>%
  mutate(Metric = dplyr::recode(Order.q, `0`="Richness", `1`="Hill_Shannon", `2`="Hill_Simpson")) %>%
  dplyr::select(Plot, Type, year, plot_year, Metric, Value = qD, common_coverage = SC) %>%
  left_join(precipitation_yearly %>% dplyr::select(year, same, previous1), by = "year") %>%
  mutate(same_c = same - mean_same, previous1_c = previous1 - mean_prev)

dd_rich <- long_cov %>%
  filter(Metric == "Richness") %>%
  drop_na(Value, Plot, Type, same, previous1, same_c, previous1_c) %>%
  mutate(Type = forcats::fct_relevel(Type, "Reserve", "Fragment"))

# Primary model, now WITH year as a random effect (matches R-05's
# m_rich_full_year_reml specification)
m_rich_full <- lmer(Value ~ Type * previous1_c + Type * same_c + (1 | Plot) + (1 | year),
                    data = dd_rich, REML = TRUE)
cat("Singular fit?", lme4::isSingular(m_rich_full), "\n")

print(car::Anova(m_rich_full, type = "III"))

#### 2) Extract site-level random-effect intercepts (BLUPs) -----------------

blups <- ranef(m_rich_full)$Plot %>%
  as.data.frame() %>%
  rownames_to_column("Plot") %>%
  rename(Intercept_BLUP = `(Intercept)`)

blups_geo <- blups %>%
  left_join(site_meta, by = "Plot") %>%
  drop_na(Latitude, Longitude)

cat("\nSites with BLUPs matched to coordinates:", nrow(blups_geo),
    "/", nrow(blups), "\n")
if (nrow(blups_geo) < nrow(blups)) {
  cat("WARNING: some Plot names did not match site_meta$Plot. Check:\n")
  print(setdiff(blups$Plot, site_meta$Plot))
}

#### 3) Moran's I on site-level BLUPs, using great-circle distance ----------

# Pairwise great-circle distance (meters) -> inverse-distance weights
coords <- blups_geo %>% dplyr::select(Longitude, Latitude) %>% as.matrix()
dist_m <- geosphere::distm(coords, fun = distHaversine)

# Inverse-distance weight matrix; zero diagonal; guard against div-by-zero
w <- 1 / dist_m
diag(w) <- 0
w[!is.finite(w)] <- 0

moran_result <- ape::Moran.I(blups_geo$Intercept_BLUP, w)
cat("\n--- Moran's I on site-level richness-model random intercepts ---\n")
print(moran_result)

#### 3b) Moran's I on Type itself: are reserves spatially clustered? --------
# This is the test that speaks *directly* to the reviewer's concern -- not
# "is there leftover unexplained spatial structure after modeling," but
# "is habitat type (fragment vs. reserve) itself spatially clustered?"
# Uses one row per site (not site-year), coded Reserve = 1, Fragment = 0,
# with the same site set / coordinates / distance weighting as above so the
# two Moran's I results are directly comparable.

type_geo <- dd_rich %>%
  distinct(Plot) %>%
  mutate(Plot = if_else(Plot == "SWS1", "SWS1A", Plot)) %>%  # same physical site
  distinct(Plot) %>%
  left_join(site_meta, by = "Plot") %>%
  drop_na(Latitude, Longitude, Type)

type_geo <- dd_rich %>%
  distinct(Plot) %>%
  mutate(Plot = if_else(Plot == "SWS1", "SWS1A", Plot)) %>%
  distinct(Plot) %>%
  left_join(site_meta, by = "Plot") %>%
  drop_na(Latitude, Longitude, Type) %>%
  mutate(Type_num = if_else(Type == "Reserve", 1, 0))

cat("\nSites with Type matched to coordinates:", nrow(type_geo), "\n")

coords_type <- type_geo %>% dplyr::select(Longitude, Latitude) %>% as.matrix()
dist_m_type <- geosphere::distm(coords_type, fun = distHaversine)

w_type <- 1 / dist_m_type
diag(w_type) <- 0
w_type[!is.finite(w_type)] <- 0

moran_type_result <- ape::Moran.I(type_geo$Type_num, w_type)
cat("\n--- Moran's I on habitat Type (Reserve = 1, Fragment = 0) ---\n")
cat("(Tests whether reserve/fragment status is itself spatially clustered,\n")
cat(" independent of the diversity models)\n")
print(moran_type_result)

#### 4) Spline correlogram: residual correlation vs. distance ---------------
# Useful supplementary figure: shows whether/at what distance spatial
# correlation decays. lat/lon in decimal degrees works fine at this county
# scale for ncf::spline.correlog.

set.seed(123)
correlog_result <- ncf::spline.correlog(
  x = blups_geo$Longitude,
  y = blups_geo$Latitude,
  z = blups_geo$Intercept_BLUP,
  resamp = 999,
  latlon = TRUE  # treats x,y as lon,lat and computes great-circle distance
)

pdf(file.path(figures_dir, "FigureSX_spatial_correlogram.pdf"), width = 6, height = 5)
plot(correlog_result, main = "Spline correlogram: richness model site intercepts")
dev.off()
cat("\nSaved spline correlogram to", file.path(figures_dir, "FigureSX_spatial_correlogram.pdf"), "\n")

#### 5) Robustness check: refit with explicit spatial correlation structure -
# Compare a GLS with no spatial structure to one with an exponential spatial
# correlation on site coordinates. Uses SITE-level data (one row per site,
# averaged across years) -- because each site contributes exactly one row
# here, there is no year-level pseudoreplication in this specific check, so
# no (1|year) analogue is needed for the GLS comparison.

dd_rich_geo <- dd_rich %>%
  dplyr::select(-any_of("Type")) %>%   # avoid ambiguity; Type comes from site_meta below
  left_join(site_meta, by = "Plot") %>%
  drop_na(Latitude, Longitude, Type)

# corExp requires UNIQUE coordinates across rows -- but dd_rich_geo has one
# row per site-YEAR, so each site's coordinates repeat once per survey year
# (a "zero distances" error). On top of that, SWS1 and SWS1A are two Plot
# codes at the identical physical location, which would break this even
# after collapsing years. Fix: run this specific robustness check at the
# SITE level (averaging across years), with SWS1/SWS1A merged into one row
# since they are literally the same location. This is a deliberate
# simplification for the spatial-correlation test only -- the main models
# (m_rich_full, m_rich_nested) still use the full site-year data.
dd_rich_site <- dd_rich_geo %>%
  mutate(Plot = if_else(Plot == "SWS1", "SWS1A", Plot)) %>%  # merge same physical site
  group_by(Plot, Type, Longitude, Latitude) %>%
  summarise(
    Value       = mean(Value, na.rm = TRUE),
    previous1_c = mean(previous1_c, na.rm = TRUE),
    same_c      = mean(same_c, na.rm = TRUE),
    n_years     = n(),
    .groups = "drop"
  ) %>%
  mutate(Type = forcats::fct_relevel(Type, "Reserve", "Fragment"))  # match main model's reference level

cat("\n--- Site-level data for spatial GLS (averaged across years) ---\n")
cat("Rows (unique sites):", nrow(dd_rich_site), "\n")
cat("Any duplicate coordinates remaining?",
    any(duplicated(dd_rich_site[, c("Longitude","Latitude")])), "\n")

gls_nospace <- nlme::gls(
  Value ~ Type * previous1_c + Type * same_c,
  data = dd_rich_site,
  method = "ML"
)

gls_space <- nlme::gls(
  Value ~ Type * previous1_c + Type * same_c,
  data = dd_rich_site,
  correlation = nlme::corExp(form = ~ Longitude + Latitude, nugget = TRUE),
  method = "ML"
)

cat("\n--- AIC comparison: no spatial structure vs. exponential spatial correlation ---\n")
print(AIC(gls_nospace, gls_space))

cat("\n--- Type effect, no spatial structure ---\n")
print(summary(gls_nospace)$tTable)

cat("\n--- Type effect, with spatial correlation structure ---\n")
print(summary(gls_space)$tTable)

cat("\n--- Side-by-side: does the Type effect survive spatial correction? ---\n")

# Grab whichever Type MAIN EFFECT coefficient row exists (name depends on
# factor reference level) -- exclude interaction terms (TypeX:previous1_c
# etc.), which also start with "Type" and would otherwise match too.
type_row_nospace <- grep("^Type[^:]+$", rownames(summary(gls_nospace)$tTable), value = TRUE)
type_row_space   <- grep("^Type[^:]+$", rownames(summary(gls_space)$tTable),   value = TRUE)
cat("Type coefficient name in gls_nospace:", type_row_nospace, "\n")
cat("Type coefficient name in gls_space:  ", type_row_space, "\n")

type_compare <- tibble(
  Model    = c("No spatial structure", "Exponential spatial correlation"),
  Estimate = c(summary(gls_nospace)$tTable[type_row_nospace, "Value"],
               summary(gls_space)$tTable[type_row_space, "Value"]),
  SE       = c(summary(gls_nospace)$tTable[type_row_nospace, "Std.Error"],
               summary(gls_space)$tTable[type_row_space, "Std.Error"]),
  P        = c(summary(gls_nospace)$tTable[type_row_nospace, "p-value"],
               summary(gls_space)$tTable[type_row_space, "p-value"])
)
print(type_compare)

#### 6) Does clustering WITHIN habitat type (not overall) explain the Type
#### effect? Reviewer's specific concern: several reserve sites sit inside
#### the same reserve complex (e.g., 4 sites at Mission Trails), so they may
#### not be independent replicates of "reserve" the way scattered fragments
#### are independent replicates of "fragment." This is a pseudoreplication
#### question, not a smooth-distance-decay question, so it needs its own test.
####
#### 6a) Descriptive: are reserve-reserve distances systematically smaller
####     than fragment-fragment distances? (quantifies "city layout artefact")
#### 6b) Formal: does Type still predict richness once sites are nested
####     within their reserve/locality complex AND year is accounted for,
####     i.e. once within-complex non-independence is explicitly modeled as
####     its own random effect alongside the year correction used throughout
####     the rest of the paper?
##############################################################################

# Locality groupings from Table S1 (supplementary site metadata) -- reserve
# sites are grouped into their actual named reserve/protected-area complex;
# fragment sites, which are geographically dispersed through the urban
# matrix with no shared complex, are each treated as their own locality.
locality_lookup <- tribble(
  ~Plot,   ~Locality,
  "MTE2",  "Mission Trails Regional Park",
  "MTE3",  "Mission Trails Regional Park",
  "MTE4",  "Mission Trails Regional Park",
  "MTI2",  "Mission Trails Regional Park",
  "SWEA",  "San Diego National Wildlife Refuge",
  "SWI2",  "San Diego National Wildlife Refuge",
  "SWI4",  "San Diego National Wildlife Refuge",
  "ECR4",  "Elliott Chaparral Reserve",
  "ECR5",  "Elliott Chaparral Reserve"
)
# EDIT if your fragment site list differs -- any Plot not listed above is
# assigned its own singleton locality (== itself) below.

# rebuild locality-joined version off the corrected type_geo
type_geo_locality <- type_geo %>%
  left_join(locality_lookup, by = "Plot") %>%
  mutate(Locality = if_else(is.na(Locality), Plot, Locality))

pairwise_by_type <- function(df, type_filter) {
  sub <- df %>% filter(Type == type_filter)
  if (nrow(sub) < 2) return(tibble(Type = type_filter, pair_dist_km = numeric(0)))
  coords_sub <- sub %>% dplyr::select(Longitude, Latitude) %>% as.matrix()
  d <- geosphere::distm(coords_sub, fun = distHaversine) / 1000  # km
  d[upper.tri(d)] %>% tibble(Type = type_filter, pair_dist_km = .) %>% dplyr::select(Type, pair_dist_km)
}

## 6a) Pairwise distances within Reserve sites vs. within Fragment sites ----

pairwise_dists <- bind_rows(
  pairwise_by_type(type_geo_locality, "Reserve"),
  pairwise_by_type(type_geo_locality, "Fragment")
)

cat("\n--- Pairwise inter-site distance (km), by habitat type ---\n")
print(pairwise_dists %>% group_by(Type) %>%
        summarise(n_pairs = n(), mean_km = mean(pair_dist_km),
                  median_km = median(pair_dist_km), min_km = min(pair_dist_km), .groups = "drop"))

# Permutation test: is mean Reserve-Reserve distance significantly smaller
# than mean Fragment-Fragment distance, more than expected by chance given
# the overall spatial layout of all sites?
set.seed(42)
observed_diff <- pairwise_dists %>%
  group_by(Type) %>% summarise(m = mean(pair_dist_km), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = m) %>%
  mutate(diff = Fragment - Reserve) %>% pull(diff)

all_coords <- type_geo_locality %>% dplyr::select(Plot, Type, Longitude, Latitude)
n_reserve <- sum(all_coords$Type == "Reserve")

perm_diffs <- replicate(4999, {
  shuffled_type <- sample(all_coords$Type)
  tmp <- all_coords %>% mutate(Type = shuffled_type)
  d_r <- pairwise_by_type(tmp, "Reserve")$pair_dist_km
  d_f <- pairwise_by_type(tmp, "Fragment")$pair_dist_km
  mean(d_f) - mean(d_r)
})

perm_p <- mean(perm_diffs >= observed_diff)
cat("\nObserved (mean Fragment-Fragment dist) - (mean Reserve-Reserve dist):",
    round(observed_diff, 2), "km\n")
cat("Permutation P (is Reserve clustering more than chance, given site layout):",
    round(perm_p, 4), "\n")

## 6b) Nested random effect: Plot nested within Locality, PLUS year ---------
# If several reserve sites share a Locality (e.g., 4 sites at Mission
# Trails), (1 | Locality/Plot) explicitly treats Locality as the higher-level
# unit of replication, so Reserve's effective sample size shrinks to the
# number of distinct reserve complexes rather than the number of individual
# reserve plots. (1 | year) is added alongside this, matching the primary
# richness model everywhere else in the paper. Compare Type's estimate/
# significance to the original (1|Plot)+(1|year) model above.

dd_rich_locality <- dd_rich %>%
  left_join(type_geo_locality %>% dplyr::select(Plot, Locality), by = "Plot")

m_rich_nested <- lmer(
  Value ~ Type * previous1_c + Type * same_c + (1 | Locality/Plot) + (1 | year),
  data = dd_rich_locality, REML = TRUE
)
cat("\nSingular fit (nested model)?", lme4::isSingular(m_rich_nested), "\n")

cat("\n--- Richness model with Plot nested within reserve/fragment Locality, +year ---\n")
print(car::Anova(m_rich_nested, type = "III"))

cat("\n--- Variance components: how much variance sits at the Locality level? ---\n")
print(as.data.frame(VarCorr(m_rich_nested)))

cat("\n--- Compare Type effect: original (1|Plot)+(1|year) vs. nested (1|Locality/Plot)+(1|year) ---\n")
nested_compare <- tibble(
  Model = c("Original: (1|Plot)+(1|year)", "Nested: (1|Locality/Plot)+(1|year)"),
  # pull Type main-effect p-value from each ANOVA table for comparability
  Type_Chisq = c(car::Anova(m_rich_full, type = "III")["Type", "Chisq"],
                 car::Anova(m_rich_nested, type = "III")["Type", "Chisq"]),
  Type_P     = c(car::Anova(m_rich_full, type = "III")["Type", "Pr(>Chisq)"],
                 car::Anova(m_rich_nested, type = "III")["Type", "Pr(>Chisq)"])
)
print(nested_compare)

n_distinct_reserve_localities <- type_geo_locality %>%
  filter(Type == "Reserve") %>% distinct(Locality) %>% nrow()
n_distinct_fragment_localities <- type_geo_locality %>%
  filter(Type == "Fragment") %>% distinct(Locality) %>% nrow()

cat("\nEffective number of independent Reserve localities:", n_distinct_reserve_localities, "\n")
cat("Effective number of independent Fragment localities:", n_distinct_fragment_localities, "\n")

#### 7) Save all results for the rebuttal letter / supplement ---------------

saveRDS(
  list(
    moran_I_residuals    = moran_result,
    moran_I_type         = moran_type_result,
    correlog              = correlog_result,
    gls_nospace           = gls_nospace,
    gls_space             = gls_space,
    aic_comparison        = AIC(gls_nospace, gls_space),
    type_effect_compare   = type_compare,
    pairwise_dists        = pairwise_dists,
    perm_test_p           = perm_p,
    m_rich_nested         = m_rich_nested,
    nested_vs_original    = nested_compare,
    n_reserve_localities  = n_distinct_reserve_localities,
    n_fragment_localities = n_distinct_fragment_localities
  ),
  file.path(results_dir, "spatial_autocorrelation_results.rds")
)

cat("\nDone. Results saved to",
    file.path(results_dir, "spatial_autocorrelation_results.rds"), "\n")