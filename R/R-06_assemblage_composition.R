## =============================================================================
## R-06_assemblage_composition.R
##
## Purpose : Test Prediction 3 — winter precipitation and land-use change
##           alter bee assemblage composition.
##
##           (1) Redundancy analysis (RDA): constrained ordination of
##               Hellinger-transformed, effort-normalized abundances on
##               site type, WPy0, and WPy0-1. Includes variance partitioning
##               via partial RDAs and sensitivity analysis with additional
##               precipitation lags.
##
##           (2) Indicator species analysis (IndVal.g): identifies species
##               statistically associated with combinations of site type
##               (Reserve/Fragment) and precipitation conditions (Wet/Dry)
##               for both same-year and previous-year precipitation windows.
##
## Inputs  : data/cleaned/df_analysis.csv    — filtered occurrence records
##                                             with trap_days (from R-01)
##           data/cleaned/precipyrs.csv      — annual winter totals + lags
##                                             (from R-02)
##           data/raw/traits.csv             — species traits for indicator
##                                             species table annotation
##
## Outputs : figures/rda_biplot.svg                        — RDA biplot
##           figures/rda_variance_partition.svg            — variance partition bars
##           figures/indicator_euler_same.png              — Euler diagram (same-yr)
##           figures/indicator_euler_previous.png          — Euler diagram (prev-yr)
##           figures/indicator_bar_same.png                — indicator bar (same-yr)
##           figures/indicator_bar_previous.png            — indicator bar (prev-yr)
##           results/Table_S_RDA_results.xlsx              — RDA stats
##           results/Table_S_RDA_sensitivity.xlsx          — lag sensitivity
##           results/Table_S_same_indicators_traits.csv    — same-yr indicators
##           results/Table_S_previous_indicators_traits.csv — prev-yr indicators
##
## Author  : Jess Mullins
## Date    : 6-Apr-2026, edited August 24, 2026
## =============================================================================

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(tidyverse)
  library(indicspecies)
  library(eulerr)
  library(ggfittext)
  library(patchwork)
  library(scales)
  library(writexl)
  library(ggrepel)
})

## File paths
figures_dir <- "figures"
results_dir <- "results"
cleaned_dir <- "data/cleaned"
raw_dir     <- "data/raw"

set.seed(123)

# =============================================================================
# SECTION 1 — Load data
# =============================================================================

df <- read.csv(file.path(cleaned_dir, "df_analysis.csv")) %>%
  mutate(
    year      = factor(year),
    Count     = as.numeric(individualCount),
    sample_id = paste(fieldNumber, year, sep = "_")
  )

precip <- read.csv(file.path(cleaned_dir, "precipyrs.csv")) %>%
  mutate(year = factor(year))

traits_tbl        <- read.csv(file.path(raw_dir, "traits.csv"))
trait_species_col <- "scientificName"

survey_years <- c("2011","2012","2015","2016","2022","2023","2024")

# remove indicator spp 
indicators <- readRDS(file.path(cleaned_dir, "indicator_species.rds"))

# Filter indicators for sensitivity analysis 
df <- df %>%
  filter(!scientificName == "Perdita") # no genus-level indicators
# , (!scientificName %in% indicators)) # unhash to remove indicators

# =============================================================================
# SECTION 2 — Build community matrices
#
# spec_h_rate: Hellinger-transformed effort-normalized abundances
# This is the MAIN matrix for RDA per methods (Borcard et al. 2011)
# =============================================================================

# 2a) Effort-normalized rates (count / trap_days)
counts_long <- df %>%
  group_by(sample_id, scientificName, trap_days) %>%
  summarise(count = sum(Count), .groups = "drop")

comm_rate <- counts_long %>%
  mutate(
    count    = as.numeric(count),
    trap_days = as.numeric(trap_days),
    rate     = count / trap_days
  ) %>%
  group_by(sample_id, scientificName) %>%
  summarise(rate = sum(rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = scientificName, values_from = rate,
              values_fill = 0) %>%
  column_to_rownames("sample_id")

# 2b) Hellinger transformation
spec_h_rate <- decostand(comm_rate, method = "hellinger")

# =============================================================================
# SECTION 3 — Environmental matrix
# =============================================================================

env_data <- df %>%
  dplyr::select(sample_id, year, fieldNumber, Type,
                decimalLatitude, decimalLongitude, trap_days) %>%
  dplyr::arrange(sample_id) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  left_join(
    precip %>%
      dplyr::select(year, same, previous1, previous2, previous3),
    by = "year"
  ) %>%
  tibble::column_to_rownames("sample_id")

# Align env to community matrix
env_data <- env_data[rownames(spec_h_rate), , drop = FALSE]

# Standardize numeric predictors
env_stand <- env_data %>%
  mutate(Type = factor(Type)) %>%
  mutate(across(c(same, previous1, previous2, previous3, trap_days),
                ~ scale(log1p(.x))[, 1]))

env_plain <- env_data   # unscaled copy for plotting

# =============================================================================
# SECTION 4 — RDA: main model
# =============================================================================

mod_main <- rda(spec_h_rate ~ Type + same + previous1, data = env_stand)

summary(mod_main)
vif.cca(mod_main)

# --- Year-blocked permutations -----------------------------------------------
# Same rationale as the mixed models elsewhere: same/previous1 (WPy0/WPy0-1)
# are year-level variables, so sites sampled in the same year are not
# independent with respect to precipitation. Unrestricted permutation
# (shuffling site-years freely) overstates the effective sample size for
# these terms, inflating Type I error risk -- the same concern raised for
# the diversity and body-size models. Restricting permutations to within
# year (permute::how(blocks = ...)) is the RDA-appropriate analogue of
# adding (1 | year) to a mixed model.
perm_year_rda <- how(nperm = 999, blocks = env_stand$year)

# Global test
a_global <- anova(mod_main, permutations = perm_year_rda)
print(a_global)

# Axes test
a_axes <- anova(mod_main, by = "axis", permutations = perm_year_rda)
print(a_axes)

# Terms test
a_terms <- anova(mod_main, by = "terms", permutations = perm_year_rda)
print(a_terms)

RsquareAdj(mod_main)

# =============================================================================
# SECTION 5 — Variance partitioning (partial RDAs)
# =============================================================================
# not included in manuscript but interesting 

uniq_type <- anova(rda(spec_h_rate ~ Type     + Condition(same + previous1),
                       data = env_stand), permutations = perm_year_rda)$Variance[1]
uniq_same <- anova(rda(spec_h_rate ~ same     + Condition(Type + previous1),
                       data = env_stand), permutations = perm_year_rda)$Variance[1]
uniq_prev <- anova(rda(spec_h_rate ~ previous1 + Condition(Type + same),
                       data = env_stand), permutations = perm_year_rda)$Variance[1]

explained_var <- a_global$Variance[1]
residual_var  <- a_global$Variance[2]
total_var     <- explained_var + residual_var
shared        <- max(explained_var - (uniq_type + uniq_same + uniq_prev), 0)

R2  <- explained_var / total_var
R2a <- RsquareAdj(mod_main)$adj.r.squared

# Variance partition figures
df_total <- tibble(
  component = factor(c("Explained", "Residual"),
                     levels = c("Explained", "Residual")),
  value = c(explained_var, residual_var)
) %>%
  mutate(pct = value / sum(value),
         pct_lab = percent(pct, accuracy = 1))

df_expl <- tibble(
  component = factor(
    c("Unique Type", "Unique WPy0-1", "Unique WPy0", "Shared"),
    levels = c("Unique Type", "Unique WPy0-1", "Unique WPy0", "Shared")
  ),
  value = c(uniq_type, uniq_prev, uniq_same, shared)
) %>%
  mutate(pct = value / sum(value),
         pct_lab = ifelse(pct >= 0.06, percent(pct, accuracy = 1), ""))

pal_total <- c("Explained" = "grey35", "Residual" = "grey80")
pal_expl  <- c("Unique Type"    = "#D55E00",
               "Unique WPy0-1"  = "#009E73",
               "Unique WPy0"    = "#0072B2",
               "Shared"         = "#999999")

th <- theme_classic(base_size = 16) +
  theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y  = element_blank(), legend.position = "right",
        legend.title = element_blank(),
        plot.title   = element_text(face = "bold"))

p_total <- ggplot(df_total, aes(x = 1, y = pct, fill = component)) +
  geom_col(width = 0.6) +
  geom_fit_text(aes(label = pct_lab), reflow = TRUE, min.size = 8) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = pal_total) +
  labs(title    = "Total variance",
       subtitle = paste0("R² = ", round(R2, 2), "  |  R²adj = ", round(R2a, 2)),
       x = NULL, y = NULL) + th

p_expl <- ggplot(df_expl, aes(x = 1, y = pct, fill = component)) +
  geom_col(width = 0.6) +
  geom_fit_text(aes(label = pct_lab), reflow = TRUE, min.size = 8) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = pal_expl) +
  labs(title = "Explained variance", x = NULL, y = NULL) + th

p_vp <- p_total / p_expl + plot_layout(heights = c(1, 1))
print(p_vp)

ggsave(file.path(figures_dir, "rda_variance_partition.svg"),
       p_vp, width = 8, height = 5)

# =============================================================================
# SECTION 6 — RDA biplot
# =============================================================================

sc_sites <- scores(mod_main, display = "sites", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(
    Type         = env_plain[sample_id, "Type"],
    year         = env_plain[sample_id, "year"],
    previous1_mm = as.numeric(env_plain[sample_id, "previous1"])
  )

sc_bp <- scores(mod_main, display = "bp", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("var")

vector_labels <- tibble(
  var   = c("previous1", "same", "TypeReserve"),
  label = c("WPy0-1", "WPy0", "Reserve")
)

sc_bp2 <- sc_bp %>%
  left_join(vector_labels, by = "var") %>%
  mutate(label = replace_na(label, ""))

# --- Species scores: identify the 5 taxa most strongly driving the RDA
# ordination (longest vectors in RDA1-RDA2 space), for display as light
# gray arrows -- visually distinct from the bold black environmental
# predictor vectors above.
sc_species <- scores(mod_main, display = "species", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  mutate(veclen = sqrt(RDA1^2 + RDA2^2)) %>%
  arrange(desc(veclen))

cat("\n--- Top 5 taxa driving the RDA ordination (by vector length) ---\n")
print(sc_species %>% dplyr::select(species, RDA1, RDA2, veclen) %>% slice_head(n = 5))

sc_species_top5 <- sc_species %>% slice_head(n = 5)

# Species vectors are often a very different scale from site/env vectors in
# a correlation biplot; a modest constant multiplier keeps them visible
# without overwhelming the plot. Adjust species_vec_scale if arrows look
# too short/long relative to the environmental vectors once plotted.
species_vec_scale <- 1
sc_species_top5 <- sc_species_top5 %>%
  mutate(RDA1_plot = RDA1 * species_vec_scale,
         RDA2_plot = RDA2 * species_vec_scale)

eig_vals  <- eigenvals(mod_main, model = "constrained")
prop_var  <- eig_vals / sum(eig_vals)
axis1_lab <- paste0("RDA1 (", round(100 * prop_var[1], 1), "%)")
axis2_lab <- paste0("RDA2 (", round(100 * prop_var[2], 1), "%)")

p_rda <- ggplot(sc_sites, aes(RDA1, RDA2)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4) +
  # Species vectors (top 5 taxa), light gray, drawn first so they sit
  # beneath the site points and environmental vectors
  geom_segment(data = sc_species_top5,
               aes(x = 0, y = 0, xend = RDA1_plot, yend = RDA2_plot),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.02, "npc")),
               color = "gray65", linewidth = 0.5) +
  ggrepel::geom_text_repel(data = sc_species_top5,
                           aes(x = RDA1_plot, y = RDA2_plot, label = species),
                           inherit.aes = FALSE,
                           color = "gray35", fontface = "italic", size = 3,
                           bg.color = "white", bg.r = 0.15,
                           segment.color = "gray65", segment.size = 0.3,
                           min.segment.length = 0, box.padding = 0.6, point.padding = 0.3,
                           force = 8, max.iter = 10000, max.overlaps = Inf, seed = 42) +
  geom_point(aes(fill = previous1_mm, shape = Type),
             size = 6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = c(Reserve = 24, Fragment = 21),
                     name = "Site type") +
  scale_fill_distiller(palette = "YlGn", direction = 1,
                       name = "WPy0-1 (mm)") +
  geom_segment(data = sc_bp2,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.02, "npc")),
               linewidth = 0.6) +
  geom_text(data = sc_bp2 %>% filter(label != ""),
            aes(x = RDA1, y = RDA2, label = label),
            vjust = -0.4, size = 3.5) +
  coord_equal() +
  labs(x = axis1_lab, y = axis2_lab) +
  theme_classic(base_size = 13) +
  theme(legend.position = "bottom",
        legend.key.height = unit(0.8, "lines"),
        legend.key.width  = unit(0.9, "lines"))

print(p_rda)
ggsave(file.path(figures_dir, "rda_biplot.svg"),
       p_rda, width = 8, height = 8)

# =============================================================================
# SECTION 7 — PERMANOVA (adonis2)
# =============================================================================

common_ids <- intersect(rownames(spec_h_rate), rownames(env_plain))
Y <- spec_h_rate[common_ids, , drop = FALSE]
E <- env_plain[common_ids, , drop = FALSE] %>%
  mutate(Type = factor(Type),
         Plot = factor(fieldNumber),
         year = factor(year))

perm_plot <- how(nperm = 999, blocks = E$Plot)
perm_year <- how(nperm = 999, blocks = E$year)

adon_plot <- adonis2(Y ~ Type + same + previous1,
                     data = E, method = "euclidean",
                     permutations = perm_plot, by = "margin")
print(adon_plot)

adon_year <- adonis2(Y ~ Type + same + previous1,
                     data = E, method = "euclidean",
                     permutations = perm_year, by = "margin")
print(adon_year)

# =============================================================================
# SECTION 8 — Sensitivity: additional precipitation lags
# =============================================================================

# --- Candidate models for lag-sensitivity analysis --------------------------
# RECONSTRUCTED: this object was missing from the script entirely. Based on
# the stated purpose in the script header ("sensitivity analysis with
# additional precipitation lags") and the Methods description already in
# Appendix S1 ("Additional two- and three-year precipitation lags were added
# individually or in combination"), this tests whether adding WPy0-2 and/or
# WPy0-3 to the base model (Type + same + previous1) meaningfully changes
# RDA support, beyond the base WPy0/WPy0-1 model already reported as
# mod_main in Section 4.
#
# NOTE: re-run and compare output against your existing Table S11/S12
# sensitivity panel in Appendix S1 to confirm it matches whatever generated
# those original numbers -- the exact formula set may differ slightly.
sens_forms <- list(
  "Base (WPy0-1 only)" = spec_h_rate ~ Type + same + previous1,
  "+ WPy0-2"            = spec_h_rate ~ Type + same + previous1 + previous2,
  "+ WPy0-3"            = spec_h_rate ~ Type + same + previous1 + previous3,
  "+ WPy0-2 + WPy0-3"   = spec_h_rate ~ Type + same + previous1 + previous2 + previous3
)

sens_out <- lapply(names(sens_forms), function(nm) {
  m <- rda(sens_forms[[nm]], data = env_stand)
  
  global <- anova(m, permutations = perm_year_rda) %>%
    as.data.frame() %>% rownames_to_column("Row") %>%
    mutate(model = nm, .before = 1) %>%
    rename(p_value = `Pr(>F)`)
  
  terms <- anova(m, by = "terms", permutations = perm_year_rda) %>%
    as.data.frame() %>% rownames_to_column("Term") %>%
    mutate(model = nm, .before = 1) %>%
    rename(p_value = `Pr(>F)`)
  
  overview <- tibble(model = nm, adj_R2 = RsquareAdj(m)$adj.r.squared)
  
  # VIF for each predictor in this model
  vif_vals <- vegan::vif.cca(m)
  vif_tbl  <- tibble(
    model     = nm,
    predictor = names(vif_vals),
    VIF       = as.numeric(vif_vals)
  )
  
  list(overview = overview, global = global, terms = terms, vif = vif_tbl)
})
names(sens_out) <- names(sens_forms)

writexl::write_xlsx(
  list(
    overview   = bind_rows(lapply(sens_out, `[[`, "overview")),
    global     = bind_rows(lapply(sens_out, `[[`, "global")),
    term_tests = bind_rows(lapply(sens_out, `[[`, "terms")),
    VIF        = bind_rows(lapply(sens_out, `[[`, "vif"))
  ),
  file.path(results_dir, "Table_S_RDA_sensitivity.xlsx")
)

# =============================================================================
# SECTION 9 — Export RDA results
# =============================================================================

writexl::write_xlsx(
  list(
    overview   = tibble(n_sites   = nrow(spec_h_rate),
                        n_species = ncol(spec_h_rate),
                        adj_R2    = RsquareAdj(mod_main)$adj.r.squared),
    global     = a_global %>% as.data.frame() %>%
      rownames_to_column("Row") %>% rename(p_value = `Pr(>F)`),
    axes       = a_axes  %>% as.data.frame() %>%
      rownames_to_column("Axis") %>% rename(p_value = `Pr(>F)`),
    terms      = a_terms %>% as.data.frame() %>%
      rownames_to_column("Term") %>% rename(p_value = `Pr(>F)`),
    VIFs       = tibble(Predictor = names(vif.cca(mod_main)),
                        VIF       = as.numeric(vif.cca(mod_main))),
    PERMANOVA_plot = adon_plot %>% as.data.frame() %>%
      rownames_to_column("Term") %>% rename(p_value = `Pr(>F)`),
    PERMANOVA_year = adon_year %>% as.data.frame() %>%
      rownames_to_column("Term") %>% rename(p_value = `Pr(>F)`)
  ),
  file.path(results_dir, "Table_S_RDA_results.xlsx")
)

# =============================================================================
# SECTION 10 — Indicator species analysis
# =============================================================================

# Wet/dry year definitions (terciles; middle year excluded)
same_wet_years <- c("2024", "2011", "2023")
same_dry_years <- c("2015", "2022", "2012")
prev_wet_years <- c("2011", "2012", "2024")
prev_dry_years <- c("2016", "2015", "2023")

grp_cols <- c("s.Fragment_Dry", "s.Fragment_Wet",
              "s.Reserve_Dry",  "s.Reserve_Wet")

# ---- Helpers ----

standardize_species <- function(x) {
  str_trim(x) %>% str_replace_all("\\s+", "_")
}

run_indval_precip <- function(env_data, spec_mat,
                              wet_years, dry_years,
                              type_col = "Type",
                              nperm = 1200, seed = 1) {
  wet_years <- as.character(wet_years)
  dry_years <- as.character(dry_years)
  
  sub_env <- env_data %>%
    as.data.frame() %>%
    rownames_to_column("plot_year") %>%
    mutate(year_chr = as.character(year)) %>%
    filter(year_chr %in% c(wet_years, dry_years)) %>%
    mutate(
      precip_group = if_else(year_chr %in% wet_years, "Wet", "Dry"),
      group = interaction(.data[[type_col]], precip_group, sep = "_")
    ) %>%
    column_to_rownames("plot_year") %>%
    droplevels()
  
  keep_ids <- intersect(rownames(spec_mat), rownames(sub_env))
  sub_env  <- sub_env[keep_ids, , drop = FALSE]
  sub_comm <- spec_mat[keep_ids, , drop = FALSE]
  stopifnot(identical(rownames(sub_env), rownames(sub_comm)))
  
  nz_cols <- colSums(sub_comm, na.rm = TRUE) > 0
  sub_mat <- as.matrix(sub_comm[, nz_cols, drop = FALSE])
  storage.mode(sub_mat) <- "double"
  
  set.seed(seed)
  ind_res <- multipatt(sub_mat, sub_env$group,
                       func = "IndVal.g", duleg = FALSE,
                       control = how(nperm = nperm))
  
  list(ind_res     = ind_res,
       sub_env     = sub_env,
       group_sizes = table(sub_env$group),
       years       = list(wet = wet_years, dry = dry_years))
}

extract_sig_df <- function(ind_res, timing_label,
                           grp_cols, p_cutoff = 0.05) {
  stopifnot(all(grp_cols %in% colnames(ind_res$sign)))
  ind_res$sign %>%
    as.data.frame() %>%
    rownames_to_column("species_raw") %>%
    mutate(species = standardize_species(species_raw),
           Timing  = timing_label) %>%
    filter(!is.na(p.value) & p.value <= p_cutoff) %>%
    transmute(species, Timing,
              across(all_of(grp_cols)),
              IndVal  = stat,
              P_value = p.value) %>%
    rename_with(~ str_replace(.x, "^s\\.", ""), starts_with("s.")) %>%
    arrange(P_value)
}

add_traits <- function(sig_df, traits_tbl, trait_species_col) {
  traits_clean <- traits_tbl %>%
    mutate(species = standardize_species(.data[[trait_species_col]])) %>%
    select(-all_of(trait_species_col)) %>%
    relocate(species)
  left_join(sig_df, traits_clean, by = "species")
}

make_euler_from_sig <- function(sig_df,
                                sets = c("Fragment_Dry","Fragment_Wet","Reserve_Dry","Reserve_Wet")) {
  stopifnot(all(sets %in% colnames(sig_df)))
  set_list <- lapply(sets, function(g)
    unique(sig_df$species[sig_df[[g]] > 0]))
  names(set_list) <- sets
  euler(set_list)
}

save_euler_png <- function(venn_fit, filename,
                           width_px = 4800, height_px = 3000, res = 600,
                           col_map = NULL, lab_map = NULL) {
  set_names <- rownames(venn_fit$ellipses)
  if (is.null(col_map))
    col_map <- c("Fragment_Dry" = "#F2E5D7", "Fragment_Wet" = "grey70",
                 "Reserve_Dry"  = "#A6D96A", "Reserve_Wet"  = "#1B7837")
  if (is.null(lab_map))
    lab_map <- c("Fragment_Dry" = "F Dry", "Fragment_Wet" = "F Wet",
                 "Reserve_Dry"  = "R Dry", "Reserve_Wet"  = "R Wet")
  png(filename, width = width_px, height = height_px, res = res)
  plot(venn_fit,
       fill       = col_map[set_names],
       alpha      = 0.65,
       quantities = list(type = "counts", cex = 1.1),
       labels     = list(labels = lab_map[set_names],
                         font = 2, cex = 1.2, col = "black"),
       shape      = "ellipse")
  dev.off()
}

# ---- Run analyses ----

same_out <- run_indval_precip(env_data, spec_h_rate,
                              same_wet_years, same_dry_years, seed = 1)
prev_out <- run_indval_precip(env_data, spec_h_rate,
                              prev_wet_years, prev_dry_years, seed = 1)

print(same_out$group_sizes)
print(prev_out$group_sizes)

# ---- Extract significant species ----

same_sig <- extract_sig_df(same_out$ind_res, "Same",     grp_cols)
prev_sig <- extract_sig_df(prev_out$ind_res, "Previous", grp_cols)

same_sig <- same_sig %>%
  mutate(Also_previous = if_else(species %in% prev_sig$species, "Yes", "No"))
prev_sig <- prev_sig %>%
  mutate(Also_same = if_else(species %in% same_sig$species, "Yes", "No"))

cat("SAME: n species =", n_distinct(same_sig$species), "\n")
cat("PREV: n species =", n_distinct(prev_sig$species), "\n")
cat("Shared =", length(intersect(same_sig$species, prev_sig$species)), "\n")

# ---- Add traits + export ----

col_display <- c("Fragment_Dry","Fragment_Wet","Reserve_Dry","Reserve_Wet")

same_sig_traits <- add_traits(same_sig, traits_tbl, trait_species_col) %>%
  mutate(Species = str_replace_all(species, "_", " ")) %>%
  relocate(Species, .before = species)

prev_sig_traits <- add_traits(prev_sig, traits_tbl, trait_species_col) %>%
  mutate(Species = str_replace_all(species, "_", " ")) %>%
  relocate(Species, .before = species)

write_csv(same_sig_traits,
          file.path(results_dir, "Table_S_same_indicators_traits.csv"))
write_csv(prev_sig_traits,
          file.path(results_dir, "Table_S_previous_indicators_traits.csv"))

# ---- Euler diagrams ----
# ended up making this manually in Google Slides
venn_same <- make_euler_from_sig(same_sig)
venn_prev <- make_euler_from_sig(prev_sig)

save_euler_png(venn_same,
               file.path(figures_dir, "indicator_euler_same.png"))
save_euler_png(venn_prev,
               file.path(figures_dir, "indicator_euler_previous.png"))

# get counts for all combinations
# ── Species sets ──────────────────────────────────────────────────────────────
same_sp <- unique(same_sig$species)
prev_sp <- unique(prev_sig$species)

# Groups within same_sig
same_frag_dry  <- same_sig$species[same_sig$Fragment_Dry  == 1]
same_frag_wet  <- same_sig$species[same_sig$Fragment_Wet  == 1]
same_res_dry   <- same_sig$species[same_sig$Reserve_Dry   == 1]
same_res_wet   <- same_sig$species[same_sig$Reserve_Wet   == 1]

# Groups within prev_sig
prev_frag_dry  <- prev_sig$species[prev_sig$Fragment_Dry  == 1]
prev_frag_wet  <- prev_sig$species[prev_sig$Fragment_Wet  == 1]
prev_res_dry   <- prev_sig$species[prev_sig$Reserve_Dry   == 1]
prev_res_wet   <- prev_sig$species[prev_sig$Reserve_Wet   == 1]

# ── Counts ────────────────────────────────────────────────────────────────────
count_summary <- tribble(
  ~Region,                              ~N,
  # Top-level overlap
  "Same only",                          length(setdiff(same_sp, prev_sp)),
  "Previous only",                      length(setdiff(prev_sp, same_sp)),
  "Same & Previous",                    length(intersect(same_sp, prev_sp)),
  
  # Within Same — by group
  "Same: Fragment Dry only",            length(setdiff(same_frag_dry,  c(same_frag_wet, same_res_dry, same_res_wet))),
  "Same: Fragment Wet only",            length(setdiff(same_frag_wet,  c(same_frag_dry, same_res_dry, same_res_wet))),
  "Same: Reserve Dry only",             length(setdiff(same_res_dry,   c(same_frag_dry, same_frag_wet, same_res_wet))),
  "Same: Reserve Wet only",             length(setdiff(same_res_wet,   c(same_frag_dry, same_frag_wet, same_res_dry))),
  "Same: Frag Dry & Frag Wet",          length(intersect(same_frag_dry,  same_frag_wet)),
  "Same: Frag Dry & Res Dry",           length(intersect(same_frag_dry,  same_res_dry)),
  "Same: Frag Dry & Res Wet",           length(intersect(same_frag_dry,  same_res_wet)),
  "Same: Frag Wet & Res Dry",           length(intersect(same_frag_wet,  same_res_dry)),
  "Same: Frag Wet & Res Wet",           length(intersect(same_frag_wet,  same_res_wet)),
  "Same: Res Dry & Res Wet",            length(intersect(same_res_dry,   same_res_wet)),
  "Same: all 4 groups",                 length(Reduce(intersect, list(same_frag_dry, same_frag_wet, same_res_dry, same_res_wet))),
  
  # Within Previous — by group
  "Prev: Fragment Dry only",            length(setdiff(prev_frag_dry,  c(prev_frag_wet, prev_res_dry, prev_res_wet))),
  "Prev: Fragment Wet only",            length(setdiff(prev_frag_wet,  c(prev_frag_dry, prev_res_dry, prev_res_wet))),
  "Prev: Reserve Dry only",             length(setdiff(prev_res_dry,   c(prev_frag_dry, prev_frag_wet, prev_res_wet))),
  "Prev: Reserve Wet only",             length(setdiff(prev_res_wet,   c(prev_frag_dry, prev_frag_wet, prev_res_dry))),
  "Prev: Frag Dry & Frag Wet",          length(intersect(prev_frag_dry,  prev_frag_wet)),
  "Prev: Frag Dry & Res Dry",           length(intersect(prev_frag_dry,  prev_res_dry)),
  "Prev: Frag Dry & Res Wet",           length(intersect(prev_frag_dry,  prev_res_wet)),
  "Prev: Frag Wet & Res Dry",           length(intersect(prev_frag_wet,  prev_res_dry)),
  "Prev: Frag Wet & Res Wet",           length(intersect(prev_frag_wet,  prev_res_wet)),
  "Prev: Res Dry & Res Wet",            length(intersect(prev_res_dry,   prev_res_wet)),
  "Prev: all 4 groups",                 length(Reduce(intersect, list(prev_frag_dry, prev_frag_wet, prev_res_dry, prev_res_wet)))
)

print(count_summary, n = Inf)

# let's save the list of indicator spp for downstream analyses
# ── Previous-year indicator species list ─────────────────────────────────────

prev_indicators <- sort(unique(gsub("_", " ", prev_sig$species)))

cat("Previous-year indicators (n =", length(prev_indicators), "):\n")
cat(paste0('"', prev_indicators, '"', collapse = ",\n"), "\n\n")

saveRDS(prev_indicators, file.path(cleaned_dir, "indicator_species.rds"))

# ── % of occurrences that were indicator species ─────────────────────────────

total_occurrences <- nrow(df_analysis)
indicator_occurrences <- sum(df_analysis$scientificName %in% indicators)

pct_indicators <- (indicator_occurrences / total_occurrences) * 100

cat("Total occurrences:          ", total_occurrences, "\n")
cat("Indicator spp occurrences:  ", indicator_occurrences, "\n")
cat("% indicator occurrences:    ", round(pct_indicators, 1), "%\n")

# =============================================================================
# SECTION 11 — Variance decomposition: wet vs. dry years (WPy0-1)
#
# Reviewer request: quantify how much of the spread in WPy0-1 is explained
# by the wet/dry grouping used in the indicator species analysis, vs. left
# over as within-group variance. Uses the same prev_wet_years/prev_dry_years
# grouping already defined above (Section 10).
# =============================================================================

wetdry_precip <- precip %>%
  mutate(year_chr = as.character(year)) %>%
  filter(year_chr %in% c(prev_wet_years, prev_dry_years)) %>%
  mutate(group = if_else(year_chr %in% prev_wet_years, "Wet", "Dry")) %>%
  dplyr::select(year_chr, group, previous1)

print(wetdry_precip)

aov_wetdry <- aov(previous1 ~ group, data = wetdry_precip)
print(summary(aov_wetdry))

# % variance explained: between-group (Wet vs. Dry) vs. within-group (residual)
ss <- summary(aov_wetdry)[[1]][["Sum Sq"]]
variance_tbl <- tibble(
  source  = c("Between groups (Wet vs. Dry)", "Within groups (residual)"),
  sum_sq  = ss,
  pct_var = round(100 * ss / sum(ss), 2)
)
print(variance_tbl)

cat("\nWPy0-1 by group (mean \u00b1 SD):\n")
wetdry_precip %>%
  group_by(group) %>%
  summarise(mean_previous1 = mean(previous1), sd_previous1 = sd(previous1),
            .groups = "drop") %>%
  print()





# Run in your R-06 session (needs prev_sig, and the same_frag_dry/wet,
# same_res_dry/wet etc. group-membership vectors already built for the
# Euler diagram / count_summary).

# ---- 1) Full species list with group memberships (WPy0-1 / prev_sig) -----
cat("=== Full indicator species list (WPy0-1), n =", n_distinct(prev_sig$species), "===\n")
print(prev_sig %>%
        select(species, Fragment_Dry, Fragment_Wet, Reserve_Dry, Reserve_Wet, IndVal, P_value) %>%
        arrange(species), n = Inf)

# ---- 2) Species names for each Euler region (for eyeballing against the figure) --
cat("\n--- Reserve-Wet ONLY ---\n")
print(setdiff(prev_res_wet, c(prev_res_dry, prev_frag_dry, prev_frag_wet)))

cat("\n--- Reserve-Dry ONLY ---\n")
print(setdiff(prev_res_dry, c(prev_res_wet, prev_frag_dry, prev_frag_wet)))

cat("\n--- Fragment-Wet ONLY ---\n")
print(setdiff(prev_frag_wet, c(prev_res_wet, prev_res_dry, prev_frag_dry)))

cat("\n--- Fragment-Dry ONLY ---\n")
print(setdiff(prev_frag_dry, c(prev_res_wet, prev_res_dry, prev_frag_wet)))

cat("\n--- Reserve-Wet & Reserve-Dry (both, no fragment) ---\n")
print(setdiff(intersect(prev_res_wet, prev_res_dry), c(prev_frag_dry, prev_frag_wet)))

cat("\n--- Reserve-Wet & Fragment-Dry (no other) ---\n")
print(setdiff(intersect(prev_res_wet, prev_frag_dry), c(prev_res_dry, prev_frag_wet)))

cat("\n--- Reserve-Dry & Fragment-Dry (no other) ---\n")
print(setdiff(intersect(prev_res_dry, prev_frag_dry), c(prev_res_wet, prev_frag_wet)))

cat("\n--- Fragment-Dry & Fragment-Wet (no reserve) ---\n")
print(setdiff(intersect(prev_frag_dry, prev_frag_wet), c(prev_res_wet, prev_res_dry)))

cat("\n--- Reserve-Dry & Fragment-Wet (no other) ---\n")
print(setdiff(intersect(prev_res_dry, prev_frag_wet), c(prev_res_wet, prev_frag_dry)))

cat("\n--- Reserve-Wet & Fragment-Wet (no other) ---\n")
print(setdiff(intersect(prev_res_wet, prev_frag_wet), c(prev_res_dry, prev_frag_dry)))

cat("\n--- All 4 groups ---\n")
print(Reduce(intersect, list(prev_frag_dry, prev_frag_wet, prev_res_dry, prev_res_wet)))

# ---- 3) Sanity check: do the region counts sum to the total unique species? ----
region_counts <- c(
  reserve_wet_only  = length(setdiff(prev_res_wet, c(prev_res_dry, prev_frag_dry, prev_frag_wet))),
  reserve_dry_only  = length(setdiff(prev_res_dry, c(prev_res_wet, prev_frag_dry, prev_frag_wet))),
  frag_wet_only     = length(setdiff(prev_frag_wet, c(prev_res_wet, prev_res_dry, prev_frag_dry))),
  frag_dry_only     = length(setdiff(prev_frag_dry, c(prev_res_wet, prev_res_dry, prev_frag_wet))),
  res_wet_res_dry   = length(setdiff(intersect(prev_res_wet, prev_res_dry), c(prev_frag_dry, prev_frag_wet))),
  res_wet_frag_dry  = length(setdiff(intersect(prev_res_wet, prev_frag_dry), c(prev_res_dry, prev_frag_wet))),
  res_dry_frag_dry  = length(setdiff(intersect(prev_res_dry, prev_frag_dry), c(prev_res_wet, prev_frag_wet))),
  frag_dry_frag_wet = length(setdiff(intersect(prev_frag_dry, prev_frag_wet), c(prev_res_wet, prev_res_dry))),
  res_dry_frag_wet  = length(setdiff(intersect(prev_res_dry, prev_frag_wet), c(prev_res_wet, prev_frag_dry))),
  res_wet_frag_wet  = length(setdiff(intersect(prev_res_wet, prev_frag_wet), c(prev_res_dry, prev_frag_dry))),
  all_four          = length(Reduce(intersect, list(prev_frag_dry, prev_frag_wet, prev_res_dry, prev_res_wet)))
)

cat("\n=== Region counts ===\n")
print(region_counts)
cat("\nSum of all regions:", sum(region_counts), "\n")
cat("Total unique species (n_distinct(prev_sig$species)):", n_distinct(prev_sig$species), "\n")
cat("Match?", sum(region_counts) == n_distinct(prev_sig$species), "\n")


# Corrected: includes the 4 triple-overlap regions missing from the first check

cat("--- Reserve-Wet & Reserve-Dry & Fragment-Dry (not Fragment-Wet) ---\n")
r1 <- setdiff(Reduce(intersect, list(prev_res_wet, prev_res_dry, prev_frag_dry)), prev_frag_wet)
print(r1)

cat("\n--- Reserve-Wet & Reserve-Dry & Fragment-Wet (not Fragment-Dry) ---\n")
r2 <- setdiff(Reduce(intersect, list(prev_res_wet, prev_res_dry, prev_frag_wet)), prev_frag_dry)
print(r2)

cat("\n--- Reserve-Wet & Fragment-Dry & Fragment-Wet (not Reserve-Dry) ---\n")
r3 <- setdiff(Reduce(intersect, list(prev_res_wet, prev_frag_dry, prev_frag_wet)), prev_res_dry)
print(r3)

cat("\n--- Reserve-Dry & Fragment-Dry & Fragment-Wet (not Reserve-Wet) ---\n")
r4 <- setdiff(Reduce(intersect, list(prev_res_dry, prev_frag_dry, prev_frag_wet)), prev_res_wet)
print(r4)

triple_counts <- c(
  res_wet_res_dry_frag_dry = length(r1),
  res_wet_res_dry_frag_wet = length(r2),
  res_wet_frag_dry_frag_wet = length(r3),
  res_dry_frag_dry_frag_wet = length(r4)
)

cat("\n=== Triple-overlap counts ===\n")
print(triple_counts)

cat("\nSum of original 11 regions: 39\n")
cat("Sum of triple overlaps:", sum(triple_counts), "\n")
cat("New total:", 39 + sum(triple_counts), "\n")
cat("Target (n_distinct(prev_sig$species)): 47\n")
cat("Match?", (39 + sum(triple_counts)) == 47, "\n")
