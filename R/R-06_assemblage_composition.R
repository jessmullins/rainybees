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
## Date    : 6-Apr-2026
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

indicators <- c("Halictus farinosus", "Lasioglossum sisymbrii",
                "Augochlorella pomoniella", "Eucera dorsata", "Eucera tricinctella",
                "Lasioglossum punctatoventre", "Lasioglossum robustum",
                "Megachile subnigra", "Melissodes plumosus", "Anthidium jocosum",
                "Ceratina nanula", "Lasioglossum MSN turgiventre",
                "Lasioglossum perparvum", "Melissodes tessellatus",
                "Micralictoides ruficaudus", "Andrena osmioides",
                "Perdita interrupta interrupta", "Dianthidium dubium",
                "Dufourea sandhouseae", "Andrena ablegata",
                "Dianthidium pudicum consimile", "Megachile coquilletti",
                "Osmia gabrielis", "Anthophorula torticornis", "Bombus vosnesenskii",
                "Lasioglossum titusi", "Lasioglossum actinosum",
                "Lasioglossum stictaspis complex", "Osmia granulosa",
                "Halictus ligatus", "Bombus californicus", "Osmia clarescens",
                "Hylaeus mesillae cressoni", "Diadasia opuntiae",
                "Calliopsis obscurella", "Hesperapis fuchsi", "Melissodes stearnsi",
                "Ashmeadiella californica", "Dufourea rhamni", "Andrena atypica",
                "Lasioglossum perichlarum", "Megachile onobrychidis",
                "Dufourea brevicornis", "Melissodes velutinus", "Diadasia ochracea",
                "Perdita claypolei australior")

# Filter indicators for sensitivity analysis — comment out to include all species
# df <- df %>%
 # filter(!scientificName %in% indicators)

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

# Global test
a_global <- anova(mod_main, permutations = 999)
print(a_global)

# Axes test
a_axes <- anova(mod_main, by = "axis", permutations = 999)
print(a_axes)

# Terms test
a_terms <- anova(mod_main, by = "terms", permutations = 999)
print(a_terms)

RsquareAdj(mod_main)

# =============================================================================
# SECTION 5 — Variance partitioning (partial RDAs)
# =============================================================================
# not included in manuscript but interesting 

uniq_type <- anova(rda(spec_h_rate ~ Type     + Condition(same + previous1),
                       data = env_stand), permutations = 999)$Variance[1]
uniq_same <- anova(rda(spec_h_rate ~ same     + Condition(Type + previous1),
                       data = env_stand), permutations = 999)$Variance[1]
uniq_prev <- anova(rda(spec_h_rate ~ previous1 + Condition(Type + same),
                       data = env_stand), permutations = 999)$Variance[1]

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

eig_vals  <- eigenvals(mod_main, model = "constrained")
prop_var  <- eig_vals / sum(eig_vals)
axis1_lab <- paste0("RDA1 (", round(100 * prop_var[1], 1), "%)")
axis2_lab <- paste0("RDA2 (", round(100 * prop_var[2], 1), "%)")

p_rda <- ggplot(sc_sites, aes(RDA1, RDA2)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_point(aes(fill = previous1_mm, shape = Type),
             size = 6, color = "black", stroke = 0.6) +
  scale_shape_manual(values = c(Reserve = 22, Fragment = 21),
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

sens_forms <- list(
  base          = spec_h_rate ~ Type + same + previous1,
  add_prev2     = spec_h_rate ~ Type + same + previous1 + previous2,
  add_prev3     = spec_h_rate ~ Type + same + previous1 + previous3,
  prev1_prev2   = spec_h_rate ~ Type + same + previous1 + previous2,
  prev2_prev3   = spec_h_rate ~ Type + same + previous1 + previous2 + previous3
)

sens_out <- lapply(names(sens_forms), function(nm) {
  m <- rda(sens_forms[[nm]], data = env_stand)
  global <- anova(m, permutations = 999) %>%
    as.data.frame() %>% rownames_to_column("Row") %>%
    mutate(model = nm, .before = 1) %>%
    rename(p_value = `Pr(>F)`)
  terms <- anova(m, by = "terms", permutations = 999) %>%
    as.data.frame() %>% rownames_to_column("Term") %>%
    mutate(model = nm, .before = 1) %>%
    rename(p_value = `Pr(>F)`)
  overview <- tibble(model = nm, adj_R2 = RsquareAdj(m)$adj.r.squared)
  list(overview = overview, global = global, terms = terms)
})
names(sens_out) <- names(sens_forms)

writexl::write_xlsx(
  list(
    overview   = bind_rows(lapply(sens_out, `[[`, "overview")),
    global     = bind_rows(lapply(sens_out, `[[`, "global")),
    term_tests = bind_rows(lapply(sens_out, `[[`, "terms"))
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

