###Activity"

#packages to be used: 

library(projpred)
library(lsr)
library(tidybayes)
#install.packages("brms")
#install.packages("cmdstanr")  # Optional: much faster backend
library("brms")
brms::brm  # will call Stan under the hood
library(dplyr)
library(tidyr)
library(ggplot2)
library(bayesplot)
#smoothing plot
# Load required packages
library(brms)
library(ggplot2)
library(patchwork)
library(emmeans)
library(forcats)
library(loo)          # for loo_compare
library(performance)    
library(dplyr)
library(corrr)
library(ggcorrplot)
library(car)  # for VIF

#install.packages("ggcorrplot")  # only once
library(ggcorrplot)
library(ape)
library(ggplot2)
#install.packages("ggtree")
library(ggtree)
library(ggtreeExtra)
library(rstatix)

#setting working directory:
setwd("/Users/dillanhoyt/Documents/PhD/Meta Analysis/Statistics/Meta Analysis")

meta_ACT <- read.csv("effect_size_totalactivity.csv")
n_distinct(meta_ACT$Paper_ID) #8 papers



# ---------------- Prepare factors (idempotent) ----------------
meta_ACT <- meta_ACT %>%
  mutate(
    Forest_Group        = as.factor(Forest_Group),
    Matrix_mean         = as.factor(Matrix_mean),
    Conservation_status = as.factor(Conservation_status),
    Diet                = as.factor(Diet),
    Strata              = as.factor(Strata),
    Foraging_strat      = as.factor(Foraging_strat),
    Habitat_affinity_new= as.factor(Habitat_affinity_new)
  )


#scaling all numerical predictors: 

meta_ACT <- meta_ACT %>%
  mutate(
    # numeric transforms
    Elevation_s   = as.numeric(scale(Elevation)),
    Rainfall_mm_s = as.numeric(scale(Rainfall_mm)),
    Latitude_s    = as.numeric(scale(Latitude)),
    Longitude_s   = as.numeric(scale(Longitude)),
    FA_length_s   = as.numeric(scale(FA_length)),
    Wing_loading_s  = as.numeric(scale(Wing_loading)),
    Aspect_ratio_s  = as.numeric(scale(Aspect_ratio)),
    Averagemass_g_s  = as.numeric(scale(Averagemass_g)))

#removing things that aren't species: 

# strict but reasonably permissive binomial pattern:
binomial_rgx <- "^[A-Z][A-Za-z'`-]* [a-z][A-Za-z'`-]*$"

# prepare species column (trim whitespace)
meta_ACT <- meta_ACT %>%
  mutate(Species_clean = str_trim(Species))

# logical mask: TRUE when looks like a single binomial
is_binomial <- str_detect(meta_ACT$Species_clean, binomial_rgx)

# subset to only binomials
meta_ACT_binomial <- meta_ACT %>% filter(is_binomial)

unique(meta_ACT_binomial$Species) #58

meta_ACT<- meta_ACT_binomial

#*********** MATCHING PHYLOGENY ***************
# ---- phylogeny mapping for acoustic activity dataset (meta_ACT) ----
library(dplyr)
library(stringr)
library(tibble)
library(ape)

# --- 0. checks & helpers ----------------------------------------------------
if (!exists("meta_ACT")) stop("meta_ACT not found in workspace; load it first.")

# Reuse / re-declare canonicalizers similar to your species-level script
canon <- function(x) {
  x2 <- ifelse(is.na(x), NA_character_, as.character(x))
  x2 <- trimws(x2)
  x2 <- gsub("_", " ", x2)
  x2 <- gsub("[^a-zA-Z ]+", " ", x2)
  x2 <- gsub("\\s+", " ", x2)
  tolower(trimws(x2))
}
clean_species_name <- function(x) {
  y <- as.character(x)
  y <- trimws(y)
  y <- gsub("_ott\\d+$", "", y)
  y <- gsub("_", " ", y)
  y <- iconv(y, from = "", to = "ASCII//TRANSLIT")
  y <- gsub("\\b(sp|spp|cf)\\b\\.?","", y, ignore.case = TRUE)
  y <- gsub("\\(.*?\\)", "", y)
  y <- gsub("[\"'`]", "", y)
  y <- gsub("[^A-Za-z0-9\\s-]", " ", y)
  y <- gsub("\\s+", " ", y)
  tolower(trimws(y))
}

# --- 1. prepare unique species names from meta_ACT --------------------------
species_raw <- meta_ACT %>%
  mutate(Species_orig = as.character(Species)) %>%
  distinct(Species_orig) %>%
  pull(Species_orig) %>%
  as.character()

species_clean <- unique(vapply(species_raw, clean_species_name, character(1), USE.NAMES = FALSE))
cat("Prepared", length(species_clean), "cleaned unique species names from meta_ACT\n")

# --- 2. TNRS matching (requires tnrs_match_names in env) --------------------
if (!exists("tnrs_match_names")) stop("Function tnrs_match_names() not found; load TNRS helper first.")
tnrs_matches <- tnrs_match_names(species_clean)
tnrs_df <- as_tibble(tnrs_matches) %>%
  mutate(search_string = as.character(search_string),
         unique_name   = as.character(unique_name),
         ott_id        = suppressWarnings(as.integer(ott_id)))

bad_rows <- tnrs_df %>% filter(is.na(unique_name) | unique_name == "" | is.na(ott_id))
cat("TNRS: unmatched rows (missing unique_name/ott_id):", nrow(bad_rows), "\n")

# retry pass for problematic names (optional)
still_bad <- character(0)
if (nrow(bad_rows) > 0) {
  problem_names <- unique(bad_rows$search_string)
  double_clean <- function(z) {
    z2 <- gsub("[0-9\\-]+", " ", z)
    z2 <- gsub("\\s+", " ", z2)
    trimws(z2)
  }
  cleaned2 <- unique(vapply(problem_names, double_clean, character(1), USE.NAMES = FALSE))
  cat("Re-running TNRS on", length(cleaned2), "re-cleaned problem names ...\n")
  tnrs_retry <- tnrs_match_names(cleaned2)
  retry_df <- as_tibble(tnrs_retry) %>%
    mutate(search_string = as.character(search_string),
           unique_name   = as.character(unique_name),
           ott_id        = suppressWarnings(as.integer(ott_id)))
  tnrs_df <- tnrs_df %>%
    bind_rows(retry_df) %>%
    distinct(search_string, .keep_all = TRUE)
  bad_rows <- tnrs_df %>% filter(is.na(unique_name) | unique_name == "" | is.na(ott_id))
  still_bad <- unique(bad_rows$search_string)
  cat("After retry, still unmatched:", length(still_bad), "\n")
}

# write unmatched for manual inspection
write.csv(bad_rows, "meta_ACT_tnrs_unmatched.csv", row.names = FALSE)
cat("Wrote unmatched TNRS inputs to meta_ACT_tnrs_unmatched.csv\n")

# --- 3. build ott_lookup ----------------------------------------------------
ott_lookup <- tnrs_df %>%
  mutate(original_name = trimws(tolower(search_string)),
         Species_phylo  = ifelse(!is.na(unique_name) & !is.na(ott_id),
                                 paste0(gsub(" ", "_", unique_name), "_ott", ott_id),
                                 NA_character_)) %>%
  select(original_name, Species_phylo, ott_id) %>%
  distinct(original_name, .keep_all = TRUE)

cat("Unique ott ids extracted:", length(na.omit(unique(ott_lookup$ott_id))), "\n")

# --- 4. build tree & phylogeny matrix (requires tol_induced_subtree) ------
ott_ids <- na.omit(unique(ott_lookup$ott_id))
if (length(ott_ids) == 0) stop("No ott_ids found from TNRS output — inspect tnrs_df.")
if (!exists("tol_induced_subtree")) stop("Function tol_induced_subtree() not found; ensure TOL helper is available.")
cat("Building tree from", length(ott_ids), "ott ids ...\n")
tree <- tol_induced_subtree(ott_ids = ott_ids)
tree_bl <- ape::compute.brlen(tree, method = "Grafen")
phylo_cor_full <- ape::vcv.phylo(tree_bl, corr = TRUE)
cat("phylogeny matrix dimension:", paste(dim(phylo_cor_full), collapse = " x "), "\n")

# --- 5. join ott_lookup back into meta_ACT (clean keys) ---------------------
df <- meta_ACT %>%
  mutate(Species_lower = clean_species_name(as.character(Species))) %>%
  left_join(ott_lookup %>% rename(Species_lower = original_name),
            by = "Species_lower") %>%
  mutate(Species_phylo = ifelse(is.na(Species_phylo), NA_character_, as.character(Species_phylo)),
         ott_id = as.integer(ott_id))

# --- 6. precedence and Species_phylo_use -----------------------------------
df <- df %>%
  mutate(
    .phy_existing = if ("Species_phylo_use" %in% names(.)) as.character(.data$Species_phylo_use) else NA_character_,
    .phy_explicit = if ("Species_phylo" %in% names(.)) as.character(.data$Species_phylo) else NA_character_,
    .phy_ott      = as.character(Species_phylo)
  ) %>%
  mutate(Species_phylo_use = dplyr::coalesce(.phy_existing, .phy_explicit, .phy_ott)) %>%
  select(-starts_with(".phy_"))

# --- 7. canonical fallback: try exact matches vs phylogeny names -------------
phylo_raw  <- rownames(phylo_cor_full)
phylo_base <- gsub("_ott\\d+$", "", phylo_raw)
phylo_df <- tibble(
  Species_phylo       = phylo_raw,
  Species_phylo_base  = phylo_base,
  Species_phylo_canon = canon(phylo_base)
)

missing_idx <- which(is.na(df$Species_phylo_use))
if (length(missing_idx) > 0) {
  cj <- df %>%
    mutate(row_id = dplyr::row_number(),
           Species_canon = canon(as.character(Species))) %>%
    left_join(phylo_df %>% select(Species_phylo, Species_phylo_canon),
              by = c("Species_canon" = "Species_phylo_canon"))
  df$Species_phylo_use[missing_idx] <- cj$Species_phylo[match(missing_idx, cj$row_id)]
}

# --- 8. genus-level proxy for still-missing --------------------------------
still_missing_rows <- df %>% filter(is.na(Species_phylo_use) | Species_phylo_use == "")
if (nrow(still_missing_rows) > 0) {
  cat("Attempting genus-level proxy for", nrow(still_missing_rows), "rows\n")
  phylo_base_vec <- phylo_base
  phylo_genus <- tolower(sub("_.*", "", phylo_base_vec))
  phylo_lookup <- tibble(phylo_label = phylo_raw, phylo_base = phylo_base_vec, phylo_genus = phylo_genus)
  df <- df %>%
    mutate(
      Species_clean = gsub('["“”]', '', as.character(Species)),
      genus = tolower(sub(" .*", "", Species_clean))
    )
  genus_proxy <- phylo_lookup %>%
    group_by(phylo_genus) %>%
    slice(1) %>%
    ungroup() %>%
    select(phylo_genus, proxy_phylo = phylo_label)
  df <- df %>%
    left_join(genus_proxy, by = c("genus" = "phylo_genus")) %>%
    mutate(
      Species_phylo_use = dplyr::coalesce(as.character(Species_phylo_use), as.character(proxy_phylo))
    ) %>%
    select(-proxy_phylo)
}

# --- 9. final prune: keep only rows with Species_phylo_use present in tree ----
present <- rownames(phylo_cor_full)
df_final <- df %>%
  filter(!is.na(Species_phylo_use) & Species_phylo_use %in% present) %>%
  mutate(Species = factor(Species_phylo_use, levels = unique(Species_phylo_use))) %>%
  droplevels()

# pruned phylogeny correlation matrix
species_order <- levels(df_final$Species)
phylo_cor <- phylo_cor_full[species_order, species_order, drop = FALSE]

# --- 10. diagnostics & outputs ----------------------------------------------
cat("Rows before mapping:", nrow(meta_ACT), "\n")
cat("Rows after mapping/prune:", nrow(df_final), "\n")
cat("Distinct species (phylo labels) kept:", length(unique(as.character(df_final$Species))), "\n")
cat("phylo_cor dim:", dim(phylo_cor)[1], "x", dim(phylo_cor)[2], "\n")

# write outputs for inspection
meta_ACT_phylo <- df_final
write.csv(meta_ACT_phylo, "meta_ACT_phylo_mapped.csv", row.names = FALSE)
cat("Wrote meta_ACT_phylo_mapped.csv (relabelled + pruned meta_ACT)\n")

# Save phylo_cor as RDS for modelling
saveRDS(phylo_cor, file = "phylo_cor_meta_ACT.rds")
cat("Saved pruned phylogenetic correlation matrix to phylo_cor_meta_ACT.rds\n")

# Save list of still-bad TNRS inputs for manual mapping
if (exists("still_bad") && length(still_bad) > 0) {
  write.csv(data.frame(bad = still_bad), "meta_ACT_tnrs_still_bad.csv", row.names = FALSE)
  cat("Wrote still-bad TNRS inputs to meta_ACT_tnrs_still_bad.csv\n")
}

#**********************************************
meta_ACT <- meta_ACT_phylo

######### adding weight *****************
vi_floor <- 0.2

meta_ACT <- meta_ACT %>%
  mutate(
    # focal replication
    nt = as.numeric(Sample_size_sum),
    
    # partner (continuous pooled) replication
    nc = as.numeric(Cont_Sample_size_sum),
    
    # safe fallbacks
    nt_fix = ifelse(is.na(nt) | nt <= 0, 0.5, nt),
    nc_fix = ifelse(is.na(nc) | nc <= 0, 0.5, nc),
    
    # sampling variance
    vi_raw = 1/nt_fix + 1/nc_fix,
    vi     = pmax(vi_raw, vi_floor),
    
    # measurement error SD
    se_obs = sqrt(vi)
  )

#*******************

#Habitat_affinity_new 
meta_ACT %>%
  filter(!is.na(Habitat_affinity_new) & Habitat_affinity_new != "") %>%
  distinct(Species) %>%
  nrow() #37

#which are the species with missing: Habitat_affinity_new

meta_ACT %>%
  filter(is.na(Habitat_affinity_new) | Habitat_affinity_new == "") %>%
  distinct(Species) %>%
  arrange(Species)

#                             Species
# 1    Hesperoptenus_tickelli_ott445742
# 2       Hipposideros_pomona_ott905428
# 3  Miniopterus_fuliginosus_ott1011476
# 4       Miniopterus_pusillus_ott61867
# 5   Pipistrellus_ceylonicus_ott743848
# 6      Rhinolophus_beddomei_ott994640
# 7   Rhinolophus_indorouxii_ott7067837
# 8       Rhinolophus_lepidus_ott217411
# 9       Rhinolophus_rouxii_ott1047994
# 10     Noctilio_albiventris_ott604407
# 11       Noctilio_leporinus_ott604404
# 12    Peropteryx_trinitatis_ott432533
# 13    Saccopteryx_canescens_ott232381
# 14     Perimyotis_subflavus_ott977400
# 15       Emballonura_atrata_ott156176
# 16  Hipposideros_commersoni_ott156174
# 17       Scotomanes_ornatus_ott759645

unique(meta_ACT$Habitat_affinity_new)

meta_ACT <- meta_ACT %>%
  mutate(
    Habitat_affinity_new = case_when(
      
      Species %in% c("Noctilio_leporinus_ott604404",
                     "Noctilio_albiventris_ott604407") ~ "Water-adapted",
      TRUE ~ Habitat_affinity_new
    )
  )

meta_ACT <- meta_ACT %>%
  mutate(
    Habitat_affinity_new = case_when(
      
      Species %in% c("Rhinolophus_beddomei_ott994640",
                     "Rhinolophus_indorouxii_ott7067837",
                     "Rhinolophus_lepidus_ott217411",
                     "Rhinolophus_rouxii_ott1047994",
                     "Hipposideros_pomona_ott905428",
                     "Hipposideros_commersoni_ott156174",
                     "Scotomanes_ornatus_ott759645") ~ "Clutter-adapted",
      TRUE ~ Habitat_affinity_new
    )
  )
      
meta_ACT <- meta_ACT %>%
  mutate(
    Habitat_affinity_new = case_when(
      
      Species %in% c("Miniopterus_fuliginosus_ott1011476",
                     "Miniopterus_pusillus_ott61867") ~ "Open-adapted",
      TRUE ~ Habitat_affinity_new
    )
  )

meta_ACT <- meta_ACT %>%
  mutate(
    Habitat_affinity_new = case_when(
      
      Species %in% c("Pipistrellus_ceylonicus_ott743848",
                     "Hesperoptenus_tickelli_ott445742",
                     "Peropteryx_trinitatis_ott432533",
                     "Saccopteryx_canescens_ott232381",
                     "Perimyotis_subflavus_ott977400",
                     "Emballonura_atrata_ott156176") ~ "Edge-adapted",
      TRUE ~ Habitat_affinity_new
    )
  )

meta_ACT %>%
  filter(!is.na(Habitat_affinity_new) & Habitat_affinity_new != "") %>%
  distinct(Species) %>%
  nrow() #54

#********* call type ********
##call type
meta_ACT %>%
  filter(!is.na(Call_type) & Call_type != "") %>%
  distinct(Species) %>%
  nrow() #29



unique(meta_ACT$Call_type)

meta_ACT <- meta_ACT %>%
  mutate(
    Call_type_clean = case_when(
      is.na(Call_type) ~ NA_character_,
      
      # CF-dominated (any string starting with CF)
      str_detect(Call_type, "^CF") ~ "CF",
      
      # FM.QCF, FM_QCF, FM-QCF -> treat as FM-dominated
      str_detect(Call_type, "^FM") ~ "FM",
      
      # Pure QCF
      Call_type == "QCF" ~ "QCF",
      
      TRUE ~ NA_character_
    )
  )

# Check result
table(meta_ACT$Call_type, meta_ACT$Call_type_clean, useNA = "ifany")

unique(meta_ACT$Call_type_clean)

meta_ACT %>%
  filter(is.na(Call_type_clean) | Call_type_clean == "") %>%
  distinct(Species) %>%
  arrange(Species)


meta_ACT <- meta_ACT %>%
  mutate(
    Call_type_clean = case_when(
      
      # ---------------- CF specialists ----------------
      Species %in% c(
        "Rhinolophus_beddomei_ott994640",
        "Rhinolophus_indorouxii_ott7067837",
        "Rhinolophus_lepidus_ott217411",
        "Rhinolophus_rouxii_ott1047994",
        "Hipposideros_commersoni_ott156174",
        "Hipposideros_pomona_ott905428"
      ) ~ "CF",
      
      # ---------------- QCF (narrowband open-space molossid style) ----------------
      Species %in% c(
        "Molossus_currentium_ott3614007",
        "Molossus_rufus_ott267993",
        "Promops_centralis_ott301613",
        "Nyctinomops_macrotis_ott14996",        # if present in your dataset
        "Nyctinomops_laticaudatus_ott14995"     # optional if present
      ) ~ "QCF",
      
      # ---------------- FM-dominated (broadband vespertilionid / emballonurid) ----------------
      Species %in% c(
        "Hesperoptenus_tickelli_ott445742",
        "Miniopterus_fuliginosus_ott1011476",
        "Miniopterus_pusillus_ott61867",
        "Miniopterus_manavi_ott18883",
        "Pipistrellus_ceylonicus_ott743848",
        "Peropteryx_trinitatis_ott432533",
        "Saccopteryx_canescens_ott232381",
        "Emballonura_atrata_ott156176",
        "Perimyotis_subflavus_ott977400",
        "Myotis_albescens_ott353889",
        "Myotis_nigricans_ott307133",
        "Myotis_goudoti_ott31985",
        "Eptesicus_fuscus_ott10737",
        "Lasionycteris_noctivagans_ott401283",
        "Lasiurus_borealis_ott61860",
        "Aeorestes_cinereus_ott369537",
        "Cyttarops_alecto_ott599540",
        "Diclidurus_albus_ott1060477",
        "Dasypterus_ega_ott635119"
      ) ~ "FM",
      
      TRUE ~ Call_type_clean
    )
  )



#Foraging_strat
meta_ACT %>%
  filter(!is.na(Foraging_strat) & Foraging_strat != "") %>%
  distinct(Species) %>%
  nrow() #21


# looking at just F-C contrasts: 



meta_ACT_FC <- meta_ACT %>%
  # keep only fragment–continuous contrasts
  filter(Comparison == "Frag_vs_Cont") %>%
  
  # ensure lnRR is numeric
  mutate(lnRR = as.numeric(lnRR)) %>%
  
  # keep only valid effect sizes and variances
  filter(
    !is.na(lnRR), is.finite(lnRR),
    !is.na(vi),   is.finite(vi), vi > 0
  ) %>%
  
  # ensure modelling IDs are factors
  mutate(
    Paper_ID = factor(Paper_ID),
    Species  = factor(Species)
  ) %>%
  droplevels()


#how many species wihtin this contrast? 
n_distinct(meta_ACT_FC$Species) #41 species! 

#how many papers? 
n_distinct(meta_ACT_FC$Paper_ID) #5 papers.. hmm 

#total bats 

#before running models
options(mc.cores = parallel::detectCores())


# Log response ratio model of species richness vs fragment size 
# Your species richness model in brms

#global model: 
colnames(meta_ACT_FC)
meta_ACT_FC$Forest_Group <- as.factor(meta_ACT_FC$Forest_Group)
meta_ACT_FC$Matrix_mean <- as.factor(meta_ACT_FC$Matrix_mean)


sapply(meta_ACT_FC[, c("Matrix_mean", "Forest_Group")], nlevels)
str(meta_ACT_FC)

unique(meta_ACT_FC$Paper_ID) #5 papers included. 

#look at variable names within the column Matrix_mean
unique(meta_ACT_FC$Matrix_mean)

#41 species in total 
#Aspect_ratio
meta_ACT_FC %>%
  filter(!is.na(Aspect_ratio) & Aspect_ratio != "") %>%
  distinct(Species) %>%
  nrow() #24 Species 

#Wing_loading
meta_ACT_FC %>%
  filter(!is.na(Wing_loading) & Wing_loading != "") %>%
  distinct(Species) %>%
  nrow() #22

#Averagemass_g 
meta_ACT_FC %>%
  filter(!is.na(Averagemass_g) & Averagemass_g != "") %>%
  distinct(Species) %>%
  nrow() #38

#FA 
meta_ACT_FC %>%
  filter(!is.na(FA_length_s) & FA_length_s != "") %>%
  distinct(Species) %>%
  nrow() #41


colnames(meta_ACT_FC) # Call_type #Foraging_strategy


#trying to figure out how to maximise species to include in model with variable limitations: 
#creating function to count retained species for any variable set: 
count_species_for_vars <- function(df, vars) {
  
  required <- c("lnRR", "se_obs", vars)
  
  df_complete <- df %>%
    filter(if_all(all_of(required), ~ !is.na(.x)))
  
  tibble(
    n_rows = nrow(df_complete),
    n_species = n_distinct(df_complete$Species),
    n_papers = n_distinct(df_complete$Paper_ID)
  )
}

vars_current <- c(
  "Averagemass_g_s",
  "Aspect_ratio_s",
  "Conservation_status",
  "Habitat_affinity_new",
  "Call_type_clean"
  
)

count_species_for_vars(meta_ACT_FC, vars_current) #20 species and 4 papers 

#trying alternative combinations: 
count_species_for_vars(meta_ACT_FC,
                       c("Averagemass_g_s",
                         "Call_type_clean", 
                         "Habitat_affinity_new")) #38 species 

#**************

#running a model for traits: 
m_batTraits_species.gaussian_CF_ACT  <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~
    Averagemass_g_s  +
    Habitat_affinity_new + Call_type_clean +
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_FC,
  data2 = list(phylo_cor = phylo_cor),
  family = gaussian(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  seed = 2025
)

summary(m_batTraits_species.gaussian_CF_ACT); pp_check(m_batTraits_species.gaussian_CF_ACT)
saveRDS(m_batTraits_species.gaussian_CF_ACT, file = "m_batTraits_species.gaussian_CF_ACT.rds")

#saving: 

#checking model 
pp_check(m_batTraits_species.gaussian_CF_ACT, type = "hist")
pp_check(m_batTraits_species.gaussian_CF_ACT, type = "scatter_avg")


mcmc_intervals(as_draws_df(m_batTraits_species.gaussian_CF_ACT), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")

#running a student model :
m_batTraits_species.student_CF_ACT  <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~
    Averagemass_g_s  +
    Habitat_affinity_new + Call_type_clean +
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_FC,
  data2 = list(phylo_cor = phylo_cor),
  family = student(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 2025
)


#saving: 
summary(m_batTraits_species.student_CF_ACT); pp_check(m_batTraits_species.student_CF_ACT)
saveRDS(m_batTraits_species.student_CF_ACT, file = "m_batTraits_species.student_CF_ACT.rds")


mcmc_intervals(as_draws_df(m_batTraits_species.student_CF_ACT), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")

plot(conditional_effects(m_batTraits_species.student_CF_ACT, effects = "Habitat_affinity_new"), ask = FALSE)
plot(conditional_effects(m_batTraits_species.student_CF_ACT, effects = "Call_type_clean"), ask = FALSE)

#loo to compare gaussian vs student: 
loo_trait_gaussian <- loo(m_batTraits_species.gaussian_CF_ACT)
loo_trait_student <- loo(m_batTraits_species.student_CF_ACT)
loo_compare(loo_trait_gaussian, loo_trait_student) #student is significantly better. continue with student


#landscape (there is just one kind of forest group, so I will not include this within the model)

m_species.student_CF_ACT_landscape  <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~
    Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean  +
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_FC,
  data2 = list(phylo_cor = phylo_cor),
  family = student(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 2025
)

#saving: 
summary(m_species.student_CF_ACT_landscape); pp_check(m_species.student_CF_ACT_landscape)
saveRDS(m_species.student_CF_ACT_landscape, file = "m_species.student_CF_ACT_landscape.rds")

mcmc_intervals(as_draws_df(m_species.student_CF_ACT_landscape), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")


#************************* fragment fragment contrasts *************************


meta_ACT_FF <- meta_ACT %>%
  filter(Comparison == "Frag_vs_Frag") %>%
  mutate(
    size_small = Fragment_size_ha,   # the smaller of the pair (from your construction)
    log_size_small = log(size_small)
  ) %>%
  filter(is.finite(log_size_small)) #creating a logged number that is the smallest fragment of the two being compared, because the size of the smaller fragment is the variable that changes between comparisons — the large fragment just acts as a reference.

#scaling log_size_small: 
meta_ACT_FF <- meta_ACT_FF %>%
  mutate(log_size_small_s = as.numeric(scale(log_size_small)))



m_batTraits_species.student_FF_ACT <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~
    Averagemass_g_s  + Call_type_clean +
    Habitat_affinity_new + 
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_FF,
  data2 = list(phylo_cor = phylo_cor),
  family = student(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.99, max_treedepth = 12)
)

saveRDS(m_batTraits_species.student_FF_ACT, file = "m_batTraits_species.student_FF_ACT.rds")

summary(m_batTraits_species.student_FF_ACT); pp_check(m_batTraits_species.student_FF_ACT)

mcmc_intervals(as_draws_df(m_batTraits_species.student_FF_ACT), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")

#only 5 papers, 39 species, 4167 observations
#we're seeing similar relationships with aerial insectivores, whereby larger bodied bats do worse in small fragments

#************** combining them together ************


# control: set behaviour for continuous forest replacement if you want other behaviour
force_cf <- FALSE   # or TRUE, or "missing_or_zero"

meta_ACT_All <- meta_ACT %>%
  # 1. canonicalise contrasts & flip lnRR for Cont_vs_Frag
  mutate(
    .orig_comp = Comparison,
    lnRR = if_else(.orig_comp == "Cont_vs_Frag", -lnRR, lnRR),
    Comparison = case_when(
      .orig_comp %in% c("Frag_vs_Frag", "Frag_vs_Cont") ~ .orig_comp,
      .orig_comp == "Cont_vs_Frag"                       ~ "Frag_vs_Cont",
      TRUE                                               ~ NA_character_
    ),
    Fragment_size_ha         = as.numeric(Fragment_size_ha),
    Partner_Fragment_size_ha = as.numeric(Partner_Fragment_size_ha)
  ) %>%
  
  # keep only recognised contrasts
  filter(!is.na(Comparison)) %>%
  
  # 2. For Frag_vs_Cont rows where the recorded fragment (numerator) has NO size,
  #    treat the continuous forest as numerator: flip lnRR and assign CF = 10000.
  mutate(
    frag_missing = (Comparison == "Frag_vs_Cont" & Treatment == "Frag" & is.na(Fragment_size_ha)),
    
    # flip lnRR and Treatment where fragment side is missing (so CF becomes numerator)
    lnRR = if_else(frag_missing, -lnRR, lnRR),
    Treatment = if_else(frag_missing, "Cont", Treatment),
    
    # set Partner_Fragment_size_ha to 10000 for the continuous partner that will be used as numerator
    Partner_Fragment_size_ha = case_when(
      frag_missing ~ 10000,
      force_cf == TRUE & Comparison == "Frag_vs_Cont" & Treatment == "Cont" ~ 10000,
      identical(force_cf, "missing_or_zero") & Comparison == "Frag_vs_Cont" & Treatment == "Cont" & (is.na(Partner_Fragment_size_ha) | Partner_Fragment_size_ha == 0) ~ 10000,
      TRUE ~ Partner_Fragment_size_ha
    )
  ) %>%
  
  # 3. compute numerator area depending on canonicalised Comparison/Treatment
  mutate(
    area_num_ha = case_when(
      Comparison == "Frag_vs_Frag"                            ~ Fragment_size_ha,                      # smaller fragment stored in Fragment_size_ha
      Comparison == "Frag_vs_Cont" & Treatment == "Frag"      ~ Fragment_size_ha,                      # fragment is numerator
      Comparison == "Frag_vs_Cont" & Treatment == "Cont"      ~ Partner_Fragment_size_ha,              # partner (CF) is numerator (possibly set to 10000 above)
      TRUE                                                   ~ NA_real_
    )
  ) %>%
  
  # 4. log-transform and scale
  mutate(
    log_area_num   = as.numeric(log(area_num_ha)),
    log_area_num_s = as.numeric(scale(log_area_num))
  ) %>%
  
  # 5. drop rows w/ missing area and remove helpers
  filter(is.finite(log_area_num)) %>%
  select(-.orig_comp, -frag_missing)

length(unique(meta_ACT_All$Species))   # should be 54 species 

#quick checks to see if it worked: 
table(meta_ACT_All$Comparison)

# Frag_vs_Cont Frag_vs_Frag 
# 437         4437

summary(meta_ACT_All$area_num_ha)

sum(meta_ACT_All$area_num_ha == 10000, na.rm = TRUE) #41

#making sure that the smaller frag is the numerator
meta_ACT_All <- meta_ACT_All %>%
  mutate(
    # ensure for Frag_vs_Frag the numerator area is the smaller of the two recorded fragment sizes
    area_num_ha = case_when(
      Comparison == "Frag_vs_Frag" ~ pmin(Fragment_size_ha, Partner_Fragment_size_ha, na.rm = TRUE),
      TRUE                          ~ area_num_ha
    ),
    log_area_num   = as.numeric(log(area_num_ha)),
    log_area_num_s = as.numeric(scale(log_area_num))
  ) %>%
  filter(is.finite(log_area_num))


# 1) Create s_obs (prefer existing se_obs, otherwise compute from vi)
meta_ACT_All <- meta_ACT_All %>%
  mutate(
    se_obs = if ("se_obs" %in% names(.)) se_obs else NA_real_,   # keep if exists
    vi     = if ("vi"     %in% names(.)) vi     else NA_real_,   # keep if exists
    # create s_obs: prefer se_obs, otherwise sqrt(vi), otherwise NA
    s_obs = case_when(
      !is.na(se_obs) ~ as.numeric(se_obs),
      !is.na(vi)     ~ as.numeric(sqrt(vi)),
      TRUE           ~ NA_real_
    )
  )


#just a summary fo rthe paper in terms of how many weights I actually used: 
# 2) How many variances are exactly at the floor?
vi_floor <- 0.2
sum(meta_ACT_All$vi == vi_floor, na.rm = TRUE)
round(100 * sum(meta_ACT_All$vi == vi_floor, na.rm = TRUE) / nrow(meta_ACT_All), 3)


# 2) Make sure Comparison is a factor and relevel to Frag_vs_Cont
meta_ACT_All$Comparison <- factor(meta_ACT_All$Comparison, levels = c("Frag_vs_Cont", "Frag_vs_Frag"))
meta_ACT_All$Comparison <- relevel(meta_ACT_All$Comparison, ref = "Frag_vs_Cont")


meta_ACT_All %>%
  distinct(Species, Family) %>%   # count each species once
  count(Family, name = "n_species") %>%
  arrange(desc(n_species))

m_batTraits_species.student_ACT_ALL <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~ Comparison +
    Averagemass_g_s  +
    Habitat_affinity_new + Call_type_clean  +
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_All,
  data2 = list(phylo_cor = phylo_cor),
  family = student(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

summary(m_batTraits_species.student_ACT_ALL); pp_check(m_batTraits_species.student_ACT_ALL)
#8 papers, 4601 observations, 50 species.

mcmc_intervals(as_draws_df(m_batTraits_species.student_ACT_ALL), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")

saveRDS(m_batTraits_species.student_ACT_ALL, file = "m_batTraits_species.student_ACT_ALL.rds")

#and looking at landscpae variables in the same context: 


m_species.student_ACT_ALL_landscape <- brm(
  lnRR | se(se_obs, sigma = TRUE) ~ Comparison +
    Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean  +
    (1 | Paper_ID) + (1 | gr(Species, cov = phylo_cor)),
  data  = meta_ACT_All,
  data2 = list(phylo_cor = phylo_cor),
  family = student(),
  prior  = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 4000, warmup = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

summary(m_species.student_ACT_ALL_landscape); pp_check(m_species.student_ACT_ALL_landscape)
#54 species, 8 papers and 4874 observations
saveRDS(m_species.student_ACT_ALL_landscape, file = "m_species.student_ACT_ALL_landscape.rds")


mcmc_intervals(as_draws_df(m_species.student_ACT_ALL_landscape), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")















#landscape predictors between continuous forest and fragments 
#this is looking at forest size: 

#gaussian or student?
m_ACT_landscape_CF <- brm(
  lnRR | se(sqrt(vi), sigma = TRUE) ~  Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean + Forest_Group + (1 | Paper_ID),
  data = act_total_FC,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4, iter = 2000, warmup = 1000
)

# count divergences
library(rstan)
sampler_params <- rstan::get_sampler_params(m_ACT_landscape_CF$fit, inc_warmup = FALSE)
sum(sapply(sampler_params, function(x) sum(x[,"divergent__"])))

# how many Paper_ID levels?
length(unique(act_total_FC$Paper_ID))

# which rows were dropped because of NAs in model vars?
vars <- all.vars(formula(m_ACT_landscape_CF))
nrow(act_total_FC); nrow(model.frame(m_ACT_landscape_CF))
colSums(is.na(act_total_FC[, vars]))
#there are very few papers for paper id to be a random effect, so I am using it as a fixed effect instead: 

m_ACT_landscape_CF.fixed <- brm(
  lnRR | se(sqrt(vi), sigma = TRUE) ~ 
    Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean + Forest_Group + Paper_ID,
  data = act_total_FC,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sigma")
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 4000, warmup = 2000
)
saveRDS(m_ACT_landscape_CF.fixed, file = "m_ACT_landscape_CF.fixed.rds")

summary(m_ACT_landscape_CF.fixed)

#checking model 
pp_check(m_ACT_landscape_CF.fixed, type = "hist")
pp_check(m_ACT_landscape_CF.fixed, type = "scatter_avg")

mcmc_intervals(as_draws_df(m_ACT_landscape_CF.fixed), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")

#student: 

m_ACT_landscape_CF_student <- brm(
  lnRR | se(sqrt(vi), sigma = TRUE) ~ 
    Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean + Forest_Group + Paper_ID,
  data = act_total_FC,
  family = student(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sigma")
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 4000, warmup = 2000
)

#checking model 
pp_check(m_ACT_landscape_CF_student, type = "hist")
pp_check(m_ACT_landscape_CF_student, type = "scatter_avg")

mcmc_intervals(as_draws_df(m_ACT_landscape_CF_student), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")


loo_act_gaussian <- loo(m_ACT_landscape_CF)
loo_act_student <- loo(m_ACT_landscape_CF_student)
loo_compare(loo_act_gaussian, loo_act_student) #student is worse. stick with gaussian. 

#**************************************************************************************
#**************************************************************************************
#Fragment vs fragment

#fragment vs fragment: 
act_total_FF <- act_total %>%
  filter(Comparison == "Frag_vs_Frag") %>%
  mutate(
    size_small = Fragment_size_ha,   # the smaller of the pair (from your construction)
    log_size_small = log(size_small)
  ) %>%
  filter(is.finite(log_size_small)) #creating a logged number that is the smallest fragment of the two being compared, because the size of the smaller fragment is the variable that changes between comparisons — the large fragment just acts as a reference.

#scaling log_size_small: 
act_total_FF <- act_total_FF %>%
  mutate(log_size_small_s = as.numeric(scale(log_size_small)))

#I want to exclude Papers 34 (which has fragments being >9800 ha),
act_total_FF <-subset(act_total_FF,Paper_ID!= 34)

model_ACT_FF <-  brm(
  lnRR | se(sqrt(vi), sigma = TRUE) ~ log_size_small_s +
    Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean + Forest_Group + Paper_ID,
  data = act_total_FF,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sigma")
  ),
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  chains = 4, iter = 4000, warmup = 2000)
)
saveRDS(model_ACT_FF, file = "m_activty_FF.rds")

summary(model_ACT_FF)
  
#checking model 
pp_check(m_ACT_landscape_CF.fixed, type = "hist")
pp_check(m_ACT_landscape_CF.fixed, type = "scatter_avg")
  
mcmc_intervals(as_draws_df(m_ACT_landscape_CF.fixed), pars = vars(starts_with("b_")),
                 prob = 0.95) +
    vline_at(0, linetype = "dashed")

#*********
#*#all combined 
#*


#applying a 10,000 ha cap on continuous forests
ACT_ALL <- act_total %>%
  mutate(
    # area of the site in the NUMERATOR of lnRR
    area_num_ha = case_when(
      Comparison == "Frag_vs_Frag" ~ Fragment_size_ha,
      Comparison == "Frag_vs_Cont" & Treatment == "Frag" ~ Fragment_size_ha,
      Comparison == "Cont_vs_Frag" & Treatment == "Cont" ~ 10000,
      TRUE ~ NA_real_
    ),
    
    # cap values above 10,000 ha
    area_num_ha = if_else(area_num_ha > 10000, 10000, area_num_ha),
    
    # log-transform + scale
    log_area_num   = log(area_num_ha),
    log_area_num_s = as.numeric(scale(log_area_num))
  ) %>%
  filter(is.finite(log_area_num))

#I want to exclude Papers 34 (which has fragments being >9800 ha),
ACT_ALL <-subset(ACT_ALL,Paper_ID!= 34)

#how many papers are in this dataset?
n_distinct(ACT_ALL$Paper_ID)

m_activity_all <-brm(
  lnRR | se(sqrt(vi), sigma = TRUE) ~  Comparison + log_area_num_s
  + Elevation_s + Rainfall_mm_s + Latitude_s + Longitude_s +
    Matrix_mean + Forest_Group + (1 | Paper_ID),
  data = ACT_ALL,
  family = gaussian(),
  prior = c(
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  ),
  chains = 4, cores = 4,
  iter = 3000, warmup = 1500,
  control = list(
    adapt_delta = 0.99,         # ← KEY FIX
    max_treedepth = 15          # optional but helps
  ),
  seed = 2025
)

summary(m_activity_all)
saveRDS(m_activity_all, file = "m_activty_ALL.rds")


#checking model 
pp_check(m_activity_all, type = "hist")
pp_check(m_activity_all, type = "scatter_avg")

mcmc_intervals(as_draws_df(m_activity_all), pars = vars(starts_with("b_")),
               prob = 0.95) +
  vline_at(0, linetype = "dashed")


























