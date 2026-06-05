# Global meta-analysis of bats' responses to forest fragmentation

This repository contains the data, fitted Bayesian models, and analysis scripts used in a global meta-analysis examining how bats respond to forest fragmentation across tropical and temperate regions.

Analyses are organised into four major components:

- Species richness (assemblage-level responses)
- Total abundance (assemblage-level responses)
- Species-level abundance (capture-based data)
- Species activity (acoustic data)

Within each component, analyses are separated into:

- Fragment–continuous forest comparisons (F–C), comparing fragmented forests with nearby continuous forest.
- Fragment–fragment comparisons (F–F), comparing fragments of different sizes.
- Combined models, which incorporate both comparison types and include comparison type as a fixed effect.

The repository includes all data files, fitted model objects (.rds), and Quarto documents used to reproduce the figures and analyses presented in the manuscript.

**Model information: **
## Model files

This repository contains fitted Bayesian hierarchical models used in the manuscript analyses. Models are grouped by response variable and comparison type.

### Species richness

#### Fragment–continuous forest (F-C)

- m_sr_all_new.rds — primary species richness model.

#### Fragment–fragment forest (F–F)

- m_SR_FF_new.rds — primary fragment–fragment species richness model.

#### Combined models

- m_sr_all_randomeffect_andinteraction_new.rds — species richness model including landscape covariates, random effects, and interaction terms.

#### Sensitivity analyses

- model_sr_full_gaussian.rds — Gaussian error model.
- model_sr_full_gaussian_hlfeps.rds — half-epsilon continuity-correction sensitivity analysis.
- model_sr_full_gaussian_dbleps.rds — double-epsilon continuity-correction sensitivity analysis.

---

### Total abundance

#### Fragment–continuous forest (C–F)

- m_abundance_CF_student_new.rds — primary abundance model (Student-t error).
- m_abundance_CF_gaussian_new.rds — Gaussian sensitivity model.

#### Fragment–fragment forest (F–F)

- m_abundance_FF_new.rds — primary fragment–fragment abundance model.

#### Combined models

- m_abundance_ALL_new.rds — combined abundance model.
- m_abundance_interaction_ALL_new.rds — abundance model including landscape interactions.
- m_abundance_ALL_FragSize.rds — abundance model including fragment-size effects.

#### Sensitivity analyses

- m_abundance_ALL_6K.rds — reduced-draw sensitivity analysis.
- m_abundance_ALL_20K.rds — increased-draw sensitivity analysis.
- model_AB_full_doubleeps.rds — double-epsilon continuity-correction sensitivity analysis.

---

### Activity

#### Fragment–continuous forest (C–F)

- m_batTraits_species.student_CF_ACT.rds — primary activity model.
- m_batTraits_species.gaussian_CF_ACT.rds — Gaussian sensitivity model.

#### Fragment–fragment forest (F–F)

- m_batTraits_species.student_FF_ACT.rds — primary fragment–fragment activity model.

#### Combined models

- m_batTraits_species.student_ACT_ALL.rds — combined activity model.


### Species-level responses

#### Fragment–continuous forest (C–F)

- m_batTraits_species.student_CF.rds — primary species-level model used in manuscript figures.
- m_batTraits_species.gaussian_CF.rds — Gaussian sensitivity model.
- m_batTraits_species.forag_habitat_FC.rds — foraging-habitat model.
- m_batTraits_species.strata_wl_FC.rds — wing-loading and forest-strata model.
- m_batTraits_species.student_CF_AspectRatio_Strata.rds — aspect-ratio and strata model.
- m_batTraits_species.student_CF_bodymassDiet_new.rds — body-mass and diet interaction model.
- m_batTraits_species.student_CF_StrataDiet.rds — strata and diet interaction model.

#### Fragment–fragment forest (F–F)

- m_batTraits_species.student_FF_new.rds — primary fragment–fragment species-level model.

#### Combined models

- m_batTraits_species.student_ALL_new_NEW.rds — combined species-level model.

#### Sensitivity analyses

- m_batTraits_species.uncorr1_student.rds —     species-level model with phylogenetic random effects and comparison-specific species variation.
- m_batTraits_species.gaussian_CF_halfeps.rds — half-epsilon continuity-correction sensitivity analysis.
- m_batTraits_species.gaussian_CF_doubleeps.rds — double-epsilon continuity-correction sensitivity analysis.

---

### Notes

Primary manuscript figures were generated from the models referenced within the corresponding Quarto documents. Additional model files are provided to facilitate reproducibility of sensitivity analyses, alternative model structures, and supplementary analyses.

