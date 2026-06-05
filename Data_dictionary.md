## Data dictionary

### Study and site identifiers

| Column | Description |
|----------|-------------|
| SiteUID | Unique identifier for a sampling site within a study. |
| Paper_ID | Unique identifier for the source publication. |
| SiteUID_new | Updated site identifier generated during data harmonisation. |
| UniqueFragID | Unique identifier for a forest fragment. |
| UniqueFragID_site | Fragment identifier associated with a site-level observation. |
| Unique_name | Original fragment/site name reported in the publication. |

---

### Treatment and comparison information

| Column | Description |
|----------|-------------|
| Treatment | Habitat category used in the comparison (e.g., fragment or continuous forest). |
| Treatment_site | Site-level treatment category. |
| Comparison | Comparison type used in the meta-analysis (Fragment–Continuous or Fragment–Fragment). |
| Type | Response category used in the analysis. |
| Partner_Paper_ID | Publication associated with the comparison partner site. |
| Partner_SiteUID | Site identifier of the comparison partner. |
| Partner_UniqueFragID | Fragment identifier of the comparison partner. |

---

### Fragment characteristics

| Column | Description |
|----------|-------------|
| Fragment_size_ha | Area of the focal forest fragment (ha). |
| Fragment_size_ha_site | Fragment area associated with a site-level observation. |
| Cont_Fragment_size_ha | Area assigned to the continuous-forest comparison site. |
| Partner_Fragment_size_ha | Area of the comparison partner fragment. |
| forest_area_ref_ha | Reference forest area used in analyses. |
| log_forest_area_ref | Natural log-transformed reference forest area. |
| log_forest_area_ref_s | Standardised log-transformed reference forest area. |
| log_size_small | Log-transformed area of the smaller fragment in fragment–fragment comparisons. |
| log_size_small_s | Standardised smaller-fragment area. |
| Fragment_size_ha_s | Standardised fragment area. |

---

### Response variables

| Column | Description |
|----------|-------------|
| Response | Original response value reported in the study. |
| Response_site | Site-level response value. |
| site_value | Site-specific value used for comparison calculations. |
| Partner_site_value | Site value of the comparison partner. |
| Cont_Value | Response value for the continuous-forest comparison site. |
| Value | Final response value used in the meta-analysis. |
| lnRR | Log response ratio effect size. |
| lnRR_input | Raw values used to calculate lnRR. |
| lnRR_for_plot | Effect size used for plotting. |
| lnRR_plot | Plotting version of lnRR after processing. |

---

### Sampling effort

| Column | Description |
|----------|-------------|
| Mist_net_hours | Total mist-net effort (net-hours). |
| Mist_net_hours_site | Site-level mist-net effort. |
| Sample_size | Number of sampling units. |
| Sample_size_sum | Total sample size for the focal site. |
| Cont_Sample_size_sum | Total sample size for the continuous-forest comparison site. |
| Partner_Sample_size_sum | Total sample size for the comparison partner. |
| Rate | Response standardised by sampling effort. |

---

### Continuous forest matching variables

| Column | Description |
|----------|-------------|
| Cont_Paper_ID | Publication identifier for the matched continuous-forest site. |
| Cont_SiteUID | Original identifier for the matched continuous-forest site. |
| Cont_SiteUID_new | Updated identifier for the matched continuous-forest site. |
| Cont_SiteUID_generated_new | Generated identifier for continuous-forest matching. |
| Cont_UniqueFragID | Fragment identifier assigned to the continuous-forest comparison site. |
| Cont_loc_key_new | Location key used during continuous-site matching. |
| pooled_cont_new | Indicator used when pooling continuous-forest sites. |

---

### Meta-analysis calculations

| Column | Description |
|----------|-------------|
| nt | Sample size of the treatment group. |
| nc | Sample size of the comparison/control group. |
| nt_fix | Corrected treatment sample size used when original values were missing. |
| nc_fix | Corrected comparison sample size used when original values were missing. |
| vi_raw | Raw sampling variance estimate. |
| vi | Sampling variance used in the meta-analysis. |
| se_obs | Standard error associated with the effect size estimate. |

---

### Geographic variables

| Column | Description |
|----------|-------------|
| Latitude | Site latitude (decimal degrees). |
| Longitude | Site longitude (decimal degrees). |
| Country | Country in which the study was conducted. |

---

### Environmental variables

| Column | Description |
|----------|-------------|
| Matrix_mean | Dominant matrix type surrounding the focal forest fragment. |
| Forest_Group | Broad forest biome classification. |
| Forest_type | Original forest classification reported by the study. |
| Elevation | Elevation above sea level (m). |
| Rainfall_mm | Mean annual rainfall (mm). |
| Elevation_s | Standardised elevation. |
| Rainfall_mm_s | Standardised rainfall. |
| Latitude_s | Standardised latitude. |
| Longitude_s | Standardised longitude. |

---

### Publication metadata

| Column | Description |
|----------|-------------|
| Publication.Year | Year of publication of the source study. |
