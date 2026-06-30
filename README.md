# FAO_CC_Risks_to_Fed_Mariculture

Repository for work being done under the FAO project "Climate change risks to interconnected food production systems with a focus on fed mariculture".

For the Shiny app repository, see [Fish-MIP/FAO_CCRFM_shiny](https://github.com/Fish-MIP/FAO_CCRFM_shiny).

## Prepped data folders

* `data/prepped_data/mapped_farms_meanyear`: A subset of locations mapped in @clawson_mapping_2022 (those that might be assigned to a species in this analysis) are extracted from the "meanyear" temperature rasters prepared in `data-prep-scripts/historical-sst.R`. These files (one per 365-day meanyear per farm_ID) are categorised by country, FAO fishing region, and species group. These files were used to assign mapped farms to species in this analysis.
* `data/prepped_data/assigned_farms_meanyear`: The "meanyear" temperatures of all assigned farms used in this analysis. "Meanyear" temperatures were extracted from rasters prepared in `data-prep-scripts/historical-sst.R` and farm locations were assinged in `documentation-qmds/04_assigning_locations.qmd`. All farm locations have a unique ID. IDs up to 95443 correspond to the mapped farm IDs from @clawson_mapping_2022, while farm IDs > 95443 were created during the assignment process. 
* `data/prepped_data/assigned_farms_alltemps`: All temperature timeseries (2025-2099) for all assigned farms used in this analysis. Farms were assigned in `documentation-qmds/04_assigning_locations.qmd` and temperatures were extracted in `documentation-qmds/05_extracting_future_temperatures.qmd`.
