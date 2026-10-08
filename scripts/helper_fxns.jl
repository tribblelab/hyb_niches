using CSV, DataFrames
using RCall

R"""
library(gatoRs)
library(ggplot2)
library(sf)
library(ggspatial)
library(gridExtra)
library(CoordinateCleaner)
library(readxl)
library(dplyr) #needs to be loaded in last, so it defaults to correct `filter` fxn
"""

# ─────────────────────────────────────────────────────────────────────────────
# Calceolaria-specific helpers
#
# the general functions (file lookup, pulling, cleaning, plotting) live in
# niche-modeling-utils/scripts — include those alongside this file.
# ─────────────────────────────────────────────────────────────────────────────

#resolve the best available CSV path for a taxon
# Shared helper: resolve the best available CSV path for a taxon, printing/erroring appropriately
function resolve_taxon_path(taxon::String)
    clean_dir = "data/occurrence_data/pt_occs_clean"
    filename = replace(taxon, " " => "_")
    
    path = resolve_clean_file(filename, clean_dir)
    if isempty(path)
        error("No cleaned or georef_merged file found for taxon: $taxon")
    end
    
    if endswith(path, "_georef_merged.csv")
        println("  Using georef_merged file: $path")
    else
        println("  Using cleaned file (no georef_merged found): $path")
    end
    
    return path
end


# List of safe taxa based on user request
const safe_taxa = ["uniflora", "triandra", "tripartita", "purpurea", "cana", "hyssopifolia", "fothergillii", "boliviana", "lanigera", "tenella", "scapiflora", "martinessi", "glacialis", "nitida", "picta"]

# Helper to check if a taxon is "safe"
function is_safe_taxon(taxon::String)
    t_lower = lowercase(taxon)
    for s in safe_taxa
        if occursin(s, t_lower)
            return true
        end
    end
    return false
end

function filter_countries(taxon::String, countries::Vector{String}; save::Bool=true)
    clean_path = resolve_taxon_path(taxon)
    filtered_path = joinpath("data/occurrence_data/pt_occs_filtered", basename(clean_path))
    filter_countries(clean_path, filtered_path, countries; save=save)
end

function filter_coords(taxon::String; save::Bool=false, kwargs...)
    clean_path = resolve_taxon_path(taxon)
    filtered_path = joinpath("data/occurrence_data/pt_occs_filtered", basename(clean_path))
    filter_coords(clean_path, filtered_path; taxon_label=taxon, save=save, kwargs...)
end

