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

# Helper to check if a taxon is "safe": matches on the whole specific epithet, so that
# ex. "cana" doesn't also match "talcana". Takes names with spaces or underscores.
function is_safe_taxon(taxon::String)
    parts = split(lowercase(taxon), r"[ _]+")
    return length(parts) >= 2 && parts[2] in safe_taxa
end

# Resolve the paths a per-taxon filter reads from & writes to. Reads from the filtered file if one
# has already been written (so successive filters on a taxon stack), otherwise from the clean file.
function resolve_filter_paths(taxon::String)
    clean_path = resolve_taxon_path(taxon)
    filtered_path = joinpath("data/occurrence_data/pt_occs_filtered", basename(clean_path))
    input_path = isfile(filtered_path) ? filtered_path : clean_path
    return input_path, filtered_path
end

function filter_countries(taxon::String, countries::Vector{String}; save::Bool=true)
    input_path, filtered_path = resolve_filter_paths(taxon)
    filter_countries(input_path, filtered_path, countries; save=save)
end

function filter_coords(taxon::String; save::Bool=false, kwargs...)
    input_path, filtered_path = resolve_filter_paths(taxon)
    filter_coords(input_path, filtered_path; taxon_label=taxon, save=save, kwargs...)
end

function filter_scientific_names(taxon::String, names_to_remove::Vector{String}; save::Bool=true)
    input_path, filtered_path = resolve_filter_paths(taxon)
    filter_scientific_names(input_path, filtered_path, names_to_remove; save=save)
end
