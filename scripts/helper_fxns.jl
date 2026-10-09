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

# List of safe taxa based on user request
const safe_taxa = ["uniflora", "triandra", "tripartita", "purpurea", "cana", "hyssopifolia", "fothergillii", "boliviana", "lanigera", "tenella", "scapiflora", "martinessi", "glacialis", "nitida", "picta"]

# Helper to check if a taxon is "safe": matches on the whole specific epithet, so that
# ex. "cana" doesn't also match "talcana". Takes names with spaces or underscores.
function is_safe_taxon(taxon::AbstractString)
    parts = split(lowercase(taxon), r"[ _]+")
    return length(parts) >= 2 && parts[2] in safe_taxa
end

# true for the rows to keep: everything for safe taxa; for the rest, iNat observations
# (HUMAN_OBSERVATION) are only kept if their occurrenceID is in `safe_urls`
function inat_mask(df::AbstractDataFrame, taxon::AbstractString, safe_urls)
    is_safe_taxon(taxon) && return fill(true, nrow(df))
    return Bool[
        ismissing(row.basisOfRecord) || row.basisOfRecord != "HUMAN_OBSERVATION" ||
        (!ismissing(row.occurrenceID) && row.occurrenceID in safe_urls)
        for row in eachrow(df)
    ]
end

# Preview a lat/lon filter on a taxon's final points (green = kept, red × = removed),
# ex. to pick the bounds for a new row of `geo_filters` in 3_occ_cleaning.qmd
function preview_coords(taxon::AbstractString; final_dir::String="data/occurrence_data/pt_occs_final", kwargs...)
    df = read_occs(joinpath(final_dir, replace(taxon, " " => "_") * ".csv"))
    keep = coords_mask(df; kwargs...)
    println("  $(count(!, keep)) point(s) would be removed, $(count(keep)) kept")
    plot_kept_removed(df[keep, :], df[.!keep, :]; title=taxon)
end
