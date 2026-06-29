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
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_taxa_data(traits_path) -> (taxa_traits, taxa_to_nativerange_dict)

Load taxa traits CSV and build a scientificName → botanicalCountries lookup.
"""
function load_taxa_data(traits_path::String)
    taxa_traits = DataFrame(CSV.File(traits_path))
    taxa_to_nativerange_dict = Dict(
        row.scientificName => row.botanicalCountries for row in eachrow(taxa_traits)
    )
    return taxa_traits, taxa_to_nativerange_dict
end

"""
    load_georef_df(taxa_filename, georef_dir) -> DataFrame or nothing

Look for a georeferenced points file (.csv or .xlsx) for a taxon, filter to
viable rows, and return a DataFrame with Float64 lat/lon columns — or `nothing`
if no usable file is found.
"""
function load_georef_df(taxa_filename::String, georef_dir::String)
    for ext in [".csv", ".xlsx"]
        candidates = filter(
            f -> startswith(f, taxa_filename) && endswith(f, ext),
            readdir(georef_dir)
        )
        isempty(candidates) && continue

        georef_path = joinpath(georef_dir, first(candidates))
        println("  Found georeferenced file: $georef_path")

        try
            df = if ext == ".csv"
                DataFrame(CSV.File(georef_path))
            else
                @rput georef_path
                R"georef_r <- readxl::read_excel(georef_path)"
                @rget georef_r
                georef_r
            end

            viable_col = findfirst(c -> uppercase(string(c)) == "VIABLE", names(df))
            if isnothing(viable_col)
                println("  WARNING: No 'Viable' column found in $georef_path — skipping georef points")
                return nothing
            end

            rename!(df, names(df)[viable_col] => :Viable)
            df[!, :Viable] = map(v -> uppercase(strip(string(v))), df.Viable)
            filter!(row -> row.Viable == "TRUE", df)

            df[!, :latitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.latitude)
            df[!, :longitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.longitude)
            filter!(row -> !ismissing(row.latitude) && !ismissing(row.longitude), df)

            println("  Georeferenced viable points: $(nrow(df))")
            return nrow(df) > 0 ? df : nothing

        catch e
            println("  WARNING: Could not read georef file $georef_path: $e")
            return nothing
        end
    end
    return nothing
end

"""
    taxon_stem(fname::String) -> String

Strip suffixes from a filename (like .csv, _georef_merged, or date strings)
to get the bare taxon stem name.
"""
function taxon_stem(fname::String)::String
    base = replace(fname, r"\.csv$" => "")
    base = replace(base, r"_georef_merged$" => "")
    base = replace(base, r"-\d{4}_\d{2}_\d{2}_cleaned$" => "")
    return base
end

"""
    resolve_clean_file(filename, clean_dir) -> String

Return the best available cleaned CSV path for a taxon in `clean_dir`:
prefers `*_georef_merged.csv`, falls back to any `*_cleaned.csv`.
Returns an empty string if nothing is found.
"""
function resolve_clean_file(filename::String, clean_dir::String)
    stem = taxon_stem(filename)
    
    georef_candidates = filter(
        f -> startswith(f, stem) && endswith(f, "_georef_merged.csv"),
        readdir(clean_dir)
    )
    cleaned_candidates = filter(
        f -> startswith(f, stem) && occursin("cleaned", f) && endswith(f, ".csv"),
        readdir(clean_dir)
    )
    if !isempty(georef_candidates)
        return joinpath(clean_dir, first(georef_candidates))
    elseif !isempty(cleaned_candidates)
        return joinpath(clean_dir, first(cleaned_candidates))
    else
        return ""
    end
end

"""
    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)

Transfer all data needed by the R plotting block to R's environment.
"""
function push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean=false)
    @rput taxa raw_file clean_file country_codes only_clean
    has_georef = !isnothing(georef_df) && nrow(georef_df) > 0
    @rput has_georef
    if has_georef
        georef_ids = string.(georef_df.ID)
        @rput georef_ids
    end
end

# Core R plotting block (rendered into a grid.arrange object in R)

const RPLOT_BLOCK = """
    raw_df  <- read.csv(raw_file)
    clean_df <- read.csv(clean_file)
    raw_df  <- raw_df[!is.na(raw_df\$latitude) & !is.na(raw_df\$longitude), ]

    if (nrow(raw_df) == 0 || nrow(clean_df) == 0) stop("empty data")

    raw_lon_min <- min(raw_df\$longitude, na.rm=TRUE) - 3
    raw_lon_max <- max(raw_df\$longitude, na.rm=TRUE) + 3
    raw_lat_min <- min(raw_df\$latitude,  na.rm=TRUE) - 3
    raw_lat_max <- max(raw_df\$latitude,  na.rm=TRUE) + 3

    all_clean_lon <- clean_df\$longitude
    all_clean_lat <- clean_df\$latitude
    clean_lon_min <- min(all_clean_lon, na.rm=TRUE) - 2
    clean_lon_max <- max(all_clean_lon, na.rm=TRUE) + 2
    clean_lat_min <- min(all_clean_lat, na.rm=TRUE) - 2
    clean_lat_max <- max(all_clean_lat, na.rm=TRUE) + 2

    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    countries <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)

    bot_underlay <- bot_regions[bot_regions\$LEVEL3_COD %in% country_codes, ]
    underlay_layer <- if (nrow(bot_underlay) > 0) {
        geom_sf(data=bot_underlay, fill="khaki1", color="goldenrod4",
                alpha=0.25, linewidth=0.3)
    } else { NULL }

    p1 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=raw_df, aes(x=longitude, y=latitude),
                   color="blue", size=1.5, alpha=0.6) +
        coord_sf(xlim=c(raw_lon_min, raw_lon_max),
                 ylim=c(raw_lat_min, raw_lat_max)) +
        labs(title=paste0("Raw (n=", nrow(raw_df), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))

    clean_title <- paste0("Cleaned (n=", nrow(clean_df), ")")

    clean_df\$plot_category <- "PRESERVED_SPECIMEN"
    if ("basisOfRecord" %in% names(clean_df)) {
        clean_df\$plot_category[clean_df\$basisOfRecord == "HUMAN_OBSERVATION"] <- "HUMAN_OBSERVATION"
    }
    if (has_georef) {
        is_georefed <- as.character(clean_df\$ID) %in% as.character(georef_ids)
        clean_df\$plot_category[is_georefed & clean_df\$plot_category == "PRESERVED_SPECIMEN"] <- "PRESERVED_SPECIMEN (georeferenced)"
    }
    
    plot_colors <- c("HUMAN_OBSERVATION" = "blue",
                     "PRESERVED_SPECIMEN" = "darkgreen",
                     "PRESERVED_SPECIMEN (georeferenced)" = "darkorange")

    p2 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=clean_df, aes(x=longitude, y=latitude, color=plot_category),
                   size=1.5, alpha=0.8) +
        scale_color_manual(values=plot_colors) +
        coord_sf(xlim=c(clean_lon_min, clean_lon_max),
                 ylim=c(clean_lat_min, clean_lat_max)) +
        labs(title=clean_title, x="Longitude", y="Latitude", color="Record Type") +
        annotation_scale(location="bl") +
        annotation_north_arrow(location="tl",
                               height=unit(0.8,"cm"), width=unit(0.8,"cm")) +
        theme_minimal() +
        theme(plot.title=element_text(size=10),
              legend.position="bottom",
              legend.title=element_text(size=9),
              legend.text=element_text(size=8))

    if (only_clean) {
        combined <- grid.arrange(p2, ncol=1, top=taxa)
    } else {
        combined <- grid.arrange(p1, p2, ncol=2, top=taxa)
    }
    print(combined)
"""

# useful functions

"""
    plot_all_occurrence_maps(;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
        pdf_file       = "data/occurrence_data/summarytables_plots/occurrence_maps_before_after.pdf",
        date_suffix    = "2026_04_23")

Iterate over every taxon in `traits_path`, build a before/after occurrence map
for each, and write all pages to `pdf_file`.  Returns the number of taxa
successfully plotted.
"""
function plot_all_occurrence_maps(;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    filtered_dir="data/occurrence_data/pt_occs_filtered",
    shapefile_path="data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
    pdf_file="data/occurrence_data/summarytables_plots/occurrence_maps_before_after.pdf",
    date_suffix="2026_04_23",
    only_clean=false,
    min_points=15
)
    taxa_traits, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    @rput pdf_file
    R"pdf(pdf_file, width=14, height=7)"

    plotted_count = 0

    for (idx, taxa) in enumerate(taxa_traits.scientificName)
        println("Processing $idx/$(length(taxa_traits.scientificName)): $taxa")

        filename = replace(taxa, " " => "_")
        raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
        clean_file = resolve_clean_file(filename, clean_dir)

        if !isfile(raw_file) || isempty(clean_file)
            println("  Skipping - files not found")
            continue
        end

        filtered_file = joinpath(filtered_dir, basename(clean_file))
        if isfile(filtered_file)
            clean_file = filtered_file
        end

        raw_df = filter(
            row -> !ismissing(row.latitude) && !ismissing(row.longitude),
            DataFrame(CSV.File(raw_file))
        )
        clean_df = DataFrame(CSV.File(clean_file))

        if nrow(raw_df) == 0 || nrow(clean_df) == 0
            println("  Skipping - no valid coordinates")
            continue
        end

        if nrow(clean_df) < min_points
            println("  Skipping - fewer than $min_points points ($(nrow(clean_df)))")
            continue
        end

        georef_df = load_georef_df(filename, georef_dir)
        country_codes = get(taxa_to_nativerange_dict, taxa, String[])

        push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
        reval(RPLOT_BLOCK)
        plotted_count += 1
    end

    R"dev.off()"
    println("\nDone. Plotted $plotted_count taxa → $pdf_file")
    return plotted_count
end


"""
    plot_species(taxa;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
        date_suffix    = "2026_04_23")

Display a before/after occurrence map for a single `taxa` string in the R
graphics viewer (no PDF output).

# Example
```julia
plot_species("Calceolaria alba")
```
"""
function plot_species(taxa::String;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    filtered_dir="data/occurrence_data/pt_occs_filtered",
    shapefile_path="data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
    date_suffix="2026_04_23",
    only_clean=false,
    min_points=0
)
    _, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    filename = replace(taxa, " " => "_")
    raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
    clean_file = resolve_clean_file(filename, clean_dir)

    isfile(raw_file) || error("Raw file not found: $raw_file")
    isempty(clean_file) && error("No cleaned or georef_merged file found for: $taxa")

    filtered_file = joinpath(filtered_dir, basename(clean_file))
    if isfile(filtered_file)
        clean_file = filtered_file
    end

    raw_df = filter(
        row -> !ismissing(row.latitude) && !ismissing(row.longitude),
        DataFrame(CSV.File(raw_file))
    )
    clean_df = DataFrame(CSV.File(clean_file))

    nrow(raw_df) > 0 || error("No valid coordinates in raw file for $taxa")
    nrow(clean_df) > 0 || error("No rows in cleaned file for $taxa")
    nrow(clean_df) >= min_points || error("Only $(nrow(clean_df)) points in cleaned file, which is < min_points=$min_points")

    georef_df = load_georef_df(filename, georef_dir)
    country_codes = get(taxa_to_nativerange_dict, taxa, String[])

    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
    reval(RPLOT_BLOCK)

    println("Plot displayed for: $taxa")
    return nothing
end

"""
    prepare_geolocate_files(;
        input_dir  = "data/occurrence_data/pt_occs_to_georeference",
        output_dir = "data/occurrence_data/pt_occs_to_georeference_geolocate_format")
 
Read every CSV in `input_dir`, reformat columns to GEOLocate batch input
format, and write the result to `output_dir` (same filename).
 
Expected input columns: locality, country, stateProvince, county, latitude,
longitude, ID, scientificName, basisOfRecord.
 
Output columns (GEOLocate format): "locality string", country, state, county,
latitude, longitude, "correction status", precision, "error polygon",
"multiple results", ID, name, basis.
"""
function prepare_geolocate_files(;
    input_dir="data/occurrence_data/pt_occs_to_georeference",
    output_dir="data/occurrence_data/pt_occs_to_georeference_geolocate_format"
)
    isdir(output_dir) || mkpath(output_dir)

    files = DataFrames.filter(f -> endswith(f, ".csv"), readdir(input_dir))

    if isempty(files)
        println("No CSV files found in $input_dir")
        return 0
    end

    converted_count = 0

    for (idx, file) in enumerate(files)
        println("Processing $idx/$(length(files)): $file")

        input_file = joinpath(input_dir, file)
        output_file = joinpath(output_dir, file)

        @rput input_file output_file

        R"""
        rawdf_GeoRef <- read.csv(input_file)

        if (nrow(rawdf_GeoRef) == 0) {
            cat("  Skipping - empty file\n")
            next
        }

        rawdf_GeoRef <- rawdf_GeoRef %>%
            dplyr::select("locality string" = locality,
                          country,
                          state = stateProvince,
                          county,
                          latitude,
                          longitude,
                          ID,
                          name = scientificName,
                          basis = basisOfRecord)

        rawdf_GeoRef$'correction status' <- ""
        rawdf_GeoRef$precision           <- ""
        rawdf_GeoRef$'error polygon'     <- ""
        rawdf_GeoRef$'multiple results'  <- ""

        rawdf_GeoRef2 <- rawdf_GeoRef[, c("locality string", "country",
                                           "state", "county", "latitude",
                                           "longitude", "correction status",
                                           "precision", "error polygon",
                                           "multiple results", "ID",
                                           "name", "basis")]

        write.csv(rawdf_GeoRef2, output_file, row.names = FALSE)
        """

        converted_count += 1
    end

    println("\nDone. Converted $converted_count files → $output_dir")
    return converted_count
end

"""
    load_pt_occs_df(;
        clean_dir = "data/occurrence_data/pt_occs_clean") -> DataFrame

Stack all cleaned occurrence CSVs into a single DataFrame, selecting the best
available file per taxon (`_georef_merged.csv` preferred over `_cleaned.csv`,
via `resolve_clean_file`).  Two extra columns are prepended:

    taxon       – bare taxon name (spaces → underscores, no date/suffix)
    source_file – basename of the file that was loaded

All other columns are passed through as-is.
"""
function load_pt_occs_df(;
    clean_dir="data/occurrence_data/pt_occs_clean"
)
    all_files = filter(f -> endswith(f, ".csv"), readdir(clean_dir))
    taxa = sort(unique(taxon_stem.(all_files)))

    chunks = DataFrame[]
    for t in taxa
        fpath = resolve_clean_file(t, clean_dir)
        isempty(fpath) && continue

        df = CSV.read(fpath, DataFrame; missingstring=["", "NA"])
        insertcols!(df, 1,
            :taxon => t,
            :source_file => basename(fpath)
        )
        push!(chunks, df)
    end

    return vcat(chunks...; cols=:union)
end

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

function filter_countries(taxon::String, countries::Vector{String})
    path = resolve_taxon_path(taxon)
    cleandf = DataFrame(CSV.File(path))

    # check that at least one country is present in df
    dfcountries = cleandf.country
    if length(intersect(countries, dfcountries)) == 0
        error("no input country provided in `countries` list matches countries existing in dataframe")
    end

    # filter every country out of df
    for country in countries
        filter!(:country => x -> x != country, cleandf)
    end
    CSV.write(path, cleandf)
end

# Remove rows outside the specified lat/lon bounding box.
# All bounds are optional — specify only the axes you want to constrain.
# Points with missing coordinates are dropped when a bound is supplied for that axis.
# Always plots a preview map (green = kept, red × = removed).
# Set save=true to write the filtered result back to the file.
function filter_coords(taxon::String;
    lat_min=nothing, lat_max=nothing,
    lon_min=nothing, lon_max=nothing,
    save::Bool=false)

    if all(isnothing, (lat_min, lat_max, lon_min, lon_max))
        error("At least one of lat_min, lat_max, lon_min, lon_max must be specified")
    end

    path = resolve_taxon_path(taxon)
    df = DataFrame(CSV.File(path))

    # Build a row-wise keep predicate
    function keep(row)
        lat = row.latitude
        lon = row.longitude
        if !isnothing(lat_min) || !isnothing(lat_max)
            (ismissing(lat) || !isa(lat, Number)) && return false
            !isnothing(lat_min) && lat < lat_min && return false
            !isnothing(lat_max) && lat > lat_max && return false
        end
        if !isnothing(lon_min) || !isnothing(lon_max)
            (ismissing(lon) || !isa(lon, Number)) && return false
            !isnothing(lon_min) && lon < lon_min && return false
            !isnothing(lon_max) && lon > lon_max && return false
        end
        return true
    end

    keep_mask = [keep(row) for row in eachrow(df)]
    kept_df = df[keep_mask, :]
    removed_df = df[.!keep_mask, :]
    n_kept = nrow(kept_df)
    n_removed = nrow(removed_df)

    println("  Preview: $n_removed point(s) will be removed, $n_kept will be kept")

    # Plot preview: green = kept, red × = removed
    taxon_label = taxon
    kept_lat = Vector{Union{Missing,Float64}}(kept_df.latitude)
    kept_lon = Vector{Union{Missing,Float64}}(kept_df.longitude)
    removed_lat = Vector{Union{Missing,Float64}}(removed_df.latitude)
    removed_lon = Vector{Union{Missing,Float64}}(removed_df.longitude)
    @rput taxon_label kept_lat kept_lon removed_lat removed_lon
    R"""
    kept_pts    <- data.frame(latitude  = as.numeric(kept_lat),
                              longitude = as.numeric(kept_lon))
    removed_pts <- data.frame(latitude  = as.numeric(removed_lat),
                              longitude = as.numeric(removed_lon))
    all_lat <- c(kept_pts$latitude,  removed_pts$latitude)
    all_lon <- c(kept_pts$longitude, removed_pts$longitude)
    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    borders  <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)
    p <- ggplot() + world + borders +
        geom_point(data=kept_pts,
                   aes(x=longitude, y=latitude),
                   color="darkgreen", size=2, alpha=0.7) +
        { if (nrow(removed_pts) > 0)
              geom_point(data=removed_pts,
                         aes(x=longitude, y=latitude),
                         color="red", shape=4, size=3.5, stroke=1.2, alpha=0.9)
          else NULL } +
        coord_sf(xlim=c(min(all_lon, na.rm=TRUE)-2, max(all_lon, na.rm=TRUE)+2),
                 ylim=c(min(all_lat, na.rm=TRUE)-2, max(all_lat, na.rm=TRUE)+2)) +
        labs(title=paste0(taxon_label, "  —  green: kept (", nrow(kept_pts),
                          ")   red ×: removed (", nrow(removed_pts), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))
    print(p)
    """

    if save
        CSV.write(path, kept_df)
        println("  Written: $path")
    else
        println("  Preview only — re-run with save=true to apply.")
    end
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

# ─────────────────────────────────────────────────────────────────────────────
# Plotting Utilities
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_bot_regions(shapefile_path)

Load TDWG botanical country polygons into R (once, globally).
"""
function load_bot_regions(shapefile_path::String)
    @rput shapefile_path
    R"bot_regions = st_read(shapefile_path, quiet = TRUE)"
end

"""
    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)

Transfer all data needed by the R plotting block to R's environment.
"""
function push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean=false)
    @rput taxa raw_file clean_file country_codes only_clean
    has_georef = !isnothing(georef_df) && nrow(georef_df) > 0
    @rput has_georef
    if has_georef
        georef_ids = string.(georef_df.ID)
        @rput georef_ids
    end
end

# Core R plotting block (rendered into a grid.arrange object in R)

const RPLOT_BLOCK = """
    raw_df  <- read.csv(raw_file)
    clean_df <- read.csv(clean_file)
    raw_df  <- raw_df[!is.na(raw_df\$latitude) & !is.na(raw_df\$longitude), ]

    if (nrow(raw_df) == 0 || nrow(clean_df) == 0) stop("empty data")

    raw_lon_min <- min(raw_df\$longitude, na.rm=TRUE) - 3
    raw_lon_max <- max(raw_df\$longitude, na.rm=TRUE) + 3
    raw_lat_min <- min(raw_df\$latitude,  na.rm=TRUE) - 3
    raw_lat_max <- max(raw_df\$latitude,  na.rm=TRUE) + 3

    all_clean_lon <- clean_df\$longitude
    all_clean_lat <- clean_df\$latitude
    clean_lon_min <- min(all_clean_lon, na.rm=TRUE) - 2
    clean_lon_max <- max(all_clean_lon, na.rm=TRUE) + 2
    clean_lat_min <- min(all_clean_lat, na.rm=TRUE) - 2
    clean_lat_max <- max(all_clean_lat, na.rm=TRUE) + 2

    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    countries <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)

    bot_underlay <- bot_regions[bot_regions\$LEVEL3_COD %in% country_codes, ]
    underlay_layer <- if (nrow(bot_underlay) > 0) {
        geom_sf(data=bot_underlay, fill="khaki1", color="goldenrod4",
                alpha=0.25, linewidth=0.3)
    } else { NULL }

    p1 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=raw_df, aes(x=longitude, y=latitude),
                   color="blue", size=1.5, alpha=0.6) +
        coord_sf(xlim=c(raw_lon_min, raw_lon_max),
                 ylim=c(raw_lat_min, raw_lat_max)) +
        labs(title=paste0("Raw (n=", nrow(raw_df), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))

    clean_title <- paste0("Cleaned (n=", nrow(clean_df), ")")

    clean_df\$plot_category <- "PRESERVED_SPECIMEN"
    if ("basisOfRecord" %in% names(clean_df)) {
        clean_df\$plot_category[clean_df\$basisOfRecord == "HUMAN_OBSERVATION"] <- "HUMAN_OBSERVATION"
    }
    if (has_georef) {
        is_georefed <- as.character(clean_df\$ID) %in% as.character(georef_ids)
        clean_df\$plot_category[is_georefed & clean_df\$plot_category == "PRESERVED_SPECIMEN"] <- "PRESERVED_SPECIMEN (georeferenced)"
    }
    
    plot_colors <- c("HUMAN_OBSERVATION" = "blue",
                     "PRESERVED_SPECIMEN" = "darkgreen",
                     "PRESERVED_SPECIMEN (georeferenced)" = "darkorange")

    p2 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=clean_df, aes(x=longitude, y=latitude, color=plot_category),
                   size=1.5, alpha=0.8) +
        scale_color_manual(values=plot_colors) +
        coord_sf(xlim=c(clean_lon_min, clean_lon_max),
                 ylim=c(clean_lat_min, clean_lat_max)) +
        labs(title=clean_title, x="Longitude", y="Latitude", color="Record Type") +
        annotation_scale(location="bl") +
        annotation_north_arrow(location="tl",
                               height=unit(0.8,"cm"), width=unit(0.8,"cm")) +
        theme_minimal() +
        theme(plot.title=element_text(size=10),
              legend.position="bottom",
              legend.title=element_text(size=9),
              legend.text=element_text(size=8))

    if (only_clean) {
        combined <- grid.arrange(p2, ncol=1, top=taxa)
    } else {
        combined <- grid.arrange(p1, p2, ncol=2, top=taxa)
    }
    print(combined)
"""

# ─────────────────────────────────────────────────────────────────────────────
# Plotting Functions
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_all_occurrence_maps(;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
        pdf_file       = "data/occurrence_data/summarytables_plots/occurrence_maps_before_after.pdf",
        date_suffix    = "2026_04_23")

Iterate over every taxon in `traits_path`, build a before/after occurrence map
for each, and write all pages to `pdf_file`.  Returns the number of taxa
successfully plotted.
"""
function plot_all_occurrence_maps(;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    filtered_dir="data/occurrence_data/pt_occs_filtered",
    shapefile_path="data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
    pdf_file="data/occurrence_data/summarytables_plots/occurrence_maps_before_after.pdf",
    date_suffix="2026_04_23",
    only_clean=false,
    min_points=15
)
    taxa_traits, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    @rput pdf_file
    R"pdf(pdf_file, width=14, height=7)"

    plotted_count = 0

    for (idx, taxa) in enumerate(taxa_traits.scientificName)
        println("Processing $idx/$(length(taxa_traits.scientificName)): $taxa")

        filename = replace(taxa, " " => "_")
        raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
        clean_file = resolve_clean_file(filename, clean_dir)

        if !isfile(raw_file) || isempty(clean_file)
            println("  Skipping - files not found")
            continue
        end

        filtered_file = joinpath(filtered_dir, basename(clean_file))
        if isfile(filtered_file)
            clean_file = filtered_file
        end

        raw_df = filter(
            row -> !ismissing(row.latitude) && !ismissing(row.longitude),
            DataFrame(CSV.File(raw_file))
        )
        clean_df = DataFrame(CSV.File(clean_file))

        if nrow(raw_df) == 0 || nrow(clean_df) == 0
            println("  Skipping - no valid coordinates")
            continue
        end

        if nrow(clean_df) < min_points
            println("  Skipping - fewer than $min_points points ($(nrow(clean_df)))")
            continue
        end

        georef_df = load_georef_df(filename, georef_dir)
        country_codes = get(taxa_to_nativerange_dict, taxa, String[])

        push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
        reval(RPLOT_BLOCK)
        plotted_count += 1
    end

    R"dev.off()"
    println("\nDone. Plotted $plotted_count taxa → $pdf_file")
    return plotted_count
end


"""
    plot_species(taxa;
        traits_path    = "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
        raw_dir        = "data/occurrence_data/pt_occs_raw",
        clean_dir      = "data/occurrence_data/pt_occs_clean",
        georef_dir     = "data/occurrence_data/pt_occs_georeferenced",
        shapefile_path = "data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
        date_suffix    = "2026_04_23")

Display a before/after occurrence map for a single `taxa` string in the R
graphics viewer (no PDF output).

# Example
```julia
plot_species("Calceolaria alba")
```
"""
function plot_species(taxa::String;
    traits_path="data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv",
    raw_dir="data/occurrence_data/pt_occs_raw",
    clean_dir="data/occurrence_data/pt_occs_clean",
    georef_dir="data/occurrence_data/pt_occs_georeferenced",
    filtered_dir="data/occurrence_data/pt_occs_filtered",
    shapefile_path="data/occurrence_data/supp_data/bot_country_shapefiles/level3.shp",
    date_suffix="2026_04_23",
    only_clean=false,
    min_points=0
)
    _, taxa_to_nativerange_dict = load_taxa_data(traits_path)
    load_bot_regions(shapefile_path)

    filename = replace(taxa, " " => "_")
    raw_file = joinpath(raw_dir, "$(filename)-$(date_suffix).csv")
    clean_file = resolve_clean_file(filename, clean_dir)

    isfile(raw_file) || error("Raw file not found: $raw_file")
    isempty(clean_file) && error("No cleaned or georef_merged file found for: $taxa")

    filtered_file = joinpath(filtered_dir, basename(clean_file))
    if isfile(filtered_file)
        clean_file = filtered_file
    end

    raw_df = filter(
        row -> !ismissing(row.latitude) && !ismissing(row.longitude),
        DataFrame(CSV.File(raw_file))
    )
    clean_df = DataFrame(CSV.File(clean_file))

    nrow(raw_df) > 0 || error("No valid coordinates in raw file for $taxa")
    nrow(clean_df) > 0 || error("No rows in cleaned file for $taxa")
    nrow(clean_df) >= min_points || error("Only $(nrow(clean_df)) points in cleaned file, which is < min_points=$min_points")

    georef_df = load_georef_df(filename, georef_dir)
    country_codes = get(taxa_to_nativerange_dict, taxa, String[])

    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
    reval(RPLOT_BLOCK)

    println("Plot displayed for: $taxa")
    return nothing
end
