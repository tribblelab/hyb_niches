"""
Fetch non-research-grade iNaturalist observations by aacocucci and aliciasersic
(both observed-by and identified-by), filtered to only taxa in
traits_species_cleaned-2026_02_25.csv, then format them to match the
pt_occs_clean CSV schema (iNaturalist rows).

Output: data/occurrence_data/supp_data/non_research_obs_raw.json        (full API details, filtered taxa)
        data/occurrence_data/supp_data/non_research_inat_formatted.csv  (ready to append to pt_occs_clean files)
"""

import urllib.request
import json
import time
import csv
import os
from datetime import datetime

HEADERS    = {"User-Agent": "calceolaria-scraper/1.0"}
REPO_ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TRAITS_CSV = os.path.join(REPO_ROOT, "data/taxonomy_trait_data/traits_species_cleaned-2026_02_25.csv")
OUT_JSON   = os.path.join(REPO_ROOT, "data/occurrence_data/supp_data/non_research_obs_raw.json")
OUT_CSV    = os.path.join(REPO_ROOT, "data/occurrence_data/supp_data/non_research_inat_formatted.csv")

# ── 0. Load taxa list from traits CSV ────────────────────────────────────────

def load_traits_taxa(path):
    """Read scientificName column from traits CSV.
    Returns:
      - traits_names: set of full names as in the CSV (e.g. "Calceolaria X subsp. Y")
      - inat_to_traits: dict mapping iNat-style names to traits-style names
        iNat uses "Calceolaria X Y" for subspecies (no "subsp." token)
    """
    traits_names = set()
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["scientificName"].strip()
            if name:
                traits_names.add(name)

    # Build mapping: strip "subsp." so we can match iNat names
    inat_to_traits = {}
    for name in traits_names:
        # "Calceolaria X subsp. Y"  -> iNat "Calceolaria X Y"
        inat_name = name.replace(" subsp.", "").replace(" var.", "").replace(" f.", "")
        inat_to_traits[inat_name] = name
        # also map the exact name in case iNat does include it
        inat_to_traits[name] = name

    return traits_names, inat_to_traits

print(f"Loading traits taxa from: {TRAITS_CSV}")
traits_names, inat_to_traits = load_traits_taxa(TRAITS_CSV)
print(f"  {len(traits_names)} taxa in traits CSV")

# ── 1. Pull non-research-grade observations ──────────────────────────────────

def api_get(url, params):
    from urllib.parse import urlencode
    full_url = url + "?" + urlencode(params)
    req = urllib.request.Request(full_url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)

users     = ["aacocucci", "aliciasersic"]
users_str = ",".join(users)
all_obs   = {}   # obs_id -> full obs dict
skipped_taxa = set()

for role_param in ["user_id", "ident_user_id"]:
    page = 1
    while True:
        data = api_get(
            "https://api.inaturalist.org/v1/observations",
            {
                role_param:      users_str,
                "taxon_id":      129854,          # Calceolaria genus
                "quality_grade": "needs_id",      # NON-research grade only
                "rank":          "species",
                "per_page":      200,
                "page":          page,
                "order":         "desc",
                "order_by":      "created_at",
            }
        )
        results = data.get("results", [])
        total   = data.get("total_results", 0)
        print(f"{role_param} | page {page}: got {len(results)} of {total}")

        for obs in results:
            taxon = obs.get("taxon")
            if not taxon:
                continue
            inat_name = taxon.get("name", "")
            if " " not in inat_name:      # skip genus-only IDs
                continue

            # ── filter: only keep taxa in our traits list ──
            if inat_name not in inat_to_traits:
                skipped_taxa.add(inat_name)
                continue

            obs_id = obs["id"]
            if obs_id not in all_obs:
                all_obs[obs_id] = obs

        if len(results) < 200 or page * 200 >= total:
            break
        page += 1
        time.sleep(1)

print(f"\nTotal unique non-research obs (in-taxa): {len(all_obs)}")
if skipped_taxa:
    print(f"Skipped (not in traits CSV): {sorted(skipped_taxa)}")

# Save raw
with open(OUT_JSON, "w") as f:
    json.dump(list(all_obs.values()), f, indent=2)
print(f"Saved raw JSON: {OUT_JSON}")

# ── 2. Format to match pt_occs_clean CSV schema ──────────────────────────────

CSV_COLS = [
    "scientificName", "genus", "specificEpithet", "infraspecificEpithet",
    "ID", "occurrenceID", "basisOfRecord", "eventDate",
    "year", "month", "day",
    "institutionCode", "recordedBy",
    "country", "county", "stateProvince",
    "locality",
    "latitude", "longitude", "coordinateUncertaintyInMeters",
    "informationWithheld", "habitat",
    "aggregator", "accepted_name",
]

def parse_inat_name(inat_name, inat_to_traits):
    """Resolve iNat name to traits-canonical name, then parse parts."""
    canonical = inat_to_traits.get(inat_name, inat_name)
    parts = canonical.split()
    genus   = parts[0] if len(parts) >= 1 else ""
    epithet = parts[1] if len(parts) >= 2 else ""
    rank_tokens = {"subsp.", "var.", "f.", "ssp."}
    infra_parts = [p for p in parts[2:] if p not in rank_tokens]
    infra = infra_parts[0] if infra_parts else None
    return canonical, canonical, genus, epithet, infra

def build_locality(obs):
    desc = obs.get("description") or "NA"
    place_guess = obs.get("place_guess") or "NA"
    desc_clean = desc.replace("\n", " ").replace("\r", " ").strip() if desc != "NA" else "NA"
    return f"locality:  NA, occurrenceRemarks: {desc_clean}, verbatimLocality: {place_guess}"

rows = []
for obs in all_obs.values():
    taxon     = obs.get("taxon", {})
    inat_name = taxon.get("name", "")

    sci_name, accepted_name, genus, epithet, infra = parse_inat_name(inat_name, inat_to_traits)

    obs_id = obs["id"]
    occ_id = f"https://www.inaturalist.org/observations/{obs_id}"

    observed_on = obs.get("observed_on") or ""
    if observed_on:
        try:
            dt = datetime.strptime(observed_on, "%Y-%m-%d")
            yr, mo, dy = dt.year, dt.month, dt.day
            event_date = observed_on
        except Exception:
            yr = mo = dy = None
            event_date = observed_on
    else:
        yr = mo = dy = None
        event_date = None

    location = obs.get("location") or ""
    if location and "," in location:
        lat_s, lon_s = location.split(",", 1)
        try:
            lat = round(float(lat_s), 6)
            lon = round(float(lon_s), 6)
        except Exception:
            lat = lon = None
    else:
        lat = lon = None

    coord_uncertainty = obs.get("positional_accuracy") or None

    info_withheld = None
    if obs.get("obscured") or obs.get("geoprivacy") not in (None, "open"):
        info_withheld = "coordinates obscured on iNaturalist"

    recorded_by = obs.get("user", {}).get("login") or None

    place_guess = obs.get("place_guess") or ""
    country = state_province = county = None
    if place_guess:
        parts_pg = [p.strip() for p in place_guess.split(",")]
        if len(parts_pg) >= 2:
            country        = parts_pg[-1]
            state_province = parts_pg[-2]
        elif len(parts_pg) == 1:
            country = parts_pg[0]

    locality_str = build_locality(obs)

    row = {
        "scientificName":                sci_name,
        "genus":                         genus,
        "specificEpithet":               epithet,
        "infraspecificEpithet":          infra,
        "ID":                            obs_id,
        "occurrenceID":                  occ_id,
        "basisOfRecord":                 "HUMAN_OBSERVATION",
        "eventDate":                     event_date,
        "year":                          yr,
        "month":                         mo,
        "day":                           dy,
        "institutionCode":               "iNaturalist",
        "recordedBy":                    recorded_by,
        "country":                       country,
        "county":                        county,
        "stateProvince":                 state_province,
        "locality":                      locality_str,
        "latitude":                      lat,
        "longitude":                     lon,
        "coordinateUncertaintyInMeters": coord_uncertainty,
        "informationWithheld":           info_withheld,
        "habitat":                       None,
        "aggregator":                    "iNaturalist",
        "accepted_name":                 accepted_name,
    }
    rows.append(row)
    print(f"  {obs_id} | {accepted_name} | {recorded_by} | {obs.get('quality_grade')} | lat={lat}, lon={lon}")

with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=CSV_COLS, extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)

print(f"\nSaved formatted CSV ({len(rows)} rows): {OUT_CSV}")

taxon_counts = {}
for r in rows:
    t = r["accepted_name"]
    taxon_counts[t] = taxon_counts.get(t, 0) + 1
print("\nBy taxon:")
for k, v in sorted(taxon_counts.items()):
    print(f"  {k}: {v}")
