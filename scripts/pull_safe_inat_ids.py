"""
Fetch every species-level Calceolaria iNaturalist observation that was observed
by or identified by aacocucci or aliciasersic. These are the "safe" iNat
observations that are kept in 3_occ_cleaning.qmd for taxa where we don't
trust iNat IDs in general.

Output: data/occurrence_data/supp_data/calceolaria_ids.json
        (or the path given as the first argument)
"""

import urllib.request
import json
import time
import os
import sys
from urllib.parse import urlencode

HEADERS   = {"User-Agent": "calceolaria-scraper/1.0"}
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO_ROOT, "data/occurrence_data/supp_data/calceolaria_ids.json")

def api_get(url, params):
    full_url = url + "?" + urlencode(params)
    req = urllib.request.Request(full_url, headers=HEADERS)
    with urllib.request.urlopen(req, timeout=20) as r:
        return json.load(r)

users     = ["aacocucci", "aliciasersic"]
users_str = ",".join(users)
all_obs   = {}   # obs_id -> summary dict

# Query 1: observations BY these users
# Query 2: observations IDENTIFIED BY these users
for role_param in ["user_id", "ident_user_id"]:
    page = 1
    while True:
        data = api_get(
            "https://api.inaturalist.org/v1/observations",
            {
                role_param: users_str,
                "taxon_id": 129854,          # Calceolaria genus
                "rank":     "species",       # only species-level IDs
                "per_page": 200,
                "page":     page,
                "order":    "desc",
                "order_by": "created_at",
            },
        )
        results = data.get("results", [])
        total   = data.get("total_results", 0)
        print(f"{role_param} | page {page}: got {len(results)} of {total}")

        for obs in results:
            taxon = obs.get("taxon")
            if not taxon:
                continue
            name = taxon.get("name", "")
            # ensure it's a species (has a space = genus + epithet)
            if " " not in name:
                continue
            obs_id = obs["id"]
            if obs_id not in all_obs:
                all_obs[obs_id] = {
                    "obs_id":        obs_id,
                    "observer":      obs["user"]["login"],
                    "taxon_name":    name,
                    "taxon_id":      taxon.get("id"),
                    "date":          str(obs.get("created_at", ""))[:10],
                    "quality_grade": obs.get("quality_grade", "none"),
                    "url":           f"https://www.inaturalist.org/observations/{obs_id}",
                }

        if len(results) < 200 or page * 200 >= total:
            break
        page += 1
        time.sleep(1)

with open(OUT_JSON, "w") as f:
    json.dump(list(all_obs.values()), f, indent=2)

print(f"\n{len(all_obs)} species-level Calceolaria observations saved to {OUT_JSON}")
