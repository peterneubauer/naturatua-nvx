# Ustorp Nature Value Index (NVX)

A test of the Ocean Health Index method ([ohicore](https://github.com/OHI-Science/ohicore),
[Blekinge Health Index](https://github.com/OHI-Baltic/bhi-blekinge)) applied on land, to the nature value
biotopes of the Swedish-standard inventories (NVI, SS 199000:2023) at Eksjö Ustorp 1:6, Sweden.

**App:** https://peterneubauer.github.io/naturatua-nvx/

- Status: indicator → rescaled score (`normalize_data`) → weighted means per sub-goal and category
  (`weighted_power_mean`) → NVX 0–100 → nature value class per SS 199000, calibrated against the expert classes.
- Dimensions (`ohicore::CalculateGoalIndex`): pressure (per-category rank matrix), resilience (per-category
  weights), trend, likely future, score.
- Spatial (`ohicore::CalculateAll`, BHI `SpatialIntersections`): area-weighted aggregation from biotope to BDU
  plot, BDU area, land parcel and property, plus groupings by habitat type and expert class. Point data
  (acoustic recorders) is attached to polygons by search radius.

This is a Python port of the R logic, not the R code itself. The formulas are the same; see the method
section on the page and the comments in `nvindex.py` for the mapping to the R functions.

## Files

- `parse_nvi.py` reads the object sheets of NVI 2024 (text export of the report, not in the repo) plus NVI 2026
  → `data/biotopes_raw.json`.
- `nvindex.py` builds the index tree (`STRUCTURE`, `PRESSURES`, `RESILIENCE`), gapfills, computes all units and
  writes `data/index_scores.json` and the app `nvx.html` from `nvx_template.html`.
- `data/layers/` holds the map layers (same as on lab.naturatua.com). `sweref.py` converts to SWEREF 99 TM.

```
python -m venv venv && venv/bin/pip install shapely
venv/bin/python -I nvindex.py
```

## What is real and what is test data

Real layers: NVI 2024 (Ekologigruppen AB) and NVI 2026 (Salmonic AB, draft), the geometry, the acoustic
point positions and the management agreement. Hypothetical layers, seeded random (seed 2026): structural
variation, connectivity, soil carbon and acoustic species richness. eDNA is missing (delivery 2027).
Objects 27–30 lack object sheets in the text export and are gapfilled. Weights, reference points and class
thresholds are a first proposal.
