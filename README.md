# Ustorp Naturvärdesindex (NVX)

Ett test av Ocean Health Index-metoden ([ohicore](https://github.com/OHI-Science/ohicore),
[Blekinge Hälsoindex](https://github.com/OHI-Science/bhi-blekinge)) på naturvärdesbiotoperna i Eksjö Ustorp 1:6.

**Appen:** https://peterneubauer.github.io/naturatua-nvx/

- Status: indikator → normaliserad poäng (`normalize_data`) → viktade medelvärden per delmål och kategori
  (`weighted_power_mean`) → NVX 0–100 → naturvärdesklass enligt SS 199000, kalibrerad mot expertklasserna.
- Dimensioner (`ohicore::CalculateGoalIndex`): tryck, resiliens, trend, trolig framtid, poäng.
- Rumsligt (`ohicore::CalculateAll`, BHI `SpatialIntersections`): areaviktad aggregering från biotop till
  BDU-ruta, BDU-område, skifte och fastighet. Punktdata (ljud) aggregeras till polygon via sökradie.

## Filer

- `parse_nvi.py` läser objektsbladen i NVI 2024 (textexport av rapporten, inte i repot) plus NVI 2026
  → `data/biotopes_raw.json`.
- `nvindex.py` bygger indexträdet (`STRUCTURE`, `PRESSURES`, `RESILIENCE`), gapfyller, räknar alla
  enheter och skriver `data/index_scores.json` samt appen `nvx.html` från `nvx_template.html`.
- `data/layers/` är kartlagren (samma som på lab.naturatua.com). `sweref.py` räknar om till SWEREF 99 TM.

```
python -m venv venv && venv/bin/pip install shapely
venv/bin/python -I nvindex.py
```

## Vad som är verkligt och vad som är test

Verkliga lager: NVI 2024 (Ekologigruppen) och NVI 2026 (Salmonic, utkast), geometrin, skötselavtalet.
Hypotetiska lager, slumpade med frö 2026: strukturell variation, konnektivitet, markkol och akustisk
artrikedom. eDNA saknas (levereras 2027). Objekt 27–30 saknar objektsblad i textexporten och gapfylls.
Vikter, trösklar och klassgränser är ett första förslag.
