# Naturvärdesindex (NVX) för Ustorp: Ocean Health Index-ramverket (ohicore + BHI Blekinge) tillämpat på
# naturvärdesbiotoper.
#
# Status:     indikator -> normaliserad poäng 0-1 (normalize_data) -> viktade medelvärden per delmål och
#             kategori (weighted_power_mean, BHI) -> NVX 0-100 -> naturvärdesklass (SS 199000).
# Dimensioner (ohicore CalculateGoalIndex): tryck (pressures), resiliens, trend, trolig framtid, poäng
#             xF = DISCOUNT*(1 + BETA*trend + (1-BETA)*(resiliens - tryck)) * status,  score = (status + xF)/2
# Rumsligt (ohicore CalculateAll + BHI SpatialIntersections): biotoppolygonerna är rapportzoner. Varje
#             bedömningsenhet (BDU-ruta, BDU-område, skifte, fastighet) får weighted.mean(score, area) per
#             gren och dimension, med fraction_of_assess_unit som vikt. Punktdata aggregeras till polygon.
#
# Användning:  venv/bin/python -I index/nvindex.py   -> index/data/index_scores.json + index/nvx.html
import json, math, random, sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from sweref import to_tm
from shapely.geometry import shape, Point, mapping
from shapely.ops import unary_union
from shapely import transform as sh_transform

LAYERS = HERE / 'data/layers'
YEAR = 2026
BETA, DISCOUNT = 0.67, 1.0  # ohicore goal_beta / goal_discount


def tm(geom):
    return sh_transform(shape(geom), lambda a: [list(to_tm(x, y)) for x, y in a[:, :2]]).buffer(0)


# ── 1. Indexstruktur (IndexStructure.csv i BHI) ─────────────────────────────────────────────────────
# fn='norm': threshold = värde som ger 1.0, zero = värde som ger 0.0. fn='mean': power mean med exponent.
STRUCTURE = [
    dict(branch='NVX', main=None, w=1, fn='mean', exponent=0.5, name='Nature Value Index'),
    dict(branch='BV', main='NVX', w=1, fn='mean', exponent=1, name='Biotope value', short='Biotope'),
    dict(branch='AV', main='NVX', w=1, fn='mean', exponent=1, name='Species value', short='Species'),
    dict(branch='LM', main='NVX', w=1, fn='mean', exponent=1, name='Landscape & soil', short='Landscape'),
    dict(branch='BV_TILL', main='BV', w=1, fn='mean', exponent=1, name='Condition'),
    dict(branch='BV_SALL', main='BV', w=1, fn='mean', exponent=1, name='Rarity'),
    dict(branch='BV_FUNK', main='BV', w=1, fn='mean', exponent=1, name='Ecological function'),
    dict(branch='AV_VARD', main='AV', w=1, fn='mean', exponent=1, name='Indicator species'),
    dict(branch='AV_HOT', main='AV', w=1, fn='mean', exponent=1, name='Threatened & protected'),
    dict(branch='AV_SENS', main='AV', w=1, fn='mean', exponent=1, name='Sensor data'),
    dict(branch='LM_STRUK', main='LM', w=1, fn='mean', exponent=1, name='Structure & shape'),
    dict(branch='LM_FUNK', main='LM', w=1, fn='mean', exponent=1, name='Landscape function'),
    # indikatorer. source NVI/GEO/AVTAL = verkliga, DRONE/SAT/SOIL/AUDIO/EDNA = hypotetiska testlager
    dict(branch='naturlighet', main='BV_TILL', w=2, fn='norm', threshold=3, zero=0, units='SIS level 0–3', source='NVI', name='Naturalness and impact'),
    dict(branch='dodved', main='BV_TILL', w=1, fn='norm', threshold=20, zero=0, units='m³/ha', source='NVI', name='Dead wood'),
    dict(branch='naturvardestrad', main='BV_TILL', w=1, fn='norm', threshold=3, zero=0, units='trees/ha', source='NVI', name='Habitat trees'),
    dict(branch='sallsynthet', main='BV_SALL', w=2, fn='norm', threshold=3, zero=0, units='SIS level 0–3', source='NVI', name='Biotope rarity'),
    dict(branch='natura', main='BV_SALL', w=1, fn='norm', threshold=1, zero=0, units='0/1', source='NVI', name='Natura 2000 habitat type'),
    dict(branch='betydelse', main='BV_FUNK', w=2, fn='norm', threshold=3, zero=0, units='SIS level 0–3', source='NVI', name='Importance for populations'),
    dict(branch='nyckelbiotop', main='BV_FUNK', w=1, fn='norm', threshold=1, zero=0, units='0/1', source='NVI', name='Woodland key habitat quality'),
    dict(branch='signalsumma', main='AV_VARD', w=2, fn='norm', threshold=20, zero=0, units='weighted sum', source='NVI', name='Indicator value, sum'),
    dict(branch='artantal', main='AV_VARD', w=1, fn='norm', threshold=15, zero=0, units='species', source='NVI', name='Conservation species count'),
    dict(branch='rodlistade', main='AV_HOT', w=1, fn='norm', threshold=3, zero=0, units='species', source='NVI', name='Red-listed species'),
    dict(branch='skyddade', main='AV_HOT', w=1, fn='norm', threshold=3, zero=0, units='species', source='NVI', name='Protected species (ASF)'),
    dict(branch='akustik', main='AV_SENS', w=1, fn='norm', threshold=25, zero=5, units='bird species/point', source='AUDIO', name='Acoustic species richness', hypothetical=True, footprint=True),
    dict(branch='edna', main='AV_SENS', w=1, fn='norm', threshold=400, zero=50, units='OTUs/sample', source='EDNA', name='eDNA richness', hypothetical=True),
    dict(branch='struktur', main='LM_STRUK', w=1, fn='norm', threshold=8, zero=0, units='m (SD canopy height)', source='DRONE', name='Structural variation', hypothetical=True, footprint=True),
    dict(branch='kanteffekt', main='LM_STRUK', w=1, fn='norm', threshold=0.3, zero=1, units='share within 30 m of edge', source='GEO', name='Core area (inverted edge effect)'),
    dict(branch='konnektivitet', main='LM_FUNK', w=1, fn='norm', threshold=80, zero=20, units='% natural land within 500 m', source='SAT', name='Connectivity', hypothetical=True),
    dict(branch='markkol', main='LM_FUNK', w=1, fn='norm', threshold=8, zero=1, units='% SOC 0–20 cm', source='SOIL', name='Soil carbon', hypothetical=True),
]
BY = {s['branch']: s for s in STRUCTURE}
CHILDREN = {}
for s in STRUCTURE:
    CHILDREN.setdefault(s['main'], []).append(s['branch'])
INDICATORS = [s['branch'] for s in STRUCTURE if s['fn'] == 'norm']

# Tryckmatris (OHI pressures matrix, CalculatePressuresAll): rank 0-3 per kategori anger hur hårt ett tryck
# slår mot just den kategorin. tryck_c = min(1, sum(rank_c * värde) / 3). 0 = trycket påverkar inte kategorin.
PRESSURES = [
    dict(layer='p_invasiva', name='Invasive species', ranks=dict(BV=2, AV=3, LM=1), source='NVI'),
    dict(layer='p_avverkning', name='Felling / planting', ranks=dict(BV=3, AV=2, LM=2), source='NVI'),
    dict(layer='p_dikning', name='Drainage / machine damage', ranks=dict(BV=3, AV=2, LM=2), source='NVI'),
    dict(layer='p_godsling', name='Fertiliser impact', ranks=dict(BV=2, AV=3, LM=1), source='NVI'),
    dict(layer='p_igenvaxning', name='Overgrowth', ranks=dict(BV=1, AV=2, LM=1), source='NVI'),
    dict(layer='p_exploatering', name='Buildings / roads', ranks=dict(BV=3, AV=3, LM=3), source='NVI'),
]
# Resiliensmatris (CalculateResilienceAll): vikt per kategori, resiliens_c = viktat medel. 0 = ej relevant.
RESILIENCE = [
    dict(layer='r_skotsel', name='Active conservation management (BDU area 2026–2031)', weights=dict(BV=2, AV=2, LM=1), source='CONTRACT'),
    dict(layer='r_skydd', name='Key habitat / Natura 2000 status', weights=dict(BV=1, AV=1, LM=1), source='NVI'),
    dict(layer='r_konnektivitet', name='Connectivity', weights=dict(BV=0, AV=1, LM=2), source='SAT', hypothetical=True),
]
CATS = [s['branch'] for s in STRUCTURE if s['main'] == 'NVX']
DIMS = ('status', 'pressure', 'resilience', 'trend', 'future', 'score')

KLASS_NAMN = {1: 'Highest', 2: 'High', 3: 'Considerable', 4: 'Some', 5: 'No nature value'}


# ── 2. Branch-funktioner (BranchFunctions.R) ────────────────────────────────────────────────────────
def normalize(value, threshold, zero):
    return min(1, max(0, (value - zero) / (threshold - zero)))


def power_mean(values, weights, p):
    w = sum(weights)
    if p == 0:
        return math.exp(sum(wi * math.log(max(v, 1e-6)) for v, wi in zip(values, weights)) / w)
    return (sum(wi * v ** p for v, wi in zip(values, weights)) / w) ** (1 / p)


def clamp(v, lo=0, hi=1):
    return min(hi, max(lo, v))


def calibrate_cuts(pairs):
    """Välj klassgränser (NVX 0-100 -> klass) som maximerar överensstämmelsen med expertklassen.
    Girig sökning uppifrån: gränsen för klass k väljs bland mittpunkterna mellan sorterade poäng."""
    pairs = sorted(pairs, key=lambda t: -t[0])
    scores = [s for s, _ in pairs]
    cands = sorted({round((a + b) / 2, 1) for a, b in zip(scores, scores[1:])} | {100, 0}, reverse=True)
    cuts, lo = [], 101
    for k in (1, 2, 3, 4):
        best = None
        for c in [x for x in cands if x < lo]:
            agree = sum(1 for s, e in pairs if (e == k) == (c <= s < lo))
            if best is None or agree > best[0]:
                best = (agree, c)
        cuts.append((best[1], k)); lo = best[1]
    return cuts


def to_class(score, cuts):
    if score is None:
        return None
    for cut, k in cuts:
        if score >= cut:
            return k
    return 5


# ── 3. Rumsliga lager ───────────────────────────────────────────────────────────────────────────────
def load_zones():
    """Rapportzoner = naturvärdesbiotoper. NVI 2026 (Salmonic) ersätter 2024 där de överlappar."""
    feats = json.load(open(LAYERS / 'nvi_biotopes.json'))['features']
    feats += json.load(open(LAYERS / 'nvi_2026_salmonic.geojson'))['features']
    newer = [tm(f['geometry']) for f in feats if f['properties'].get('Utf_org') == 'Salmonic AB']
    zones, n_ej = [], 0
    for f in feats:
        p = f['properties']
        oid = p['Objektid']
        if oid == 'Ej naturvärde':
            n_ej += 1
            oid = f'E{n_ej:02d}'
        g = tm(f['geometry'])
        if p.get('Utf_org') != 'Salmonic AB':
            for n in newer:
                g = g.difference(n)
        if g.area > 100:
            zones.append(dict(id=oid, nvklass=int(p['Nvklass']) or 5, biotop=p.get('Biotop') or '',
                              naturtyp=p.get('Naturtyp') or '', natura=p.get('Natura_kod'), org=p.get('Utf_org'), geom=g))
    return zones


def load_units(zones):
    """Bedömningsenheter (SpatialAssessmentUnits): fyra rumsliga nivåer ovanför biotopen, plus två
    grupperingar (naturtyp, expertklass) som också aggregeras areaviktat."""
    units = []
    for key, level, name in (('naturtyp', 'naturtyp', lambda k: k), ('nvklass', 'klass', lambda k: f'Expert class {k}' if k < 5 else 'No nature value')):
        groups = {}
        for z in zones:
            k = z[key].split(',')[0].strip().capitalize() if key == 'naturtyp' else z[key]
            groups.setdefault(k or 'Unknown', []).append(z['geom'])
        for k, gs in sorted(groups.items(), key=lambda t: str(t[0])):
            units.append(dict(id=f"{level}:{k}", level=level, name=name(k), geom=unary_union(gs)))
    for f in json.load(open(LAYERS / 'bdu_rutor_ustorp.geojson'))['features']:
        units.append(dict(id=f['properties']['bdu_id'], level='ruta', name=f['properties']['bdu_id'], geom=tm(f['geometry'])))
    area = tm(json.load(open(LAYERS / 'bdu_omrade_ustorp.geojson'))['features'][0]['geometry'])
    units.append(dict(id='BDU', level='omrade', name='BDU area', geom=area))
    skiften = []
    for f in json.load(open(LAYERS / 'ustorp_property_borders.json'))['features']:
        g = tm(f['geometry'])
        skiften.append(g)
        units.append(dict(id=f['properties']['name'].replace(' ', '_'), level='skifte', name=f['properties']['name'].title().replace('Skifte', 'Parcel'), geom=g))
    units.append(dict(id='USTORP', level='fastighet', name='Property Eksjö Ustorp 1:6', geom=unary_union(skiften)))
    return units, area


def intersections(zones, units):
    """SpatialIntersections: fraction_of_zone och fraction_of_assess_unit för varje (zon, enhet)."""
    rows = []
    for u in units:
        for z in zones:
            a = z['geom'].intersection(u['geom']).area
            if a > 1:
                rows.append(dict(zone_id=z['id'], assess_unit_id=u['id'], area=a,
                                 fraction_of_zone=a / z['geom'].area, fraction_of_assess_unit=a / u['geom'].area))
    return rows


# ── 4. Data per zon: verkliga NVI-lager + geometri + hypotetiska testlager ───────────────────────────
def field_indicators(raw, zone, area_ha):
    """Indikatorvärden ur objektsbladen. None = saknas (gapfylls)."""
    if raw is None:
        return {k: None for k in INDICATORS if BY[k]['source'] == 'NVI'} | {'natura': 1 if zone['natura'] else 0}
    sp = raw['species']
    return dict(
        naturlighet=raw['bv_naturlighet'], dodved=raw['deadwood_m3ha'], naturvardestrad=raw['old_trees'] / max(area_ha, 0.1),
        sallsynthet=raw['bv_sallsynthet'], natura=1 if zone['natura'] else 0, betydelse=raw['bv_betydelse'],
        nyckelbiotop=1 if raw['nyckelbiotop'] else 0,
        signalsumma=sum(s['signal'] for s in sp.values()), artantal=len(sp),
        rodlistade=sum(1 for s in sp.values() if s['redlist']), skyddade=sum(1 for s in sp.values() if s['protected']),
    )


def pressure_layers(raw, zone):
    """Tryckvärden 0-1 ur 'Påverkan:' i objektsbladen och biotoptypen. None = saknas."""
    bio = zone['biotop'].lower()
    expl = 1 if any(k in bio for k in ('tomt', 'väg', 'byggnad', 'gräsmatta')) else 0
    if raw is None:
        if zone['nvklass'] == 5:  # hyggen, tomtmark, väg: inga objektsblad men känd påverkan
            return dict(p_invasiva=0, p_avverkning=1 if any(k in bio for k in ('hygge', 'plantering', 'produktion')) else 0,
                        p_dikning=0, p_godsling=0, p_igenvaxning=0, p_exploatering=expl)
        return dict(p_invasiva=None, p_avverkning=None, p_dikning=None, p_godsling=None, p_igenvaxning=None, p_exploatering=expl)
    pv = ' '.join(raw['paverkan']).lower()
    return dict(
        p_invasiva=1 if raw['invasive'] else 0,
        p_avverkning=1 if any(k in pv for k in ('kalavverkat', 'planterat')) else 0.5 if any(k in pv for k in ('plockhugget', 'första generationen', 'röjt')) else 0,
        p_dikning=1 if any(k in pv for k in ('dikat', 'dränerande', 'körskador')) else 0,
        p_godsling=1 if 'tydlig gödsel' in pv else 0.5 if 'svag gödsel' in pv else 0,
        p_igenvaxning=1 if 'igenväxning' in pv else 0,
        p_exploatering=expl,
    )


def hypothetical_layers(zones):
    """HYPOTETISKA testlager, slumpade med frö så att de är reproducerbara. Byts ut mot riktiga
    leveranser (Conifer Vision hösten 2026, eDNA våren 2027, Sentinel-2, jordprover, ljudanalysen)."""
    rnd = random.Random(2026)
    forest = lambda z: 'skog' in (z['naturtyp'] + z['biotop']).lower()
    wet = lambda z: any(k in z['biotop'].lower() for k in ('sump', 'kärr', 'bäck', 'alkärr'))
    hyp = {}
    for z in zones:
        k = z['nvklass']
        base = {1: 7, 2: 6, 3: 4, 4: 2.5, 5: 1.5}[k]
        hyp[z['id']] = dict(
            konnektivitet=round(min(100, max(0, rnd.gauss(70 if k <= 3 else 50, 12)))),
            markkol=round(max(0.5, rnd.gauss(12 if wet(z) else (4 if forest(z) else 3), 1.0) * (1.1 if k <= 2 else 1)), 1))
    # Ljudpunkter: sex AudioMoth-punkter från fältarbetet i juli 2025. Lägena är VERKLIGA: medelpunkt av de
    # tre panoramaplaceringarna per punkt i fältmarkeringarna (OsmAnd-GPX i Drive), kontrollerade mot
    # drönarblockens KML. Artantalen är däremot påhittade. N88 och PSU hade trasiga enheter -> NA.
    pts = []
    for name, lat, lon, val in [('DN9', 57.54182, 15.18240, 22), ('JRE', 57.54295, 15.17790, 17), ('PSU', 57.53997, 15.17640, None),
                                ('UE8', 57.53888, 15.18045, 31), ('NMN', 57.53508, 15.17803, 28), ('N88', 57.53259, 15.18528, None)]:
        x, y = to_tm(lon, lat)
        pts.append(dict(id=name, layer='akustik', value=val, x=x, y=y, lat=lat, lon=lon, radius=150,
                        geom=Point(x, y).buffer(150)))
    footprints = pts + drone_footprints(rnd)
    return hyp, footprints


def footprints_to_zones(footprints, zones, layer):
    """Fotavtryck (punkt med sökradie, drönarblock) -> zon. Motsvarar BHI:s SpatialIntersections för lager som
    rapporteras på en annan geometri än zonerna: värdet är medel över fotavtrycken, viktat med skärningsarean
    (= fraction_of_assess_unit * zonarea). Ett fotavtryck som täcker flera zoner informerar var och en i
    proportion till hur mycket av zonen det täcker. Aggregering uppåt sker sedan bara på zonareor, så samma
    mätning räknas aldrig dubbelt på BDU- eller fastighetsnivå. Täckningen (andel av zonen under något
    fotavtryck) sparas som cov_<lager>."""
    vals, cov = {}, {}
    for z in zones:
        parts = [(f['value'], f['geom'].intersection(z['geom']).area) for f in footprints if f['layer'] == layer]
        parts = [(v, a) for v, a in parts if a > 1]
        covered = unary_union([f['geom'] for f in footprints if f['layer'] == layer]).intersection(z['geom']).area
        cov[z['id']] = round(covered / z['geom'].area, 3)
        with_val = [(v, a) for v, a in parts if v is not None]
        vals[z['id']] = round(sum(v * a for v, a in with_val) / sum(a for _, a in with_val), 1) if with_val else None
    return vals, cov


# Drönarblock (1 ha) från fältarbetet i juli 2025, hörn ur KML-filerna i Drive. VERKLIGA lägen; värdet
# (strukturell variation) är HYPOTETISKT tills Conifer Visions analys finns. NMN saknar drönarblock.
DRONE_BLOCKS = {
    'DN9': [(15.1832289, 57.5422703), (15.1832244, 57.5413720), (15.1815539, 57.5413744), (15.1815584, 57.5422727)],
    'JRE': [(15.1788084, 57.5434023), (15.1788040, 57.5425041), (15.1771335, 57.5425064), (15.1771378, 57.5434047)],
    'N88': [(15.1863314, 57.5331697), (15.1863268, 57.5322715), (15.1846568, 57.5322739), (15.1846613, 57.5331722)],
    'PSU': [(15.1774898, 57.5403165), (15.1774854, 57.5394183), (15.1758151, 57.5394206), (15.1758194, 57.5403189)],
    'UE8': [(15.1811928, 57.5391383), (15.1811884, 57.5382400), (15.1795180, 57.5382424), (15.1795224, 57.5391407)],
}


def drone_footprints(rnd):
    from shapely.geometry import Polygon
    out = []
    for name, ring in DRONE_BLOCKS.items():
        g = Polygon([to_tm(lon, lat) for lon, lat in ring])
        out.append(dict(id=name, layer='struktur', value=round(max(0, rnd.gauss(4.5, 2.0)), 1), geom=g))
    return out


# ── 5. Gapfilling (rumslig variant av gapfill_data) ─────────────────────────────────────────────────
def gapfill(values, zones, keys):
    """GF-S01: saknat värde ersätts med medianen för zoner med samma expertklass. Returnerar flaggor."""
    byid = {z['id']: z for z in zones}
    flags = {}
    for name in keys:
        for k in range(1, 6):
            vals = sorted(v[name] for zid, v in values.items() if v.get(name) is not None and byid[zid]['nvklass'] == k)
            med = vals[len(vals) // 2] if vals else None
            for zid, v in values.items():
                if byid[zid]['nvklass'] == k and v.get(name) is None:
                    if med is None:
                        flags.setdefault(zid, {})[name] = 'missing'
                    else:
                        v[name] = med
                        flags.setdefault(zid, {})[name] = 'gapfilled'
    return flags


# ── 6. Status per zon (calculate_status / weighted_power_mean) ──────────────────────────────────────
def score_tree(ind, flags):
    scores, fill, miss = {}, {}, {}

    def rec(b):
        s = BY[b]
        if s['fn'] == 'norm':
            if ind.get(b) is None:
                scores[b], fill[b], miss[b] = None, 0, 1
            else:
                scores[b] = normalize(ind[b], s['threshold'], s['zero'])
                fill[b], miss[b] = (1 if flags.get(b) == 'gapfilled' else 0), 0
            return
        kids = CHILDREN[b]
        for k in kids:
            rec(k)
        avail = [(scores[k], BY[k]['w']) for k in kids if scores[k] is not None]
        wsum = sum(BY[k]['w'] for k in kids)
        scores[b] = power_mean([v for v, _ in avail], [w for _, w in avail], s['exponent']) if avail else None
        fill[b] = sum(fill[k] * BY[k]['w'] for k in kids) / wsum
        miss[b] = sum(miss[k] * BY[k]['w'] for k in kids) / wsum

    rec('NVX')
    return scores, fill, miss


def pressure_score(pl, cat):
    """CalculatePressuresAll: sum(rank_c*värde)/3, tak 1. Saknade lager räknas som 0 (na.rm)."""
    return clamp(sum(p['ranks'][cat] * (pl.get(p['layer']) or 0) for p in PRESSURES) / 3)


def resilience_score(rl, cat):
    avail = [(rl[r['layer']], r['weights'][cat]) for r in RESILIENCE if rl.get(r['layer']) is not None and r['weights'][cat]]
    return sum(v * w for v, w in avail) / sum(w for _, w in avail) if avail else 0


def goal_index(status, trend, pressure, resilience):
    """ohicore CalculateGoalIndex: trolig framtid och poäng för ett mål."""
    r = min(resilience, pressure)  # resiliens kan inte överstiga trycket
    xF = clamp(DISCOUNT * (1 + BETA * trend + (1 - BETA) * (r - pressure)) * status)
    return xF, (status + xF) / 2


def all_dims(sc, v, trend):
    """Dimensioner per kategori (mål) och för NVX som helhet. NVX-nivån är kategoriernas viktade medel,
    som index = weighted.mean(goal score) i CalculateAll; status är dock power mean (p=0.5) som i trädet."""
    cat = {}
    for c in CATS:
        st = sc[c]
        if st is None:
            continue
        p, r = pressure_score(v, c), resilience_score(v, c)
        xF, score = goal_index(st, trend, p, r)
        cat[c] = dict(status=st, pressure=p, resilience=r, trend=trend, future=xF, score=score)
    w = {c: BY[c]['w'] for c in cat}
    top = {d: sum(cat[c][d] * w[c] for c in cat) / sum(w.values()) for d in DIMS}
    top['status'] = sc['NVX']
    return top, cat


# ── 7. Rumslig aggregering (CalculateAll: weighted.mean(score, area) per mål och dimension) ──────────
def aggregate(unit, zone_rows, inter):
    rows = [(r, zone_rows[r['zone_id']]) for r in inter if r['assess_unit_id'] == unit['id']]
    tot = sum(r['area'] for r, _ in rows)
    if not tot:
        return None

    def wmean(get):
        vals = [(get(z), r['area']) for r, z in rows if get(z) is not None]
        return round(sum(v * a for v, a in vals) / sum(a for _, a in vals), 2) if vals else None

    layer_keys = INDICATORS + [p['layer'] for p in PRESSURES] + [r['layer'] for r in RESILIENCE] + ['cov_akustik', 'cov_struktur']
    return dict(
        id=unit['id'], level=unit['level'], name=unit['name'], area_ha=round(unit['geom'].area / 1e4, 2),
        coverage=round(tot / unit['geom'].area, 3),
        scores={b: wmean(lambda z: z['scores'][b]) for b in BY},
        dims={d: wmean(lambda z: z['dims'][d]) for d in DIMS},
        cat_dims={c: {d: wmean(lambda z: z['cat_dims'].get(c, {}).get(d)) for d in DIMS} for c in CATS},
        layers={l: wmean(lambda z: z['layers'][l]) for l in layer_keys},
        expert_klass=wmean(lambda z: z['expert_klass']),
        fraction_filled=wmean(lambda z: z['fraction_filled']), fraction_missing=wmean(lambda z: z['fraction_missing']),
        parts=sorted([dict(id=z['id'], area_ha=round(r['area'] / 1e4, 3), fraction=round(r['area'] / tot, 4), fraction_of_zone=round(r['fraction_of_zone'], 3),
                           status=z['dims']['status'], klass=z['nvx_klass']) for r, z in rows], key=lambda p: -p['fraction']),
        geom=mapping(unit['geom']))


# ── 8. Kör ───────────────────────────────────────────────────────────────────────────────────────────
zones = load_zones()
units, area = load_units(zones)
inter = intersections(zones, units)
raw = {r['id']: r for r in json.load(open(HERE / 'data/biotopes_raw.json'))}
hyp, footprints = hypothetical_layers(zones)
points = [f for f in footprints if f['layer'] == 'akustik']
akustik, cov_akustik = footprints_to_zones(footprints, zones, 'akustik')
struktur, cov_struktur = footprints_to_zones(footprints, zones, 'struktur')
in_bdu = {r['zone_id']: r['fraction_of_zone'] for r in inter if r['assess_unit_id'] == 'BDU'}

values = {}
for z in zones:
    ha = z['geom'].area / 1e4
    r = raw.get(z['id'])
    if z['nvklass'] == 5 and r is None:
        # Ej naturvärde (hyggen, tomtmark, väg): inga objektsblad. SIS-minimum, inte gapfyllt.
        v = {k: 0 for k in INDICATORS if BY[k]['source'] == 'NVI'}
    else:
        v = field_indicators(r, z, ha)
    v |= hyp[z['id']]
    v['akustik'] = akustik[z['id']]
    v['struktur'] = struktur[z['id']]
    v['cov_akustik'], v['cov_struktur'] = cov_akustik[z['id']], cov_struktur[z['id']]
    v['edna'] = None  # levereras våren 2027
    v['kanteffekt'] = round(1 - z['geom'].buffer(-30).area / z['geom'].area, 2)  # ur geometrin
    v |= pressure_layers(r, z)
    v['r_skotsel'] = round(in_bdu.get(z['id'], 0), 2)
    v['r_skydd'] = 1 if (r and r['nyckelbiotop']) else 0.5 if z['natura'] else 0
    v['r_konnektivitet'] = v['konnektivitet'] / 100
    values[z['id']] = v

flags = gapfill(values, zones, INDICATORS + [p['layer'] for p in PRESSURES])

zone_rows = {}
for z in zones:
    v = values[z['id']]
    sc, fill, miss = score_tree(v, flags.get(z['id'], {}))
    top, cat = all_dims(sc, v, 0.0)  # trend 0: ingen tidsserie ännu
    zone_rows[z['id']] = dict(
        id=z['id'], level='biotop', name=z['id'], biotop=z['biotop'], naturtyp=z['naturtyp'], org=z['org'],
        area_ha=round(z['geom'].area / 1e4, 2), expert_klass=z['nvklass'],
        scores={k: (round(s * 100, 1) if s is not None else None) for k, s in sc.items()},
        dims={d: round(top[d] * 100, 1) for d in DIMS},
        cat_dims={c: {d: round(cat[c][d] * 100, 1) for d in DIMS} for c in cat},
        layers=v, flags=flags.get(z['id'], {}),
        fraction_filled=round(fill['NVX'], 3), fraction_missing=round(miss['NVX'], 3), geom=mapping(z['geom']))

cuts = calibrate_cuts([(r['dims']['status'], r['expert_klass']) for r in zone_rows.values()])
for r in zone_rows.values():
    r['nvx_klass'] = to_class(r['dims']['status'], cuts)

agg = [a for a in (aggregate(u, zone_rows, inter) for u in units) if a]
for a in agg:
    a['nvx_klass'] = to_class(a['dims']['status'], cuts)

LEVELS = [dict(id='biotop', name='Biotope', spatial=True), dict(id='ruta', name='BDU plot', spatial=True),
          dict(id='omrade', name='BDU area', spatial=True), dict(id='skifte', name='Parcel', spatial=True),
          dict(id='fastighet', name='Property', spatial=True),
          dict(id='naturtyp', name='Habitat type', spatial=False), dict(id='klass', name='Expert class', spatial=False)]
out = dict(year=YEAR, beta=BETA, discount=DISCOUNT, cats=CATS, structure=STRUCTURE, pressures=PRESSURES, resilience=RESILIENCE, class_cuts=cuts,
           levels=LEVELS, units=list(zone_rows.values()) + agg,
           points=[{k: v for k, v in q.items() if k != 'geom'} for q in points],
           footprints=[dict(id=f['id'], layer=f['layer'], value=f['value'], geom=mapping(f['geom'])) for f in footprints if f['layer'] == 'struktur'],
           intersections=[{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in inter])
(HERE / 'data/index_scores.json').write_text(json.dumps(out, ensure_ascii=False))

# ── 9. Rapport ───────────────────────────────────────────────────────────────────────────────────────
rows = list(zone_rows.values())
agree = sum(1 for r in rows if r['expert_klass'] == r['nvx_klass'])
print(f"klassgränser {cuts}; överensstämmelse {agree}/{len(rows)}; ±1: {sum(1 for r in rows if abs(r['expert_klass'] - r['nvx_klass']) <= 1)}/{len(rows)}\n")
print(f"{'enhet':>10} {'ha':>6} {'täck':>5} {'exp':>4} {'status':>6} {'kl':>2} {'tryck':>5} {'resil':>5} {'framt':>5} {'poäng':>5}")
for r in rows + agg:
    if (r['level'] == 'biotop' and r['id'].startswith('E')) or r['level'] in ('naturtyp', 'klass'):
        continue
    d = r['dims']
    print(f"{r['id']:>10} {r['area_ha']:>6} {r.get('coverage', 1):>5} {r['expert_klass']:>4.1f} {d['status']:>6} {r['nvx_klass']:>2} "
          f"{d['pressure']:>5.0f} {d['resilience']:>5.0f} {d['future']:>5.0f} {d['score']:>5.0f}")

for q in points:
    inside = [u['name'] for u in units if u['level'] == 'skifte' and u['geom'].contains(Point(q['x'], q['y']))]
    print(f"{q['id']}: {inside or 'utanför fastigheten'}")

# ── 10. Visualisering ───────────────────────────────────────────────────────────────────────────────
tpl = (HERE / 'nvx_template.html').read_text()
js = json.dumps(out, ensure_ascii=False).replace('</', '<\\/')
(HERE / 'nvx.html').write_text(tpl.replace('/*DATA*/', js))
print(f"\nskrev index/nvx.html ({len(js) // 1024} kB data)")
