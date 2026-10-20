# Läser objektsbladen i NVI 2024 (Ekologigruppen, textexport från Drive) och NVI 2026 (Salmonic, utkast)
# och skriver rådata per naturvärdesbiotop till data/biotopes_raw.json.
# Användning: python -I index/parse_nvi.py
import json, re
from pathlib import Path

HERE = Path(__file__).parent
TXT = (HERE / 'data/nvi2024_ekologigruppen.txt').read_text()

LEVELS = ['obetydligt', 'visst', 'påtagligt', 'högt', 'mycket högt']
SIGNAL = {'Mycket högt': 8, 'Högt': 4, 'Påtagligt': 2, 'Visst': 1, 'Obetydligt': 0.25}
REDLIST = {'NT': 1, 'VU': 2, 'EN': 3, 'CR': 4}

# Standardfraserna i SS 199000 för de tre biotopvärdeskriterierna, ordnade lågt -> högt (0..3).
BV_NATURLIGHET = ['Kraftigt negativt påverkad', 'Negativt påverkad biotop med viss naturlighet',
                  'endast liten negativ påverkan', 'referensförhållande eller idealtillstånd och']
BV_SALLSYNTHET = ['vanligt förekommande', 'mindre vanlig', 'ovanlig i ett', 'sällsynt eller påtagligt minskande']
BV_BETYDELSE = ['inte har någon särskild betydelse', 'viss särskild betydelse',
                'påtaglig särskild betydelse', 'hög särskild betydelse']


def ordinal(text, phrases):
    hits = [i for i, p in enumerate(phrases) if p.lower() in text.lower()]
    return max(hits) if hits else None


def sheets():
    starts = [(m.start(), int(m.group(1)), m.group(2).strip())
              for m in re.finditer(r'^(\d+)\\\. ([^\n]*?)(?=Naturvärdesklass|\n)', TXT, re.M)]
    starts = [s for s in starts if 'Naturvärdesklass:' in TXT[s[0]:s[0] + 300]]
    for i, (pos, oid, title) in enumerate(starts):
        end = starts[i + 1][0] if i + 1 < len(starts) else len(TXT)
        yield oid, title, TXT[pos:end]


def species(block):
    # Artraderna i objektsbladet: "Namn (Latin) [frekvens] Signalvärde typer"
    sec = block.split('Naturvårdsarter', 1)[-1]
    sec = sec.split('Invasiva arter', 1)[0] + ' ' + sec.split('Invasiva arter', 1)[-1][:2000]
    out = {}
    # Latinska namn är ibland omkastade i exporten: "Gullgröppa aureus) (Pseudomerulius"
    ms = list(re.finditer(r'([A-ZÅÄÖ][a-zåäö\-]+(?: [a-zåäö\-]+)?) (?:\([A-Z][a-z]+[^)]*\)|[a-z\-]+\) \([A-Z][a-z]+)', sec))
    for i, m in enumerate(ms):
        name = m.group(1).strip()
        rest = sec[m.end():ms[i + 1].start() if i + 1 < len(ms) else m.end() + 200]
        sig = next((v for k, v in SIGNAL.items() if re.search(r'\b' + k + r'\b', rest)), 0)
        red = re.search(r'\((NT|VU|EN|CR)\)', rest)
        out[name] = {'signal': sig, 'redlist': red.group(1) if red else None,
                     'protected': 'AFS' in rest or '§' in rest}
    return out


def deadwood(block):
    # Volymklasser för död ved: sällsynt <1, enstaka-sparsamt 1-5, tämligen allmänt 5-15, allmänt 15-30 m3/ha
    vol = {'Sällsynt (\\<1 m3': 0.5, '(1-5 m3': 3, '(5-15 m3': 10, '(15-30 m3': 22, '(\\>30 m3': 35}
    return round(sum(block.count(k) * v for k, v in vol.items()), 1)


def oldtrees(block):
    return len(re.findall(r'Levande \w+ (?:Mycket grov|Grov|Gammal|Nästan gammal|Mycket gammal)', block))


def paverkan(flat):
    # "Påverkan: Dikat, Kalavverkat" -> lista; används som tryckindikatorer
    m = re.search(r'Påverkan: (.*?)(?: §|Inventerad av|Värdeelement|Kontinuitet)', flat)
    return [x.strip() for x in m.group(1).split(',') if x.strip()] if m else []


def parse_2024():
    out = []
    for oid, title, b in sheets():
        flat = re.sub(r'\s+', ' ', b)
        cls = int(re.search(r'naturvärdesklass (\d)', flat).group(1))
        lv = re.search(r'bedöms ha ett ((?:mycket )?\w+) artvärde och ((?:mycket )?\w+) biotopvärde', flat)
        mbv = re.search(r'Motiv för biotopvärde:(.*?)Motiv för artvärde:', flat)
        mav = re.search(r'Motiv för artvärde:(.*?)(?:Preliminär|Naturvårdsförutsättningar|Värdeelement)', flat)
        sp = species(b)
        inv = re.search(r'Invasiva arter, negativa indikatorarter och övriga arter (.{0,60})', flat)
        out.append({
            'id': str(oid), 'title': title, 'source': 'NVI 2024 Ekologigruppen', 'nvklass': cls,
            'expert_av': lv.group(1) if lv else None, 'expert_bv': lv.group(2) if lv else None,
            'bv_naturlighet': ordinal(mbv.group(1), BV_NATURLIGHET) if mbv else None,
            'bv_sallsynthet': ordinal(mbv.group(1), BV_SALLSYNTHET) if mbv else None,
            'bv_betydelse': ordinal(mbv.group(1), BV_BETYDELSE) if mbv else None,
            'av_motiv': mav.group(1).strip() if mav else None,
            'deadwood_m3ha': deadwood(b), 'old_trees': oldtrees(b),
            'species': sp,
            'invasive': bool(inv and 'Ej noterat' not in inv.group(1) and 'lupin' in flat.lower()),
            'wet': bool(re.search(r'Markfuktighet: (Våt|Blöt|Fuktig)', flat)),
            'paverkan': paverkan(flat), 'nyckelbiotop': 'Uppfyller kvalitetskrav på nyckelbiotop' in flat,
        })
    for r in out:
        r.update(OVERRIDES.get(r['id'], {}))
    return out


# Textexporten är trasig på några ställen; värdena är lästa ur objektsbladen för hand.
OVERRIDES = {
    '19': {'bv_naturlighet': 2},  # "god naturlighet ... liten negativ påverkan", ordföljden förvanskad i exporten
}


def parse_2026():
    # LS01 och LS02 ur Salmonics utkast 2026-08-16 (Tabell 3 och objektbeskrivningen), inmatade för hand.
    ls01_species = {
        'Spillkråka': (1, None, True), 'Vanlig padda': (0.25, None, True), 'Blodrot': (0.5, None, False),
        'Blåsippa': (1, None, True), 'Blåsuga': (1.5, None, False), 'Ekbräken': (1, None, False),
        'Gökärt': (0.5, None, False), 'Hultbräken': (1, None, False), 'Kärrviol': (1, None, False),
        'Liljekonvalj': (1, None, True), 'Ormbär': (2, None, False), 'Stjärnstarr': (1, None, False),
        'Stor blåklocka': (1, None, False), 'Kransmossa': (1, None, False), 'Långfliksmossa': (0.5, None, False),
        'Skuggsprötmossa': (1, None, False), 'Stor tujamossa': (0.5, None, False)}
    return [
        {'id': 'LS01', 'title': 'Mosaikartad, delvis örtrik igenväxningsskog', 'source': 'NVI 2026 Salmonic (utkast)',
         'nvklass': 4, 'expert_av': 'visst', 'expert_bv': 'visst',
         # "Mellan bra och dåligt tillstånd", "Vanlig", "Viss ekologisk funktion"
         'bv_naturlighet': 1, 'bv_sallsynthet': 0, 'bv_betydelse': 1,
         'av_motiv': 'Måttliga förekomster av arter med visst signalvärde och sparsamma förekomster av arter med påtagligt signalvärde.',
         'deadwood_m3ha': 3, 'old_trees': 2,
         'species': {k: {'signal': s, 'redlist': r, 'protected': p} for k, (s, r, p) in ls01_species.items()},
         'invasive': False, 'wet': True, 'paverkan': ['Kalavverkat', 'Igenväxning med träd'], 'nyckelbiotop': False},
        {'id': 'LS02', 'title': 'Tät igenväxningsskog', 'source': 'NVI 2026 Salmonic (utkast)',
         'nvklass': 0, 'expert_av': 'obetydligt', 'expert_bv': 'obetydligt',
         'bv_naturlighet': 0, 'bv_sallsynthet': 0, 'bv_betydelse': 0,
         'av_motiv': 'Enstaka värdearter.', 'deadwood_m3ha': 0, 'old_trees': 0,
         'species': {}, 'invasive': False, 'wet': False, 'paverkan': ['Kalavverkat', 'Igenväxning med träd'], 'nyckelbiotop': False},
    ]


if __name__ == '__main__':
    rows = parse_2024() + parse_2026()
    (HERE / 'data/biotopes_raw.json').write_text(json.dumps(rows, ensure_ascii=False, indent=1))
    for r in rows:
        sp = r['species']
        print(f"{r['id']:>4} k{r['nvklass']} AV={r['expert_av']!s:<10} BV={r['expert_bv']!s:<10} "
              f"nat={r['bv_naturlighet']} sal={r['bv_sallsynthet']} bet={r['bv_betydelse']} "
              f"dw={r['deadwood_m3ha']:>4} old={r['old_trees']} n_sp={len(sp):>2} "
              f"sig={sum(s['signal'] for s in sp.values()):>5} red={sum(1 for s in sp.values() if s['redlist'])} "
              f"inv={r['invasive']}")
