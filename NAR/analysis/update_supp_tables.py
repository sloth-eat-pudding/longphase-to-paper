"""Write the revision (a032938) LongPhase-TO numbers into the Supplementary tables of NAR/supplementary.tex.

Only LongPhase-TO values change. Raw-caller rows, ASCAT and PURPLE columns keep the published values.
  - Somatic tables (ONT SNV/indel, PacBio SNV/indel): rows 'ClairS-TO+LongPhase-TO' and
    'DeepSomatic-TO+LongPhase-TO' (= DeepSomatic-TO + ClairS-TO PoN on both platforms).
  - Per-mixture TDF tables: columns LP-TO (C) and LP-TO (D) (D = DeepSomatic-TO + ClairS-TO PoN).
  - tab:tdf-mae-by-fraction and tab:tdf-accuracy-summary: the LongPhase-TO rows are recomputed from the
    per-mixture tables (expected fraction = the table's Expected column), with the same definitions.
Source: ../../Revision_2026-11/. Run: python update_supp_tables.py  (edits ../supplementary.tex in place).
"""
import csv, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.join(HERE, '..', '..', 'Revision_2026-11')
SUPP = os.path.join(HERE, '..', 'supplementary.tex')
LP_C, LP_D = 'LongPhase-TO(ClairS-TO)', 'LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)'


def tsv(name):
    return list(csv.DictReader((l for l in open(os.path.join(REV, name)) if not l.startswith('#')), delimiter='\t'))


def ds_name(platform, d):
    if platform == 'ONT':
        return 'HCC1937_UCSC_with_duplicate_reads' if d == 'HCC1937_UCSC' else d
    return d.replace('_CASTLE', '_CMH')


s = open(SUPP).read()


def block(label):
    i = s.index('\\label{%s}' % label)
    j = s.index('\\end{longtable}', i) if 'longtable' in s[i - 400:i] or label.startswith('tab:pacbio_') or label in (
        'tab:somatic_purity', 'tab:somatic_indel_purity', 'tab:tumor-purity', 'tab:pacbio-tumor-purity') else s.index('\\end{tabular}', i)
    return i, j


# ---- somatic tables ----
M = {(r['platform'], r['dataset'], float(r['target_tdf']), r['variant_type'], r['method']): r
     for r in tsv('metrics_snv_indel_ONT_PacBio.tsv')}
n_som = 0
for label, plat, vt in [('tab:somatic_purity', 'ONT', 'SNV'), ('tab:somatic_indel_purity', 'ONT', 'INDEL'),
                        ('tab:pacbio_somatic_snv', 'PacBio', 'SNV'), ('tab:pacbio_somatic_indel', 'PacBio', 'INDEL')]:
    i, j = block(label)
    lines = s[i:j].split('\n')
    for k, l in enumerate(lines):
        c = [x.strip() for x in l.rstrip().rstrip('\\').split('&')]
        if len(c) == 6 and c[2] in ('ClairS-TO+LongPhase-TO', 'DeepSomatic-TO+LongPhase-TO'):
            r = M[(plat, ds_name(plat, c[0].replace('\\_', '_')), float(c[1]), vt, LP_C if c[2].startswith('ClairS') else LP_D)]
            lines[k] = '  ' + ' & '.join(c[:3] + ['%.3f' % float(r[x]) for x in ('precision', 'recall', 'f1')]) + ' \\\\'
            n_som += 1
    s = s[:i] + '\n'.join(lines) + s[j:]

# ---- per-mixture TDF tables ----
TD = {'ONT': tsv('tdf_prediction_all_methods_ont_plot_values.tsv'), 'PacBio': tsv('tdf_prediction_all_methods_pacbio_plot_values.tsv')}
TMAP = {'ONT': {'COLO829_NYGC': 'COLO829_R10', 'H1437_UCSC': 'H1437', 'H2009_UCSC': 'H2009', 'HCC1937_UCSC': 'HCC1937', 'HCC1954_UCSC': 'HCC1954'},
        'PacBio': {'HCC1395_PacBio': 'HCC1395', 'COLO829_PacBio': 'COLO829', 'HCC1954_CASTLE': 'HCC1954', 'HCC1937_CASTLE': 'HCC1937',
                   'H2009_CASTLE': 'H2009', 'H1437_CASTLE': 'H1437'}}
mix = {'ONT': [], 'PacBio': []}   # (dataset, expected, {method: estimate})
for label, plat in [('tab:tumor-purity', 'ONT'), ('tab:pacbio-tumor-purity', 'PacBio')]:
    i, j = block(label)
    lines = s[i:j].split('\n')
    for k, l in enumerate(lines):
        c = [x.strip() for x in l.rstrip().rstrip('\\').split('&')]
        if len(c) == 7 and re.match(r'^[0-9.]+$', c[1]):
            d = c[0].replace('\\_', '_'); exp = float(c[1])
            cand = [x for x in TD[plat] if x['dataset'] == TMAP[plat].get(d, d)]
            r = min(cand, key=lambda x: abs(float(x['simulated_tdf']) - exp))
            c[5] = '%.3f' % float(r['clairs_to_longphase_to_tdf']); c[6] = '%.3f' % float(r['deepsomatic_to_pon_longphase_to_tdf'])
            lines[k] = '  ' + ' & '.join(c) + ' \\\\'
            mix[plat].append((d, exp, {'C': float(c[5]), 'D': float(c[6])}))
    s = s[:i] + '\n'.join(lines) + s[j:]


# ---- summary statistics for the LongPhase-TO rows ----
def stats(rows, key):
    x = [e for _, e, _ in rows]; y = [m[key] for _, _, m in rows]; n = len(x)
    err = [abs(a - b) for a, b in zip(x, y)]
    mx, my = sum(x) / n, sum(y) / n
    vx = sum((a - mx) ** 2 for a in x) / n; vy = sum((b - my) ** 2 for b in y) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(x, y)) / n
    ccc = 2 * cov / (vx + vy + (mx - my) ** 2); slope = cov / vx; icpt = my - slope * mx
    ds = sorted({d for d, _, _ in rows})
    jk = [sum(abs(e - m[key]) for d, e, m in rows if d != o) / sum(1 for d, _, _ in rows if d != o) for o in ds]
    lv = {}
    for (d, e, m), er in zip(rows, err):
        lv.setdefault(round(e * 5) / 5, []).append(er)
    return dict(mae=sum(err) / n, jk=(min(jk), max(jk)), ccc=ccc, slope=slope, icpt=icpt,
                lv=[sum(lv[v]) / len(lv[v]) for v in (0.2, 0.4, 0.6, 0.8, 1.0)])


def sgn(v):
    v = round(v, 2)
    return '0.00' if v == 0 else ('$+$%.2f' % v if v > 0 else '$-$%.2f' % -v)


for label in ('tab:tdf-mae-by-fraction', 'tab:tdf-accuracy-summary'):
    i, j = block(label)
    body = s[i:j]
    parts = body.split('\\multicolumn')
    for p_idx, plat in [(1, 'ONT'), (2, 'PacBio')]:
        for key, name in [('C', 'LongPhase-TO (ClairS-TO)'), ('D', 'LongPhase-TO (DeepSomatic-TO)')]:
            st = stats(mix[plat], key)
            if label == 'tab:tdf-mae-by-fraction':
                new = name + ' & ' + ' & '.join('%.2f' % v for v in st['lv']) + ' & %.2f \\\\' % st['mae']
            else:
                new = (name + ' & %.2f & [%.2f, %.2f] & %.2f & %.2f & %s & 0 \\\\'
                       % (st['mae'], st['jk'][0], st['jk'][1], st['ccc'], st['slope'], sgn(st['icpt'])))
            parts[p_idx] = re.sub(re.escape(name) + r' &[^\n]*\\\\', lambda _: new, parts[p_idx], count=1)
    s = s[:i] + '\\multicolumn'.join(parts) + s[j:]

open(SUPP, 'w').write(s)
print('somatic rows updated:', n_som, '| TDF mixtures:', {k: len(v) for k, v in mix.items()})
