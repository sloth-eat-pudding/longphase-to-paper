"""Mean precision/recall/F1 per platform, variant type and method, with and without the HCC1395 datasets
(the triplet-graph models were fitted on HCC1395_HKU at 0.6). Source: Revision_2026-11/metrics_snv_indel_ONT_PacBio.tsv.
HCC1937_UCSC_dedup is left out (the paper uses the with-duplicate-reads set). Writes hcc1395_excluded.tsv."""
import csv, os
from collections import defaultdict
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'Revision_2026-11', 'metrics_snv_indel_ONT_PacBio.tsv')
rows = [r for r in csv.DictReader(open(SRC), delimiter='\t') if r['dataset'] != 'HCC1937_UCSC_dedup']
BASE = {'LongPhase-TO(ClairS-TO)': 'ClairS-TO raw PASS', 'LongPhase-TO(DeepSomatic-TO)': 'DeepSomatic-TO raw PASS',
        'LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)': 'DeepSomatic-TO raw PASS'}
out = []
for subset, keep in [('all', lambda d: True), ('excluding HCC1395', lambda d: not d.startswith('HCC1395'))]:
    agg = defaultdict(list)
    for r in rows:
        if keep(r['dataset']):
            agg[(r['platform'], r['variant_type'], r['method'])].append(r)
    for (p, v, m), rs in sorted(agg.items()):
        mean = lambda k: sum(float(x[k]) for x in rs) / len(rs)
        rec = dict(subset=subset, platform=p, variant_type=v, method=m, n=len(rs),
                   precision=round(mean('precision'), 4), recall=round(mean('recall'), 4), f1=round(mean('f1'), 4))
        if m in BASE:
            b = {(x['dataset'], x['mixture']): float(x['f1']) for x in rows
                 if keep(x['dataset']) and x['platform'] == p and x['variant_type'] == v and x['method'] == BASE[m]}
            d = [float(x['f1']) - b[(x['dataset'], x['mixture'])] for x in rs]
            rec.update(delta_f1=round(sum(d) / len(d), 4), n_improved=sum(x > 0 for x in d), n_drop_gt_0_01=sum(x < -0.01 for x in d))
        out.append(rec)
cols = ['subset', 'platform', 'variant_type', 'method', 'n', 'precision', 'recall', 'f1', 'delta_f1', 'n_improved', 'n_drop_gt_0_01']
with open(os.path.join(HERE, 'hcc1395_excluded.tsv'), 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', restval=''); w.writeheader(); w.writerows(out)
