"""Draft of Figure 4 (NAR revision): tumor DNA fraction, ploidy and HG008 copy number.

Sources (raw numbers only):
  ../../Revision_2026-11/tdf_prediction_all_methods_{ont,pacbio}_plot_values.tsv  (PURPLE, LongPhase-TO, PacBio ASCAT)
  ../supplementary.tex, tab:tumor-purity  (ONT ASCAT T+N / TO and expected fraction; the student's ONT ASCAT
                                            columns differ from the published table, so the published values are used)
  ../../Revision_2026-11/ploidy_accuracy_ONT_PacBio.tsv
  ../../Revision_2026-11/hg008_ONT_cn_loh.tsv
LongPhase-TO = DeepSomatic-TO + ClairS-TO PoN configuration. HCC1937 ONT = the with-duplicate-reads set.
Run with a Python that has numpy + matplotlib. Writes drafts/fig4_draft.png.
"""
import csv, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.join(HERE, '..', '..', 'Revision_2026-11')
SUPP = os.path.join(HERE, '..', 'supplementary.tex')

INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
LP, CMP = '#2a78d6', '#8a8985'          # LongPhase-TO = series 1 blue; comparators = neutral gray
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6,
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'legend.frameon': False})


def tsv(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith('#')), delimiter='\t'))


def supp_table(label):
    s = open(SUPP).read()
    i = s.index('\\label{%s}' % label); j = s.index('\\end{longtable}', i)
    return [[c.strip().replace('\\_', '_') for c in l.strip().rstrip('\\').split('&')]
            for l in s[i:j].splitlines() if l.count('&') >= 6 and not l.strip().startswith('Dataset')]


# ---------- TDF points: (method, platform, truth, estimate) ----------
ONT_MAP = {'COLO829_NYGC': 'COLO829_R10', 'H1437_UCSC': 'H1437', 'H2009_UCSC': 'H2009',
           'HCC1937_UCSC': 'HCC1937', 'HCC1954_UCSC': 'HCC1954'}
pts = []
ont = tsv(os.path.join(REV, 'tdf_prediction_all_methods_ont_plot_values.tsv'))
for r in supp_table('tab:tumor-purity'):
    d, exp, a_tn, a_to = r[0], float(r[1]), r[2], r[3]
    cand = [x for x in ont if x['dataset'] == ONT_MAP.get(d, d)]
    s = min(cand, key=lambda x: abs(float(x['simulated_tdf']) - exp))
    for m, v in [('ASCAT (T+N)', a_tn), ('ASCAT (TO)', a_to), ('PURPLE (TO)', s['purple_to_tdf']),
                 ('LongPhase-TO', s['deepsomatic_to_pon_longphase_to_tdf'])]:
        if v not in ('n.a.', '', 'NA'):
            pts.append((m, 'ONT', exp, float(v)))
for s in tsv(os.path.join(REV, 'tdf_prediction_all_methods_pacbio_plot_values.tsv')):
    for m, k in [('ASCAT (T+N)', 'ascat_tn_tdf'), ('ASCAT (TO)', 'ascat_to_tdf'), ('PURPLE (TO)', 'purple_to_tdf'),
                 ('LongPhase-TO', 'deepsomatic_to_pon_longphase_to_tdf')]:
        if s[k] not in ('', 'NA', 'nan'):
            pts.append((m, 'PacBio', float(s['simulated_tdf']), float(s[k])))

# ---------- ploidy ----------
pl = [r for r in tsv(os.path.join(REV, 'ploidy_accuracy_ONT_PacBio.tsv'))
      if r['method'] == 'LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)' and 'dedup' not in r['dataset']]
LINES = ['COLO829', 'H1437', 'H2009', 'HCC1395', 'HCC1937', 'HCC1954']
hg = tsv(os.path.join(REV, 'hg008_ONT_cn_loh.tsv'))

fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.05, 0.9], hspace=0.62, wspace=0.38)

# (a) small multiples, one per method
row_first = []
METHODS = ['ASCAT (T+N)', 'ASCAT (TO)', 'PURPLE (TO)', 'LongPhase-TO']
for i, m in enumerate(METHODS):
    ax = fig.add_subplot(gs[0, i])
    col = LP if m == 'LongPhase-TO' else CMP
    ax.plot([0, 1.05], [0, 1.05], color=GRID, lw=1, zorder=0)
    for plat, mk, fc in [('ONT', 'o', col), ('PacBio', 's', 'white')]:
        xy = np.array([(t, e) for mm, p, t, e in pts if mm == m and p == plat])
        ax.scatter(xy[:, 0], xy[:, 1], s=14, marker=mk, facecolor=fc, edgecolor=col, lw=0.8, zorder=2,
                   label=plat if i == 0 else None)
    err = np.array([abs(e - t) for mm, p, t, e in pts if mm == m])
    ax.text(0.98, 0.03, f'MAE {err.mean():.3f}\n{(err <= 0.05).sum()}/{len(err)} within ±0.05',
            transform=ax.transAxes, va='bottom', ha='right', fontsize=6.3, color=INK)
    ax.set_title(m, fontsize=7.5, color=INK, pad=4)
    ax.set_xlim(0.1, 1.05); ax.set_ylim(0, 1.05); ax.set_aspect('equal')
    ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel('Expected DNA fraction')
    if i == 0:
        row_first.append(ax)
        ax.set_ylabel('Estimated DNA fraction')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# (b) ploidy per cell line
ax = fig.add_subplot(gs[1, :]); row_first.append(ax)
rng = np.random.default_rng(1)
for j, cl in enumerate(LINES):
    rs = [r for r in pl if r['dataset'].startswith(cl)]
    lit = float(rs[0]['literature_ploidy'])
    ax.add_patch(plt.Rectangle((j - 0.42, lit - 0.2), 0.84, 0.4, color='#f0efec', lw=0, zorder=0))
    ax.plot([j - 0.42, j + 0.42], [lit, lit], color=INK2, lw=0.8, zorder=1)
    for off, key, col in [(-0.17, 'purple_java_ploidy', CMP), (0.17, 'ploidy', LP)]:
        for plat, mk, fc in [('ONT', 'o', col), ('PacBio', 's', 'white')]:
            v = [float(r[key]) for r in rs if r['platform'] == plat and r[key] not in ('NA', '')]
            ax.scatter(j + off + rng.uniform(-0.07, 0.07, len(v)), v, s=12, marker=mk, facecolor=fc,
                       edgecolor=col, lw=0.8, zorder=2)
ax.set_xticks(range(len(LINES))); ax.set_xticklabels(LINES)
ax.set_ylabel('Tumor ploidy'); ax.set_xlim(-0.6, len(LINES) - 0.4)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
nl = sum(r['within_0.2'] == 'True' for r in pl); npj = sum(r['purple_java_within_0.2'] == 'True' for r in pl)
npt = sum(r['purple_java_ploidy'] not in ('NA', '') for r in pl)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], ls='', marker='o', mfc=CMP, mec=CMP, ms=4, label=f'PURPLE (TO)  {npj}/{npt} within ±0.2'),
                   Line2D([], [], ls='', marker='o', mfc=LP, mec=LP, ms=4, label=f'LongPhase-TO  {nl}/{len(pl)} within ±0.2'),
                   plt.Rectangle((0, 0), 1, 1, color='#f0efec', label='Literature ploidy ±0.2'),
                   Line2D([], [], ls='', marker='o', mfc=INK2, mec=INK2, ms=4, label='ONT'),
                   Line2D([], [], ls='', marker='s', mfc='white', mec=INK2, ms=4, label='PacBio HiFi')],
          loc='upper left', fontsize=6.3, ncol=5, handletextpad=0.3, columnspacing=1.2)
ax.set_ylim(1, 5.6)


# (c) HG008: TDF error, copy-number agreement, LOH precision against expected fraction
specs = [('estimated_tdf', 'Estimated DNA fraction', True), ('ploidy', 'Ploidy', False),
         ('cn_within_0.5', 'Genome within 0.5 of\nNIST copy number', False), ('loh_precision', 'LOH precision (by length)', False)]
for i, (k, lab, diag) in enumerate(specs):
    ax = fig.add_subplot(gs[2, i])
    for meth, col in [('PURPLE (Hartwig Java, tumor-only)', CMP), ('LongPhase-TO', LP)]:
        rs = sorted([r for r in hg if r['method'] == meth], key=lambda r: float(r['truth_tdf']))
        x = [float(r['truth_tdf']) for r in rs]; y = [float(r[k]) for r in rs]
        ax.plot(x, y, color=col, lw=1.5, marker='o', ms=3.5, mfc=col, mec='white', mew=0.6,
                label=('PURPLE (TO)' if col == CMP else 'LongPhase-TO') if i == 0 else None)
    if i == 0:
        row_first.append(ax)
    if k == 'ploidy':
        ax.axhline(1.47, color=INK2, lw=0.8, ls=(0, (3, 2)), zorder=0)
        ax.text(0.15, 1.76, 'dashed: NIST mean\ncopy number 1.47', fontsize=5.8, color=INK2, va='top')
    if diag:
        ax.plot([0, 1.05], [0, 1.05], color=GRID, lw=1, zorder=0); ax.set_ylim(0, 1.05)
    ax.set_xticks([0.2, 0.6, 1.0]); ax.set_xlim(0.1, 1.05)
    ax.set_xlabel('Expected DNA fraction'); ax.set_title(lab, fontsize=7, color=INK, pad=4)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    if i == 0:
        ax.legend(loc='lower right', fontsize=6, handlelength=1.5)
for lab, ax0 in zip('abc', row_first):
    ax0.text(-0.42 if lab != 'b' else -0.07, 1.22 if lab != 'b' else 1.06, lab, transform=ax0.transAxes, fontsize=10, fontweight='bold', va='bottom')
row_first[2].text(-0.2, 1.22, '    HG008 (external cell line; truth = NIST v0.6 copy-number benchmark)', transform=row_first[2].transAxes, fontsize=7, color=INK2, va='bottom')

os.makedirs(os.path.join(HERE, 'drafts'), exist_ok=True)
fig.savefig(os.path.join(HERE, 'drafts', 'fig4_draft.png'), dpi=300, bbox_inches='tight', facecolor='white')
