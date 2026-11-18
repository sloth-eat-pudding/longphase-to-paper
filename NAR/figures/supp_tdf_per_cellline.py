"""Supplementary Figure (draft): estimated against expected DNA fraction, one panel per dataset (option 3 style).

Data from fig4ab_merged.py (same sources as Figure 4 and Supplementary Figure S6).
LongPhase-TO = LongPhase-TO(DeepSomatic-TO), no ClairS-TO PoN (the configuration of Figure 4).
Compact 3 x 5 grid: left block = 9 nanopore datasets (3 x 3), right block = 6 PacBio HiFi datasets (3 x 2).
Dataset names sit inside each panel (bottom right, an area that no data point uses). Writes drafts/supp_tdf_per_cellline_draft.png.
"""
import importlib.util, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ab', os.path.join(HERE, 'fig4ab_merged.py'))
ab = importlib.util.module_from_spec(spec); spec.loader.exec_module(ab)
pts, COL, MK = ab.pts, ab.COL, ab.MK
INK, INK2, MUTED = ab.INK, ab.INK2, ab.MUTED
ORDER = ['LongPhase-TO', 'PURPLE (TO)', 'ASCAT (T+N)', 'ASCAT (TO)']
Z = {'LongPhase-TO': 5, 'PURPLE (TO)': 4, 'ASCAT (T+N)': 3, 'ASCAT (TO)': 2}
ONT = ['COLO829_ONT', 'COLO829_NYGC', 'H1437_UCSC', 'H2009_UCSC', 'HCC1395_HKU', 'HCC1395_NYGC',
       'HCC1937_UCSC', 'HCC1954_UCSC', 'HG008']
PB = ['COLO829_PacBio', 'H1437_CASTLE', 'H2009_CASTLE', 'HCC1395_PacBio', 'HCC1937_CASTLE', 'HCC1954_CASTLE']


def panel(ax, plat, ds):
    ax.plot([0, 1.1], [0, 1.1], color=MUTED, lw=0.5, zorder=1)
    for m in ORDER[::-1]:
        q = sorted([p for p in pts if p['platform'] == plat and p['dataset'] == ds and p['method'] == m],
                   key=lambda p: p['expected'])
        if not q:
            continue
        x = [p['expected'] for p in q]; y = [p['estimate'] for p in q]
        ax.plot(x, y, color=COL[m], lw=1.1 if m == 'LongPhase-TO' else 0.75, zorder=Z[m])
        ax.scatter(x, y, s=9, marker=MK[m], facecolor=COL[m], edgecolor='white', lw=0.35, zorder=Z[m] + 0.5)
    ax.set_xlim(0.1, 1.08); ax.set_ylim(0.1, 1.08)
    ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([0.2, 0.6, 1.0])
    ax.tick_params(labelsize=5.4, pad=1.5)
    ax.spines['top'].set_visible(True); ax.spines['right'].set_visible(True)
    for sp in ax.spines.values():
        sp.set_color('#c9c8c3'); sp.set_linewidth(0.5)
    ax.text(0.96, 0.05, ds.replace('_', ' '), transform=ax.transAxes, ha='right', va='bottom', fontsize=6, color=INK)


fig = plt.figure(figsize=(7.1, 4.35))
gs = fig.add_gridspec(3, 6, width_ratios=[1, 1, 1, 0.16, 1, 1], wspace=0.07, hspace=0.07)
ref = None
for blk, (plat, names, cols) in enumerate([('ONT', ONT, [0, 1, 2]), ('PacBio HiFi', PB, [4, 5])]):
    for k, ds in enumerate(names):
        r, c = divmod(k, len(cols))
        ax = fig.add_subplot(gs[r, cols[c]], sharex=ref, sharey=ref); ref = ref or ax
        panel(ax, plat, ds)
        ax.tick_params(labelbottom=(r == 2), labelleft=(cols[c] in (0, 4)))
        if r == 0 and c == 0:
            ax.text(0.0, 1.05, 'Nanopore' if plat == 'ONT' else plat, transform=ax.transAxes, fontsize=7,
                    fontweight='bold', color=INK, va='bottom')
fig.supxlabel('Expected DNA fraction', fontsize=6.5, y=0.035)
fig.supylabel('Estimated DNA fraction', fontsize=6.5, x=0.075)
h = [Line2D([], [], color=COL[m], lw=1.0, marker=MK[m], ms=3.6, mfc=COL[m], mec='white', mew=0.4, label=m) for m in ORDER]
fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 0.995), ncol=4, fontsize=6.3,
           handlelength=1.8, columnspacing=1.6, handletextpad=0.4)
fig.text(0.5, -0.005, 'ASCAT was not run on HG008.', ha='center', fontsize=5.6, color=INK2)
fig.savefig(os.path.join(HERE, 'drafts', 'supp_tdf_per_cellline_draft.png'), bbox_inches='tight', facecolor='white', dpi=220)
fig.savefig(os.path.join(HERE, 'supp_tdf_per_dataset.pdf'), bbox_inches='tight', facecolor='white')  # used by ../supplementary.tex
