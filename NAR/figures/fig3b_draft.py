"""Draft of Figure 3(b) (NAR revision): mean F1 against tumor DNA fraction, ONT and PacBio HiFi, SNV and indel.

Source: ../../Revision_2026-11/metrics_snv_indel_ONT_PacBio.tsv (HCC1937_UCSC_dedup left out).
DeepSomatic-TO recalibration = the DeepSomatic-TO + ClairS-TO PoN configuration on both platforms.
Each point = mean over datasets (ONT 8, PacBio 6). Writes drafts/fig3b_draft.png.
"""
import csv, os
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'Revision_2026-11', 'metrics_snv_indel_ONT_PacBio.tsv')
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
C_CS, C_DS = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.6, 'legend.frameon': False})

rows = [r for r in csv.DictReader(open(SRC), delimiter='\t') if r['dataset'] != 'HCC1937_UCSC_dedup']
agg = defaultdict(list)
for r in rows:
    agg[(r['platform'], r['variant_type'], r['method'], float(r['target_tdf']))].append(float(r['f1']))

SERIES = [('ClairS-TO raw PASS', 'ClairS-TO', C_CS, (0, (3, 2)), 'o', 'white'),
          ('LongPhase-TO(ClairS-TO)', 'ClairS-TO + LongPhase-TO', C_CS, '-', 'o', C_CS),
          ('DeepSomatic-TO raw PASS', 'DeepSomatic-TO', C_DS, (0, (3, 2)), 's', 'white'),
          ('LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)', 'DeepSomatic-TO + LongPhase-TO', C_DS, '-', 's', C_DS)]
fig, axs = plt.subplots(2, 2, figsize=(4.6, 3.8), sharex=True)
for i, vt in enumerate(['SNV', 'INDEL']):
    for j, plat in enumerate(['ONT', 'PacBio']):
        ax = axs[i, j]
        for m, lab, col, ls, mk, fc in SERIES:
            xs = sorted({k[3] for k in agg if k[0] == plat and k[1] == vt and k[2] == m})
            ys = [sum(agg[(plat, vt, m, x)]) / len(agg[(plat, vt, m, x)]) for x in xs]
            ax.plot(xs, ys, color=col, ls=ls, lw=1.5, marker=mk, ms=4, mfc=fc, mec=col, mew=0.9,
                    label=lab if (i, j) == (0, 0) else None)
        ax.set_title(f"{'ONT' if plat == 'ONT' else 'PacBio HiFi'}, {'SNV' if vt == 'SNV' else 'indel'}",
                     fontsize=7.5, color=INK, pad=4)
        ax.grid(axis='y', color=GRID, lw=0.6); ax.set_axisbelow(True)
        ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_ylim(0, 0.8 if vt == 'SNV' else 0.5)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
        if j == 0:
            ax.set_ylabel('Mean F1-score')
        if i == 1:
            ax.set_xlabel('Tumor DNA fraction')
fig.legend(loc='upper center', ncol=2, fontsize=6.5, bbox_to_anchor=(0.5, 1.07), handlelength=2.6)
fig.tight_layout()
os.makedirs(os.path.join(HERE, 'drafts'), exist_ok=True)
fig.savefig(os.path.join(HERE, 'drafts', 'fig3b_draft.png'), dpi=300, bbox_inches='tight', facecolor='white')
