"""Supplementary Figures (drafts): PacBio HiFi somatic SNV and indel recalibration, one column per dataset.

Source: ../../Revision_2026-11/metrics_snv_indel_ONT_PacBio.tsv (platform PacBio).
Rows: F1, precision, recall. Style follows Figure 3(b): ClairS-TO blue, DeepSomatic-TO orange,
caller alone = dashed line with open markers, caller + LongPhase-TO = solid line with filled markers.
DeepSomatic-TO + LongPhase-TO = LongPhase-TO(DeepSomatic-TO), no ClairS-TO PoN (as in the main text).
Dataset suffix _CMH in the TSV is shown as _CASTLE. Writes drafts/supp_pacbio_{snv,indel}_recalibration_draft.png.
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'Revision_2026-11', 'metrics_snv_indel_ONT_PacBio.tsv')
INK, INK2, MUTED, GRID = '#1a1a19', '#52514e', '#8a8985', '#e9e8e4'
C_CS, C_DS = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.edgecolor': INK2, 'axes.labelcolor': INK,
                     'xtick.color': INK2, 'ytick.color': INK2, 'axes.linewidth': 0.5, 'legend.frameon': False,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                     'savefig.dpi': 400})
SERIES = [('ClairS-TO raw PASS', 'ClairS-TO', C_CS, (0, (3, 2)), 'o', 'white'),
          ('LongPhase-TO(ClairS-TO)', 'ClairS-TO + LongPhase-TO', C_CS, '-', 'o', C_CS),
          ('DeepSomatic-TO raw PASS', 'DeepSomatic-TO', C_DS, (0, (3, 2)), 's', 'white'),
          ('LongPhase-TO(DeepSomatic-TO)', 'DeepSomatic-TO + LongPhase-TO', C_DS, '-', 's', C_DS)]
DS = ['COLO829_PacBio', 'H1437_CMH', 'H2009_CMH', 'HCC1395_PacBio', 'HCC1937_CMH', 'HCC1954_CMH']
ROWS = [('f1', 'F1-score'), ('precision', 'Precision'), ('recall', 'Recall')]
T = [0.2, 0.4, 0.6, 0.8, 1.0]
rows = [r for r in csv.DictReader(open(SRC), delimiter='\t') if r['platform'] == 'PacBio']
D = {(r['dataset'], float(r['target_tdf']), r['method'], r['variant_type']): r for r in rows}


def draw(vt, out):
    fig, axs = plt.subplots(3, 6, figsize=(7.1, 3.9), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.07, hspace=0.07))
    for c, ds in enumerate(DS):
        for r, (key, ylab) in enumerate(ROWS):
            ax = axs[r, c]
            for m, lab, col, ls, mk, fc in SERIES:
                y = [float(D[(ds, t, m, vt)][key]) for t in T]
                solid = ls == '-'
                ax.plot(T, y, color=col, ls=ls, lw=1.1 if solid else 0.8, zorder=3 if solid else 2)
                ax.scatter(T, y, s=9, marker=mk, facecolor=fc, edgecolor=col, lw=0.7, zorder=4 if solid else 3)
            ax.set_xlim(0.1, 1.1); ax.set_ylim(0, 1.02)
            ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([0, 0.5, 1.0])
            ax.grid(axis='y', color=GRID, lw=0.4); ax.set_axisbelow(True)
            ax.tick_params(labelsize=5.6, pad=1.5)
            for side in ('top', 'right'):
                ax.spines[side].set_visible(True)
            for sp in ax.spines.values():
                sp.set_color('#c9c8c3'); sp.set_linewidth(0.5)
            if r == 0:
                ax.set_title(ds.replace('_CMH', '_CASTLE').replace('_', ' '), fontsize=6.3, color=INK, pad=3)
            if c == 0:
                ax.set_ylabel(ylab, fontsize=6.5, color=INK)
    fig.supxlabel('Tumor DNA fraction', fontsize=6.5, y=0.02)
    h = [Line2D([], [], color=col, ls=ls, lw=1.1 if ls == '-' else 0.8, marker=mk, ms=3.6, mfc=fc, mec=col, mew=0.7, label=lab)
         for m, lab, col, ls, mk, fc in SERIES]
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=6.3,
               handlelength=2.2, columnspacing=1.6, handletextpad=0.4)
    fig.subplots_adjust(top=0.88)
    fig.savefig(os.path.join(HERE, 'drafts', out), bbox_inches='tight', facecolor='white', dpi=220)


draw('SNV', 'supp_pacbio_snv_recalibration_draft.png')
draw('INDEL', 'supp_pacbio_indel_recalibration_draft.png')
