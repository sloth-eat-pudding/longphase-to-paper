"""Supplementary Figure: LongPhase-TO against QUAL-threshold filtering of its input caller.

Precision-recall curves of ClairS-TO and DeepSomatic-TO when the QUAL threshold of their FILTER=PASS calls is raised
(0th to 98th percentile, sweep A), with the LongPhase-TO point of each caller, for all 14 cell-line datasets
at one tumor DNA fraction and variant type. ONT HCC1937 = HCC1937_UCSC_dedup.
Data: issue #1 delivery, branch revision-2026-11-results, Revision_2026-11/qual_threshold_baseline/
(copied to ../analysis/qual_threshold_baseline/). Writes supp_qual_threshold_pr.pdf.
"""
import os, sys, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, '..', 'analysis', 'qual_threshold_baseline')
TDF = float(sys.argv[1]) if len(sys.argv) > 1 else 0.4; VT = sys.argv[2] if len(sys.argv) > 2 else 'SNV'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.edgecolor': '#898781', 'xtick.color': '#52514e', 'ytick.color': '#52514e',
                     'legend.frameon': False, 'pdf.fonttype': 42})
q = pd.read_csv(os.path.join(D, 'qual_sweep.tsv'), sep='\t'); l = pd.read_csv(os.path.join(D, 'lp_vs_curve.tsv'), sep='\t')
keep = lambda d: d[~d.dataset.str.contains('with_duplicate') & (d.sweep == 'A')]
q, l = keep(q), keep(l)
PAIR = {'ClairS-TO': 'LongPhase-TO(ClairS-TO)', 'DeepSomatic-TO': 'LongPhase-TO(DeepSomatic-TO)'}
CC = {'ClairS-TO': '#2a78d6', 'DeepSomatic-TO': '#eb6834'}
l = l[l.method.isin(PAIR.values())]
DS = [d for p in ['ONT', 'PacBio'] for d in sorted(l[l.platform == p].dataset.unique())]
name = lambda d: d.replace('_dedup', '').replace('_', '\n', 1)
fig, ax = plt.subplots(2, 8, figsize=(7.2, 2.6), sharex=True, sharey=True)
pos = [(0, j) for j in range(8)] + [(1, j) for j in range(6)]
for (r, c), ds in zip(pos, DS):
    a = ax[r, c]
    for cal in PAIR:
        s = q[(q.dataset == ds) & (q.caller == cal) & (q.target_tdf == TDF) & (q.variant_type == VT) & (q.threshold != 'default')]
        s = s.sort_values('percentile'); a.plot(s.recall, s.precision, color=CC[cal], lw=1.0, alpha=0.55)
        p = l[(l.dataset == ds) & (l.caller == cal) & (l.target_tdf == TDF) & (l.variant_type == VT)]
        a.plot(p.lp_recall, p.lp_precision, 'D', color=CC[cal], ms=4.5, mec='white', mew=0.6)
    a.set_title(name(ds), fontsize=6); a.set_xlim(0, 1); a.set_ylim(0, 1); a.grid(color='#e1e0d9', lw=0.4)
    a.set_xticks([0, 1]); a.set_yticks([0, 0.5, 1]); a.set_aspect('equal')
for a in ax[:, 0]: a.set_ylabel('Precision')
for a in list(ax[1, :6]) + [ax[0, 6], ax[0, 7]]: a.set_xlabel('Recall'); a.tick_params(labelbottom=True)
fig.text(0.0, 0.74, 'Nanopore', rotation=90, fontsize=7, fontweight='bold', va='center')
fig.text(0.0, 0.27, 'PacBio HiFi', rotation=90, fontsize=7, fontweight='bold', va='center')
ax[1, 6].axis('off'); ax[1, 7].axis('off')
h = [plt.Line2D([], [], color=CC[c], lw=1, label=f'{c}, QUAL threshold') for c in PAIR] + \
    [plt.Line2D([], [], ls='', marker='D', color=CC[c], label=PAIR[c]) for c in PAIR]
ax[1, 6].legend(handles=h, loc='center left', fontsize=6, bbox_to_anchor=(-0.1, 0.5))
fig.subplots_adjust(left=0.07, right=0.99, wspace=0.25, hspace=0.55)
fig.savefig(os.path.join(HERE, 'supp_qual_threshold_pr.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(HERE, 'drafts', 'supp_qual_threshold_pr.png'), dpi=200, bbox_inches='tight')
