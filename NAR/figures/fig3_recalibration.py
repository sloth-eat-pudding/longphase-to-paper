"""Figure 3 (F1 of each dataset: (a) nanopore SNV and indel, (b) PacBio HiFi SNV and indel) and the four per-dataset Supplementary grids (F1, precision, recall).

Source: the somatic tables of ../supplementary.tex, so that figures and tables show the same values:
  tab:somatic_purity, tab:somatic_indel_purity (nanopore), tab:pacbio_somatic_snv, tab:pacbio_somatic_indel (PacBio HiFi).
DeepSomatic-TO + LongPhase-TO = LongPhase-TO(DeepSomatic-TO), no ClairS-TO PoN (as in the tables).
Style (author choice, 2026-11-17, "gray before, color after"): caller alone = gray dashed line with open markers
(circle = ClairS-TO, square = DeepSomatic-TO); caller + LongPhase-TO = solid line, ClairS-TO navy #3B4992,
DeepSomatic-TO red #EE0000 (AAAS). Open axes (left and bottom only), no grid.
Writes ../Figure_3.pdf, supp_recal_{ont,pacbio}_{snv,indel}.pdf and PNG copies in drafts/.
"""
import os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
SUPP = os.path.join(HERE, '..', 'supplementary.tex')
INK, AXIS, GRAY = '#222222', '#7a7a7a', '#a3a3a3'
C_CS, C_DS = '#3B4992', '#EE0000'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.edgecolor': AXIS, 'axes.labelcolor': INK,
                     'xtick.color': '#52514e', 'ytick.color': '#52514e', 'axes.linewidth': 0.5, 'legend.frameon': False,
                     'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400})
SERIES = [  # label, color, line style, line width, marker, marker face, zorder
    ('ClairS-TO', GRAY, (0, (2.5, 1.5)), 0.8, 'o', 'white', 2),
    ('DeepSomatic-TO', GRAY, (0, (2.5, 1.5)), 0.8, 's', 'white', 2),
    ('ClairS-TO+LongPhase-TO', C_CS, '-', 1.1, 'o', C_CS, 3),
    ('DeepSomatic-TO+LongPhase-TO', C_DS, '-', 1.1, 's', C_DS, 3)]
T = [0.2, 0.4, 0.6, 0.8, 1.0]
TABLES = {('ONT', 'SNV'): 'tab:somatic_purity', ('ONT', 'indel'): 'tab:somatic_indel_purity',
          ('PacBio', 'SNV'): 'tab:pacbio_somatic_snv', ('PacBio', 'indel'): 'tab:pacbio_somatic_indel'}
ORDER = ['COLO829', 'H1437', 'H2009', 'HCC1395', 'HCC1937', 'HCC1954']
src = open(SUPP).read()


def table(label):
    i = src.index('\\label{%s}' % label); j = src.index('\\end{longtable}', i); out = {}
    for l in src[i:j].splitlines():
        c = [x.strip().replace('\\_', '_') for x in l.strip().rstrip('\\').split('&')]
        if len(c) == 6 and re.match(r'^[0-9.]+$', c[1]):
            out[(c[0].replace('_CMH', '_CASTLE'), float(c[1]), c[2])] = dict(precision=float(c[3]), recall=float(c[4]), f1=float(c[5]))
    return out


DATA = {k: table(v) for k, v in TABLES.items()}


FIRST = ('_ONT', '_HKU')   # replicate order as in the DNA-fraction Supplementary figure


def names(d):
    return sorted({k[0] for k in d}, key=lambda x: (ORDER.index(x.split('_')[0]), not x.endswith(FIRST), x))


def label(ax, letter, x, y=1.12):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom', ha='left', color=INK)


def handles():
    order = [SERIES[0], SERIES[2], SERIES[1], SERIES[3]]
    return [Line2D([], [], color=col, ls=ls, lw=lw + 0.2, marker=mk, ms=3.2, mfc=fc, mec=col, mew=0.6,
                   label=lab.replace('+', ' + ')) for lab, col, ls, lw, mk, fc, z in order]


def style_axes(ax):
    for sp in ('left', 'bottom'):
        ax.spines[sp].set_color(AXIS); ax.spines[sp].set_linewidth(0.5)
    ax.tick_params(labelsize=5.5, pad=1.5, length=2, color=AXIS)


def series(ax, getter):
    for lab, col, ls, lw, mk, fc, z in SERIES:
        ax.plot(T, [getter(t, lab) for t in T], color=col, ls=ls, lw=lw, zorder=z, marker=mk,
                ms=2.6 if z == 3 else 2.4, mfc=fc, mec=col, mew=0.6)


# ---------------- Figure 3: F1 of each dataset ----------------
def figure3(by='platform', draft=None):
    """by='variant': (a) SNV, (b) indel, rows = platforms. by='platform': (a) nanopore, (b) PacBio HiFi, rows = SNV, indel."""
    fig = plt.figure(figsize=(7.1, 5.2))
    outer = fig.add_gridspec(2, 1, hspace=0.3)
    blocks = ([[('SNV', 'ONT'), ('SNV', 'PacBio')], [('indel', 'ONT'), ('indel', 'PacBio')]] if by == 'variant' else
              [[('SNV', 'ONT'), ('indel', 'ONT')], [('SNV', 'PacBio'), ('indel', 'PacBio')]])
    for b, rows_ in enumerate(blocks):
        blk = outer[b].subgridspec(2, 1, hspace=0.45)
        for r, (vt, plat) in enumerate(rows_):
            top = 1.0 if vt == 'SNV' else 0.8
            d = DATA[(plat, vt)]; ds = names(d)
            row = blk[r].subgridspec(1, len(ds), wspace=0.12)
            first = None
            for c, name in enumerate(ds):
                ax = fig.add_subplot(row[0, c], sharey=first); first = first or ax
                series(ax, lambda t, lab: d[(name, t, lab)]['f1'])
                ax.set_xlim(0.12, 1.08); ax.set_ylim(0, top); ax.set_xticks([0.2, 0.6, 1.0])
                ax.set_yticks([0, 0.5, 1.0] if vt == 'SNV' else [0, 0.4, 0.8])
                style_axes(ax)
                ax.set_title(name.replace('_', ' '), fontsize=5.8, color=INK, pad=2)
                if c > 0:
                    ax.tick_params(labelleft=False)
                else:
                    ax.set_ylabel(f"{'Nanopore' if plat == 'ONT' else 'PacBio HiFi'}\n{vt} F1", fontsize=6.3, color=INK)
                    if r == 0:
                        label(ax, 'ab'[b], x=-0.75 if plat == 'ONT' else -0.55)
    fig.supxlabel('Tumor DNA fraction', fontsize=6.5, y=0.045)
    fig.legend(handles=handles(), loc='upper center', bbox_to_anchor=(0.5, 0.965), ncol=4, fontsize=6.3,
               handlelength=2.4, columnspacing=1.6, handletextpad=0.4)
    fig.subplots_adjust(top=0.9)
    if draft:
        fig.savefig(os.path.join(HERE, 'drafts', draft), bbox_inches='tight', facecolor='white', dpi=220); return
    fig.savefig(os.path.join(HERE, '..', 'Figure_3.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(HERE, 'drafts', 'Figure_3_new.png'), bbox_inches='tight', facecolor='white', dpi=220)


# ---------------- per-dataset grids (Supplement) ----------------
def grid(plat, vt):
    d = DATA[(plat, vt)]; ds = names(d)
    fig, axs = plt.subplots(3, len(ds), figsize=(7.1, 3.9), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.12, hspace=0.12))
    for c, name in enumerate(ds):
        for r, (key, ylab) in enumerate([('f1', 'F1-score'), ('precision', 'Precision'), ('recall', 'Recall')]):
            ax = axs[r, c]
            series(ax, lambda t, lab: d[(name, t, lab)][key])
            ax.set_xlim(0.12, 1.08); ax.set_ylim(0, 1.02); ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([0, 0.5, 1.0])
            style_axes(ax)
            if r == 0:
                ax.set_title(name.replace('_', ' '), fontsize=6 if len(ds) > 6 else 6.3, color=INK, pad=3)
            if c == 0:
                ax.set_ylabel(ylab, fontsize=6.5, color=INK)
    fig.supxlabel('Tumor DNA fraction', fontsize=6.5, y=0.02)
    fig.legend(handles=handles(), loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=6.3,
               handlelength=2.4, columnspacing=1.6, handletextpad=0.4)
    fig.subplots_adjust(top=0.88)
    stem = f"supp_recal_{'ont' if plat == 'ONT' else 'pacbio'}_{vt.lower()}"
    fig.savefig(os.path.join(HERE, stem + '.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(os.path.join(HERE, 'drafts', stem + '.png'), bbox_inches='tight', facecolor='white', dpi=220)


if __name__ == '__main__':
    figure3()
    for plat in ('ONT', 'PacBio'):
        for vt in ('SNV', 'indel'):
            grid(plat, vt)
