"""Figure 4(c): HG008 total copy number against the NIST v0.6 benchmark, from the real segment tables (issue #5).

Data: Revision_2026-11/hg008_segments/ (NIST truth segments, LongPhase-TO at lambda = 0.8 and Hartwig PURPLE
tumor-only segments for the five ONT mixtures). Each method segment is cut to the NIST truth records on chr1-22
(the evaluated regions). Agreement = length-weighted fraction with |estimated - NIST total CN| <= tolerance;
at tolerance 0.5 this reproduces cn_within_0.5 of hg008_ONT_cn_loh.tsv.

Left: pure tumor (t44_n00) along the genome. Right, by variant:
  A  tolerance curve, pure tumor only (layout of the traced draft)
  B  tolerance curves, pure tumor bold + the four mixtures thin
  C  genome within +-0.25 copies against the expected DNA fraction (five mixtures)
Writes drafts/fig4c_real_{A,B,C}.png.
"""
import csv, os
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
SEGDIR = os.path.join(HERE, '..', '..', 'Revision_2026-11', 'hg008_segments')
OUT = os.path.join(HERE, 'drafts')
INK, INK2, MUTED, GRID, BAND, TRUTH = '#1a1a19', '#52514e', '#8a8985', '#e9e8e4', '#f6f5f2', '#c9c8c3'
LP, PU = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400})

CHROMS = [f'chr{i}' for i in range(1, 23)]
L = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
     135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
     46709983, 50818468]
OFF = dict(zip(CHROMS, np.concatenate([[0], np.cumsum(L)])[:-1] / 1e6))   # Mb offsets
OFF_END = np.concatenate([[0], np.cumsum(L)]) / 1e6
METHODS = [('PURPLE (Hartwig Java, tumor-only)', 'pu', PU, 'PURPLE (TO)'), ('LongPhase-TO', 'lp', LP, 'LongPhase-TO')]
YMAX = 4.5

truth = [r for r in csv.DictReader(open(os.path.join(SEGDIR, 'hg008_nist_truth_segments.tsv')), delimiter='\t')
         if r['chrom'] in OFF]
segs = defaultdict(list)
TDF = {}
for r in csv.DictReader(open(os.path.join(SEGDIR, 'hg008_method_segments.tsv')), delimiter='\t'):
    segs[(r['mixture'], r['method'])].append(r)
    TDF[r['mixture']] = float(r['truth_tdf'])
MIX = sorted(TDF, key=TDF.get)
TABLE = {(r['mixture'], r['method']): float(r['cn_within_0.5'])
         for r in csv.DictReader(open(os.path.join(SEGDIR, 'hg008_ONT_cn_loh.tsv')), delimiter='\t')}


def pieces(mix, method):
    """(chrom, start, end, nist_cn, est_cn) for each overlap of a method segment with a NIST truth record."""
    by = defaultdict(list)
    for s in segs[(mix, method)]:
        by[s['chrom']].append((int(s['start']), int(s['end']), float(s['total_cn'])))
    out = []
    for t in truth:
        ts, te, tc = int(t['start']), int(t['end']), float(t['total_cn'])
        for a, b, v in by[t['chrom']]:
            lo, hi = max(ts, a), min(te, b)
            if hi > lo:
                out.append((t['chrom'], lo, hi, tc, v))
    return out


P = {(m, meth): pieces(m, meth) for m in MIX for meth, *_ in METHODS}


def within(p, tol):
    w = np.array([b - a for _, a, b, _, _ in p], float)
    e = np.array([abs(v - t) for *_, t, v in p])
    return 100 * w[e <= tol + 1e-9].sum() / w.sum()


# 0.5 must reproduce the student's table
for k, p in P.items():
    assert abs(within(p, 0.5) / 100 - TABLE[k]) < 0.0015, (k, within(p, 0.5), TABLE[k])

PURE = MIX[-1]
TOL = np.linspace(0, 1, 201)


def chrom_axis(ax, labels=True):
    for i in range(22):
        if i % 2 == 0:
            ax.axvspan(OFF_END[i], OFF_END[i + 1], color=BAND, lw=0, zorder=0)
    ax.set_xlim(0, OFF_END[-1]); ax.set_xticks([])
    if labels:
        for i in range(22):
            ax.text((OFF_END[i] + OFF_END[i + 1]) / 2, -0.04, str(i + 1) if i < 18 or i % 2 == 1 else '',
                    transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=5.3, color=INK2)
    ax.spines['bottom'].set_visible(False)


def merged(rows):
    """Merge adjacent rows (chrom, start, end, value) with the same value, in genome Mb."""
    out = []
    for c, a, b, v in rows:
        x0, x1 = OFF[c] + a / 1e6, OFF[c] + b / 1e6
        if out and abs(out[-1][1] - x0) < 1e-6 and abs(out[-1][2] - v) < 1e-9:
            out[-1][1] = x1
        else:
            out.append([x0, x1, v])
    return out


TRUTH_RUNS = merged(sorted(((t['chrom'], int(t['start']), int(t['end']), float(t['total_cn'])) for t in truth),
                           key=lambda r: (CHROMS.index(r[0]), r[1])))


def tracks(fig, gs_col):
    for r, (meth, k, col, name) in enumerate(METHODS):
        ax = fig.add_subplot(gs_col[r, 0]); chrom_axis(ax, labels=(r == 1))
        for x0, x1, v in TRUTH_RUNS:
            ax.add_patch(plt.Rectangle((x0, v - 0.17), x1 - x0, 0.34, color=TRUTH, lw=0, zorder=1))
        for x0, x1, v in merged([(c, a, b, v) for c, a, b, _, v in P[(PURE, meth)]]):
            ax.plot([x0, x1], [min(v, YMAX - 0.08)] * 2, color=col, lw=1.4, solid_capstyle='butt', zorder=2)
            if v > YMAX:
                ax.plot((x0 + x1) / 2, YMAX - 0.08, marker='^', ms=2.2, color=col, mew=0, zorder=3)
        ax.set_ylim(0, YMAX); ax.set_yticks([0, 1, 2, 3, 4]); ax.set_ylabel('Total CN')
        ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
        ax.text(0.003, 0.95, name, transform=ax.transAxes, color=col, fontsize=6.5, fontweight='bold', va='top')
        if r == 0:
            ax.text(-0.05, 1.12, 'c', transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')
            ax.set_title('HG008 pure tumor: total copy number', loc='left', pad=6)
            h = [plt.Rectangle((0, 0), 1, 1, color=TRUTH, label='NIST v0.6 truth')]
            ax.legend(handles=h, loc='lower right', bbox_to_anchor=(1.0, 0.98), fontsize=5.8, handlelength=1.6)


def tol_axis(ax):
    ax.axvline(0.5, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.set_xlim(0, 1); ax.set_ylim(0, 102); ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
    ax.set_xlabel('Tolerance (copies)'); ax.set_ylabel('Genome within tolerance (%)')
    ax.set_title('Agreement with NIST', loc='left')


def label_25(ax, k, col, v):
    ax.plot([0.25], [v], marker='o', ms=3, color=col, mec='white', mew=0.4, zorder=4)
    ax.text(0.28, v + (2 if k == 'lp' else -2), f'{v:.1f}%', color=col, fontsize=5.8,
            va='bottom' if k == 'lp' else 'top')


def right_A(ax):
    for meth, k, col, name in METHODS:
        ax.plot(TOL, [within(P[(PURE, meth)], t) for t in TOL], color=col, lw=1.2)
        label_25(ax, k, col, within(P[(PURE, meth)], 0.25))
    tol_axis(ax)
    ax.text(0.98, 15, 'LongPhase-TO', color=LP, fontsize=5.8, ha='right', va='center')
    ax.text(0.98, 6, 'PURPLE (TO)', color=PU, fontsize=5.8, ha='right', va='center')


def right_B(ax):
    for meth, k, col, name in METHODS:
        for m in MIX[:-1]:
            ax.plot(TOL, [within(P[(m, meth)], t) for t in TOL], color=col, lw=0.5, alpha=0.35)
        ax.plot(TOL, [within(P[(PURE, meth)], t) for t in TOL], color=col, lw=1.3)
        label_25(ax, k, col, within(P[(PURE, meth)], 0.25))
    tol_axis(ax)
    ax.text(0.98, 24, 'thick: pure tumor', color=INK2, fontsize=5.3, ha='right', va='center', bbox=KEYBOX)
    ax.text(0.98, 17, 'thin: DNA fraction 0.2–0.8', color=INK2, fontsize=5.3, ha='right', va='center', bbox=KEYBOX)
    ax.text(0.98, 9, 'LongPhase-TO', color=LP, fontsize=5.8, ha='right', va='center', bbox=KEYBOX)
    ax.text(0.98, 2, 'PURPLE (TO)', color=PU, fontsize=5.8, ha='right', va='center', bbox=KEYBOX)


def right_C(ax):
    x = np.array([TDF[m] for m in MIX])
    for meth, k, col, name in METHODS:
        y = [within(P[(m, meth)], 0.25) for m in MIX]
        ax.plot(x, y, color=col, lw=1.0, marker='o' if k == 'lp' else 's', ms=3.2, mec='white', mew=0.4)
        ax.text(1.04, y[-1], f'{y[-1]:.1f}%', color=col, fontsize=5.8, ha='left', va='center')
    ax.set_xlim(0.1, 1.22); ax.set_ylim(0, 102); ax.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_xlabel('Expected DNA fraction'); ax.set_ylabel('Genome within ±0.25 copies (%)')
    ax.set_title('Agreement with NIST', loc='left')
    ax.text(0.15, 20, 'LongPhase-TO', color=LP, fontsize=5.8, va='center')
    ax.text(0.15, 11, 'PURPLE (TO)', color=PU, fontsize=5.8, va='center')


KEYBOX = dict(facecolor='white', edgecolor='none', pad=0.6)
RIGHT = {'A': right_A, 'B': right_B, 'C': right_C}


def panel_c(fig=None, variant='A'):
    own = fig is None
    if own:
        fig = plt.figure(figsize=(7.1, 2.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[4.3, 1], hspace=0.25, wspace=0.18)
    tracks(fig, gs)
    RIGHT[variant](fig.add_subplot(gs[:, 1]))
    if own:
        fig.savefig(os.path.join(OUT, f'fig4c_real_{variant}.png'), bbox_inches='tight', facecolor='white', dpi=200)
    return fig


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    for v in RIGHT:
        panel_c(variant=v)
    print('mixture  method        within0.5(table)  within0.25  within0.1')
    for (m, meth), p in P.items():
        print(m, meth[:12], f'{within(p, .5):.1f} ({100 * TABLE[(m, meth)]:.1f})', f'{within(p, .25):.1f}',
              f'{within(p, .1):.1f}')
