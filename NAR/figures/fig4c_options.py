"""Three layout options for Figure 4(c): HG008 segment-level copy number against the NIST v0.6 benchmark.

MOCK-UP DATA. The segment values come from drafts/MOCKUP_hg008_segments_traced_from_png.json, which was traced
from the student's bitmap Revision_2026-11/hg008_cn_track_1.0.png (about 0.9 Mb per pixel). They only show the
layout. The final figure needs the segment tables (NIST truth, LongPhase-TO, PURPLE) for all five mixtures.
The right-hand summary (genome within 0.5 of NIST, five mixtures) uses the real values of hg008_ONT_cn_loh.tsv.
Writes drafts/fig4c_option_{1,2,3}.png.
"""
import csv, json, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.join(HERE, '..', '..', 'Revision_2026-11')
OUT = os.path.join(HERE, 'drafts')
INK, INK2, MUTED, GRID, BAND, TRUTH = '#1a1a19', '#52514e', '#8a8985', '#e9e8e4', '#f6f5f2', '#c9c8c3'
LP, PU = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400})

L = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422,
     135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167,
     46709983, 50818468]
OFF = np.concatenate([[0], np.cumsum(L)]) / 1e6          # Mb offsets
seg = json.load(open(os.path.join(OUT, 'MOCKUP_hg008_segments_traced_from_png.json')))

# 1-Mb grid: truth, LongPhase-TO, PURPLE
G = int(OFF[-1]) + 1
grid = {k: np.full(G, np.nan) for k in ('nist', 'lp', 'pu')}
for k in grid:
    for ch, a, b, v in seg[k]:
        grid[k][int(OFF[ch - 1] + a / 1e6):int(OFF[ch - 1] + b / 1e6) + 1] = v
# in the bitmap, PURPLE is drawn over LongPhase-TO; where both agree, only PURPLE is visible
fill = np.isnan(grid['lp']) & ~np.isnan(grid['pu']) & ~np.isnan(grid['nist'])
grid['lp'][fill] = grid['pu'][fill]
for k in ('lp', 'pu'):
    grid[k][np.isnan(grid['nist'])] = np.nan              # evaluate only where NIST has truth
xg = np.arange(G) + 0.5

hg = [r for r in csv.DictReader(open(os.path.join(REV, 'hg008_ONT_cn_loh.tsv')), delimiter='\t')]
SUM = {m: sorted([(float(r['truth_tdf']), float(r['cn_within_0.5'])) for r in hg if r['method'].startswith(m)])
       for m in ('PURPLE', 'LongPhase-TO')}


def chrom_axis(ax, labels=True):
    for i in range(22):
        if i % 2 == 0:
            ax.axvspan(OFF[i], OFF[i + 1], color=BAND, lw=0, zorder=0)
    ax.set_xlim(0, OFF[-1]); ax.set_xticks([])
    if labels:
        for i in range(22):
            ax.text((OFF[i] + OFF[i + 1]) / 2, -0.04, str(i + 1) if i < 18 or i % 2 == 1 else '', transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=5.3, color=INK2)
    ax.spines['bottom'].set_visible(False)


def runs(y):
    """Yield (x0, x1, value) for constant runs of a 1-Mb grid (NaN breaks)."""
    i = 0
    while i < len(y):
        if np.isnan(y[i]):
            i += 1; continue
        j = i
        while j + 1 < len(y) and not np.isnan(y[j + 1]) and abs(y[j + 1] - y[i]) < 0.03:
            j += 1
        yield i, j + 1, y[i]; i = j + 1


def summary(ax):
    """Fraction of the evaluated genome whose |estimated - NIST| copy number is <= t (TDF 1.0, traced)."""
    t = np.linspace(0, 1, 201)
    for k, col, name in [('pu', PU, 'PURPLE (TO)'), ('lp', LP, 'LongPhase-TO')]:
        ok = ~np.isnan(grid[k]); e = np.abs(grid[k][ok] - grid['nist'][ok])
        ax.plot(t, [100 * np.mean(e <= v) for v in t], color=col, lw=1.2)
        v25 = 100 * np.mean(e <= 0.25)
        ax.plot([0.25], [v25], marker='o', ms=3, color=col, mec='white', mew=0.4)
        ax.text(0.28, v25 + (2 if k == 'lp' else -2), f'{v25:.0f}%', color=col, fontsize=5.8, va='bottom' if k == 'lp' else 'top')
    ax.axvline(0.5, color=MUTED, lw=0.5, ls=(0, (3, 2)))
    ax.text(0.52, 55, '0.5', fontsize=5.2, color=MUTED, va='bottom')
    ax.set_xlim(0, 1); ax.set_ylim(0, 102); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1'])
    ax.set_xlabel('Tolerance (copies)'); ax.set_ylabel('Genome within tolerance (%)')
    ax.text(0.98, 15, 'LongPhase-TO', color=LP, fontsize=5.8, ha='right', va='center')
    ax.text(0.98, 6, 'PURPLE (TO)', color=PU, fontsize=5.8, ha='right', va='center')
    ax.set_title('Agreement with NIST', loc='left')


def letter(ax, x=-0.06, y=1.1):
    ax.text(x, y, 'c', transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')


MOCK = 'MOCK-UP: all values traced from the student PNG (pure tumor, about 1 Mb resolution); not for publication'


# ---------- Option 1: split genome tracks, truth band under each method ----------
def option1(fig=None, mock=True):
    own = fig is None
    if own:
        fig = plt.figure(figsize=(7.1, 2.3))
    gs = fig.add_gridspec(2, 2, width_ratios=[4.3, 1], hspace=0.25, wspace=0.18)
    for r, (k, col, name) in enumerate([('pu', PU, 'PURPLE (TO)'), ('lp', LP, 'LongPhase-TO')]):
        ax = fig.add_subplot(gs[r, 0]); chrom_axis(ax, labels=(r == 1))
        for a, b, v in runs(grid['nist']):
            ax.add_patch(plt.Rectangle((a, v - 0.17), b - a, 0.34, color=TRUTH, lw=0, zorder=1))
        for a, b, v in runs(grid[k]):
            ax.plot([a, b], [v, v], color=col, lw=1.4, solid_capstyle='butt', zorder=2)
        ax.set_ylim(0, 3.6); ax.set_yticks([0, 1, 2, 3]); ax.set_ylabel('Total CN')
        ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)
        w = np.nanmean(np.abs(grid[k] - grid['nist']) <= 0.5) * 100
        ax.text(0.003, 0.93, name, transform=ax.transAxes, color=col, fontsize=6.5, fontweight='bold', va='top')
        if r == 0:
            letter(ax, -0.05, 1.12)
            ax.set_title('HG008 pure tumor: total copy number', loc='left', pad=6)
            h = [plt.Rectangle((0, 0), 1, 1, color=TRUTH, label='NIST v0.6 truth')]
            ax.legend(handles=h, loc='lower right', bbox_to_anchor=(1.0, 0.98), fontsize=5.8, handlelength=1.6)
    summary(fig.add_subplot(gs[:, 1]))
    if mock:
        fig.text(0.5, -0.06 if own else 0.0, MOCK, ha='center', fontsize=5.5, color='#c2402f')
    if own:
        fig.savefig(os.path.join(OUT, 'fig4c_option_1.png'), bbox_inches='tight', facecolor='white')


# ---------- Option 2: signed error along the genome ----------
def option2():
    fig = plt.figure(figsize=(7.1, 1.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[4.3, 1], wspace=0.18)
    ax = fig.add_subplot(gs[0, 0]); chrom_axis(ax)
    ax.axhspan(-0.5, 0.5, color='#ecebe7', lw=0, zorder=0.5)
    ax.axhline(0, color=MUTED, lw=0.5, zorder=1)
    for k, col, name, dy in [('pu', PU, 'PURPLE (TO)', 0), ('lp', LP, 'LongPhase-TO', 0)]:
        e = grid[k] - grid['nist']
        ax.fill_between(xg, 0, np.where(np.isnan(e), np.nan, e), color=col, alpha=0.85 if k == 'lp' else 0.75, lw=0,
                        step='mid', zorder=3 if k == 'lp' else 2)
    ax.set_ylim(-1.0, 1.0); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_ylabel('Estimated − NIST CN')
    ax.text(0.003, 0.97, 'PURPLE (TO)', transform=ax.transAxes, color=PU, fontsize=6.5, fontweight='bold', va='top')
    ax.text(0.19, 0.97, 'LongPhase-TO', transform=ax.transAxes, color=LP, fontsize=6.5, fontweight='bold', va='top')
    ax.text(0.997, 0.97, 'gray band: ±0.5', transform=ax.transAxes, color=INK2, fontsize=5.8, ha='right', va='top')
    ax.set_title('HG008 pure tumor: copy-number error along the genome', loc='left', pad=6)
    letter(ax, -0.05, 1.08)
    summary(fig.add_subplot(gs[0, 1]))
    fig.text(0.5, -0.08, MOCK, ha='center', fontsize=5.5, color='#c2402f')
    fig.savefig(os.path.join(OUT, 'fig4c_option_2.png'), bbox_inches='tight', facecolor='white')


# ---------- Option 3: length-weighted agreement by NIST copy-number state ----------
def option3():
    fig = plt.figure(figsize=(7.1, 2.1))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.35, 1], wspace=0.35)
    states = [1, 2, 3]
    for c, (k, col, name) in enumerate([('pu', PU, 'PURPLE (TO)'), ('lp', LP, 'LongPhase-TO')]):
        ax = fig.add_subplot(gs[0, c])
        ax.plot([0.5, 3.5], [0.5, 3.5], color=GRID, lw=0.8, zorder=0)
        for s in states:
            ax.add_patch(plt.Rectangle((s - 0.4, s - 0.5), 0.8, 1.0, color=BAND, lw=0, zorder=0))
            m = np.abs(grid['nist'] - s) < 0.25
            v = grid[k][m]; v = v[~np.isnan(v)]
            if not len(v):
                continue
            b = np.round(v * 10) / 10
            u, n = np.unique(b, return_counts=True)
            ax.scatter(np.full(len(u), s), u, s=n * 0.35, color=col, alpha=0.85, lw=0, zorder=2)
            ax.text(s + 0.38, s - 0.47, f'{n.sum():.0f} Mb', fontsize=5.2, color=INK2, ha='right', va='bottom')
        w = np.nanmean(np.abs(grid[k] - grid['nist']) <= 0.5) * 100
        ax.set_xlim(0.5, 3.5); ax.set_ylim(0, 3.6); ax.set_xticks(states); ax.set_yticks([0, 1, 2, 3])
        ax.set_xlabel('NIST total copy number'); ax.set_ylabel('Estimated copy number' if c == 0 else '')
        ax.set_title(name, loc='left', color=col, fontweight='bold')
        if c == 0:
            letter(ax, -0.25, 1.06)
    ax = fig.add_subplot(gs[0, 2]); summary(ax)
    fig.text(0.39, 0.99, 'HG008 pure tumor; circle area = genome length (Mb); shaded box = within 0.5 of truth',
             ha='center', fontsize=5.8, color=INK2)
    fig.text(0.5, -0.08, MOCK, ha='center', fontsize=5.5, color='#c2402f')
    fig.savefig(os.path.join(OUT, 'fig4c_option_3.png'), bbox_inches='tight', facecolor='white')


if __name__ == '__main__':
  os.makedirs(OUT, exist_ok=True)
  option1(); option2(); option3()
  for k in ('lp', 'pu'):
    ok = ~np.isnan(grid[k]); print(k, "traced within 0.5:", round(float(np.mean(np.abs(grid[k][ok] - grid["nist"][ok]) <= 0.5)), 3))
