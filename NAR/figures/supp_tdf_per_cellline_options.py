"""Supplementary Figure options (drafts): DNA fraction for each dataset. Data from fig4ab_merged.py.
LongPhase-TO = DeepSomatic-TO + ClairS-TO PoN configuration. Writes drafts/supp_percell_option_{1,2,3}.png.
  1: signed error, grid of cell line x platform
  2: heatmap of signed error, one row per dataset, one block per method
  3: estimate vs expected, grid of cell line x platform, quiet style
"""
import importlib.util, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ab', os.path.join(HERE, 'fig4ab_merged.py'))
ab = importlib.util.module_from_spec(spec); spec.loader.exec_module(ab)
pts, COL, MK = ab.pts, ab.COL, ab.MK
INK, INK2, MUTED, GRID, BAND = ab.INK, ab.INK2, ab.MUTED, ab.GRID, ab.BAND
OUT = os.path.join(HERE, 'drafts')
ORDER = ['LongPhase-TO', 'PURPLE (TO)', 'ASCAT (T+N)', 'ASCAT (TO)']   # front to back in label order
LINES = ['COLO829', 'H1437', 'H2009', 'HCC1395', 'HCC1937', 'HCC1954', 'HG008']
PLATS = [('ONT', 'Nanopore'), ('PacBio HiFi', 'PacBio HiFi')]
ONT = ['COLO829_ONT', 'COLO829_NYGC', 'H1437_UCSC', 'H2009_UCSC', 'HCC1395_HKU', 'HCC1395_NYGC',
       'HCC1937_UCSC', 'HCC1954_UCSC', 'HG008']
PB = ['COLO829_PacBio', 'H1437_CASTLE', 'H2009_CASTLE', 'HCC1395_PacBio', 'HCC1937_CASTLE', 'HCC1954_CASTLE']


def get(plat, ds, m):
    return sorted([p for p in pts if p['platform'] == plat and p['dataset'] == ds and p['method'] == m],
                  key=lambda p: p['expected'])


def datasets(plat, line):
    return [d for d in (ONT if plat == 'ONT' else PB) if d.split('_')[0] == line]


def legend(fig, y, extra=()):
    h = [Line2D([], [], color=COL[m], lw=1.0, marker=MK[m], ms=3.8, mfc=COL[m], mec='white', mew=0.4, label=m)
         for m in ORDER] + list(extra)
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, y), ncol=len(h), fontsize=6.3,
               handlelength=1.8, columnspacing=1.5, handletextpad=0.4)


def grid_fig(draw, ylabel, fname, note, h=2.9):
    fig, axs = plt.subplots(2, 7, figsize=(7.1, h), sharex=True, sharey=True,
                            gridspec_kw=dict(wspace=0.12, hspace=0.18))
    for r, (plat, pname) in enumerate(PLATS):
        for c, line in enumerate(LINES):
            ax = axs[r, c]
            dss = datasets(plat, line)
            if not dss:
                ax.axis('off'); ax.text(0.5, 0.5, 'not sequenced', transform=ax.transAxes, ha='center', va='center',
                                        fontsize=5.6, color=MUTED); continue
            draw(ax, plat, dss)
            if r == 0:
                ax.set_title(line, fontsize=6.8, pad=3, color=INK)
            if c == 0:
                ax.set_ylabel(pname, fontsize=6.8, color=INK, labelpad=3)
            ax.tick_params(labelsize=5.6)
    fig.supxlabel('Expected DNA fraction', fontsize=6.5, y=-0.01)
    fig.text(0.055, 0.5, ylabel, rotation=90, ha='center', va='center', fontsize=6.5)
    second = Line2D([], [], color=MUTED, lw=0.8, ls=(0, (2, 1.5)), marker='o', ms=3.6, mfc='white', mec=MUTED, mew=0.7,
                    label='second dataset')
    legend(fig, 1.04, [second])
    fig.text(0.5, -0.07, note, ha='center', fontsize=5.6, color=INK2)
    fig.savefig(os.path.join(OUT, fname), bbox_inches='tight', facecolor='white', dpi=220)


def series(ax, plat, dss, yfun):
    for m in ORDER[::-1]:
        z = {'LongPhase-TO': 5, 'PURPLE (TO)': 4, 'ASCAT (T+N)': 3, 'ASCAT (TO)': 2}[m]
        for k, ds in enumerate(dss):
            q = get(plat, ds, m)
            if not q:
                continue
            x = np.array([p['expected'] for p in q]) + (0.02 if k else 0)
            y = [yfun(p) for p in q]
            lw = 1.1 if m == 'LongPhase-TO' else 0.75
            ax.plot(x, y, color=COL[m], lw=lw, ls='-' if k == 0 else (0, (2, 1.5)), zorder=z, alpha=1 if m != 'ASCAT (TO)' else 0.8)
            ax.scatter(x, y, s=9, marker=MK[m], facecolor=COL[m] if k == 0 else 'white', edgecolor='white' if k == 0 else COL[m],
                       lw=0.4 if k == 0 else 0.7, zorder=z + 0.5)


# ---------- option 1: signed error ----------
def opt1_draw(ax, plat, dss):
    ax.axhspan(-0.05, 0.05, color=BAND, lw=0, zorder=0); ax.axhline(0, color=MUTED, lw=0.5, zorder=1)
    series(ax, plat, dss, lambda p: p['estimate'] - p['expected'])
    ax.set_xlim(0.1, 1.1); ax.set_ylim(-0.65, 0.88); ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([-0.5, 0, 0.5])
    ax.grid(axis='y', color=GRID, lw=0.4); ax.set_axisbelow(True)


# ---------- option 3: identity, quiet ----------
def opt3_draw(ax, plat, dss):
    ax.plot([0, 1.1], [0, 1.1], color=MUTED, lw=0.5, zorder=1)
    series(ax, plat, dss, lambda p: p['estimate'])
    ax.set_xlim(0.1, 1.08); ax.set_ylim(0.1, 1.08); ax.set_xticks([0.2, 0.6, 1.0]); ax.set_yticks([0.2, 0.6, 1.0])
    ax.set_aspect('equal', adjustable='box')


# ---------- option 2: heatmap ----------
def option2():
    rows = [('ONT', d) for d in ONT] + [('PacBio HiFi', d) for d in PB]
    lev = [0.2, 0.4, 0.6, 0.8, 1.0]
    cmap = LinearSegmentedColormap.from_list('div', ['#3d5fa8', '#9db3d9', '#f3f2ee', '#e7a98f', '#b8402a'])
    norm = TwoSlopeNorm(vcenter=0, vmin=-0.6, vmax=0.6)
    fig = plt.figure(figsize=(7.1, 3.7))
    ax = fig.add_axes([0.14, 0.08, 0.76, 0.8])
    gapx, gapy = 0.6, 0.6
    for b, m in enumerate(ORDER):
        x0 = b * (5 + gapx)
        ax.text(x0 + 2.5, -1.25, m, ha='center', va='bottom', fontsize=6.8, color=INK, fontweight='bold')
        for j, lv in enumerate(lev):
            ax.text(x0 + j + 0.5, -0.25, f'{lv:.1f}', ha='center', va='bottom', fontsize=5.6, color=INK2)
        for i, (plat, ds) in enumerate(rows):
            y = i + (gapy if plat != 'ONT' else 0)
            q = {p['level']: p for p in get(plat, ds, m)}
            for j, lv in enumerate(lev):
                if lv in q:
                    e = q[lv]['estimate'] - q[lv]['expected']
                    ax.add_patch(Rectangle((x0 + j + 0.04, y + 0.04), 0.92, 0.92, color=cmap(norm(np.clip(e, -0.6, 0.6))), lw=0))
                    if abs(e) >= 0.1:
                        ax.text(x0 + j + 0.5, y + 0.52, f'{e:+.2f}'.replace('+0.', '+.').replace('-0.', '−.'),
                                ha='center', va='center', fontsize=4.6, color='white' if abs(e) > 0.35 else INK)
                else:
                    ax.add_patch(Rectangle((x0 + j + 0.04, y + 0.04), 0.92, 0.92, facecolor='white', edgecolor=GRID,
                                           hatch='//////', lw=0.3))
    for i, (plat, ds) in enumerate(rows):
        y = i + (gapy if plat != 'ONT' else 0)
        ax.text(-0.25, y + 0.5, ds.replace('_', ' '), ha='right', va='center', fontsize=5.8, color=INK)
    ax.text(-5.2, 4.5, 'Nanopore', rotation=90, ha='center', va='center', fontsize=6.8, fontweight='bold', color=INK)
    ax.text(-5.2, 9 + gapy + 3, 'PacBio HiFi', rotation=90, ha='center', va='center', fontsize=6.8, fontweight='bold', color=INK)
    W = 4 * 5 + 3 * gapx
    ax.set_xlim(-0.1, W + 0.1); ax.set_ylim(len(rows) + gapy + 0.1, -1.4); ax.axis('off')
    ax.text(W / 2, len(rows) + gapy + 0.55, 'Expected DNA fraction (columns)', ha='center', va='top', fontsize=6.3, color=INK2)
    cax = fig.add_axes([0.925, 0.3, 0.012, 0.42])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[-0.6, -0.3, 0, 0.3, 0.6])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.6, length=2)
    cb.set_label('Estimated minus expected', fontsize=6)
    fig.text(0.14, 0.015, 'Numbers: errors of 0.10 or more.  Hatched: no result (ASCAT failed or was not run).',
             fontsize=5.6, color=INK2)
    fig.savefig(os.path.join(OUT, 'supp_percell_option_2.png'), bbox_inches='tight', facecolor='white', dpi=220)


note = 'COLO829 and HCC1395 (nanopore): filled = ONT / HKU dataset, open dashed = NYGC dataset.  ASCAT was not run on HG008.'
grid_fig(opt1_draw, 'Estimated minus expected DNA fraction', 'supp_percell_option_1.png', note + '  Band: ±0.05.')
option2()
grid_fig(opt3_draw, 'Estimated DNA fraction', 'supp_percell_option_3.png', note, h=2.45)
