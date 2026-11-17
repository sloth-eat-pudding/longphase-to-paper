"""Three layout options for Figure 4 (panels a: DNA fraction, b: ploidy, c: HG008). Draft PNGs for the author to choose.

Data loading is shared with fig4_draft.py (same sources and rules). Writes drafts/fig4_option_{A,B,C}.png.
Colors: LongPhase-TO blue, PURPLE orange, ASCAT (T+N) aqua (slots 1-3, validated all-pairs), ASCAT (TO) neutral gray.
"""
import csv, os, re, importlib.util
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.join(HERE, '..', '..', 'Revision_2026-11')
SUPP = os.path.join(HERE, '..', 'supplementary.tex')
OUT = os.path.join(HERE, 'drafts')

INK, INK2, MUTED, GRID, BAND = '#1a1a19', '#52514e', '#8a8985', '#e9e8e4', '#f1f0ec'
COL = {'LongPhase-TO': '#2a78d6', 'PURPLE (TO)': '#eb6834', 'ASCAT (T+N)': '#1baf7a', 'ASCAT (TO)': '#a3a29d'}
MK = {'LongPhase-TO': 'o', 'PURPLE (TO)': 's', 'ASCAT (T+N)': '^', 'ASCAT (TO)': 'v'}
METHODS = ['ASCAT (T+N)', 'ASCAT (TO)', 'PURPLE (TO)', 'LongPhase-TO']
plt.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5,
                     'axes.edgecolor': INK2, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
                     'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 400})


def tsv(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith('#')), delimiter='\t'))


def supp_table(label):
    s = open(SUPP).read()
    i = s.index('\\label{%s}' % label); j = s.index('\\end{longtable}', i)
    return [[c.strip().replace('\\_', '_') for c in l.strip().rstrip('\\').split('&')]
            for l in s[i:j].splitlines() if l.count('&') == 6 and re.match(r'^\s*[0-9.]+\s*$', l.split('&')[1])]


# ---------- data ----------
pts = []   # dict(method, platform, dataset, expected, estimate)
for label, plat in [('tab:tumor-purity', 'ONT'), ('tab:pacbio-tumor-purity', 'PacBio HiFi')]:
    for r in supp_table(label):
        for m, k in zip(['ASCAT (T+N)', 'ASCAT (TO)', 'PURPLE (TO)', 'LongPhase-TO'], [2, 3, 4, 6]):
            if r[k] != 'n.a.':
                pts.append(dict(method=m, platform=plat, dataset=r[0], expected=float(r[1]), estimate=float(r[k])))
LEVELS = [0.2, 0.4, 0.6, 0.8, 1.0]
for p in pts:
    p['level'] = min(LEVELS, key=lambda v: abs(v - p['expected'])); p['err'] = abs(p['estimate'] - p['expected'])

pl = [r for r in tsv(os.path.join(REV, 'ploidy_accuracy_ONT_PacBio.tsv'))
      if r['method'] == 'LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)' and 'dedup' not in r['dataset']]
LINES = ['COLO829', 'H1437', 'H2009', 'HCC1395', 'HCC1937', 'HCC1954']
hg = tsv(os.path.join(REV, 'hg008_ONT_cn_loh.tsv'))
HG = {m: sorted([r for r in hg if r['method'].startswith(m[:6])], key=lambda r: float(r['truth_tdf']))
      for m in ['PURPLE (TO)', 'LongPhase-TO']}


def label(ax, letter, x=-0.18, y=1.08):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom', ha='left', color=INK)


def method_legend(fig, methods, y=1.0, ncol=4):
    h = [Line2D([], [], ls='', marker=MK[m], ms=4, mfc=COL[m], mec='white', mew=0.4, label=m) for m in methods]
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, y), ncol=ncol, fontsize=6.5, handletextpad=0.2, columnspacing=1.4)


def hg008_panels(fig, ax1, ax2):
    for m in ['PURPLE (TO)', 'LongPhase-TO']:
        rs = HG[m]; x = [float(r['truth_tdf']) for r in rs]
        ax1.plot(x, [float(r['estimated_tdf']) for r in rs], color=COL[m], lw=1, marker=MK[m], ms=3.2, mfc=COL[m], mec='white', mew=0.4)
        ax2.plot(x, [float(r['ploidy']) for r in rs], color=COL[m], lw=1, marker=MK[m], ms=3.2, mfc=COL[m], mec='white', mew=0.4)
    ax1.plot([0, 1], [0, 1], color=GRID, lw=0.8, zorder=0); ax1.set_xlim(0.1, 1.05); ax1.set_ylim(0.1, 1.05)
    ax1.set_xlabel('Expected DNA fraction'); ax1.set_ylabel('Estimated DNA fraction'); ax1.set_title('HG008, DNA fraction', loc='left')
    ax2.axhline(1.47, color=MUTED, lw=0.6, ls=(0, (3, 2)), zorder=0)
    ax2.text(0.13, 1.52, 'NIST mean CN 1.47', fontsize=5.5, color=MUTED, ha='left', va='bottom')
    ax2.set_xlim(0.1, 1.05); ax2.set_ylim(1.3, 1.9); ax2.set_xlabel('Expected DNA fraction'); ax2.set_ylabel('Ploidy')
    ax2.set_title('HG008, ploidy', loc='left')
    for a in (ax1, ax2):
        a.set_xticks([0.2, 0.6, 1.0])


# ================= Option A: error by dilution level (strip) + ploidy deviation =================
def option_a():
    fig = plt.figure(figsize=(7.1, 4.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.2, 1, 1], height_ratios=[1, 1], hspace=0.55, wspace=0.42)
    ax = fig.add_subplot(gs[0, :])
    rng = np.random.default_rng(3); w = 0.17
    for li, lv in enumerate(LEVELS):
        if li % 2 == 0:
            ax.axvspan(li - 0.5, li + 0.5, color='#f7f6f3', lw=0, zorder=0)
        for mi, m in enumerate(METHODS):
            e = [p['err'] for p in pts if p['method'] == m and p['level'] == lv]
            x0 = li + (mi - 1.5) * w
            ax.scatter(x0 + rng.uniform(-0.05, 0.05, len(e)), e, s=6, marker=MK[m], color=COL[m], lw=0, alpha=0.85, zorder=2)
            ax.plot([x0 - 0.07, x0 + 0.07], [np.median(e)] * 2, color=INK, lw=0.9, zorder=3)
    ax.set_xticks(range(5)); ax.set_xticklabels(['0.2', '0.4', '0.6', '0.8', '1.0']); ax.set_xlim(-0.5, 4.5)
    ax.set_xlabel('Expected DNA fraction'); ax.set_ylabel('Absolute error'); ax.set_ylim(-0.02, 0.85)
    ax.grid(axis='y', color=GRID, lw=0.5); ax.set_axisbelow(True)
    mae = {m: np.mean([p['err'] for p in pts if p['method'] == m]) for m in METHODS}
    method_legend(fig, METHODS, y=1.0)
    ax.text(1.0, 1.02, '   '.join(f'{m} {mae[m]:.3f}' for m in METHODS), transform=ax.transAxes, ha='right', va='bottom',
            fontsize=5.8, color=INK2)
    ax.text(0.0, 1.02, 'Mean absolute error:', transform=ax.transAxes, ha='left', va='bottom', fontsize=5.8, color=INK2)
    label(ax, 'a', x=-0.07)
    # b: ploidy deviation
    axb = fig.add_subplot(gs[1, 0])
    axb.axhspan(-0.2, 0.2, color=BAND, lw=0, zorder=0); axb.axhline(0, color=MUTED, lw=0.5, zorder=1)
    for j, cl in enumerate(LINES):
        rs = [r for r in pl if r['dataset'].startswith(cl)]
        for off, key, m in [(-0.15, 'purple_java_ploidy', 'PURPLE (TO)'), (0.15, 'ploidy', 'LongPhase-TO')]:
            v = [float(r[key]) - float(r['literature_ploidy']) for r in rs if r[key] not in ('NA', '')]
            axb.scatter(j + off + rng.uniform(-0.05, 0.05, len(v)), v, s=7, marker=MK[m], color=COL[m], lw=0, alpha=0.9, zorder=2)
    axb.set_xticks(range(6)); axb.set_xticklabels(LINES, rotation=30, ha='right'); axb.set_ylim(-2.2, 2.2)
    axb.set_ylabel('Ploidy minus published'); axb.set_title('Ploidy (band: ±0.2)', loc='left')
    label(axb, 'b', x=-0.16)
    ax1 = fig.add_subplot(gs[1, 1]); ax2 = fig.add_subplot(gs[1, 2]); hg008_panels(fig, ax1, ax2); label(ax1, 'c', x=-0.38)
    fig.savefig(os.path.join(OUT, 'fig4_option_A.png'), bbox_inches='tight', facecolor='white')


# ================= Option B: identity scatter per platform + ploidy scatter =================
def option_b():
    fig = plt.figure(figsize=(7.1, 4.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1], hspace=0.55, wspace=0.55)
    for k, plat in enumerate(['ONT', 'PacBio HiFi']):
        ax = fig.add_subplot(gs[0, 2 * k:2 * k + 2])
        ax.plot([0, 1.05], [0, 1.05], color=GRID, lw=0.8, zorder=0)
        for m in METHODS:
            q = [p for p in pts if p['method'] == m and p['platform'] == plat]
            jit = {'ASCAT (T+N)': -0.018, 'ASCAT (TO)': -0.006, 'PURPLE (TO)': 0.006, 'LongPhase-TO': 0.018}[m]
            ax.scatter([p['expected'] + jit for p in q], [p['estimate'] for p in q], s=9, marker=MK[m], color=COL[m],
                       edgecolor='white', lw=0.3, zorder=3 if m == 'LongPhase-TO' else 2)
        ax.set_xlim(0.1, 1.06); ax.set_ylim(0, 1.06); ax.set_aspect('equal')
        ax.set_xticks(LEVELS); ax.set_xlabel('Expected DNA fraction'); ax.set_ylabel('Estimated DNA fraction')
        ax.set_title(plat, loc='left')
        y = 0.97
        for m in METHODS:
            e = np.mean([p['err'] for p in pts if p['method'] == m and p['platform'] == plat])
            ax.text(1.03, y, f'{e:.3f}', transform=ax.transAxes, color=COL[m] if m != 'ASCAT (TO)' else INK2, fontsize=5.8,
                    ha='left', va='top', fontweight='bold' if m == 'LongPhase-TO' else 'normal'); y -= 0.07
        ax.text(1.03, 1.04, 'MAE', transform=ax.transAxes, fontsize=5.8, color=INK2, ha='left', va='top')
        if k == 0:
            label(ax, 'a', x=-0.25, y=1.02)
    method_legend(fig, METHODS, y=1.0)
    axb = fig.add_subplot(gs[1, 0:2])
    xs = np.linspace(2, 4.7, 2); axb.fill_between(xs, xs - 0.2, xs + 0.2, color=BAND, lw=0, zorder=0)
    axb.plot(xs, xs, color=MUTED, lw=0.5, zorder=1)
    rng = np.random.default_rng(5)
    for key, m in [('purple_java_ploidy', 'PURPLE (TO)'), ('ploidy', 'LongPhase-TO')]:
        q = [r for r in pl if r[key] not in ('NA', '')]
        axb.scatter([float(r['literature_ploidy']) + rng.uniform(-0.04, 0.04) + (0.05 if m == 'LongPhase-TO' else -0.05) for r in q],
                    [float(r[key]) for r in q], s=8, marker=MK[m], color=COL[m], edgecolor='white', lw=0.3, zorder=2)
        n = sum(round(abs(float(r[key]) - float(r['literature_ploidy'])), 2) <= 0.2 for r in q)
        axb.text(0.02, 0.95 - (0.09 if m == 'LongPhase-TO' else 0), f'{m}: {n}/{len(q)} within ±0.2', transform=axb.transAxes,
                 color=COL[m], fontsize=5.8, va='top')
    for cl in LINES:
        lit = float([r for r in pl if r['dataset'].startswith(cl)][0]['literature_ploidy'])
        axb.text(lit, 1.15, cl, rotation=90, fontsize=5, color=MUTED, ha='center', va='bottom')
    axb.set_xlim(2.1, 4.7); axb.set_ylim(1.1, 5.0); axb.set_xlabel('Published ploidy'); axb.set_ylabel('Estimated ploidy')
    label(axb, 'b', x=-0.12)
    ax1 = fig.add_subplot(gs[1, 2]); ax2 = fig.add_subplot(gs[1, 3]); hg008_panels(fig, ax1, ax2); label(ax1, 'c', x=-0.4)
    fig.savefig(os.path.join(OUT, 'fig4_option_B.png'), bbox_inches='tight', facecolor='white')


# ================= Option C: error heatmaps =================
def option_c():
    ds = [r[0] for r in supp_table('tab:tumor-purity')][::5] + [r[0] for r in supp_table('tab:pacbio-tumor-purity')][::5]
    plat_of = {p['dataset']: p['platform'] for p in pts}
    fig = plt.figure(figsize=(7.1, 4.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[3.2, 1.6, 1.25], height_ratios=[1, 1], hspace=0.45, wspace=0.35)
    ax = fig.add_subplot(gs[:, 0])
    seq = LinearSegmentedColormap.from_list('seq', ['#f4f8fd', '#9ec5f0', '#2a78d6', '#0d3a73'])
    M = np.full((len(ds), 20), np.nan)
    for i, d in enumerate(ds):
        for mi, m in enumerate(METHODS):
            for li, lv in enumerate(LEVELS):
                q = [p['err'] for p in pts if p['dataset'] == d and p['method'] == m and p['level'] == lv]
                if q: M[i, mi * 5 + li] = q[0]
    im = ax.imshow(M, cmap=seq, vmin=0, vmax=0.6, aspect='auto', interpolation='nearest')
    for i, j in zip(*np.where(np.isnan(M))):
        ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor='white', edgecolor=GRID, hatch='////', lw=0))
    for b in range(1, 4): ax.axvline(b * 5 - 0.5, color='white', lw=2.5)
    ax.axhline(7.5, color='white', lw=2.5)
    ax.set_yticks(range(len(ds))); ax.set_yticklabels([d.replace('_CASTLE', '').replace('_PacBio', '') for d in ds], fontsize=5.8)
    ax.set_xticks(range(20)); ax.set_xticklabels([f'{v:.1f}' for v in LEVELS] * 4, fontsize=5)
    for mi, m in enumerate(METHODS):
        mae = np.nanmean(M[:, mi * 5:mi * 5 + 5])
        ax.text(mi * 5 + 2, -1.0, m, ha='center', va='bottom', fontsize=6.3, color=INK, fontweight='bold' if m == 'LongPhase-TO' else 'normal')
        ax.text(mi * 5 + 2, -0.55, f'MAE {mae:.3f}', ha='center', va='bottom', fontsize=5.5, color=INK2)
    ax.set_xlabel('Expected DNA fraction'); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.text(-0.33, 3.5, 'Nanopore', rotation=90, transform=ax.get_yaxis_transform(), ha='center', va='center', fontsize=6, color=INK2)
    ax.text(-0.33, 10.5, 'PacBio HiFi', rotation=90, transform=ax.get_yaxis_transform(), ha='center', va='center', fontsize=6, color=INK2)
    cb = fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.035, pad=0.11, aspect=30)
    cb.set_label('Absolute error of DNA fraction', fontsize=6); cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.5, length=2)
    label(ax, 'a', x=-0.2, y=1.07)
    # b: signed ploidy deviation heatmap
    axb = fig.add_subplot(gs[:, 1])
    div = LinearSegmentedColormap.from_list('div', ['#c2402f', '#f2c4b8', '#f1f0ec', '#b9d3f2', '#2a5fb0'])
    P = np.full((len(ds), 10), np.nan)
    order = {'ONT': 'ONT', 'PacBio HiFi': 'PacBio'}
    for i, d in enumerate(ds):
        dd = d.replace('_CASTLE', '_CMH'); dd = dd + '_with_duplicate_reads' if dd == 'HCC1937_UCSC' else dd
        rs = sorted([r for r in pl if r['dataset'] == dd], key=lambda r: float(r['target_tdf']))
        for k, key in enumerate(['purple_java_ploidy', 'ploidy']):
            for li, r in enumerate(rs):
                if r[key] not in ('NA', ''): P[i, k * 5 + li] = round(float(r[key]) - float(r['literature_ploidy']), 2)
    im2 = axb.imshow(P, cmap=div, norm=TwoSlopeNorm(0, -2, 2), aspect='auto', interpolation='nearest')
    axb.axvline(4.5, color='white', lw=2.5); axb.axhline(7.5, color='white', lw=2.5)
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if not np.isnan(P[i, j]) and abs(P[i, j]) > 0.2 + 1e-9:
                axb.text(j, i, '×', ha='center', va='center', fontsize=5, color=INK)
    axb.set_yticks([]); axb.set_xticks(range(10)); axb.set_xticklabels([f'{v:.1f}' for v in LEVELS] * 2, fontsize=5)
    for k, m in enumerate(['PURPLE (TO)', 'LongPhase-TO']):
        n = int(np.sum(np.abs(P[:, k * 5:k * 5 + 5]) <= 0.2 + 1e-9)); tot = int(np.sum(~np.isnan(P[:, k * 5:k * 5 + 5])))
        axb.text(k * 5 + 2, -1.0, m, ha='center', va='bottom', fontsize=6.3, color=INK, fontweight='bold' if k else 'normal')
        axb.text(k * 5 + 2, -0.55, f'{n}/{tot} within ±0.2', ha='center', va='bottom', fontsize=5.5, color=INK2)
    axb.set_xlabel('Expected DNA fraction'); axb.tick_params(length=0)
    for sp in axb.spines.values(): sp.set_visible(False)
    cb2 = fig.colorbar(im2, ax=axb, orientation='horizontal', fraction=0.035, pad=0.11, aspect=20)
    cb2.set_label('Ploidy minus published  (× : off by >0.2)', fontsize=6); cb2.outline.set_visible(False); cb2.ax.tick_params(labelsize=5.5, length=2)
    label(axb, 'b', x=-0.08, y=1.07)
    ax1 = fig.add_subplot(gs[0, 2]); ax2 = fig.add_subplot(gs[1, 2]); hg008_panels(fig, ax1, ax2); label(ax1, 'c', x=-0.35)
    h = [Line2D([], [], color=COL[m], marker=MK[m], ms=3, lw=1, mec='white', mew=0.3, label=m) for m in ['PURPLE (TO)', 'LongPhase-TO']]
    ax1.legend(handles=h, loc='lower right', fontsize=5.5)
    fig.savefig(os.path.join(OUT, 'fig4_option_C.png'), bbox_inches='tight', facecolor='white')


os.makedirs(OUT, exist_ok=True)
option_a(); option_b(); option_c()
