"""Figure 4(a, b) with HG008 merged into the cell-line panels (layouts A and B). Draft PNGs for the author to choose.

Data loading is taken from fig4_options.py (same sources and rules). HG008 adds:
  (a) 5 nanopore mixtures, LongPhase-TO and PURPLE (TO) only (ASCAT was not run on HG008);
  (b) one more cell line, reference ploidy = NIST v0.6 mean total copy number (1.47).
Source for HG008: ../../Revision_2026-11/hg008_ONT_cn_loh.tsv. Writes drafts/fig4ab_merged_{A,B}.png.
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'fig4_options.py')).read()
exec(src[:src.index('# ================= Option A')])          # data, style, helpers

HG_REF = 1.47
for m in ['PURPLE (TO)', 'LongPhase-TO']:
    for r in HG[m]:
        e = float(r['truth_tdf'])
        pts.append(dict(method=m, platform='ONT', dataset='HG008', expected=e, estimate=float(r['estimated_tdf']),
                        level=min(LEVELS, key=lambda v: abs(v - e)), err=abs(float(r['estimated_tdf']) - e)))
for a, b in zip(HG['LongPhase-TO'], HG['PURPLE (TO)']):
    pl.append({'dataset': 'HG008', 'platform': 'ONT', 'target_tdf': a['truth_tdf'], 'literature_ploidy': str(HG_REF),
               'ploidy': a['ploidy'], 'purple_java_ploidy': b['ploidy']})
LINES7 = LINES + ['HG008']
NOTE = 'HG008 added: open symbols with dark outline'


def hg_style(p):
    return p.get('dataset') == 'HG008'


def within(rows, key):
    q = [r for r in rows if r[key] not in ('NA', '')]
    return sum(round(abs(float(r[key]) - float(r['literature_ploidy'])), 2) <= 0.2 for r in q), len(q)


# ================= A =================
def merged_a(fig=None):
    own = fig is None
    if own:
        fig = plt.figure(figsize=(7.1, 3.9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 0.9], hspace=0.6)
    ax = fig.add_subplot(gs[0]); rng = np.random.default_rng(3); w = 0.17
    for li, lv in enumerate(LEVELS):
        if li % 2 == 0:
            ax.axvspan(li - 0.5, li + 0.5, color='#f7f6f3', lw=0, zorder=0)
        for mi, m in enumerate(METHODS):
            q = [p for p in pts if p['method'] == m and p['level'] == lv]
            x0 = li + (mi - 1.5) * w
            for p in q:
                hgp = hg_style(p)
                ax.scatter(x0 + rng.uniform(-0.05, 0.05), p['err'], s=13 if hgp else 6, marker=MK[m],
                           facecolor='white' if hgp else COL[m], edgecolor=COL[m] if hgp else 'none', lw=0.8 if hgp else 0,
                           alpha=1 if hgp else 0.85, zorder=4 if hgp else 2)
            ax.plot([x0 - 0.07, x0 + 0.07], [np.median([p['err'] for p in q])] * 2, color=INK, lw=0.9, zorder=3)
    ax.set_xticks(range(5)); ax.set_xticklabels(['0.2', '0.4', '0.6', '0.8', '1.0']); ax.set_xlim(-0.5, 4.5)
    ax.set_xlabel('Expected DNA fraction'); ax.set_ylabel('Absolute error'); ax.set_ylim(-0.02, 0.85)
    ax.grid(axis='y', color=GRID, lw=0.5); ax.set_axisbelow(True)
    mae = {m: np.mean([p['err'] for p in pts if p['method'] == m]) for m in METHODS}
    method_legend(fig, METHODS, y=1.0)
    ax.text(0.0, 1.02, 'Mean absolute error:', transform=ax.transAxes, ha='left', va='bottom', fontsize=5.8, color=INK2)
    ax.text(1.0, 1.02, '   '.join(f'{m} {mae[m]:.3f}' for m in METHODS), transform=ax.transAxes, ha='right', va='bottom',
            fontsize=5.8, color=INK2)
    label(ax, 'a', x=-0.07)
    axb = fig.add_subplot(gs[1])
    axb.axhspan(-0.2, 0.2, color=BAND, lw=0, zorder=0); axb.axhline(0, color=MUTED, lw=0.5, zorder=1)
    for j, cl in enumerate(LINES7):
        rs = [r for r in pl if r['dataset'].startswith(cl)]
        for off, key, m in [(-0.15, 'purple_java_ploidy', 'PURPLE (TO)'), (0.15, 'ploidy', 'LongPhase-TO')]:
            v = [float(r[key]) - float(r['literature_ploidy']) for r in rs if r[key] not in ('NA', '')]
            hgp = cl == 'HG008'
            axb.scatter(j + off + rng.uniform(-0.05, 0.05, len(v)), v, s=11 if hgp else 7, marker=MK[m],
                        facecolor='white' if hgp else COL[m], edgecolor=COL[m] if hgp else 'none', lw=0.7, alpha=0.9, zorder=2)
    axb.axvline(5.5, color=GRID, lw=0.6)
    axb.set_xticks(range(7)); axb.set_xticklabels(LINES7); axb.set_ylim(-2.2, 2.2); axb.set_xlim(-0.5, 6.5)
    axb.set_ylabel('Ploidy minus reference')
    nP, tP = within(pl, 'purple_java_ploidy'); nL, tL = within(pl, 'ploidy')
    axb.set_title(f'Ploidy (band ±0.2).  Within ±0.2: LongPhase-TO {nL}/{tL}, PURPLE (TO) {nP}/{tP}', loc='left')
    axb.text(6, -2.1, 'reference:\nNIST mean CN', fontsize=5.2, color=MUTED, ha='center', va='bottom')
    label(axb, 'b', x=-0.07)
    if own:
        fig.savefig(os.path.join(OUT, 'fig4ab_merged_A.png'), bbox_inches='tight', facecolor='white')


# ================= B =================
def merged_b(fig=None):
    own = fig is None
    if own:
        fig = plt.figure(figsize=(7.1, 2.55))
    gs = fig.add_gridspec(1, 3, wspace=0.62)
    for k, plat in enumerate(['ONT', 'PacBio HiFi']):
        ax = fig.add_subplot(gs[0, k])
        ax.plot([0, 1.05], [0, 1.05], color=GRID, lw=0.8, zorder=0)
        for m in METHODS:
            jit = {'ASCAT (T+N)': -0.018, 'ASCAT (TO)': -0.006, 'PURPLE (TO)': 0.006, 'LongPhase-TO': 0.018}[m]
            for p in [p for p in pts if p['method'] == m and p['platform'] == plat]:
                hgp = hg_style(p)
                ax.scatter(p['expected'] + jit, p['estimate'], s=13 if hgp else 9, marker=MK[m],
                           facecolor='white' if hgp else COL[m], edgecolor=COL[m] if hgp else 'white', lw=0.7 if hgp else 0.3,
                           zorder=3 if m == 'LongPhase-TO' else 2)
        ax.set_xlim(0.1, 1.06); ax.set_ylim(0, 1.06); ax.set_aspect('equal')
        ax.set_xticks([0.2, 0.6, 1.0]); ax.set_xlabel('Expected DNA fraction'); ax.set_ylabel('Estimated DNA fraction')
        ax.set_title('Nanopore (incl. HG008)' if plat == 'ONT' else plat, loc='left')
        y = 0.40
        ax.text(0.97, y + 0.06, 'MAE', transform=ax.transAxes, fontsize=5.6, color=INK2, ha='right', va='top')
        for m in METHODS:
            e = np.mean([p['err'] for p in pts if p['method'] == m and p['platform'] == plat])
            ax.text(0.97, y, f'{e:.3f}', transform=ax.transAxes, color=COL[m] if m != 'ASCAT (TO)' else INK2, fontsize=5.6,
                    ha='right', va='top', fontweight='bold' if m == 'LongPhase-TO' else 'normal'); y -= 0.065
        if k == 0:
            label(ax, 'a', x=-0.32, y=1.04)
    method_legend(fig, METHODS, y=0.98)
    axb = fig.add_subplot(gs[0, 2])
    xs = np.array([1.2, 4.7]); axb.fill_between(xs, xs - 0.2, xs + 0.2, color=BAND, lw=0, zorder=0)
    axb.plot(xs, xs, color=MUTED, lw=0.5, zorder=1)
    rng = np.random.default_rng(5)
    for key, m in [('purple_java_ploidy', 'PURPLE (TO)'), ('ploidy', 'LongPhase-TO')]:
        for r in [r for r in pl if r[key] not in ('NA', '')]:
            hgp = r['dataset'] == 'HG008'
            axb.scatter(float(r['literature_ploidy']) + rng.uniform(-0.04, 0.04) + (0.05 if m == 'LongPhase-TO' else -0.05),
                        float(r[key]), s=11 if hgp else 8, marker=MK[m], facecolor='white' if hgp else COL[m],
                        edgecolor=COL[m] if hgp else 'white', lw=0.7 if hgp else 0.3, zorder=2)
        n, t = within(pl, key)
        axb.text(0.0, 1.03 + (0.0 if m == 'LongPhase-TO' else 0.08), f'{n}/{t} within ±0.2', transform=axb.transAxes,
                 color=COL[m], fontsize=5.6, va='bottom', ha='left', fontweight='bold' if m == 'LongPhase-TO' else 'normal')
    for cl in LINES7:
        ref = float([r for r in pl if r['dataset'].startswith(cl)][0]['literature_ploidy'])
        dx = {'H1437': -0.09, 'COLO829': 0.09}.get(cl, 0)
        axb.text(ref + dx, 0.62, cl, rotation=90, fontsize=4.8, color=MUTED, ha='center', va='bottom')
    axb.set_xlim(1.2, 4.7); axb.set_ylim(0.6, 5.0); axb.set_aspect('equal')
    axb.set_xticks([2, 3, 4]); axb.set_yticks([1, 2, 3, 4, 5])
    axb.set_xlabel('Reference ploidy'); axb.set_ylabel('Estimated ploidy'); 
    label(axb, 'b', x=-0.32, y=1.12)
    fig.text(0.5, -0.04, 'Open symbols: HG008 (reference ploidy = NIST mean copy number 1.47)', ha='center', fontsize=5.6, color=INK2)
    if own:
        fig.savefig(os.path.join(OUT, 'fig4ab_merged_B.png'), bbox_inches='tight', facecolor='white')
        # Supplementary Figure (per-platform DNA fraction and ploidy), used by ../supplementary.tex
        fig.savefig(os.path.join(HERE, 'supp_tdf_ploidy_by_platform.pdf'), bbox_inches='tight', facecolor='white')


if __name__ == '__main__':
    merged_a(); merged_b()
