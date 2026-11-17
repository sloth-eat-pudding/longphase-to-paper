"""Figure 3 style options (drafts): the per-dataset F1 layout of fig3_per_dataset_draft.py in three journal styles.

1 Nature (NPG): navy #3C5488 / red #E64B35. Caller alone = same hue, thin dashed, no markers.
2 Okabe-Ito (colorblind-safe): blue #0072B2 / vermillion #D55E00. Caller alone = light tint of the same hue, solid.
3 Before gray, after color (AAAS navy #3B4992 / red #EE0000): both callers alone in gray, marker shape = caller.
All: open axes (left and bottom only), thin lines, small markers. Writes drafts/fig3_style_{1,2,3}.png.
"""
import importlib.util, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('f3', os.path.join(HERE, 'fig3_recalibration.py'))
f3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f3)
FIRST = ('_ONT', '_HKU')
LABS = ['ClairS-TO', 'ClairS-TO+LongPhase-TO', 'DeepSomatic-TO', 'DeepSomatic-TO+LongPhase-TO']


def S(color, ls='-', lw=1.0, alpha=1.0, marker=None, ms=2.6, mfc=None, z=3):
    return dict(color=color, ls=ls, lw=lw, alpha=alpha, marker=marker, ms=ms, mfc=mfc or color, z=z)


STYLES = {
    1: {'ClairS-TO': S('#3C5488', (0, (2.5, 1.5)), 0.8, 0.6, z=2),
        'ClairS-TO+LongPhase-TO': S('#3C5488', '-', 1.1, 1, 'o'),
        'DeepSomatic-TO': S('#E64B35', (0, (2.5, 1.5)), 0.8, 0.6, z=2),
        'DeepSomatic-TO+LongPhase-TO': S('#E64B35', '-', 1.1, 1, 's')},
    2: {'ClairS-TO': S('#9ECAE9', '-', 1.0, 1, z=2),
        'ClairS-TO+LongPhase-TO': S('#0072B2', '-', 1.1, 1, 'o'),
        'DeepSomatic-TO': S('#F5B88F', '-', 1.0, 1, z=2),
        'DeepSomatic-TO+LongPhase-TO': S('#D55E00', '-', 1.1, 1, 's')},
    3: {'ClairS-TO': S('#a3a3a3', (0, (2.5, 1.5)), 0.8, 1, 'o', 2.4, 'white', z=2),
        'ClairS-TO+LongPhase-TO': S('#3B4992', '-', 1.1, 1, 'o'),
        'DeepSomatic-TO': S('#a3a3a3', (0, (2.5, 1.5)), 0.8, 1, 's', 2.4, 'white', z=2),
        'DeepSomatic-TO+LongPhase-TO': S('#EE0000', '-', 1.1, 1, 's')},
}


def names(d):
    return sorted({k[0] for k in d}, key=lambda x: (f3.ORDER.index(x.split('_')[0]), not x.endswith(FIRST), x))


def draw(k):
    st = STYLES[k]
    fig = plt.figure(figsize=(7.1, 5.2))
    outer = fig.add_gridspec(2, 1, hspace=0.3)
    for b, (vt, top) in enumerate([('SNV', 1.0), ('indel', 0.8)]):
        blk = outer[b].subgridspec(2, 1, hspace=0.45)
        for r, plat in enumerate(['ONT', 'PacBio']):
            d = f3.DATA[(plat, vt)]; ds = names(d)
            row = blk[r].subgridspec(1, len(ds), wspace=0.12)
            first = None
            for c, name in enumerate(ds):
                ax = fig.add_subplot(row[0, c], sharey=first); first = first or ax
                for lab in LABS:
                    s = st[lab]; y = [d[(name, t, lab)]['f1'] for t in f3.T]
                    ax.plot(f3.T, y, color=s['color'], ls=s['ls'], lw=s['lw'], alpha=s['alpha'], zorder=s['z'],
                            marker=s['marker'], ms=s['ms'], mfc=s['mfc'], mec=s['color'], mew=0.6)
                ax.set_xlim(0.12, 1.08); ax.set_ylim(0, top); ax.set_xticks([0.2, 0.6, 1.0])
                ax.set_yticks([0, 0.5, 1.0] if vt == 'SNV' else [0, 0.4, 0.8])
                ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
                for sp in ('left', 'bottom'):
                    ax.spines[sp].set_color('#7a7a7a'); ax.spines[sp].set_linewidth(0.5)
                ax.tick_params(labelsize=5.5, pad=1.5, length=2, color='#7a7a7a')
                ax.set_title(name.replace('_', ' '), fontsize=5.8, color='#222222', pad=2)
                if c > 0:
                    ax.tick_params(labelleft=False)
                else:
                    ax.set_ylabel(f"{'Nanopore' if plat == 'ONT' else 'PacBio HiFi'}\n{vt} F1", fontsize=6.3, color='#222222')
                    if r == 0:
                        ax.text(-0.75 if plat == 'ONT' else -0.55, 1.12, 'a' if vt == 'SNV' else 'b', transform=ax.transAxes,
                                fontsize=8, fontweight='bold', va='bottom')
    fig.supxlabel('Tumor DNA fraction', fontsize=6.5, y=0.045)
    h = [Line2D([], [], color=st[l]['color'], ls=st[l]['ls'], lw=st[l]['lw'] + 0.2, alpha=st[l]['alpha'], marker=st[l]['marker'],
                ms=3.2, mfc=st[l]['mfc'], mec=st[l]['color'], mew=0.6, label=l.replace('+', ' + ')) for l in LABS]
    fig.legend(handles=h, loc='upper center', bbox_to_anchor=(0.5, 0.965), ncol=4, fontsize=6.3,
               handlelength=2.4, columnspacing=1.6, handletextpad=0.4)
    fig.subplots_adjust(top=0.9)
    fig.savefig(os.path.join(HERE, 'drafts', f'fig3_style_{k}.png'), bbox_inches='tight', facecolor='white', dpi=220)


for k in STYLES:
    draw(k)
