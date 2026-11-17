"""Full Figure 4 drafts: panels (a, b) from fig4ab_merged.py (layout A or B) + panel (c) variant B from fig4c_hg008.py.

Panel (c) uses the real HG008 segment tables (Revision_2026-11/hg008_segments/, issue #5). Writes drafts/fig4_full_{A,B}.png.
"""
import importlib.util, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


ab = load('fig4ab_merged')
c = load('fig4c_hg008')

for lay, top_h, fn in [('A', 3.6, ab.merged_a), ('B', 2.55, ab.merged_b)]:
    fig = plt.figure(figsize=(7.1, top_h + 2.5))
    top, bot = fig.subfigures(2, 1, height_ratios=[top_h, 2.5], hspace=0.04)
    fn(top)
    c.panel_c(bot, variant='B')
    fig.savefig(os.path.join(HERE, 'drafts', f'fig4_full_{lay}.png'), bbox_inches='tight', facecolor='white')

# ---- final Figure 4 (layout A, panel c variant B: tolerance curves, pure tumor bold + four mixtures thin).
fig = plt.figure(figsize=(7.1, 3.25 + 2.5))
top, bot = fig.subfigures(2, 1, height_ratios=[3.25, 2.5], hspace=0.0)
ab.merged_a(top)
c.panel_c(bot, variant='B')
fig.savefig(os.path.join(HERE, '..', 'Figure_4.pdf'), bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(HERE, 'drafts', 'Figure_4_final.png'), bbox_inches='tight', facecolor='white', dpi=200)
