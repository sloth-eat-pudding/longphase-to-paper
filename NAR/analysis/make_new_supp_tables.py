"""Write two new Supplementary tables (ploidy accuracy, HG008) into ../supplementary.tex after tab:pacbio-tumor-purity.
Sources: ../../Revision_2026-11/ploidy_accuracy_ONT_PacBio.tsv and hg008_ONT_cn_loh.tsv.
LongPhase-TO = DeepSomatic-TO + ClairS-TO PoN configuration. HCC1937 ONT = with-duplicate-reads set.
Idempotent: replaces the block between the two marker comments if present."""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.join(HERE, '..', '..', 'Revision_2026-11'); SUPP = os.path.join(HERE, '..', 'supplementary.tex')
tsv = lambda n: list(csv.DictReader((l for l in open(os.path.join(REV, n)) if not l.startswith('#')), delimiter='\t'))
NAME = lambda d: d.replace('_CMH', '_CASTLE').replace('_with_duplicate_reads', '').replace('_', '\\_')
P = [r for r in tsv('ploidy_accuracy_ONT_PacBio.tsv') if r['method'] == 'LongPhase-TO(DeepSomatic-TO+ClairS-TO PoN)' and 'dedup' not in r['dataset']]
rows, last = [], None
for r in P:
    d = NAME(r['dataset'])
    if last and d != last: rows.append('  \\addlinespace')
    last = d
    pj = 'n.a.' if r['purple_java_ploidy'] in ('NA', '') else '%.2f' % float(r['purple_java_ploidy'])
    rows.append('  %s & %s & %.1f & %.2f & %s & %.2f%s \\\\' % (d, r['platform'].replace('PacBio', 'PacBio HiFi'), float(r['target_tdf']),
                float(r['literature_ploidy']), pj, float(r['ploidy']), '$^{*}$' if r['scale_refit_2x_applied'] == 'True' else ''))
t1 = r"""\begin{small}
\begin{longtable}{@{}llccccc@{}}
  \caption[Tumor ploidy estimates]%
  {Per-mixture tumor ploidy estimated by PURPLE (TO) and LongPhase-TO on both platforms, against the published consensus ploidy of each cell line. An estimate within 0.2 of the published value counts as correct. LongPhase-TO values are from the DeepSomatic-TO configuration with the ClairS-TO normal panels. An asterisk marks runs changed by the ploidy scale correction.}
  \label{tab:ploidy-accuracy}\\
  \toprule
  Dataset & Platform & Fraction & Published & PURPLE (TO) & LongPhase-TO \\
  \midrule
  \endfirsthead
  \toprule
  Dataset & Platform & Fraction & Published & PURPLE (TO) & LongPhase-TO \\
  \midrule
  \endhead
  \bottomrule
  \endlastfoot
""" + '\n'.join(rows) + "\n\\end{longtable}\n\\end{small}\n"
H = tsv('hg008_ONT_cn_loh.tsv'); hr = []
for r in H:
    m = 'PURPLE (TO)' if r['method'].startswith('PURPLE') else 'LongPhase-TO'
    hr.append('    %.3f & %s & %.3f & %.2f & %.3f & %.3f & %.3f & %.3f \\\\' % (float(r['truth_tdf']), m, float(r['estimated_tdf']), float(r['ploidy']),
              float(r['mean_cn']), float(r['cn_within_0.5']), float(r['loh_recall']), float(r['loh_precision'])))
t2 = r"""\begin{table}[htbp]
  \centering
  \caption[HG008 external test]%
  {Tumor DNA fraction, ploidy, copy number and LOH of the HG008 nanopore mixtures, evaluated against the NIST HG008-T somatic copy-number draft benchmark v0.6 (PASS regions of chromosomes 1 to 22, chromosome 10 without truth). The NIST mean copy number is 1.47. CN within 0.5 is the fraction of evaluated length whose total copy number lies within 0.5 of the truth. LOH recall and precision are computed by length, with LOH defined as minor copy number below 0.5.}
  \label{tab:hg008}
  \small
  \begin{tabular}{@{}lccccccc@{}}
    \toprule
    Expected & Method & DNA fraction & Ploidy & Mean CN & CN within 0.5 & LOH recall & LOH precision \\
    \midrule
""" + '\n'.join(hr) + "\n    \\bottomrule\n  \\end{tabular}\n\\end{table}\n"
B, E = '% >>> revision tables (make_new_supp_tables.py)\n', '% <<< revision tables\n'
s = open(SUPP).read()
if B in s:
    s = s[:s.index(B)] + s[s.index(E) + len(E):]
i = s.index('\\label{tab:pacbio-tumor-purity}'); j = s.index('\\end{small}', i) + len('\\end{small}\n')
s = s[:j] + '\n' + B + t1 + '\n' + t2 + E + s[j:]
open(SUPP, 'w').write(s); print(len(rows), 'ploidy rows,', len(hr), 'HG008 rows')
