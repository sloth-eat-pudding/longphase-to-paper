# Re-run the improved plotting code for the grouped bar chart with a placeholder for missing data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Rectangle, Patch

# Data frame as provided
df = pd.DataFrame({
    "Dataset": [
        "COLO829_ONT", "COLO829_NYGC", "H1437_UCSC", "H2009_UCSC",
        "HCC1395_HKU", "HCC1395_NYGC", "HCC1937_UCSC", "HCC1954_UCSC"
    ],
    "LongPhase-TO (ClairS-TO)": [
        0.268455916, 0.281093435, 0.346182712, 0.359340506,
        0.508410807, 0.507369938, 0.518483294, 0.095534204
    ],
    "LongPhase-TO (DeepSomatic-TO)": [
        0.271000861, 0.290482817, 0.347934183, 0.360863194,
        0.515299790, 0.518597019, 0.522037760, 0.102781198
    ],
    "ASCAT (T+N)": [
        0.353273552, 0.291117006, 0.395085909, 0.433933056,
        0.520102097, 0.526949994, 0.561487675, np.nan
    ],
    "ASCAT (TO)": [
        0.393071152, 0.332442567, 0.447017501, 0.482634227,
        0.657210037, 0.590737696, 0.599636339, 0.587205577
    ]
})

datasets = df["Dataset"].values
methods = ["ASCAT (T+N)", "ASCAT (TO)", "LongPhase-TO (ClairS-TO)", "LongPhase-TO (DeepSomatic-TO)"]

x = np.arange(len(datasets))
bar_width = 0.18
offsets = np.linspace(-1.5*bar_width, 1.5*bar_width, num=len(methods))

# Nature-style professional color palette (colorblind-friendly)
colors = ['#332288', '#88CCEE', '#CC6677', '#DDCC77']  # Purple, Cyan, Pink, Yellow

fig, ax = plt.subplots(figsize=(9, 5))

ymax = 0.0
for i, method in enumerate(methods):
    vals = df[method].values.astype(float)
    xpos = x + offsets[i]
    # Plot real bars where values exist
    mask = ~np.isnan(vals)
    ax.bar(xpos[mask], vals[mask], width=bar_width, label=method, color=colors[i])
    ymax = max(ymax, np.nanmax(vals))

    # Add value labels (percent) for non-NaN bars
    #for j, (xm, vm) in enumerate(zip(xpos[mask], vals[mask])):
    #    ax.text(xm, vm + 0.01, f"{vm:.0%}", ha="center", va="bottom", fontsize=8, rotation=90)

# Placeholder for missing ASCAT (T+N) in HCC1954_UCSC
missing_idx = np.where(datasets == "HCC1954_UCSC")[0][0]
tn_idx = methods.index("ASCAT (T+N)")
x_missing = x[missing_idx] + offsets[tn_idx]

placeholder_height = 0.02  # 2%
rect = Rectangle((x_missing - bar_width/2, 0), bar_width, placeholder_height, fill=False, hatch="//")
ax.add_patch(rect)
ax.text(x_missing, placeholder_height + 0.005, "NA", ha="center", va="bottom", fontsize=8)

# Aesthetics
ax.set_xticks(x)
ax.set_xticklabels(datasets, rotation=35, ha="right")
ax.set_ylabel("LOH Burden (%)", labelpad=10, fontweight='bold')
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.set_ylim(0, max(0.7, ymax + 0.08))
ax.grid(axis="y", linestyle=":", linewidth=0.6)

# Legend including "Not available"
#handles, labels = ax.get_legend_handles_labels()
#handles.append(Patch(fill=False, hatch="//", label="Not available"))
#labels.append("Not available")
#ax.legend(handles, labels, ncol=2, frameon=True)
ax.legend(frameon=False, ncol=1, prop={'weight': 'bold'})

#ax.set_title("LOH burden across eight datasets and four inference configurations")

fig.tight_layout()
plt.savefig("loh_burden_grouped_bar.png", dpi=300, bbox_inches="tight")
plt.show()