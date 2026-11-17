# Create a grouped bar chart (Nature-style) comparing Block N50 (Kbp)
# for four methods across samples, using only matplotlib (no seaborn).
# Files will be saved to /mnt/data for download.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle

# Raw data from the user's table: Sample, Method, LOH N50 Size (bp)
rows = [
    ("COLO829_ONT", "ASCAT (T+N)",                     12864314),
    ("COLO829_ONT", "ASCAT (TO)",                      20362317),
    ("COLO829_ONT", "LongPhase-TO (ClairS-TO)",        50760310),
    ("COLO829_ONT", "LongPhase-TO (DeepSomatic-TO)",   50760310),
    
    ("COLO829_NYGC", "ASCAT (T+N)",                    80866668),
    ("COLO829_NYGC", "ASCAT (TO)",                     38166624),
    ("COLO829_NYGC", "LongPhase-TO (ClairS-TO)",       35105161),
    ("COLO829_NYGC", "LongPhase-TO (DeepSomatic-TO)",  45643569),
    
    ("H1437_UCSC", "ASCAT (T+N)",                      3829891),
    ("H1437_UCSC", "ASCAT (TO)",                       34097583),
    ("H1437_UCSC", "LongPhase-TO (ClairS-TO)",         43336933),
    ("H1437_UCSC", "LongPhase-TO (DeepSomatic-TO)",    34605403),
    
    ("H2009_UCSC", "ASCAT (T+N)",                      4110937),
    ("H2009_UCSC", "ASCAT (TO)",                       9209663),
    ("H2009_UCSC", "LongPhase-TO (ClairS-TO)",         31122631),
    ("H2009_UCSC", "LongPhase-TO (DeepSomatic-TO)",    34937899),
    
    ("HCC1395_HKU", "ASCAT (T+N)",                     2803775),
    ("HCC1395_HKU", "ASCAT (TO)",                      6023307),
    ("HCC1395_HKU", "LongPhase-TO (ClairS-TO)",        21772604),
    ("HCC1395_HKU", "LongPhase-TO (DeepSomatic-TO)",   21597183),
    
    ("HCC1395_NYGC", "ASCAT (T+N)",                    9642863),
    ("HCC1395_NYGC", "ASCAT (TO)",                     16795656),
    ("HCC1395_NYGC", "LongPhase-TO (ClairS-TO)",       26501707),
    ("HCC1395_NYGC", "LongPhase-TO (DeepSomatic-TO)",   29068243),
    
    ("HCC1937_UCSC", "ASCAT (T+N)",                    2214193),
    ("HCC1937_UCSC", "ASCAT (TO)",                     20879693),
    ("HCC1937_UCSC", "LongPhase-TO (ClairS-TO)",       40921298),
    ("HCC1937_UCSC", "LongPhase-TO (DeepSomatic-TO)",   43800261),
    
    # ASCAT (T+N) missing for HCC1954_UCSC; use NaN to indicate missing
    ("HCC1954_UCSC", "ASCAT (T+N)",                    np.nan),
    ("HCC1954_UCSC", "ASCAT (TO)",                     2812944),
    ("HCC1954_UCSC", "LongPhase-TO (ClairS-TO)",       16263599),
    ("HCC1954_UCSC", "LongPhase-TO (DeepSomatic-TO)",  15126962),

]

df = pd.DataFrame(rows, columns=["Sample", "Method", "N50_bp"])

# Keep only the four requested methods (already enforced above)
methods = ["ASCAT (T+N)",
           "ASCAT (TO)",
           "LongPhase-TO (ClairS-TO)",
           "LongPhase-TO (DeepSomatic-TO)"
           ]

# Convert to Kbp
df["N50_Kbp"] = df["N50_bp"] / 1000.0

# Preserve the desired sample order as they were listed
sample_order = [
    "COLO829_ONT", "COLO829_NYGC",
    "H1437_UCSC", "H2009_UCSC",
    "HCC1395_HKU", "HCC1395_NYGC",
    "HCC1937_UCSC", "HCC1954_UCSC"
]
df["Sample"] = pd.Categorical(df["Sample"], categories=sample_order, ordered=True)

# Pivot to wide for grouped bars
pivot = df.pivot_table(index="Sample", columns="Method", values="N50_Kbp", aggfunc="first").reindex(sample_order)


# Plot
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(pivot.index))
n_methods = len(methods)
bar_width = 0.18
offsets = np.linspace(-1.5*bar_width, 1.5*bar_width, n_methods)

# Nature-style professional color palette (colorblind-friendly)
colors = ['#332288', '#88CCEE', '#CC6677', '#DDCC77']  # Purple, Cyan, Pink, Yellow

ymax = 0.0
for i, method in enumerate(methods):
    y = pivot[method].values if method in pivot.columns else np.full(len(pivot.index), np.nan)
    # Plot real bars where values exist
    mask = ~np.isnan(y)
    ax.bar(x[mask] + offsets[i], y[mask], width=bar_width, label=method, color=colors[i])
    ymax = max(ymax, np.nanmax(y))

# Placeholder for missing ASCAT (T+N) in HCC1954_UCSC
missing_idx = np.where(pivot.index == "HCC1954_UCSC")[0][0]
tn_idx = methods.index("ASCAT (T+N)")
x_missing = x[missing_idx] + offsets[tn_idx]

placeholder_height = ymax * 0.02  # 2% of max height
rect = Rectangle((x_missing - bar_width/2, 0), bar_width, placeholder_height, fill=False, hatch="//")
ax.add_patch(rect)
ax.text(x_missing, placeholder_height + ymax * 0.005, "NA", ha="center", va="bottom", fontsize=8)

# Aesthetics (no explicit colors/styles)
ax.set_ylabel("LOH N50 (Kbp)", labelpad=10, fontweight='bold')
#ax.set_title("Block N50 by Method across Samples")
ax.set_xticks(x)
ax.set_xticklabels(pivot.index, rotation=35, ha="right")
ax.grid(True, axis="y", linewidth=0.2)
ax.legend(frameon=False, ncol=2, prop={'weight': 'bold'})

plt.tight_layout()

# Save outputs
png_path = "loh_n50_grouped.png"

plt.savefig(png_path, dpi=300, bbox_inches="tight")


# Display the figure inline
plt.show()
