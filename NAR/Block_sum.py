import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os

# Samples (datasets)
samples = [
    "COLO829_ONT",
    "COLO829_NYGC",
    "H1437_UCSC",
    "H2009_UCSC",
    "HCC1395_HKU",
    "HCC1395_NYGC",
    "HCC1937_UCSC",
    "HCC1954_UCSC",
]

# Block sum (bp) for each phasing tool, in the same order as samples
hapcut2 = [
    2523260676,
    2187808909,
    2625509787,
    2677986524,
    2536199864,
    1972073127,
    2430726898,
    2711566888,
]

whatshap = [
    2518385299,
    2182781406,
    2624731180,
    2676922405,
    2532455755,
    1966798050,
    2425582425,
    2710095182,
]

longphase = [
    2540883450,
    2176223422,
    2654323753,
    2706782369,
    2569265774,
    1980748528,
    2454370721,
    2736740481,
]

longphase_to = [
    2835034563,
    2664715137,
    2864794192,
    2863280796,
    2847723783,
    2726465382,
    2830053891,
    2828110080,
]

# X positions for groups
x = np.arange(len(samples))
width = 0.2  # width of each bar

fig, ax = plt.subplots(figsize=(9, 5))

# Nature-style professional color palette (colorblind-friendly)
colors = ['#4477AA', '#EE6677', '#228833', '#CCBB44']  # Blue, Red, Green, Yellow

# Convert bp to Mbp (divide by 1,000,000)
hapcut2_mbp = [val / 1e6 for val in hapcut2]
whatshap_mbp = [val / 1e6 for val in whatshap]
longphase_mbp = [val / 1e6 for val in longphase]
longphase_to_mbp = [val / 1e6 for val in longphase_to]

# Plot grouped bars
ax.bar(x - 1.5 * width, hapcut2_mbp, width, label="HapCUT2", color=colors[0])
ax.bar(x - 0.5 * width, whatshap_mbp, width, label="WhatsHap", color=colors[1])
ax.bar(x + 0.5 * width, longphase_mbp, width, label="LongPhase", color=colors[2])
ax.bar(x + 1.5 * width, longphase_to_mbp, width, label="LongPhase-TO", color=colors[3])

# Labels and title
ax.set_ylabel("Block sum (Mbp)", labelpad=10, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(samples, rotation=35, ha="right")

# Set y-axis range from 1,000 Mbp to 3,500 Mbp
ax.set_ylim(bottom=1800, top=3000)

# Thousands separator for y-axis
def thousands_formatter(x, pos):
    return f"{int(x):,}" if x >= 1 else f"{x:g}"
ax.yaxis.set_major_formatter(FuncFormatter(thousands_formatter))

ax.grid(True, axis="y", linewidth=0.2)
ax.legend(frameon=False, ncol=4, prop={'weight': 'bold'}, loc='upper center')
plt.tight_layout()

# Save to current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
png_path = os.path.join(script_dir, "block_sum_grouped_bar.png")
plt.savefig(png_path, dpi=600, bbox_inches="tight")
plt.show()