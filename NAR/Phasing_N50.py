# Grouped bar chart for Block N50 converted to Kbp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os

# Recreate the dataframe
data = [
    ("COLO829_ONT","HapCUT2",0.531188723,3569,4645832,2523260676),
    ("COLO829_ONT","WhatsHap",0.530234997,3686,4635320,2518385299),
    ("COLO829_ONT","LongPhase",0.498077451,3477,4887010,2540883450),
    ("COLO829_ONT","LongPhase-TO",0.594247857,1110,9852371,2835034563),
    ("COLO829_NYGC","HapCUT2",0.557132505,15016,584415,2187808909),
    ("COLO829_NYGC","WhatsHap",0.556136688,15261,575997,2182781406),
    ("COLO829_NYGC","LongPhase",0.521019637,15017,570504,2176223422),
    ("COLO829_NYGC","LongPhase-TO",0.606857364,6576,1403687,2664715137),
    ("H1437_UCSC","HapCUT2",0.539916273,2376,7048644,2625509787),
    ("H1437_UCSC","WhatsHap",0.538696787,2398,7048644,2624731180),
    ("H1437_UCSC","LongPhase",0.501837431,2287,7729758,2654323753),
    ("H1437_UCSC","LongPhase-TO",0.697758143,529,19323318,2864794192),
    ("H2009_UCSC","HapCUT2",0.557470906,2252,5825292,2677986524),
    ("H2009_UCSC","WhatsHap",0.556084355,2292,5723430,2676922405),
    ("H2009_UCSC","LongPhase",0.516697084,2204,7217802,2706782369),
    ("H2009_UCSC","LongPhase-TO",0.682400933,638,23169348,2863280796),
    ("HCC1395_HKU","HapCUT2",0.601509894,4553,1678506,2536199864),
    ("HCC1395_HKU","WhatsHap",0.600173280,4662,1667777,2532455755),
    ("HCC1395_HKU","LongPhase",0.560688427,4581,1754642,2569265774),
    ("HCC1395_HKU","LongPhase-TO",0.575700475,893,16768123,2847723783),
    ("HCC1395_NYGC","HapCUT2",0.498829689,21138,275106,1972073127),
    ("HCC1395_NYGC","WhatsHap",0.498078089,21510,270769,1966798050),
    ("HCC1395_NYGC","LongPhase",0.467598877,21057,279702,1980748528),
    ("HCC1395_NYGC","LongPhase-TO",0.693399557,5227,13441849,2726465382),
    ("HCC1937_UCSC","HapCUT2",0.401595056,4537,1447908,2430726898),
    ("HCC1937_UCSC","WhatsHap",0.400806954,4650,1420186,2425582425),
    ("HCC1937_UCSC","LongPhase",0.375063484,4517,1513047,2454370721),
    ("HCC1937_UCSC","LongPhase-TO",0.713396249,1135,25207324,2830053891),
    ("HCC1954_UCSC","HapCUT2",0.651938431,1969,5927639,2711566888),
    ("HCC1954_UCSC","WhatsHap",0.650486924,1999,5908073,2710095182),
    ("HCC1954_UCSC","LongPhase",0.610394471,1923,5912510,2736740481),
    ("HCC1954_UCSC","LongPhase-TO",0.659481598,1115,6477286,2828110080),
]

df = pd.DataFrame(data, columns=["sample","phasing_tool","phased_ratio","block_num","block_n50","block_sum"])
df["block_n50_kbp"] = df["block_n50"] / 1000.0

dataset_order = [
    "COLO829_ONT",
    "COLO829_NYGC",
    "H1437_UCSC",
    "H2009_UCSC",
    "HCC1395_HKU",
    "HCC1395_NYGC",
    "HCC1937_UCSC",
    "HCC1954_UCSC",
]

tools = ["HapCUT2", "WhatsHap", "LongPhase", "LongPhase-TO"]

pivot = df.pivot(index="sample", columns="phasing_tool", values="block_n50_kbp").loc[dataset_order, tools]

# Create grouped bar chart
plt.figure(figsize=(9,5))
ax = plt.gca()

n_datasets = len(pivot.index)
n_tools = len(tools)
x = range(n_datasets)
bar_width = 0.18
offsets = [(-1.5 + i) * bar_width for i in range(n_tools)]

# Nature-style professional color palette (colorblind-friendly)
colors = ['#4477AA', '#EE6677', '#228833', '#CCBB44']  # Blue, Red, Green, Yellow

for i, tool in enumerate(tools):
    ax.bar([xi + offsets[i] for xi in x], pivot[tool].values, width=bar_width, label=tool, color=colors[i])

# X/Y labels and ticks
#ax.set_xlabel("Dataset", labelpad=10)
ax.set_ylabel("Block N50 (Kbp)", labelpad=10, fontweight='bold')
ax.set_xticks(list(x))
ax.set_xticklabels(pivot.index, rotation=35, ha="right")

# Thousands separator for Kbp axis
def thousands_formatter(x, pos):
    return f"{int(x):,}" if x >= 1 else f"{x:g}"
ax.yaxis.set_major_formatter(FuncFormatter(thousands_formatter))

ax.grid(True, axis="y", linewidth=0.2)
ax.legend(frameon=False, prop={'weight': 'bold'})
plt.tight_layout()

# Save to current directory
script_dir = os.path.dirname(os.path.abspath(__file__))
png_path = os.path.join(script_dir, "phasing_n50_grouped_bar_kbp.png")
#pdf_path = os.path.join(script_dir, "phasing_n50_grouped_bar_kbp.pdf")
plt.savefig(png_path, dpi=600, bbox_inches="tight")
#plt.savefig(pdf_path, bbox_inches="tight")
plt.show()