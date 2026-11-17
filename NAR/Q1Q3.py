import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from io import StringIO

# Updated data as a tab-separated string
data = """sample	software	Expected purity	q1	q3
COLO829_ONT	LongPhase-TO (DeepSomatic-TO)	1	0.672815	0.999524
COLO829_ONT	LongPhase-TO (DeepSomatic-TO)	0.8	0.695251	0.888753
COLO829_ONT	LongPhase-TO (DeepSomatic-TO)	0.6	0.624895	0.80533
COLO829_ONT	LongPhase-TO (DeepSomatic-TO)	0.4	0.571312	0.731718
COLO829_ONT	LongPhase-TO (DeepSomatic-TO)	0.2	0.543464	0.678708
COLO829_NYGC	LongPhase-TO (DeepSomatic-TO)	1	0.707188	0.966355
COLO829_NYGC	LongPhase-TO (DeepSomatic-TO)	0.8	0.713726	0.882406
COLO829_NYGC	LongPhase-TO (DeepSomatic-TO)	0.6	0.624953	0.794514
COLO829_NYGC	LongPhase-TO (DeepSomatic-TO)	0.4	0.571208	0.718463
COLO829_NYGC	LongPhase-TO (DeepSomatic-TO)	0.2	0.540533	0.666569
H1437_UCSC	LongPhase-TO (DeepSomatic-TO)	1	0.666657	0.999371
H1437_UCSC	LongPhase-TO (DeepSomatic-TO)	0.8	0.666647	0.894623
H1437_UCSC	LongPhase-TO (DeepSomatic-TO)	0.6	0.611293	0.814938
H1437_UCSC	LongPhase-TO (DeepSomatic-TO)	0.4	0.571425	0.760228
H1437_UCSC	LongPhase-TO (DeepSomatic-TO)	0.2	0.551825	0.766483
H2009_UCSC	LongPhase-TO (DeepSomatic-TO)	1	0.666556	0.999337
H2009_UCSC	LongPhase-TO (DeepSomatic-TO)	0.8	0.666657	0.886013
H2009_UCSC	LongPhase-TO (DeepSomatic-TO)	0.6	0.605896	0.8
H2009_UCSC	LongPhase-TO (DeepSomatic-TO)	0.4	0.565217	0.733316
H2009_UCSC	LongPhase-TO (DeepSomatic-TO)	0.2	0.549997	0.713731
HCC1395_HKU	LongPhase-TO (DeepSomatic-TO)	1	0.720112	0.9995
HCC1395_HKU	LongPhase-TO (DeepSomatic-TO)	0.8	0.666664	0.922103
HCC1395_HKU	LongPhase-TO (DeepSomatic-TO)	0.6	0.621648	0.875221
HCC1395_HKU	LongPhase-TO (DeepSomatic-TO)	0.4	0.583121	0.856449
HCC1395_HKU	LongPhase-TO (DeepSomatic-TO)	0.2	0.562585	0.89988
HCC1395_NYGC	LongPhase-TO (DeepSomatic-TO)	1	0.640951	0.985027
HCC1395_NYGC	LongPhase-TO (DeepSomatic-TO)	0.8	0.647011	0.891888
HCC1395_NYGC	LongPhase-TO (DeepSomatic-TO)	0.6	0.600107	0.812461
HCC1395_NYGC	LongPhase-TO (DeepSomatic-TO)	0.4	0.567793	0.74997
HCC1395_NYGC	LongPhase-TO (DeepSomatic-TO)	0.2	0.552613	0.750295
HCC1937_UCSC	LongPhase-TO (DeepSomatic-TO)	1	0.706142	0.999995
HCC1937_UCSC	LongPhase-TO (DeepSomatic-TO)	0.8	0.666653	0.916158
HCC1937_UCSC	LongPhase-TO (DeepSomatic-TO)	0.6	0.615364	0.874896
HCC1937_UCSC	LongPhase-TO (DeepSomatic-TO)	0.4	0.583331	0.863578
HCC1937_UCSC	LongPhase-TO (DeepSomatic-TO)	0.2	0.578912	0.899944
HCC1954_UCSC	LongPhase-TO (DeepSomatic-TO)	1	0.639977	0.947228
HCC1954_UCSC	LongPhase-TO (DeepSomatic-TO)	0.8	0.625092	0.916727
HCC1954_UCSC	LongPhase-TO (DeepSomatic-TO)	0.6	0.605112	0.912904
HCC1954_UCSC	LongPhase-TO (DeepSomatic-TO)	0.4	0.590885	0.92853
HCC1954_UCSC	LongPhase-TO (DeepSomatic-TO)	0.2	0.599979	0.94319
"""

# Read into a DataFrame
df = pd.read_csv(StringIO(data), sep="\t")

# Make the purity column easier to use in code
df = df.rename(columns={"Expected purity": "Expected_purity"})

samples = df["sample"].unique()
n_samples = len(samples)

ncols = 4
nrows = int(np.ceil(n_samples / ncols))

fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), sharey=True)
axes = axes.flatten()

for idx, sample in enumerate(samples):
    ax = axes[idx]
    sub = df[df["sample"] == sample].sort_values("Expected_purity")

    x = np.arange(len(sub))
    width = 0.35

    ax.bar(x - width / 2, sub["q1"], width, label="Q1")
    ax.bar(x + width / 2, sub["q3"], width, label="Q3")

    ax.set_title(sample)
    ax.set_xlabel("Expected purity")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["Expected_purity"])
    if idx % ncols == 0:
        ax.set_ylabel("Value")

# Remove any unused subplots
for j in range(idx + 1, len(axes)):
    fig.delaxes(axes[j])

# One legend for all subplots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper right")

fig.tight_layout()
plt.show()
