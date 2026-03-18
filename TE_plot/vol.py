import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from adjustText import adjust_text
# Load file
df = pd.read_csv("GSK_vs_DMSO_sigdiff_TEonly.txt", sep="\t", index_col=False)

# Assign column names if missing
if "Name" not in df.columns:
    df.columns = ["Name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# Convert numeric columns
for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# Split genes vs TEs
df["Type"] = df["Name"].apply(lambda x: "TE" if ":" in str(x) else "Gene")

# Volcano plot
plt.figure(figsize=(10,7))

# Plot genes (blue)
#plt.scatter(
 #   df[df["Type"]=="Gene"]["log2FoldChange"],
  #  -np.log10(df[df["Type"]=="Gene"]["pvalue"]),
   # alpha=0.5, s=10, c="blue", label="Genes")

# Plot TEs (red)
plt.scatter(
    df[df["Type"]=="TE"]["log2FoldChange"],
    -np.log10(df[df["Type"]=="TE"]["pvalue"]),
    alpha=0.8, s=15, c="red", label="TEs"
)

# Threshold lines
plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
plt.axvline(-1, color="gray", linestyle="--", linewidth=0.8)
plt.axvline(1, color="gray", linestyle="--", linewidth=0.8)

# Label TEs with strong fold change
texts=[]
te_high_fc = df[(df["Type"]=="TE") & (abs(df["log2FoldChange"]) >= 1.0)]
for _, row in te_high_fc.iterrows():
    plt.text(
        row["log2FoldChange"],
        -np.log10(row["pvalue"]),
        row["Name"],
        fontsize=6,
	fontweight='bold',
        color="black",
        ha="right" if row["log2FoldChange"] < 0 else "left"
    )

adjust_text(
    texts,
    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
    only_move={'points':'y', 'texts':'y'}  # keeps x fixed, moves y
)

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot: GSK vs DMSO")
plt.legend()
plt.tight_layout()
plt.savefig("volcano_te_labeled.png", dpi=900)
plt.close()

print("volcano plot saved")

