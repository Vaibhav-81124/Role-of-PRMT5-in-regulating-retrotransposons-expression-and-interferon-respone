import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Assign column names manually
colnames = ["Name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]

# Read the file with no header, force column names
df = pd.read_csv("GSK_vs_DMSO_sigdiff_gene_TE.txt", sep="\t", header=None, names=colnames)

# Convert numeric columns
for col in ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")


# --- Improved TE Class Distribution Plot (Bar Chart) ---
import seaborn as sns

# Keep only TE rows
te_df = df[df["Name"].str.contains(":", na=False)].copy()
te_df["Class"] = te_df["Name"].apply(lambda x: x.split(":")[1] if ":" in x else "Other")

# Count TE classes
class_counts = te_df["Class"].value_counts()

# classes (bar chart)
plt.figure(figsize=(8,5))
sns.barplot(x=class_counts.values, y=class_counts.index, color="skyblue")
plt.xlabel("Count")
plt.ylabel("TE Class")
plt.title("Distribution of TE Classes (Significant DE)")
plt.tight_layout()
plt.savefig("te_class_bar.png", dpi=900)
plt.close()

#Group small classes into "Other"
threshold = 5
class_counts_grouped = class_counts.copy()
small_classes = class_counts_grouped[class_counts_grouped < threshold]
class_counts_grouped = class_counts_grouped[class_counts_grouped >= threshold]
class_counts_grouped["Other"] = small_classes.sum()

plt.figure(figsize=(6,6))
plt.pie(class_counts_grouped, labels=class_counts_grouped.index, autopct="%1.1f%%", startangle=90)
plt.title("Distribution of TE Classes (Grouped)")
plt.savefig("te_class_pie_grouped.png", dpi=900)
plt.close()

print(" Bar chart saved as te_class_bar.png")
print(" Grouped pie chart saved as te_class_pie_grouped.png")

