import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# def readScore(file : str) -> pd.DataFrame:
#   data = {}
#   with open(file, 'r') as f:
#     for line in f:
#       line = line.strip()
#       if line.startswith("Tool"):
#         tool = line.split(": ")[1]
#         data[tool] = {}
#       elif line.startswith("Median"):
#         median = float(line.split(": ")[1])
#         data[tool]["Median"] = median
#       elif line.startswith("Mean"):
#         mean = float(line.split(": ")[1])
#         data[tool]["Mean"] = mean
#   return pd.DataFrame(data)

datasetsPaths = {
  "RNA002" : {
    "H. Sapiens" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/h_sapiens/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/h_sapiens/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/h_sapiens/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/h_sapiens/comparison_w0_readQuality.csv",
    },
    "E. Coli" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/e_coli/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/e_coli/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/e_coli/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/e_coli/comparison_w0_readQuality.csv",
    },
    "SARS-CoV-2" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/sarscov2/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/sarscov2/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/sarscov2/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/sarscov2/comparison_w0_readQuality.csv",
    },
    "IVT" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/ivt/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/ivt/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/ivt/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/ivt/comparison_w0_readQuality.csv",
    },
  },
  "RNA004" : {
    "H. Sapiens" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/h_sapiens/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/h_sapiens/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/h_sapiens/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/h_sapiens/comparison_w0_readQuality.csv",
    },
    "S. Cerevisiae" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/s_cerevisiae/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/s_cerevisiae/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/s_cerevisiae/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/s_cerevisiae/comparison_w0_readQuality.csv",
    },
    "CEVd" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/cevd/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/cevd/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/cevd/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/cevd/comparison_w0_readQuality.csv",
    },
    "IVT" : {
      "Score" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/ivt/comparison_w0_score.txt",
      "Reads Segmentation Ratio" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/ivt/comparison_w0_segmentedReadsRatio.csv",
      "Read Lengths" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/ivt/comparison_w0_readLengths.csv",
      "Read Quality" : "/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/ivt/comparison_w0_readQuality.csv",
    },
  }
}

scoreData = pd.DataFrame()
readsSegmentationRatioData = pd.DataFrame()
readLengthsData = pd.DataFrame()
readQualityData = pd.DataFrame()

for pore in datasetsPaths:
  for dataset in datasetsPaths[pore]:

    scorePath = datasetsPaths[pore][dataset]["Score"]
    readsSegmentationRatioPath = datasetsPaths[pore][dataset]["Reads Segmentation Ratio"]
    readLengthsPath = datasetsPaths[pore][dataset]["Read Lengths"]
    readQualityPath = datasetsPaths[pore][dataset]["Read Quality"]

    score = pd.read_csv(scorePath, sep='\t')
    readsSegmentationRatio = pd.read_csv(readsSegmentationRatioPath, sep=",")
    readLengths = pd.read_csv(readLengthsPath, sep="\t")
    readQuality = pd.read_csv(readQualityPath, sep="\t")

    score["Pore"] = pore
    score["Dataset"] = dataset
    readsSegmentationRatio["Pore"] = pore
    readsSegmentationRatio["Dataset"] = dataset
    readLengths["Pore"] = pore
    readLengths["Dataset"] = dataset
    readQuality["Pore"] = pore
    readQuality["Dataset"] = dataset

    scoreData = pd.concat([scoreData, score], ignore_index=True)
    readsSegmentationRatioData = pd.concat([readsSegmentationRatioData, readsSegmentationRatio], ignore_index=True)
    readLengthsData = pd.concat([readLengthsData, readLengths], ignore_index=True)
    readQualityData = pd.concat([readQualityData, readQuality], ignore_index=True)


### PLotting
sns.set_theme()

dataset_order = ["H. Sapiens", "E. Coli", "SARS-CoV-2", "IVT"]

# rna002 tools vs datasets median score
score_pivot = scoreData[(scoreData["Pore"] == "RNA002") & (scoreData["Segment Quality"] == "Contrast")].pivot(index="Tool", columns="Dataset", values="Median Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma")
plt.title("RNA002 Median Contrast Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMedian.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMedian.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMedian.pdf", dpi=300)
plt.close()

score_pivot = scoreData[(scoreData["Pore"] == "RNA002") & (scoreData["Segment Quality"] == "Homogeneity")].pivot(index="Tool", columns="Dataset", values="Median Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=True).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma_r")
plt.title("RNA002 Median Homogeneity Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMedian.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMedian.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMedian.pdf", dpi=300)
plt.close()

score_pivot = scoreData[(scoreData["Pore"] == "RNA002") & (scoreData["Segment Quality"] == "Contrast")].pivot(index="Tool", columns="Dataset", values="Mean Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma")
plt.title("RNA002 Mean Contrast Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMean.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMean.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsContrastMean.pdf", dpi=300)
plt.close()

score_pivot = scoreData[(scoreData["Pore"] == "RNA002") & (scoreData["Segment Quality"] == "Homogeneity")].pivot(index="Tool", columns="Dataset", values="Mean Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=True).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma_r")
plt.title("RNA002 Mean Homogeneity Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMean.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMean.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/toolsDatasetsHomogeneityMean.pdf", dpi=300)
plt.close()

# segmented reads ratio RNA002

score_pivot = readsSegmentationRatioData[readsSegmentationRatioData["Pore"] == "RNA002"].pivot(index="Tool", columns="Dataset", values="Segmented Reads Ratio")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
# Create formatted annotation values with percentage signs
annot_values = score_pivot * 100  # Convert to percentage
annot_text = np.vectorize(lambda x: f"{x:.0f}%")(annot_values)  # Convert to string with %
sns.heatmap(annot_values, annot=annot_text, fmt="", cmap="magma")
plt.title("RNA002 Segmented Reads Ratio")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/segmentedReadsRatio.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/segmentedReadsRatio.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/segmentedReadsRatio.pdf", dpi=300)
plt.close()

########################################


dataset_order = ["H. Sapiens", "S. Cerevisiae", "CEVd", "IVT"]

# rna004 tools vs datasets median score
score_pivot = scoreData[(scoreData["Pore"] == "RNA004") & (scoreData["Segment Quality"] == "Contrast")].pivot(index="Tool", columns="Dataset", values="Median Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma")
plt.title("RNA004 Median Contrast Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMedian.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMedian.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMedian.pdf", dpi=300)
plt.close()


score_pivot = scoreData[(scoreData["Pore"] == "RNA004") & (scoreData["Segment Quality"] == "Homogeneity")].pivot(index="Tool", columns="Dataset", values="Median Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=True).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma_r")
plt.title("RNA004 Median Homogeneity Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMedian.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMedian.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMedian.pdf", dpi=300)
plt.close()

score_pivot = scoreData[(scoreData["Pore"] == "RNA004") & (scoreData["Segment Quality"] == "Contrast")].pivot(index="Tool", columns="Dataset", values="Mean Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma")
plt.title("RNA004 Mean Contrast Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMean.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMean.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsContrastMean.pdf", dpi=300)
plt.close()


score_pivot = scoreData[(scoreData["Pore"] == "RNA004") & (scoreData["Segment Quality"] == "Homogeneity")].pivot(index="Tool", columns="Dataset", values="Mean Score")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=True).index]
score_pivot = score_pivot[dataset_order]
sns.heatmap(score_pivot, annot=True, fmt=".1f", cmap="magma_r")
plt.title("RNA004 Mean Homogeneity Score")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMean.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMean.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/toolsDatasetsHomogeneityMean.pdf", dpi=300)
plt.close()

# segmented reads ratio RNA004

score_pivot = readsSegmentationRatioData[readsSegmentationRatioData["Pore"] == "RNA004"].pivot(index="Tool", columns="Dataset", values="Segmented Reads Ratio")
# Sort tools based on the highest average score across datasets (descending)
score_pivot = score_pivot.loc[score_pivot.mean(axis=1).sort_values(ascending=False).index]
score_pivot = score_pivot[dataset_order]
# Create formatted annotation values with percentage signs
annot_values = score_pivot * 100  # Convert to percentage
annot_text = np.vectorize(lambda x: f"{x:.0f}%")(annot_values)  # Convert to string with %
sns.heatmap(annot_values, annot=annot_text, fmt="", cmap="magma")
plt.title("RNA004 Segmented Reads Ratio")
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/segmentedReadsRatio.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/segmentedReadsRatio.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/segmentedReadsRatio.pdf", dpi=300)
plt.close()

########################################

# datasets Read Lengths

ax = sns.barplot(data=readLengthsData[readLengthsData["Pore"] == "RNA002"], x="Dataset", y="Value", hue="Metric")
for container in ax.containers:
    ax.bar_label(container, label_type="center", rotation=90, fontsize=10, color='white', fmt="%.0f")
plt.title("RNA002 Read Length Statistics")
plt.yscale("log")
plt.legend(title="Metric", bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
# plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit legend
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readLengths.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readLengths.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readLengths.pdf", dpi=300)
plt.close()

ax = sns.barplot(data=readLengthsData[readLengthsData["Pore"] == "RNA004"], x="Dataset", y="Value", hue="Metric")
for container in ax.containers:
    ax.bar_label(container, label_type="center", rotation=90, fontsize=10, color='white', fmt="%.0f")
plt.title("RNA004 Read Length Statistics")
plt.yscale("log")
plt.legend(title="Metric", bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
# plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit legend
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readLengths.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readLengths.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readLengths.pdf", dpi=300)
plt.close()

ax = sns.barplot(data=readQualityData[readQualityData["Pore"] == "RNA002"], x="Dataset", y="Value", hue="Metric")
for container in ax.containers:
    ax.bar_label(container, label_type="center", rotation=90, fontsize=10, color='white', fmt="%.1f")
plt.title("RNA002 Read Quality Statistics")
plt.yscale("log")
plt.legend(title="Metric", bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
# plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit legend
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readQuality.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readQuality.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna002/readQuality.pdf", dpi=300)
plt.close()

ax = sns.barplot(data=readQualityData[readQualityData["Pore"] == "RNA004"], x="Dataset", y="Value", hue="Metric")
for container in ax.containers:
    ax.bar_label(container, label_type="center", rotation=90, fontsize=10, color='white', fmt="%.1f")
plt.title("RNA004 Read Quality Statistics")
plt.yscale("log")
plt.legend(title="Metric", bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
# plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to fit legend
plt.tight_layout()
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readQuality.png", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readQuality.svg", dpi=300)
plt.savefig("/data/fass5/projects/js_dynamont/benchmark/comparison/rna004/readQuality.pdf", dpi=300)
plt.close()