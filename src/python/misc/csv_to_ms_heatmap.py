#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_csvs", nargs="+", type=str, help="List of input CSV files")
    parser.add_argument("output_heatmap", type=str, help="Path to the output heatmap image file")
    return parser.parse_args()

def process_csv(input_csv: str) -> pd.DataFrame:
    # Load the data
    df = pd.read_csv(input_csv, sep="\t")

    # Cleanup DataFrame
    df.loc[df['Tool'] == 'Control Random', 'Tool'] = 'Ctrl R.'
    df.loc[df['Tool'] == 'Control Uniform', 'Tool'] = 'Ctrl U.'
    df.loc[df['Tool'] == 'f5c Eventalign', 'Tool'] = 'f5c E.'
    df.loc[df['Tool'] == 'f5c Resquiggle', 'Tool'] = 'f5c R.'
    df['Metric'] = df['Metric'].apply(lambda x: x.replace("_", " ").lower())
    df.loc[df['Metric'] == 'total', 'Metric'] = 'total reads'
    df.loc[df['Metric'] == 'truncated', 'Metric'] = 'truncated reads'
    df.loc[df['Metric'] == 'identical', 'Metric'] = 'identical reads'
    df.loc[df['Metric'] == 'present', 'Metric'] = 'segmented reads'
    df.loc[df['Metric'] == 'missing', 'Metric'] = 'missing reads'

    # Collect meta data for controls and dorado
    # total_reads = df.loc[df['Metric'] == 'total reads', 'Value'].values[0]
    # min_length = df.loc[df['Metric'] == 'min length', 'Value'].values[0]
    # max_length = df.loc[df['Metric'] == 'max length', 'Value'].values[0]
    # mean_length = df.loc[df['Metric'] == 'mean length', 'Value'].values[0]
    # median_length = df.loc[df['Metric'] == 'median length', 'Value'].values[0]
    # n50_length = df.loc[df['Metric'] == 'n50 length', 'Value'].values[0]

    # #! Add trivial values for controls
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #             "Tool": ["Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl U.", "Ctrl U.", "Ctrl U.", "Ctrl U.", "Ctrl U."],
    #             "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed"],
    #             "Value": [total_reads, 0, 0, total_reads, 0, total_reads, 0, 0, total_reads, 0],
    #             "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #         })
    #     ], ignore_index=True
    # )

    #! Add trivial values for Dorado
    # print(df)
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #             "Tool": ["Dorado"] * 10,
    #             "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "min length", "mean length", "median length", "n50 length", "max length"],
    #             "Value": [total_reads, 0, 0, total_reads, 0, min_length, mean_length, median_length, n50_length, max_length],
    #             "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #         })
    #     ], ignore_index=True
    # )
    # print(df)

    #! include specific metrics
    df_for_agg = df[df["Metric"].isin([
        "median delta",
        "mad delta",
        "homogeneity",
        "segmented reads",
        # "missing reads",
        "truncated reads",
        # "identical reads",
        # "nt changed",
        "min length",
        # "mean length",
        # "median length",
        "n50 length",
        "max length",
        "flye total length",
        "flye n50",
        "flye mean coverage",
        "svim structural variants",
    ])]

    # Calculate metric score sum for all tools
    # ams = df.groupby("Tool")["Metric Score"].sum().reset_index()
    ams = df_for_agg.groupby('Tool')["Metric Score"].sum().reset_index()
    ams["Dataset"] = f"{input_csv.split('/')[1]} {input_csv.split('/')[2]}" # Extract dataset name from file path
    return ams

def plot_heatmap(ams: pd.DataFrame, output_file: str):
    # Pivot the data for heatmap
    heatmap_data = ams.pivot(index="Tool", columns="Dataset", values="Metric Score")

    # print(heatmap_data.columns)

    # set column order
    column_order = [
        "rna002 h_sapiens",
        "rna002 e_coli",
        "rna002 sarscov2",
        "rna002 ivt",
        "rna004 h_sapiens",
        "rna004 s_cerevisiae",
        "rna004 cevd",
        "rna004 ivt",
        "dna_r10.4.1_5kHz h_sapiens",
        "dna_r10.4.1_5kHz zymo_hmw",
        "dna_r10.4.1_5kHz s_aureus",
        "dna_r10.4.1_5kHz p_anserina",
    ]
    heatmap_data = heatmap_data[column_order]

    # Rename columns
    column_rename_map = {
        "rna002 h_sapiens": r"$H.\ sapiens$",
        "rna002 e_coli": r"$E.\ coli$",
        "rna002 sarscov2": r"SARS-CoV-2",
        "rna002 ivt": r"IVT",
        "rna004 h_sapiens": r"$H.\ sapiens$",
        "rna004 s_cerevisiae": r"$S.\ cerevisiae$",
        "rna004 cevd": r"CEVD",
        "rna004 ivt": r"IVT",
        "dna_r10.4.1_5kHz h_sapiens": r"$H.\ sapiens$",
        "dna_r10.4.1_5kHz zymo_hmw": r"Zymo HMW",
        "dna_r10.4.1_5kHz s_aureus": r"$S.\ Aureus$",
        "dna_r10.4.1_5kHz p_anserina": r"$P.\ Anserina$",
    }
    heatmap_data = heatmap_data.rename(columns=column_rename_map)

    # Add average column
    heatmap_data["tool average"] = heatmap_data.mean(axis=1)

    # Add average row (for each dataset/column, average across all tools)
    avg_row = heatmap_data.mean(axis=0)
    avg_row.name = "dataset average"
    heatmap_data = pd.concat([heatmap_data, avg_row.to_frame().T])

    # Sort tools by their mean metric score (descending order)
    print(heatmap_data)
    tool_order = heatmap_data.mean(axis=1).sort_values(ascending=False).index.tolist()
    # ensure that dorado is the top row
    if "Dorado" in tool_order:
        tool_order.remove("Dorado")
        tool_order = ["Dorado"] + [t for t in tool_order if t != "dataset average"] + (["dataset average"] if "dataset average" in tool_order else [])
    heatmap_data = heatmap_data.loc[tool_order]

    # Plot the heatmap
    plt.figure(figsize=(9, 6))  # Adjust figure size for better readability
    ax = sns.heatmap(
        heatmap_data,
        annot=True,  # Display values in cells
        fmt=".2f",  # Format values to 2 decimal places
        cmap="coolwarm",  # Use a visually appealing color palette
        cbar_kws={'label': 'Score', 'shrink': 0.8},  # Adjust color bar size
        linewidths=0.5,  # Add grey lines between cells
        linecolor="grey",  # Set the line color to grey
        annot_kws={"fontsize": 9},  # Adjust font size for annotations
        square=True,  # Make cells square
    )

    # Add superlabels above the dataset labels
    superlabels = [
        "RNA002", "RNA002", "RNA002", "RNA002",
        "RNA004", "RNA004", "RNA004", "RNA004",
        "DNA R10.4.1 5kHz", "DNA R10.4.1 5kHz", "DNA R10.4.1 5kHz", "DNA R10.4.1 5kHz"
        ""
    ]
    dataset_labels = [
        r"$H.\ sapiens$", r"$E.\ coli$", "SARS-CoV-2", "IVT",
        r"$H.\ sapiens$", r"$S.\ cerevisiae$", "CEVd", "IVT",
        r"$H.\ sapiens$", "Zymo HMW", r"$S.\ Aureus$", r"$P.\ Anserina$", "tool average"
    ]

    # Set the dataset labels
    ax.set_xticks([i + 0.5 for i in range(len(dataset_labels))])  # Center labels
    ax.set_xticklabels(dataset_labels, rotation=45, ha="right", fontsize=10)

    # Add superlabels
    for i, label in enumerate(superlabels):
        if i == 0 or superlabels[i] != superlabels[i - 1]:  # Only add label once per group
            start = i
            end = i + superlabels.count(superlabels[i]) - 1
            ax.text(
                (start + end) / 2 + 0.5, 1.25 * len(tool_order),  # Center above group
                label,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                transform=ax.transData
            )

    # Adjust layout to fit the labels
    plt.subplots_adjust(bottom=0.2, top=0.85)

    # plt.title("Aggregated Metric Score", fontsize=14)  # Add a bold title
    # plt.ylabel("Tool", fontsize=12)  # Adjust y-axis label font size
    # plt.xlabel("Dataset", fontsize=12, labelpad=25)  # Adjust x-axis label font size
    # plt.xticks(rotation=45, ha="right", fontsize=10)  # Rotate x-axis labels for better readability
    plt.xticks(rotation=25, ha="center", fontsize=9)  # Rotate x-axis labels for better readability
    plt.yticks(rotation=0, fontsize=9)  # Adjust y-axis label font size
    plt.tight_layout()  # Ensure everything fits within the figure

    # Save the heatmap
    plt.savefig(output_file, dpi=300)  # Save with high resolution for better quality
    plt.savefig(output_file.split('.')[0] + '.pdf', dpi=300)  # Save with high resolution for better quality
    plt.close()

if __name__ == "__main__":
    args = parse()

    # Process all input CSVs
    all_ams = pd.concat([process_csv(csv) for csv in args.input_csvs], ignore_index=True)

    # Plot the heatmap
    plot_heatmap(all_ams, args.output_heatmap)