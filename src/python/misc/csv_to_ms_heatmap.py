#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
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

    #! Collect meta data for controls and dorado
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

    # #! Add trivial values for Dorado
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #             "Tool": ["Dorado"] * 10,
    #             "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "min length", "mean length", "median length", "n50 length", "max length"],
    #             "Value": [total_reads, 0, 0, total_reads, 0, min_length, mean_length, median_length, n50_length, max_length],
    #             "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, df.loc[(df['Tool'] == 'Ctrl R.') & (df['Metric'] == 'min length'), 'Metric Score'].squeeze(), 1.0, 1.0, 1.0, 1.0],
    #         })
    #     ], ignore_index=True
    # )

    # Calculate metric score sum for all tools
    ams = df.groupby("Tool")["Metric Score"].sum().reset_index()
    ams["Dataset"] = f"{input_csv.split('/')[1]} {input_csv.split('/')[2]}" # Extract dataset name from file path
    return ams

def plot_heatmap(ams: pd.DataFrame, output_file: str):
    # Pivot the data for heatmap
    heatmap_data = ams.pivot(index="Tool", columns="Dataset", values="Metric Score")

    # Sort tools by their mean metric score (descending order)
    tool_order = heatmap_data.mean(axis=1).sort_values(ascending=False).index
    heatmap_data = heatmap_data.loc[tool_order]

    # Plot the heatmap
    plt.figure(figsize=(9, 7))  # Adjust figure size for better readability
    sns.heatmap(
        heatmap_data,
        annot=True,  # Display values in cells
        fmt=".2f",  # Format values to 2 decimal places
        cmap="coolwarm",  # Use a visually appealing color palette
        cbar_kws={'label': 'Metric Score'},  # Add a label to the color bar
        linewidths=0.5,  # Add grey lines between cells
        linecolor="grey",  # Set the line color to grey
        annot_kws={"fontsize": 9},  # Adjust font size for annotations
    )
    plt.title("Metric Score Heatmap", fontsize=16, fontweight="bold")  # Add a bold title
    plt.ylabel("Tool", fontsize=12)  # Adjust y-axis label font size
    plt.xlabel("Dataset", fontsize=12)  # Adjust x-axis label font size
    plt.xticks(rotation=45, ha="right", fontsize=10)  # Rotate x-axis labels for better readability
    plt.yticks(rotation=0, fontsize=10)  # Adjust y-axis label font size
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