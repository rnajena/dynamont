#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

import pandas as pd
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_csv", type=str, help="Path to the input CSV file")
    parser.add_argument("output_tex", type=str, help="Path to the output LaTeX file")
    return parser.parse_args()    

def convert_metrics_to_latex(input_csv: str, output_tex: str = None):
    # Load the data
    df = pd.read_csv(input_csv, sep="\t")

    # cleanup DataFrame
    # df.loc[df['Tool'] == 'Control Random', 'Tool'] = 'Ctrl R.'
    # df.loc[df['Tool'] == 'Control Uniform', 'Tool'] = 'Ctrl U.'
    df.loc[df['Tool'] == 'f5c Eventalign', 'Tool'] = 'f5c E.'
    df.loc[df['Tool'] == 'f5c Resquiggle', 'Tool'] = 'f5c R.'
    df['Metric'] = df['Metric'].apply(lambda x: x.replace("_", " ").lower())
    df.loc[df['Metric'] == 'total', 'Metric'] = 'total reads'
    df.loc[df['Metric'] == 'truncated', 'Metric'] = 'truncated reads'
    df.loc[df['Metric'] == 'identical', 'Metric'] = 'identical reads'
    df.loc[df['Metric'] == 'present', 'Metric'] = 'segmented reads'
    df.loc[df['Metric'] == 'missing', 'Metric'] = 'missing reads'

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

    #! collect meta data for controls and dorado
    # total_reads = df.loc[df['Metric'] == 'total reads', 'Value'].values[0]
    # min_length = df.loc[df['Metric'] == 'min length', 'Value'].values[0]
    # max_length = df.loc[df['Metric'] == 'max length', 'Value'].values[0]
    # mean_length = df.loc[df['Metric'] == 'mean length', 'Value'].values[0]
    # median_length = df.loc[df['Metric'] == 'median length', 'Value'].values[0]
    # n50_length = df.loc[df['Metric'] == 'n50 length', 'Value'].values[0]

    #! add trivial values
    ## controls
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #           "Tool": ["Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl R.", "Ctrl U.", "Ctrl U.", "Ctrl U.", "Ctrl U.", "Ctrl U."],
    #           "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed"],
    #           "Value": [total_reads, 0, 0, total_reads, 0, total_reads, 0, 0, total_reads, 0],
    #           "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #           })
    #     ], ignore_index=True)
    ## dorado
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #           "Tool": ["Dorado"] * 10,
    #           "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "min length", "mean length", "median length", "n50 length", "max length"],
    #           "Value": [total_reads, 0, 0, total_reads, 0, min_length, mean_length, median_length, n50_length, max_length],
    #         #   "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    #           })
    #     ], ignore_index=True)

    # Calculate aggregated metric score for all tools
    # metric_scores = df.groupby("Tool")["Metric Score"].sum().reset_index()
    metric_scores = df_for_agg.groupby('Tool')["Metric Score"].sum().reset_index()

    # Combine Metric Score and Value (without formatting yet)
    df["Combined"] = (
        "$" +
        df["Metric Score"].map(lambda x: f"{x:.3f}") +
        "_{~" +
        df["Value"].map(lambda x: f"{float(x):.1f}" if str(x).replace(".", "", 1).isdigit() else str(x)) +
        "}$"
    )

    # Define custom order for Metric
    custom_order = [
        "median delta", "mad delta", "homogeneity", "segmented reads", "missing reads", 
        "truncated reads", "identical reads", "nt changed", "min length", "mean length", "median length", "n50 length", "max length", "aggregated metric score",
        "flye total length", "flye n50", "flye mean coverage", "svim structural variants",
        "total reads", "time in hh:mm:ss", "memory in mb"
    ]

    # Ensure the 'Metric' column follows the custom order
    df['Metric'] = pd.Categorical(df['Metric'], categories=custom_order, ordered=True)

    # Pivot table: rows = Metric, columns = Tool, values = Combined (MetricScore_Value)
    pivot = df.pivot(index="Metric", columns="Tool", values="Combined")
    score_pivot = df.pivot(index="Metric", columns="Tool", values="Metric Score")

    # Highlight the maximum metric score per metric (row) with \cellcolor{green!15}
    for metric in pivot.index:
        if metric not in score_pivot.index:
            continue
        row_scores = score_pivot.loc[metric]
        if row_scores.dropna().empty:
            continue
        max_score = row_scores.max()
        # Find all columns (tools) with the max score (could be ties)
        max_tools = row_scores[row_scores == max_score].index
        for tool in max_tools:
            if pd.notna(pivot.at[metric, tool]):
                pivot.at[metric, tool] = f"\\cellcolor{{green!15}}{pivot.at[metric, tool]}"

    # Add the raw aggregated metric score as the last row
    agg_row = metric_scores.set_index("Tool").T
    agg_row.index = ["aggregated metric score"]  # Ensure the index is named correctly

    # Convert all values to string first to avoid dtype issues
    agg_row = agg_row.astype(str)

    # Color the max value(s) in the row
    max_val = metric_scores["Metric Score"].max()
    for tool in agg_row.columns:
        val = float(agg_row.loc["aggregated metric score", tool])
        if val == max_val:
            agg_row.loc["aggregated metric score", tool] = f"\\cellcolor{{green!15}}{val:.2f}"
        else:
            agg_row.loc["aggregated metric score", tool] = f"{val:.2f}"

    # Append the aggregated metric score row to the pivot table
    pivot = pd.concat([pivot, agg_row])

    # Map metric names to desired LaTeX names
    metric_name_map = {
        "median delta": r"median delta ($\Delta\mu$)",
        "mad delta": r"mad delta ($\Delta\sigma$)",
        "homogeneity": "homogeneity",
        "segmented reads": "segmented reads",
        "truncated reads": "truncated reads",
        "min length": "min read length",
        "n50 length": "n50 read length",
        "max length": "max read length",
        "flye total length": "flye total length",
        "flye n50": "flye n50",
        "flye mean coverage": "flye mean coverage",
        "svim structural variants": "svim structural variants",
        "aggregated metric score": "aggregated metric score"
    }

    # Only keep and rename the desired metrics
    desired_metrics = list(metric_name_map.keys())
    pivot = pivot.loc[desired_metrics]
    pivot.index = [metric_name_map[m] for m in pivot.index]

    # Reorder columns as requested
    desired_columns = ["Tool", "Dorado", "Tombo", "f5c R.", "f5c E.", "Uncalled4", "Dynamont"]
    # Only keep columns that exist in the table
    pivot = pivot[[col for col in desired_columns if col in pivot.columns]]

    latex_table = pivot.to_latex(na_rep="-", escape=False)

    if output_tex:
        with open(output_tex, "w") as f:
            f.write(latex_table)

    return latex_table

if __name__ == "__main__":
    # Example usage
    args = parse()
    convert_metrics_to_latex(args.input_csv, args.output_tex)
