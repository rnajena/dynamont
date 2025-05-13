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
    # ## dorado
    # df = pd.concat(
    #     [
    #         df, pd.DataFrame({
    #           "Tool": ["Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado"],
    #           "Metric": ["segmented reads", "missing reads", "truncated reads", "identical reads", "nt changed", "min length", "mean length", "median length", "n50 length", "max length"],
    #           "Value": [total_reads, 0, 0, total_reads, 0, min_length, mean_length, median_length, n50_length, max_length],
    #           "Metric Score": [1.0, 1.0, 1.0, 1.0, 1.0, df.loc[(df['Tool'] == 'Ctrl R.') & (df['Metric'] == 'min length'), 'Metric Score'].squeeze(), 1.0, 1.0, 1.0, 1.0],
    #           })
    #     ], ignore_index=True)

    # Calculate aggregated metric score for all tools
    metric_scores = df.groupby("Tool")["Metric Score"].sum().reset_index()
    metric_scores["Metric"] = "aggregated metric score"
    metric_scores["Value"] = "-"
    metric_scores["Combined"] = (
        "$" +
        metric_scores["Metric Score"].map(lambda x: f"{x:.3f}") +
        "_{" +
        metric_scores["Value"].map(str) +
        "}$"
    )
    
    # Append the aggregated metric scores to the DataFrame
    df = pd.concat([df, metric_scores], ignore_index=True)

    # Combine Metric Score and Value
    df["Combined"] = (
        "$" +
        df["Metric Score"].map(lambda x: f"{x:.3f}") +
        "_{" +
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
    
    # Add the raw aggregated metric score as the last row
    aggregated_metric_score_row = df[df["Metric"] == "aggregated metric score"].set_index("Metric").pivot(columns="Tool", values="Metric Score")
    aggregated_metric_score_row.index = ["aggregated metric score"]  # Ensure the index is named correctly
    # Format the aggregated metric score to 3 decimal places
    aggregated_metric_score_row = aggregated_metric_score_row.apply(
        lambda col: col.map(lambda x: f"{x:.1f}" if pd.notna(x) else x)
    )
    # Append the aggregated metric score row to the pivot table
    pivot = pd.concat([pivot, aggregated_metric_score_row])
    
    latex_table = pivot.to_latex(na_rep="-", escape=False)

    if output_tex:
        with open(output_tex, "w") as f:
            f.write(latex_table)

    return latex_table

if __name__ == "__main__":
    # Example usage
    args = parse()
    convert_metrics_to_latex(args.input_csv, args.output_tex)
