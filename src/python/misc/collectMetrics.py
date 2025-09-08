#!/usr/bin/env python
# author: Jannes Spangenberg
# e-mail: jannes.spangenberg@uni-jena.de
# github: https://github.com/JannesSP
# website: https://jannessp.github.io

from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser, Namespace
import pandas as pd
import json

def parse() -> Namespace:
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("scores", type=str, help="Path to the scores txt file")
    # read metrics after segmentation
    parser.add_argument("dynamont", type=str, help="Path to the reads metrics json file")
    parser.add_argument("uncalled4", type=str, help="Path to the reads metrics json file")
    parser.add_argument("f5c_eventalign", type=str, help="Path to the reads metrics json file")
    parser.add_argument("f5c_resquiggle", type=str, help="Path to the reads metrics json file")
    parser.add_argument("control", type=str, help="Path to the reads metrics json file")
    parser.add_argument("outfile", type=str, help="Path to the output file")
    parser.add_argument("time_dynamont", type=str, help="Path to the tools time file")
    parser.add_argument("time_uncalled4", type=str, help="Path to the tools time file")
    parser.add_argument("time_f5c_eventalign", type=str, help="Path to the tools time file")
    parser.add_argument("time_f5c_resquiggle", type=str, help="Path to the tools time file")
    parser.add_argument("subtools_dorado", type=str, help="Path to the downstream tool metrics file")
    parser.add_argument("subtools_dynamont", type=str, help="Path to the downstream tool metrics file")
    parser.add_argument("subtools_uncalled4", type=str, help="Path to the downstream tool metrics file")
    parser.add_argument("subtools_f5c_eventalign", type=str, help="Path to the downstream tool metrics file")
    parser.add_argument("subtools_f5c_resquiggle", type=str, help="Path to the downstream tool metrics file")
    parser.add_argument("--tombo", type=str, default=None, help="Path to the reads metrics json file")
    parser.add_argument("--time_tombo", type=str, default=None, help="Path to the tools time file")
    parser.add_argument("--subtools_tombo", type=str, default=None, help="Path to the downstream tool metrics time file")
    return parser.parse_args()

def main() -> None:
    args = parse()
    scores = pd.read_csv(args.scores, sep="\t")
    scores.rename(columns={"Median Score": "Value", "Segment Quality" : "Metric"}, inplace=True)

    jsons = {
        "dynamont": args.dynamont,
        "uncalled4": args.uncalled4,
        "f5c_eventalign": args.f5c_eventalign,
        "f5c_resquiggle": args.f5c_resquiggle
    }

    times = {
        "dynamont": args.time_dynamont,
        "uncalled4": args.time_uncalled4,
        "f5c_eventalign": args.time_f5c_eventalign,
        "f5c_resquiggle": args.time_f5c_resquiggle
    }

    downstream_tools = {
        "dorado" : args.subtools_dorado,
        "dynamont": args.subtools_dynamont,
        "uncalled4": args.subtools_uncalled4,
        "f5c_eventalign": args.subtools_f5c_eventalign,
        "f5c_resquiggle": args.subtools_f5c_resquiggle
    }

    # NA rna004
    if args.tombo and args.tombo != '':
        jsons["tombo"] = args.tombo

    if args.time_tombo and args.time_tombo != '':
        times["tombo"] = args.time_tombo

    if args.subtools_tombo and args.subtools_tombo != '':
        downstream_tools["tombo"] = args.subtools_tombo

    for name, json_path in jsons.items():
        with open(json_path, "r") as json_file:
            json_data = json.load(json_file)
            for metric, value in json_data.items():
                if metric == "lengths":
                    continue
                # scores.loc[name, metric] = value
                new_entry = pd.DataFrame({
                    "Tool": [name],
                    "Value": [value],
                    "Metric": [metric.lower().replace('n50', 'n50_length')]
                })
                scores = pd.concat([scores, new_entry], ignore_index=True)

    for name, time_path in times.items():
        with open(time_path, "r") as time_file:
            time = time_file.readline()[14:22]
            memory = time_file.readline().strip()[13:].split(" MB")[0]
            new_entry = pd.DataFrame({
                "Tool": [name, name],
                "Value": [time, memory],
                "Metric": ["Time in hh:mm:ss", "Memory in MB"]
            })  
            scores = pd.concat([scores, new_entry], ignore_index=True)

    for name, downstream_path in downstream_tools.items():
        with open(downstream_path, "r") as downstream_file:
            total_assembly_length = int(downstream_file.readline().strip().split(': ')[1])
            n50 = int(downstream_file.readline().strip().split(': ')[1])
            mean_cov = float(downstream_file.readline().strip().split(': ')[1])
            struct_vars = int(downstream_file.readline().strip().split(': ')[1])

            new_entry = pd.DataFrame({
                "Tool": [name, name, name, name],
                "Value": [total_assembly_length, n50, mean_cov, struct_vars],
                "Metric": ["flye total length", "flye n50", "flye mean coverage", "SVIM structural variants"]
            })
            scores = pd.concat([scores, new_entry], ignore_index=True)

    #! add default control values to dorado
    control = pd.read_csv(args.control, sep="\t")
    for _, row in control.iterrows():
        new_entry = pd.DataFrame({
            "Tool": ["Dorado"],
            "Value": [row["Value"]],
            "Metric": [row["Metric"].lower() + '_length']
        })
        scores = pd.concat([scores, new_entry], ignore_index=True)

    # print(scores.loc[scores["Metric"] == "total", "Value"].values)
    total_reads = scores.loc[scores["Metric"] == "total", "Value"].values[0]
    new_entry = pd.DataFrame({
        "Tool": ["Dorado", "Dorado", "Dorado", "Dorado", "Dorado", "Dorado"],
        "Metric": ["total", "present", "missing", "truncated", "identical", "nt changed"],
        "Value": [total_reads, total_reads, 0, 0, total_reads, 0],
    })
    scores = pd.concat([scores, new_entry], ignore_index=True)

    #! remove controls and dorado
    scores = scores[scores["Tool"] != "Control Random"]
    scores = scores[scores["Tool"] != "Control Uniform"]
    # scores = scores[scores["Tool"] != "Dorado"]

    # fix names
    scores["Tool"] = scores["Tool"].replace(
        {
            "dynamont": "Dynamont",
            "Dynamont NT": "Dynamont",
            "f5c_eventalign": "f5c Eventalign",
            "f5c_resquiggle": "f5c Resquiggle",
            "uncalled4": "Uncalled4",
            "tombo": "Tombo",
            "dorado": "Dorado",
        }
    )

    # Remove unwanted metrics
    # Exclude specific metrics (e.g., "Time in hh:mm:ss") from the Metric Score calculation
    excluded_metrics = ["missing reads", "identical reads", "Time in hh:mm:ss", "Memory in MB"]
    numeric_scores = scores[~scores["Metric"].isin(excluded_metrics)]
    numeric_scores["Value"] = pd.to_numeric(numeric_scores["Value"], errors="coerce")

    # Calculate Metric Score only for numeric values
    scores["Metric Score"] = numeric_scores.groupby("Metric")["Value"].transform(
        lambda x: x / x.max() if x.max() > 0 else 0
    )
    # print("GROUP: ", scores["Metric Score"])
    # exit(1)

    # Fill non-numeric rows with NaN for "Metric Score"
    scores["Metric Score"] = scores["Metric Score"].fillna(0)

    # calculate metric score
    # scores["Metric Score"] = scores.groupby("Metric")["Value"].transform(lambda x: x / x.max() if x.max() > 0 else 0)

    # Adjust Metric Score for specific metrics
    scores.loc[scores["Metric"].isin(["Homogeneity", "missing", "truncated", "nt_changed", "min_length"]), "Metric Score"] = 1 - scores["Metric Score"]

    # Calculate Metric Score only for numeric values
    # def metric_score(series, lower_is_better=False):
    #     if series.max() == series.min():
    #         return pd.Series([1.0] * len(series), index=series.index)
    #     if lower_is_better:
    #         return (series.max() - series) / (series.max() - series.min())
    #     else:
    #         return (series - series.min()) / (series.max() - series.min())

    # # Define which metrics are "lower is better"
    # lower_is_better_metrics = ["Homogeneity", "missing", "truncated", "nt_changed", "min_length"]

    # # Calculate scores for each metric
    # scores["Metric Score"] = 0.0
    # for metric in numeric_scores["Metric"].unique():
    #     mask = scores["Metric"] == metric
    #     lower_is_better = metric in lower_is_better_metrics
    #     values = pd.to_numeric(scores.loc[mask, "Value"], errors="coerce")
    #     scores.loc[mask, "Metric Score"] = metric_score(values, lower_is_better=lower_is_better)

    # # Fill non-numeric rows with NaN for "Metric Score"
    # scores["Metric Score"] = scores["Metric Score"].fillna(0)

    # Finalize the DataFrame
    scores = scores[["Tool", "Metric", "Value", "Metric Score"]]
    scores = scores.sort_values(by=["Metric", "Tool"])
    scores.reset_index(drop=True, inplace=True)

    print("\nWriting to", args.outfile, "\n")
    scores.to_csv(args.outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()