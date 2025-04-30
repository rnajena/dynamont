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
    parser.add_argument("--tombo", type=str, default=None, help="Path to the reads metrics json file")
    parser.add_argument("--time_tombo", type=str, default=None, help="Path to the tolls time file")
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

    # NA rna004
    if args.tombo and args.tombo != '':
        jsons["tombo"] = args.tombo

    if args.time_tombo and args.time_tombo != '':
        times["tombo"] = args.time_tombo

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

    control = pd.read_csv(args.control, sep="\t")
    for _, row in control.iterrows():
        new_entry = pd.DataFrame({
            "Tool": ["Control Random", "Control Uniform"],
            "Value": [row["Value"], row["Value"]],
            "Metric": [row["Metric"].lower() + '_length', row["Metric"].lower() + '_length']
        })
        scores = pd.concat([scores, new_entry], ignore_index=True)

    # fix names
    scores["Tool"] = scores["Tool"].replace(
        {
            "dynamont": "Dynamont",
            "Dynamont NT": "Dynamont",
            "f5c_eventalign": "f5c Eventalign",
            "f5c_resquiggle": "f5c Resquiggle",
            "uncalled4": "Uncalled4",
            "tombo": "Tombo",
        }
    )

    # Ensure Value column is numeric where needed
    # scores["Value"] = pd.to_numeric(scores["Value"], errors="coerce")

    # Exclude specific metrics (e.g., "Time in hh:mm:ss") from the Metric Score calculation
    excluded_metrics = ["Time in hh:mm:ss"]
    numeric_scores = scores[~scores["Metric"].isin(excluded_metrics)]
    numeric_scores["Value"] = pd.to_numeric(numeric_scores["Value"], errors="coerce")

    # Calculate Metric Score only for numeric values
    scores["Metric Score"] = numeric_scores.groupby("Metric")["Value"].transform(
        lambda x: x / x.max() if x.max() > 0 else 0
    )
    
    # Fill non-numeric rows with NaN for "Metric Score"
    scores["Metric Score"] = scores["Metric Score"].fillna(0)

    # calculate metric score
    # scores["Metric Score"] = scores.groupby("Metric")["Value"].transform(lambda x: x / x.max() if x.max() > 0 else 0)

    # Adjust Metric Score for specific metrics
    scores.loc[scores["Metric"].isin(["Homogeneity", "Mad Delta", "missing", "truncated", "nt_changed", "min_length"]), "Metric Score"] = 1 - scores["Metric Score"]

    # Finalize the DataFrame
    scores = scores[["Tool", "Metric", "Value", "Metric Score"]]
    scores = scores.sort_values(by=["Metric", "Tool"])
    scores.reset_index(drop=True, inplace=True)

    scores.to_csv(args.outfile, sep="\t", index=False)

if __name__ == '__main__':
    main()