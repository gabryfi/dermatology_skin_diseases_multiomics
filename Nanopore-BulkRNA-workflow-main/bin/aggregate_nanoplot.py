#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def parse_nanostats(stats_file, stage):
    """
    Parse NanoPlot NanoStats.txt file.
    """

    values = {}
    with open(stats_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Lines such as:
            # Mean read length: 858.8
            # Number of reads: 16,476,351.0
            match = re.match(r"^(.+?):\s+([0-9,.]+)\s*$", line)

            if match:
                key = match.group(1).strip()
                value = match.group(2).replace(",", "")

                try:
                    value = float(value)
                except ValueError:
                    continue

                values[key] = value

    nanoplot_dir = stats_file.parent.name # Recover barcode_condition from parent directory
    prefix = "Nanoplot_"
    sample_name = nanoplot_dir

    if sample_name.startswith(prefix):
        sample_name = sample_name[len(prefix):]

    
    sample_name = sample_name.replace("_filtered", "")
    parts = sample_name.split("_", 1)

    if len(parts) == 2:
        barcode = parts[0]
        condition = parts[1]
    else:
        barcode = sample_name
        condition = "unknown"

    return {
        "barcode": barcode,
        "condition": condition,
        "stage": stage,
        "mean_read_length": values.get("Mean read length"),
        "mean_quality": values.get("Mean read quality"),
        "median_read_length": values.get("Median read length"),
        "median_quality": values.get("Median read quality"),
        "read_length_N50": values.get("Read length N50"),
        "number_of_reads": values.get("Number of reads"),
        "total_bases": values.get("Total bases"),
    }


def find_nanostats(paths, stage):
    """
    Recursively find NanoStats.txt files.
    """

    records = []
    for path in paths:
        path = Path(path)
        if path.is_dir():
            stats_files = path.rglob("NanoStats.txt")
        elif path.name == "NanoStats.txt":
            stats_files = [path]
        else:
            continue

        for stats_file in stats_files:
            records.append(parse_nanostats(stats_file, stage))
    return records


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate NanoPlot NanoStats.txt files"
    )
    parser.add_argument(
        "--pre",
        nargs="+",
        required=True,
        help="NanoPlot pre-filter directories"
    )

    parser.add_argument(
        "--post",
        nargs="+",
        required=True,
        help="NanoPlot post-filter directories"
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # Parse NanoStats
    # ---------------------------------------------------------

    records_pre = find_nanostats(args.pre, "pre")
    records_post = find_nanostats(args.post, "post")

    records = records_pre + records_post

    if not records:
        raise RuntimeError("No NanoStats.txt files found.")

    df = pd.DataFrame(records)

    # ---------------------------------------------------------
    # Sort
    # ---------------------------------------------------------

    stage_order = {"pre": 0, "post": 1}
    df["_stage_order"] = df["stage"].map(stage_order)
    df = (
        df
        .sort_values(["barcode", "condition", "_stage_order"])
        .drop(columns="_stage_order")
    )

    # ---------------------------------------------------------
    # Save TSV
    # ---------------------------------------------------------

    summary_file = outdir / "NanoPlot_summary.tsv"
    df.to_csv(summary_file, sep="\t", index=False, float_format="%.2f")
    print(f"Summary written to: {summary_file}")

    # ---------------------------------------------------------
    # Plot helper
    # ---------------------------------------------------------

    def make_plot(column, ylabel, filename, title):

        plot_df = df.copy()
        plot_df["sample"] = (
            plot_df["barcode"] + "_" + plot_df["condition"]
        )

        pivot = plot_df.pivot(
            index="sample",
            columns="stage",
            values=column
        )

        ax = pivot.plot(
            kind="bar",
            figsize=(12, 6)
        )

        ax.set_xlabel("Sample")
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        plt.savefig(
            outdir / filename,
            dpi=300
        )

        plt.close()

    # ---------------------------------------------------------
    # Generate plots
    # ---------------------------------------------------------

    make_plot(
        "number_of_reads",
        "Number of reads",
        "NanoPlot_reads.png",
        "Number of reads: pre vs post filtering"
    )

    make_plot(
        "total_bases",
        "Total bases",
        "NanoPlot_bases.png",
        "Total bases: pre vs post filtering"
    )

    make_plot(
        "mean_quality",
        "Mean read quality",
        "NanoPlot_mean_quality.png",
        "Mean read quality: pre vs post filtering"
    )

    make_plot(
        "mean_read_length",
        "Mean read length",
        "NanoPlot_mean_read_length.png",
        "Mean read length: pre vs post filtering"
    )

    make_plot(
        "median_read_length",
        "Median read length",
        "NanoPlot_median_read_length.png",
        "Median read length: pre vs post filtering"
    )

    make_plot(
        "read_length_N50",
        "Read length N50",
        "NanoPlot_N50.png",
        "Read length N50: pre vs post filtering"
    )

    # ---------------------------------------------------------
    # HTML report
    # ---------------------------------------------------------

    html_file = outdir / "NanoPlot_summary.html"

    html = f"""
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>NanoPlot summary</title>

<style>

body {{
    font-family: Arial, sans-serif;
    margin: 40px;
}}

h1 {{
    margin-bottom: 30px;
}}

table {{
    border-collapse: collapse;
    width: 100%;
    margin-bottom: 30px;
}}

th, td {{
    border: 1px solid #cccccc;
    padding: 6px;
    text-align: right;
}}

th {{
    background-color: #eeeeee;
}}

td:first-child,
td:nth-child(2),
td:nth-child(3) {{
    text-align: left;
}}

img {{
    max-width: 100%;
    margin-bottom: 40px;
}}

</style>

</head>

<body>

<h1>NanoPlot summary</h1>

<h2>Summary table</h2>

{df.to_html(index=False, float_format=lambda x: f"{x:.2f}")}

<h2>Plots</h2>

<h3>Number of reads</h3>
<img src="NanoPlot_reads.png">

<h3>Total bases</h3>
<img src="NanoPlot_bases.png">

<h3>Mean read quality</h3>
<img src="NanoPlot_mean_quality.png">

<h3>Mean read length</h3>
<img src="NanoPlot_mean_read_length.png">

<h3>Median read length</h3>
<img src="NanoPlot_median_read_length.png">

<h3>Read length N50</h3>
<img src="NanoPlot_N50.png">

</body>
</html>
"""

    with open(html_file, "w") as f:
        f.write(html)

    print(f"HTML report written to: {html_file}")


if __name__ == "__main__":
    main()