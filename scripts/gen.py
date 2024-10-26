import os, subprocess, argparse
import pandas as pd
from pybedtools import BedTool

tsv_keep = ["chr", "TSS_start", "TSS_end", "gene_name", "gex", "strand"]

dnase_keep = [
    *tsv_keep,
    "dsc_1",
    "dnase_start",
    "dnase_end",
    "dsc_2",
    "dsc_3",
    "dsc_4",
    "dnase_val",
    "dsc_5",
    "dsc_6",
    "dsc_7",
    "dnase_dist",
]

num_peaks_keep = lambda h: [
    *tsv_keep,
    f"{h}_chr",
    f"{h}_start",
    f"{h}_end",
    f"{h}_name",
    f"{h}_signal",
    f"{h}_strand",
    f"{h}_score",
    f"{h}_p_value",
    f"{h}_q_value",
    f"{h}_read_count",
]


def tsvLoader(
    cell: str = "X1", type: str = "train", scratch: str = "./data/tmp"
) -> BedTool:
    file = pd.concat(
        [
            pd.read_csv(f"./data/CAGE-train/{cell}_{type}_info.tsv", sep="\t"),
            pd.read_csv(f"./data/CAGE-train/{cell}_{type}_y.tsv", sep="\t")[["gex"]],
        ],
        axis=1,
    )
    file.loc[:, tsv_keep].to_csv(
        scratch + "/tmp_tsv.bed", sep="\t", header=False, index=False
    )
    os.environ["LC_ALL"] = "C"
    command = [
        "sort",
        "-k1,1",
        "-k2,2n",
        scratch + "/tmp_tsv.bed",
        "-o",
        scratch + "/tmp_tsv.bed",
    ]
    subprocess.run(command, check=True)
    return BedTool(scratch + "/tmp_tsv.bed")


def bedLoader(cell: str = "X1", data: str = "DNase") -> BedTool:
    os.environ["LC_ALL"] = "C"
    command = [
        "sort",
        "-k1,1",
        "-k2,2n",
        f"./data/{data}-bed/{cell}.bed",
        "-o",
        f"./data/tmp/tmp-{data}.bed",
    ]
    subprocess.run(command, check=True)

    return BedTool(f"./data/tmp/tmp-{data}.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="specify cells and")
    parser.add_argument("-c", "--cell", type=str, help="which cell line")
    parser.add_argument("-t", "--type", type=str, help="train or val")
    args = parser.parse_args()

    tsv = tsvLoader(cell=args.cell, type=args.type)
    bed = bedLoader(cell=args.cell)
    intersect = tsv.closest(bed, d=True, t="first").to_dataframe(names=dnase_keep)
    final = intersect.drop(columns=[f"dsc_{i}" for i in range(1, 8)])
    del intersect
    histones = ["H3K4me3", "H3K27ac"]
    for h in histones:
        bed = bedLoader(cell=args.cell, data=h)
        num_peaks = tsv.window(bed, w=2000).to_dataframe(names=num_peaks_keep(h))
        closest_peak = tsv.closest(bed, d=True, t="first").to_dataframe(
            names=[
                *num_peaks_keep(h),
                f"{h}_distance",
            ]
        )
        peak_summary = (
            num_peaks.groupby("gene_name")
            .agg(
                **{
                    f"{h}_num_peaks": (f"{h}_signal", "count"),
                    f"{h}_avg_peaks": (f"{h}_signal", "mean"),
                }
            )
            .reset_index()
        )
        peak_summary = pd.merge(
            closest_peak[["gene_name", f"{h}_signal", f"{h}_distance"]],
            peak_summary,
            on="gene_name",
            how="left",
        )
        peak_summary.fillna({f"{h}_num_peaks": 0}, inplace=True)
        peak_summary.fillna({f"{h}_avg_peaks": 0}, inplace=True)
        final = pd.merge(final, peak_summary, on="gene_name", how="left")
    del num_peaks, peak_summary
    os.makedirs(f"./data/{args.cell}-{args.type}", exist_ok=True)
    final.drop(columns="gex").to_csv(
        f"./data/{args.cell}-{args.type}/features.tsv",
        sep="\t",
        index=None,
    )
    final.gex.to_csv(f"./data/{args.cell}-{args.type}/y.tsv", sep="\t", index=None)
