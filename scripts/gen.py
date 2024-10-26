import os, subprocess, argparse
import pandas as pd
from pybedtools import BedTool


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
    file.loc[:, ["chr", "TSS_start", "TSS_end", "gene_name", "gex"]].to_csv(
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


features = [0, 1, 2, 3, 5, 6, 7, 11, 15]
y = [4]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="specify cells and")
    parser.add_argument("-c", "--cell", type=str, help="which cell line")
    parser.add_argument("-t", "--type", type=str, help="train or val")
    args = parser.parse_args()

    tsv = tsvLoader(cell=args.cell, type=args.type)
    bed = bedLoader(cell=args.cell)
    intersect = tsv.closest(bed, d=True, t="first").to_dataframe(
        names=[
            "chr",
            "tss_start",
            "tss_end",
            "gene_name",
            "gex",
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
    )
    intersect = intersect.drop(columns=[f"dsc_{i}" for i in range(1, 8)])
    histones = ["H3K4me1"]
    for h in histones:
        bed = bedLoader(cell=args.cell, data=h)
        num_peaks = tsv.window(bed, w=2000, c=True).to_dataframe(
            names=["chr", "tss_start", "tss_end", "gene_name", "gex", f"num_{h}_peaks"]
        )
        intersect = pd.concat([intersect, num_peaks[[f"num_{h}_peaks"]]], axis=1)
    os.makedirs(f"./data/{args.cell}-{args.type}", exist_ok=True)
    intersect.drop(columns="gex").to_csv(
        f"./data/{args.cell}-{args.type}/features.tsv",
        sep="\t",
        index=None,
    )
    intersect.gex.to_csv(
        f"./data/{args.cell}-{args.type}/y.tsv", sep="\t", index=None
    )
