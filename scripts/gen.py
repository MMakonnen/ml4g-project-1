import os, subprocess, argparse
import pandas as pd
from pybedtools import BedTool

tsv_keep = [
    "chr",
    "TSS_start",
    "TSS_end",
    "gene_name",
    "gex",
    "strand",
    "gene_start",
    "gene_end",
]

dnase_keep = lambda h: [
    *tsv_keep,
    "dsc_1",
    f"{h}_start",
    f"{h}_end",
    "dsc_2",
    "dsc_3",
    "dsc_4",
    f"{h}_val",
    "dsc_5",
    "dsc_6",
    "dsc_7",
    # "dnase_dist",
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


def DNase(tsv: BedTool, final: pd.DataFrame, cell: str = "X1") -> pd.DataFrame:
    h = "DNase"
    bed = bedLoader(cell=cell, data=h)
    peaks = tsv.window(bed, w=140).to_dataframe(names=dnase_keep(h))
    peaks = (
        peaks.groupby("gene_name")
        .agg(
            **{
                f"{h}_num_peaks": (f"{h}_val", "count"),
                f"{h}_avg_peaks": (f"{h}_val", "mean"),
            }
        )
        .reset_index()
    )
    final = pd.merge(
        final,
        peaks,
        on="gene_name",
        how="left",
    )
    del peaks
    final.fillna({f"{h}_num_peaks": 0}, inplace=True)
    final.fillna({f"{h}_avg_peaks": 0}, inplace=True)
    return final


def H3K4me1(tsv: BedTool, final: pd.DataFrame, cell: str = "X1") -> pd.DataFrame:
    h = "H3K4me1"
    bed = bedLoader(cell=cell, data=h)
    peaks = (
        tsv.slop(b=60, genome="hg38")
        .flank(l=120, r=70, genome="hg38")
        .intersect(bed, wb=True)
        .to_dataframe(names=num_peaks_keep(h))
    )
    peaks = (
        peaks.groupby("gene_name")
        .agg(
            **{
                f"{h}_num_peaks": (f"{h}_signal", "count"),
                f"{h}_avg_peaks": (f"{h}_signal", "mean"),
            }
        )
        .reset_index()
    )
    final = pd.merge(
        final,
        peaks,
        on="gene_name",
        how="left",
    )
    del peaks
    final.fillna({f"{h}_num_peaks": 0}, inplace=True)
    final.fillna({f"{h}_avg_peaks": 0}, inplace=True)
    return final


def H3K4me3(tsv: BedTool, final: pd.DataFrame, cell: str = "X1") -> pd.DataFrame:
    h = "H3K4me3"
    bed = bedLoader(cell=cell, data=h)
    peaks = (
        tsv.slop(b=10, genome="hg38")
        .flank(l=400, r=100, genome="hg38")
        .intersect(bed, wb=True)
        .to_dataframe(names=num_peaks_keep(h))
    )
    peaks = (
        peaks.groupby("gene_name")
        .agg(
            **{
                f"{h}_num_peaks": (f"{h}_signal", "count"),
                f"{h}_avg_peaks": (f"{h}_signal", "mean"),
            }
        )
        .reset_index()
    )
    final = pd.merge(
        final,
        peaks,
        on="gene_name",
        how="left",
    )
    del peaks
    final.fillna({f"{h}_num_peaks": 0}, inplace=True)
    final.fillna({f"{h}_avg_peaks": 0}, inplace=True)
    return final


def H3K27ac(tsv: BedTool, final: pd.DataFrame, cell: str = "X1") -> pd.DataFrame:
    h = "H3K27ac"
    bed = bedLoader(cell=cell, data=h)
    peaks = (
        tsv.slop(b=10, genome="hg38")
        .flank(l=270, r=70, genome="hg38")
        .intersect(bed, wb=True)
        .to_dataframe(names=num_peaks_keep(h))
    )
    peaks = (
        peaks.groupby("gene_name")
        .agg(
            **{
                f"{h}_num_peaks": (f"{h}_signal", "count"),
                f"{h}_avg_peaks": (f"{h}_signal", "mean"),
            }
        )
        .reset_index()
    )
    final = pd.merge(
        final,
        peaks,
        on="gene_name",
        how="left",
    )
    del peaks
    final.fillna({f"{h}_num_peaks": 0}, inplace=True)
    final.fillna({f"{h}_avg_peaks": 0}, inplace=True)
    return final


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="specify cells and")
    parser.add_argument("-c", "--cell", type=str, help="which cell line")
    parser.add_argument("-t", "--type", type=str, help="train or val")
    args = parser.parse_args()

    tsv = tsvLoader(cell=args.cell, type=args.type)

    final = tsv.to_dataframe(names=tsv_keep)
    final = DNase(tsv, final, args.cell)
    final = H3K4me1(tsv, final, args.cell)
    final = H3K4me3(tsv, final, args.cell)
    final = H3K27ac(tsv, final, args.cell)

    os.makedirs(f"./data/{args.cell}-{args.type}", exist_ok=True)
    final.drop(columns="gex").to_csv(
        f"./data/{args.cell}-{args.type}/features.tsv",
        sep="\t",
        index=None,
    )
    final.gex.to_csv(f"./data/{args.cell}-{args.type}/y.tsv", sep="\t", index=None)
