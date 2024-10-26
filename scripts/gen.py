import os, subprocess
import pandas as pd
import pybedtools

def tsvLoader(path="./data/CAGE-train/X1_train_info.tsv", scratch="./data/tmp"):
    file = pd.read_csv(path, sep="\t")
    file.loc[:, ["chr", "TSS_start", "TSS_end"]].to_csv(scratch+"/tmp_tsv.bed", sep='\t', header=False, index=False)
    os.environ["LC_ALL"] = "C"
    command = ["sort", "-k1,1", "-k2,2n", scratch+"/tmp_tsv.bed", "-o", scratch+"/tmp_tsv.bed"]
    subprocess.run(command, check=True)
    d = pybedtools.BedTool(scratch+"/tmp_tsv.bed")
    return d

def bedLoader(path="./data/DNase-bed/X1.bed"):
    d = pybedtools.BedTool(path)
    return d

if __name__ == "__main__":
    tsv = tsvLoader()
    bed = bedLoader()
    tsv.closest(bed, d=True)