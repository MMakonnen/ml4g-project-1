import pandas as pd
from tqdm import tqdm

def tsvLoader(path="./data/CAGE-train/X1_train_info.tsv", scratch="./data/tmp"):
    file = pd.read_csv(path, sep="\t")
    return file.loc[:, ["chr", "TSS_start", "TSS_end"]]

def bedLoader(path="./data/DNase-bed/X1.bed"):
    file = pd.read_csv(path, sep="\t", header=None)
    return file.loc[:, [0, 1, 2, 6]]

def intersection(tsv: pd.DataFrame, bed: pd.DataFrame, window=0):
    dnase = []
    for _, tsvrow in tqdm(tsv.loc[:50].iterrows(), total=51):
        tmp = []
        for _, bedrow in bed.iterrows():
            if tsvrow["chr"] == bedrow[0] and checkWindow((tsvrow["TSS_start"], tsvrow["TSS_end"]), (bedrow[1], bedrow[2])):
                tmp.append(bedrow[6])
        dnase.append(sum(tmp))
    return tsv.loc[:50].assign(dnase=dnase)

def checkWindow(x, y, window=0):
    if (x[0] - window < y[0] and y[0] < x[1] + window) or (x[0] - window < y[1] and y[1] < x[1] + window):
        return True
    return False