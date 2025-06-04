from sys import argv
from scripts.gen import gen
from scripts.model import pipeline

if __name__ == "__main__":
    gen("X1", "train")
    gen("X1", "val")
    gen("X2", "train")
    gen("X2", "val")
    gen("X3", "test")
    pipeline(argv)
