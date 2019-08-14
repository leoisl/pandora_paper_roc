import random
import pysam
import sys
import os

SLICE = 100
SAMPLES = 5


def main():
    filepath = sys.argv[1]
    # load fasta file in and store as a string.
    fh = pysam.FastaFile(filepath)
    # get 5 random indices of 100bp (non-overlapping)
    assert fh.nreferences == 1
    ref = fh.references[0]
    length = fh.lengths[0]
    idxs = [1156]
    # idxs = random.sample(range(length - SLICE), SAMPLES)

    # make sure all idxs are more than SLICE away from each other
    while not all(idxs[i + 1] - idxs[i] > SLICE for i in range(len(idxs) - 1)):
        idxs = random.sample(range(length - SLICE), SAMPLES)

    for start in idxs:
        end = start + SLICE
        sequence = fh.fetch(ref, start=start, end=end)
        outpath = os.path.join(sys.argv[2], f"slice_{start}-{end}.fa")
        # write slices to file as "mini" fastas
        with open(outpath, "w") as fout:
            fout.write(f">{ref}\n")
            fout.write(sequence + "\n")

    fh.close()


if __name__ == "__main__":
    main()
