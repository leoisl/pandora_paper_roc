import argparse
import glob

def get_args():
    parser = argparse.ArgumentParser(description='Count number of masked bases in direcotry.')
    parser.add_argument('--dir', type=str, help='Directory with samples where **.mask.bed files are', required=True)
    args = parser.parse_args()
    return args

def get_nb_of_masked_bases(file):
    nb_of_masked_bases = 0
    with open(file) as fh:
        for line in fh:
            line_split = line.split()
            begin = int(line_split[1])
            end = int(line_split[2])
            size = end-begin
            nb_of_masked_bases += size
    return nb_of_masked_bases

def process(dir):
    files = glob.glob(dir+"/**/*.mask.bed", recursive=True)
    for file in files:
        nb_of_masked_bases = get_nb_of_masked_bases(file)
        print(f"{file}: {nb_of_masked_bases} masked bases.")

def main():
    args = get_args()
    process(args.dir)

if __name__=="__main__":
    main()