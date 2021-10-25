import re
import sys
import glob
import subprocess


def get_files(path_samples):
    files = glob.glob(f"{path_samples}/met_mat*")
    print("total files ", files)
    return files


def get_genes(files):
    genes = []
    for each_file in files:
        print(f"Start parsing file {each_file}")
        cmd = f"sed -n '2,$p' {each_file} | cut -f 1,2"
        result = subprocess.check_output(cmd, shell=True)
        pos = result.decode().strip().split("\n")
        pos = list(map(lambda k: re.sub("\t", "-", k), pos))
        genes.extend(pos)
    genes = list(set(genes))
    return genes


def write(genes, path_out):
    print(f"Start writing to {path_out}")
    with open(path_out, "w+") as out:
        for each_gene in genes:
            out.write(each_gene + "\t" + "1" + "\n")


def main():
    path_samples = sys.argv[1]
    path_out = sys.argv[2]
    files = get_files(path_samples)
    genes = get_genes(files)
    write(genes, path_out)


if __name__ == '__main__':
    main()
