import sys
import glob


def get_files(path_samples):
    files = glob.glob(f"{path_samples}/met_mat*")
    print("total files ", files)
    return files


def parse_pos(path_pos):
    print(f"Start parsing {path_pos}")
    total_pos = []
    with open(path_pos) as f:
        for line in f:
            line_list = line.strip().split("\t")
            total_pos.append(line_list[0])
    return total_pos


def parse_file(total_pos, files):
    for each_file in files:
        print(f"Start parsing {each_file}")
        title = ""
        path_out = each_file + ".update"
        this_sample_pos = dict()
        this_sample_length = []
        with open(each_file) as f:
            for line in f:
                line_list = line.strip().split("\t")
                if line.startswith("chrom"):
                    line_list.insert(2, "end")
                    title = "\t".join(line_list)
                else:
                    this_sample_length.append(len(line_list[2:]))
                    this_pos = f"{line_list[0]}-{line_list[1]}"
                    this_sample_pos[this_pos] = "\t".join(line_list[2:])
        this_sample_length = list(set(this_sample_length))
        if len(this_sample_length) > 1:
            raise ValueError(f"Invalid line length {len(this_sample_length)}")
        length = this_sample_length[0]
        print(f"Finish parsing {each_file}, start writing to {path_out}")
        with open(path_out, "w+") as out:
            out.write(title + "\n")
            for each_pos in total_pos:
                if each_pos in this_sample_pos:
                    this_pos_value = this_sample_pos[each_pos]
                    chrom, start = each_pos.split("-")
                    out_line = chrom + "\t" + start + "\t" + str(int(start) + 1) + "\t" + this_pos_value
                    out.write(out_line + "\n")
                else:
                    this_pos_value = "\t".join(["-"] * length)
                    chrom, start = each_pos.split("-")
                    out_line = chrom + "\t" + start + "\t" + str(int(start) + 1) + "\t" + this_pos_value
                    out.write(out_line + "\n")


def main():
    path_samples = sys.argv[1]
    path_pos = sys.argv[2]
    files = get_files(path_samples)
    total_pos = parse_pos(path_pos)
    parse_file(total_pos, files)


if __name__ == '__main__':
    main()
