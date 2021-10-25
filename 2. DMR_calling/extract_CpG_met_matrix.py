import re
import os
import sys
import glob
import numpy
from collections import defaultdict


def parse_subtype(path_subtype):
    print("Start reading from {}".format(path_subtype))
    normal_samples = []
    cancer_samples = []
    with open(path_subtype) as f:
        for line in f:
            if line.startswith("sample"):
                continue
            line_list = line.strip().split("\t")
            if line_list[3] == "0":
                normal_samples.append(line_list[0])
            elif line_list[3] == "1":
                cancer_samples.append(line_list[0])
            else:
                raise ValueError("Invalid sample type, expected to be only 0(normal) and 1(cancer)")
    normal_samples = list(map(lambda k: k.strip(), normal_samples))
    cancer_samples = list(map(lambda k: k.strip(), cancer_samples))
    return normal_samples, cancer_samples


def parse_ref(path_ref: str):
    print("Start reading from {}".format(path_ref))
    site_ref = defaultdict(list)
    with open(path_ref) as f:
        for line in f:
            line_list = line.strip().split()
            site_ref[line_list[0]].append(int(line_list[-1]))
    return site_ref


def parse_overlap(path_overlap: str):
    print("Start reading from {}".format(path_overlap))
    region_ref = defaultdict(dict)
    pattern = re.compile(r"^chr(.*?)$")
    with open(path_overlap) as f:
        for line in f:
            line_list = line.strip().split()
            if line_list[0] == "chrM":
                each_chr = "MT"
            else:
                each_chr = pattern.findall(line_list[0])[0]
            region_ref[each_chr][int(line_list[1])] = int(line_list[2])
    return region_ref


def compare(each_site_ref: numpy.array, each_region_ref: numpy.array, each_chr: str) -> dict:
    site_index = 0
    region_met_pos = dict()
    for each_item in each_region_ref:
        start, end = each_item
        while each_site_ref[site_index] < start:
            if site_index >= len(each_site_ref)-1:
                return region_met_pos
            site_index += 1
        while each_site_ref[site_index] <= end:
            region_met_pos[each_chr + "-" + str(each_site_ref[site_index])] = 1
            if site_index >= len(each_site_ref)-1:
                return region_met_pos
            site_index += 1
    return region_met_pos


def get_met_pos_region(site_ref, region_ref):
    print("Start searching met pos in each region...")
    total_region_met_pos = dict()
    for each_chr in region_ref:
        print("\tStart parsing chromosome {}".format(each_chr))
        each_site_ref = numpy.array(list(site_ref[each_chr]))
        each_region_ref = numpy.array(list(region_ref[each_chr].items()))
        each_chr_region_met_pos = compare(each_site_ref, each_region_ref, each_chr)
        total_region_met_pos.update(each_chr_region_met_pos)
    return total_region_met_pos


def get_sample_files(path_samples: str, path_subtype) -> defaultdict(list):
    print("Start reading from {}".format(path_samples))
    pattern = re.compile(r"^(.*?).met.cov.tmp$")
    files = glob.glob("{}/*.met.cov.tmp".format(path_samples))
    sample_files = defaultdict(list)
    normal_samples, cancer_samples = parse_subtype(path_subtype)
    for each_file in files:
        sample_name = pattern.findall(os.path.basename(each_file))[0]
        if sample_name in normal_samples:
            sample_files["normal"].append(each_file)
        elif sample_name in cancer_samples:
            sample_files["cancer"].append(each_file)
        else:
            raise ValueError("Invalid sample name, neither in normal nor cancer")
    return sample_files


def parse_met(sample_files, total_region_met_pos: dict, path_out):
    samples = []
    site_met_values = defaultdict(dict)
    site_met_values_total = defaultdict(lambda: defaultdict(dict))
    for each_type in sample_files:
        print("Start parsing type {}".format(each_type))
        for each_file in sample_files[each_type]:
            print("\tStart reading from {}".format(each_file))
            pattern = re.compile(r"^(.*?).met.cov.tmp$")
            sample_name = pattern.findall(os.path.basename(each_file))[0]
            samples.append(sample_name)
            with open(each_file) as f:
                for line in f:
                    line_list = line.strip().split()
                    pos = line_list[0] + "-" + line_list[1]
                    if pos in total_region_met_pos:
                        met_ratio = int(line_list[3]) / (int(line_list[3]) + int(line_list[4]))
                        site_met_values[pos][sample_name] = str(met_ratio)
    print("Finish parsing all files, start reshaping...")
    for each_pos in site_met_values:
        each_pos_value = site_met_values[each_pos]
        chrom, start = each_pos.split("-")
        for each_sample in samples:
            if each_sample not in each_pos_value:
                site_met_values_total[chrom][start][each_sample] = "-"
            else:
                site_met_values_total[chrom][start][each_sample] = each_pos_value[each_sample]
    write(site_met_values_total, samples, path_out)


def write(site_met_values_total, samples, path_out):
    print("Start writing to {}".format(path_out))
    with open(path_out, "w+") as out:
        out.write("chrom\tstart\t" + "\t".join(samples) + "\n")
        for each_chr in site_met_values_total:
            each_chr_value = site_met_values_total[each_chr]
            sorted_pos = sorted(each_chr_value, key=lambda k: int(k))
            for each_pos in sorted_pos:
                each_chr_pos_sample_value = [each_chr_value[each_pos][each_sample] for each_sample in samples]
                out.write(each_chr + "\t" + each_pos + "\t" + "\t".join(each_chr_pos_sample_value) + "\n")


def main():
    if len(sys.argv) < 4:
        print("Usage: python3 extractSample.py CpG.GRCh37.txt overlap.bed "
              "samples_subtype_info.xls "
              "intersection.xls")
        exit(0)
    # python extractSample.py CpG.GRCh37.txt BC_NM.sort.merged.bed
    # /data/Epi/likai/mammary_cancer/mCancer/all_sample_cov/met/met_no_scaffold BC_NM.DMR_region.met.xls
    path_ref = str(sys.argv[1])
    path_overlap = str(sys.argv[2])
    path_samples = str(sys.argv[3])
    path_subtype = str(sys.argv[4])
    path_out = str(sys.argv[5])
    site_ref = parse_ref(path_ref)
    region_ref = parse_overlap(path_overlap)
    total_region_met_pos = get_met_pos_region(site_ref, region_ref)
    sample_files = get_sample_files(path_samples, path_subtype)
    parse_met(sample_files, total_region_met_pos, path_out)


if __name__ == '__main__':
    main()
