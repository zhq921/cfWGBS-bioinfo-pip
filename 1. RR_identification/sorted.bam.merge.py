import os
import sys
import glob
import time
import subprocess
import multiprocessing


def merge(path):
    print(f"Start parsing {path}")
    path_out = path + ".merge"
    #cmd1 = f"samtools sort -@ 10 {path} -o {path_out}"
    #cmd2 = f"samtools index {path_out} -@ 10"
    cmd1 = f"bedtools merge -bed -i {path} > {path_out}"
    start_time = time.time()
    try:
        subprocess.check_output(cmd1, shell=True)
        #subprocess.check_output(cmd2, shell=True)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print(f"Finish {path}, time used: {time_used} min")
    except Exception as e:
        print(f"Error: {e}")
        exit(-1)


def main():
    if len(sys.argv) < 2:
        print("Usage: python sorted.bam.merge.py all_sample_bam")
        exit(0)
    path = sys.argv[1]
    pool = multiprocessing.Pool(processes=10)
    files = glob.glob(f"{path}/*.sorted")
    for each_file in files:
        if not os.path.exists(f"{each_file}.merge"):
            pool.apply_async(merge, (each_file,))
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
