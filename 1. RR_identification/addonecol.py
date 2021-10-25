import os
import sys
import glob
import time
import subprocess
import multiprocessing


def addone(path):
    print(f"Start parsing {path}")
    path_out = path.split("/")[-1]
    #cmd1 = f"samtools sort -@ 10 {path} -o {path_out}"
    #cmd2 = f"samtools index {path_out} -@ 10"
    cmd1 = f"sed -i s/$/\\\\t{path_out}/ {path}"
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
        print("Usage: python addonecol.py all_sample_bam")
        exit(0)
    path = sys.argv[1]
    pool = multiprocessing.Pool(processes=10)
    files = glob.glob(f"{path}/*.merge")
    for each_file in files:
        pool.apply_async(addone, (each_file,))
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
