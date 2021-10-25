import os
import sys
import glob
import time
import subprocess
import multiprocessing


def combine(path,path_out):
    print(f"Start parsing {path}")
    #cmd1 = f"samtools sort -@ 10 {path} -o {path_out}"
    #cmd2 = f"samtools index {path_out} -@ 10"
    cmd1 = f"awk '$1 ~ /(^[1-9]|^X|^Y|^MT)/' {path} >> {path_out}"
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
    os.chdir(f"{path}")
    pool = multiprocessing.Pool(processes=1)
    nm_files = subprocess.check_output("ls|grep -P '^N|^Z7|^Z18|^Z22'|grep -P 'merge$'",shell=True)
    nm_files = str(nm_files,encoding = 'utf-8').split("\n")[:-1]
    bc_files = subprocess.check_output("ls|grep -P -v '^N|^Z7|^Z18|^Z22'|grep -P 'merge$'",shell=True)
    bc_files = str(bc_files,encoding = 'utf-8').split("\n")[:-1]
    for each_file in nm_files:
        pool.apply_async(combine, (each_file,'WGBS_NM_total_cov.bed',))
    for each_file in bc_files:
        pool.apply_async(combine, (each_file,'WGBS_BC_total_cov.bed',))
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
