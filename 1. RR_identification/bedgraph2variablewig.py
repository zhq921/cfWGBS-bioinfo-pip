import sys
import time

def bg2variablewig(file):
    start_time = time.time()
    print(f"Start parsing {file}")
    f = open(file, mode = "r")
    out = open(file+".wig",mode="a")
    line = f.readline()
    chr = "1"
    print("variableStep chrom=chr{0} span=1".format(chr),file = out)
    while line:
        linels = line.strip("\n")
        linels = linels.split(sep = "\t")
        if chr == linels[0]:
            for i in range(eval(linels[1])+1,eval(linels[2])+1):
                print(i,linels[3],file=out)
        else:
            chr = linels[0]
            print("variableStep chrom=chr{0} span=1".format(chr),file = out)
            for i in range(eval(linels[1])+1,eval(linels[2])+1):
                print(i,linels[3],file=out)
        line = f.readline()
    f.close()
    out.close()
    end_time = time.time()
    time_used = (end_time - start_time) / 60
    print("Job done!!! Time used: {0: 0.2f} min".format(time_used))

def main():
    if len(sys.argv)<2:
        print("Usage: python bedgragh2variablewig.py example.bedgraph\nwith span equal to 1")
        exit(0)
    file = sys.argv[1]
    bg2variablewig(file)

if __name__ =="__main__":
    main()
