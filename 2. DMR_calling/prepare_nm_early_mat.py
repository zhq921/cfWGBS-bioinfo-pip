import time
import sys
def get_sam_lab_nm_early(file):
    print(time.asctime( time.localtime(time.time()))+" Start reading from {}".format(file))
    normal_samples = []
    early_samples = []
    with open(file) as f:
        for line in f:
            if line.startswith("sample"):
                continue
            line_list = line.strip().split("\t")
            if line_list[3] == "0":
                normal_samples.append(line_list[0])
            elif line_list[4] == "0":
                early_samples.append(line_list[0])
            else:
                continue
    normal_samples = list(map(lambda k: k.strip(), normal_samples))
    early_samples = list(map(lambda k: k.strip(), early_samples))
    return normal_samples, early_samples

def main():
    if len(sys.argv)<5:
        print("Usage: *.py sample_info out_path all_mat sam_idtolabel_out_path")
        exit(0)
    normal_samples,early_samples =get_sam_lab_nm_early(sys.argv[1])
    out = open(sys.argv[2],"a")
    sam_ID2label = open(sys.argv[4],"a")
    with open(sys.argv[3],"r") as f:
        for line in f:
            if line.startswith("chrom"):
                line = line.strip().split()
                idx = [0,1,2] #first 3 columns
                label = ["chrom","start","end"]
                samid = []
                normal_num = 0
                early_num = 0
                for i in range(3,167):
                    if line[i] in normal_samples:
                        idx.append(i)
                        normal_num+=1
                        label.append("G2_"+str(normal_num))
                        samid.append(line[i])
                    elif line[i] in early_samples:
                        idx.append(i)
                        early_num+=1
                        label.append("G1_"+str(early_num))
                        samid.append(line[i])
                    else:
                        continue
                print("\t".join(label),sep="\t",file=out)
                for everysam in range(0,len(samid)):
                    print(samid[everysam],label[everysam+3],sep="\t",file=sam_ID2label)
            else:
                line = line.strip().split()
                line_apart=[]
                for ele in idx:
                    line_apart.append(line[ele])
                print("\t".join(line_apart),sep="\t",file=out)
    out.close()
    sam_ID2label.close()
    
if __name__=="__main__":
    main()