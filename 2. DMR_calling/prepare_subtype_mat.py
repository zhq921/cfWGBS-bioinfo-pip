import time
import sys
def get_sam_lab_nm_early(file):
    print(time.asctime( time.localtime(time.time()))+" Start reading from {}".format(file))
    tnbc_samples = []
    a_samples = []
    b_her_pos_samples = []
    b_her_neg_samples = []
    her_samples = []
    with open(file) as f:
        for line in f:
            if line.startswith("sample"):
                continue
            line_list = line.strip().split("\t")
            if line_list[1] == "TNBC":
                tnbc_samples.append(line_list[0])
            elif line_list[1] == "Luminal-A":
                a_samples.append(line_list[0])
            elif line_list[1] == "Luminal-B her2+":
                b_her_pos_samples.append(line_list[0])
            elif line_list[1] == "Luminal-B her2-":
                b_her_neg_samples.append(line_list[0])
            elif line_list[1] == "HER2":
                her_samples.append(line_list[0])
            else:
                continue
    return tnbc_samples, a_samples, b_her_pos_samples, b_her_neg_samples, her_samples

def main():
    if len(sys.argv)<5:
        print("Usage: *.py sample_info out_path all_mat sam_idtolabel_out_path")
        exit(0)
    tnbc_samples, a_samples, b_her_pos_samples, b_her_neg_samples, her_samples =get_sam_lab_nm_early(sys.argv[1])
    out = open(sys.argv[2],"a")
    sam_ID2label = open(sys.argv[4],"a")
    with open(sys.argv[3],"r") as f:
        for line in f:
            if line.startswith("chrom"):
                line = line.strip().split()
                idx = [0,1,2] #first 3 columns
                label = ["chrom","start","end"]
                samid = []
                tnbc_num = 0
                a_num = 0
                b_pos_num = 0
                b_neg_num = 0
                her_num = 0
                for i in range(3,167):
                    if line[i] in tnbc_samples:
                        idx.append(i)
                        tnbc_num+=1
                        label.append("G1_"+str(tnbc_num))
                        samid.append(line[i])
                    elif line[i] in a_samples:
                        idx.append(i)
                        a_num+=1
                        label.append("G2_"+str(a_num))
                        samid.append(line[i])
                    elif line[i] in b_her_pos_samples:
                        idx.append(i)
                        b_pos_num+=1
                        label.append("G3_"+str(b_pos_num))
                        samid.append(line[i])
                    elif line[i] in b_her_neg_samples:
                        idx.append(i)
                        b_neg_num+=1
                        label.append("G4_"+str(b_neg_num))
                        samid.append(line[i])
                    elif line[i] in her_samples:
                        idx.append(i)
                        her_num+=1
                        label.append("G5_"+str(her_num))
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