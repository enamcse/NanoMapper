from sys import argv

def main():
    script, filename= argv
    print(filename)
    fp = open(filename,"r")
    fw = open(filename+"name_rs_en_as_an.txt","w+")
    cnt = 0
    for ii in fp.readlines():
        read = ii.strip('\n').split('\t')
        if (len(read)<11): continue
        name = read[0].split('_')
        ref_st = int(read[3])
        ref_en = ref_st + len(read[9])
        act_st = int(name[1])
        act_en = act_st + int(name[6])
        fw.write(read[0]+"\t"+str(ref_st)+"\t"+str(ref_en)+"\t"+str(act_st)+"\t"+str(act_en)+"\n")
    fw.close()
    fp.close()
if __name__=="__main__":main()
