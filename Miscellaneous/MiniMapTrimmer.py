def main():
    #fp = open("20K_minimap_out.mini","r");
    #tt = fp.readlines()
    #read_names = {}
    #for ii in tt:
    #    kk = ii.split('\t')
    #    if kk[0] not in read_names:
    #        read_names[kk[0]] = ii
    #fp.close()
    #fp = open("20K_minimap_out_trimmed.mini","w+");
    #for k in read_names:
    #    fp.write(read_names[k])
    #fp.close()
    fp = open("50_genseq_out.mini","r");
    tp = open("50_genseq_out_trimmed.mini","w+")
    tt = fp.readlines()
    for ii in tt:
        kk = ii.split('\t')
        tp.write("{0}\t{1}\t{2}\t{3}\n".format(kk[0],kk[7],kk[8],(int(kk[3])-int(kk[2]))))
    tp.close()
    fp.close()
if __name__ == "__main__":
    main()
