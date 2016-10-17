from sys import argv
from random import randint

def getBase(num):
    basemap = {0:'A', 1:'C', 2:'G', 3:'T'}
    return basemap[num]

def main():
    script, filename, percent = argv
    print(filename)
    error = int(percent)
    fp = open(filename,"r")
    fw = open(filename+str(2*error)+"_error_3rdbasechanged.fasta","w+")
    cnt = 0
    r3dbase = 0
    for ii in fp.readlines():
        read = ii.strip('\n').split(' ')
        l = len(read[0])
        seq = [c for c in read[0]]
        errlen = (l*error)//100
        chunk_st = randint(0,l-(errlen+errlen+errlen))
        chunk_en = chunk_st + errlen
        for x in range(chunk_st,chunk_en):
            seq[x] = getBase(randint(0,3))
        name = ">Ecloi_"+read[1]+"_aligned"+"_"+str(cnt)+"_F_"+"0_"+str(l)+"_0\n";
        for i in range(len(seq)):
            if r3dbase>=(errlen+errlen):
                break
            if ((i+1)%3==0 and i!=0) and (i>=chunk_en and i<l):
                seq[i] = getBase(randint(0,3))
                r3dbase +=1
        fw.write(name)
        fw.write(''.join(seq)+"\n")
        cnt += 1
    fw.close()
    fp.close()
if __name__=="__main__":main()
