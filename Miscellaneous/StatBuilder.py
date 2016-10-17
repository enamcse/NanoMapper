from sys import argv

def main():
    script, filename = argv
    print(filename)
    fp = open(filename,"r")
    #fw = open(filename+".stat","w+")
    flag = False
    tt = fp.readlines()
    ginbase = 0
    goutbase = 0
    linbase = 0
    loutbase = 0
    fcnt = 0
    prevname =""
    for ii in tt:
        kk = ii.split('\t')
        nmk = kk[0].split('_')
        st = int(nmk[2])
        en = int(nmk[3])
        mst = int(kk[1])
        men = int(kk[2])
        #print(st," ",en," ",mst," ",men)
        if flag==False:
            prevname = kk[0]
            flag = True
            if mst>en:
                loutbase = (men-mst)
                linbase = 0
            elif men<st:
                loutbase = (men - mst)
                linbase = 0
            elif st<=mst and en>=men:
                linbase = (men-mst)
                loutbase = 0
            elif st>=mst and en<=men:
                linbase = (en-st)
                loutbase = (st-mst) + (men - en)
            elif st<=mst and en<=men:
                linbase = (en-mst)
                loutbase = (men-en)
            elif st>=mst and en>=men:
                linbase = (st-men)
                loutbase = (st-mst)
            #print(linbase," + ",loutbase)
        else:
            if prevname==kk[0]:
                if mst>en:
                    loutbase += (men-mst)
                elif men<st:
                    loutbase += (men - mst)
                elif st<=mst and en>=men:
                    linbase += (men-mst)
                elif st>=mst and en<=men:
                    linbase += (en-st)
                    loutbase += (st-mst) + (men - en)
                elif st<=mst and en<=men:
                    linbase += (en-mst)
                    loutbase += (men-en)
                elif st>=mst and en>=men:
                    linbase += (st-men)
                    loutbase += (st-mst)
                #print(linbase," + ",loutbase)
            else:
                print(prevname)
                print("Mapped in : ",linbase," Mapped out : ",loutbase)
                prevname = kk[0]
                if mst>en:
                    loutbase = (men-mst)
                    linbase = 0
                elif men<st:
                    loutbase = (men - mst)
                    linbase = 0
                elif st<=mst and en>=men:
                    linbase = (men-mst)
                    loutbase = 0
                elif st>=mst and en<=men:
                    linbase = (en-st)
                    loutbase = (st-mst) + (men - en)
                elif st<=mst and en<=men:
                    linbase = (en-mst)
                    loutbase = (men-en)
                elif st>=mst and en>=men:
                    linbase = (st-men)
                    loutbase = (st-mst)
    fp.close()
if __name__=="__main__":main()
