import sys


def getTEcoords(TEfile,geneloc):
    with open(TEfile) as tf:
        tf.readline()
        for TEline in tf:
            TEline=TEline.rstrip().split("\t")
           #looping through the gene annotation file to find overlaps for TE coordinates
            with open(geneloc) as g:
                for geneline in g:
                    if not geneline.startswith("#") and not geneline.startswith("tig"):
                        geneline=geneline.rstrip().split("\t")
                        chrom=geneline[0].split('_')[0].strip("Chr0")
                        if geneline[0]=='Chr10_RagTag':
                            chrom=geneline[0].split('_').strip("Chr")
                        #check if the TE boundaries are within the gene boundary
                        if TEline[0] == chrom and int(TEline[6]) >= int(geneline[3]) and int(TEline[7]) <= int(geneline[4]) and geneline[2] in ("gene"):
                            if geneline[8].find("ID=") != -1:
                                #storing the whole gene ids since they have differing patterns
                                geneID = geneline[8].split(";")[0].split("=")[1]
                                #geneID = geneID.split(".")[0] + "." + geneID.split(".")[1]
                            else:
                                geneID = "-"

                        
                            newline=(TEline[0],TEline[1],TEline[2],TEline[5],geneID,geneline[3],geneline[4])
                            
                            print("\t".join(newline))



def main():
    
    TEfile=sys.argv[1]
    geneloc=sys.argv[2]
    header=["Chrom","TEID","TEpair","TEage","geneID","geneStart","geneEnd"]
    print("\t".join(header))
    getTEcoords(TEfile,geneloc)

if __name__ == '__main__':
    main()
