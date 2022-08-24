'''
This python script converts the SyRI generated special VCF to resemble the more commonly used VCF file format
Author: Krittika

Input files: 1. duplicated coordinates extracted from SyRI generated VCF; 2. SyRI generated VCF with only the SNP variants
Output files: 1.snp_repeat_coords.txt (entries with duplicated coords in the VCF)
2.snp_modified.txt (entries with modified lines based on duplicated coordinates)
3.snp_temp.txt (entries without duplicated coordinates in the original VCF)
'''


def coords(dupfile): #function to get the coordinated which are in the second column, first is the count from uniq -c
    coords=[]
    with open(dupfile,'r') as f1:
        for line in dupfile:
            line=line.lstrip().rstrip()
            line=line.split(' ')
            coords.append(line[1])
    return coords

def repeat_coords(nucvcf,coord,repout): #prints out the line from the original VCF that have matching coordinates in the coord list
    f3 = open(repout,'w+') 
    with open(nucvcf) as f2:
        for _ in range(30):
               next(f2)
        for line in f2:
            temp=line.split('\t')
            for value in coord:
                if value==temp[1]:
                    f3.write(line)
                    break
    f3.close()
    f2.close()
    return None

#function to compare two lines at a time to see if there is allele duplication in the duplicated coordinates; if different alleles found, join by a comma
def snp_modify(repout,modout): #file created from the previous step; write the output to f4
    f4 = open(modout,'w+') #change filename
    flag=[] #holds the SNP coordinates
    with open(repout,'r') as f:
        for line1,line2 in zip(f,f):
            line1=line1.rstrip("\n")
            line2=line2.rstrip("\n")
            temp1=line1.split("\t")
            temp2=line2.split("\t")
            info1 = temp1[7].split(";")
            info2 = temp2[7].split(";")

            if temp1[1]==temp2[1] and temp1[4]!=temp2[4]:

                f4.write(str(temp1[0])+"\t"+str(temp1[1])+"\t"+str(temp1[2])+"\t"+str(temp1[3])+"\t"+str(temp1[4])+","+str(temp2[4])+"\t"+str(temp1[5])+"\t"+str(temp1[6])+"\t"+info1[1]+";"+info1[4]+";"+info2[4]+";"+info1[6]+";"+"MODIFIED"+"\t"+"GT"+"\t"+"1/1")
                f4.write('\n')

            elif temp1[1] not in flag:

                f4.write(temp1[0]+"\t"+temp1[1]+"\t"+temp1[2]+"\t"+temp1[3]+"\t"+temp1[4]+"\t"+temp1[5]+"\t"+temp1[6]+"\t"+info1[1]+";"+info1[4]+";"+info2[4]+";"+info1[6]+"\t"+"GT"+"\t"+"1/1")
                f4.write("\n")

            else:
                continue

            flag.append(temp1[1])
        f.close()
        f4.close()
    return None

#function to separate out the non-repeated SNPs from the VCF and modify the lines to add GT and line name
def snp_rearrange(nucvcf,modout,tempout):
    c =[]
    f5 = open(tempout,'w+') 
    with open(modout) as f:
        for line in f:
            temp=line.split("\t")
            c.append(temp[1]) #stores coords that are processed for repeats

    coords_seen=set()


    with open(nucvcf) as f2: 
        for _ in range(31):
            next(f2)

        for line2 in f2:
            line2=line2.rstrip("\n")
            temp1=line2.split("\t")
            info = temp1[7].split(";")
            if temp1[1] not in coords_seen:
                if temp1[1] not in c:
                    f5.write(temp1[0]+"\t"+temp1[1]+"\t"+temp1[2]+"\t"+temp1[3]+"\t"+temp1[4]+"\t"+temp1[5]+"\t"+temp1[6]+"\t"+info[1]+";"+info[4]+";"+info[6]+"\t"+"GT"+"\t"+"1/1")
                    f5.write("\n")
                    coords_seen.add(temp1[1])
        f5.close()
        nucvcf.close()
        modout.close()

    return None


def main():
    dupfile = 'grassl_snp_coords_dup.txt'
    nucvcf = 'grassl_nucmer_snp.vcf'
    repout = 'snp_repeat_coords_june24.txt'
    modout = 'snp_modified_june24.txt'
    tempout = 'snp_temp.txt'

    coord = coords(dupfile)
    repeat_coords(nucvcf,coord,repout)
    snp_modify(repout,modout)
    snp_rearrange(nucvcf,modout,tempout)

if __name__ == "__main__":
    main()
