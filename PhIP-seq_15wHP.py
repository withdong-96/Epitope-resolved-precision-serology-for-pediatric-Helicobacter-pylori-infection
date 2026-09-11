# -*- coding: utf-8 -*-
#! /home/ubuntu/mambaforge/bin/python
#! /home/ubuntu/mambaforge/bin/R


### Import Library
import os
import time
import numpy as np
import pandas as pd
from copy import deepcopy
from optparse import OptionParser


### Function
def Barcode_Prepare(barcode):
    data = pd.read_csv(barcode, header = 0, index_col = 0, sep = "\t")
    barcode = np.zeros((data.shape[0]*data.shape[1],3), dtype=int)
    barcode = pd.DataFrame(barcode, columns = ["Abcode","S","N"])
    k = 0
    for i in range(data.shape[1]):
        for j in range(data.shape[0]):
            barcode.iloc[k,0] = data.iloc[j,i]
            barcode.iloc[k,1] = data.index.to_list()[j]
            barcode.iloc[k,2] = data.columns.to_list()[i]
            k = k + 1
    barcode.to_csv("Barcode_info.txt", sep = "\t", index = None)
    os.system("grep -v ^0 Barcode_info.txt |awk 'BEGIN{FS=\"\t\"}NR>1{ print $1\"\t\"$2\"\"$3;}' | awk 'BEGIN{FS=\"\t\";print \"ID\tBarcode\tClass\"}{if($0 ~ /^Blank/ || $0 ~ /^BLANK/ || $0 ~ /^control/ || $0 ~ /^blank/ || $0 ~ /^Control/ ){print $0\"\tblank\"} else{print $0\"\ttreat\"}}'>>  Barcode.txt")
    os.system("awk 'BEGIN{FS=\"\t\"}NR>1{print $2}' Barcode.txt  >> Index.txt")


def Frequency():
    path=os.getcwd() + "/2_Sample/"
    background = pd.read_csv("peptide", header = None, sep = "\t")
    data = deepcopy(background)
    data.columns = ["Peptide"]
    data = data.sort_values(by = "Peptide")
    data.index = np.arange(data.shape[0])
    
    files = os.listdir(path)
    for i in files:
        print(i)
        name = i[0:16]
        sample = path + i 
        table = pd.read_csv(sample, header = None, sep = "\t")
        table.columns = ["P", name]
        A = set(background.iloc[:,0].values.tolist())
        B = set(table.iloc[:,0].values.tolist())
        diff = list(A.difference(B))
        if len(diff) > 0:
            diff = pd.DataFrame(diff)
            diff[name] = 0
            diff.columns = ["P", name]
            temp = pd.concat([table,diff], axis = 0)
            temp = temp.sort_values(by = "P")
            temp.index = np.arange(data.shape[0])
            data = pd.concat([data,temp], axis=1)
            data = data.drop(columns = ['P'], axis = 1)
        if len(diff) ==0:
            data = table
            data.columns = ["Peptide", name]
    data.to_csv("Raw_Table.txt", sep = "\t", index=False)


def Replace(Raw):
    path=os.getcwd()
    barcode = path + "/0_Raw/Barcode.txt"
    code = pd.read_csv(barcode, header = 0, sep = "\t")
    codedict = {code.Barcode[i]:code.ID[i] for i in range(code.shape[0])}
    data = pd.read_csv(Raw, header = 0, sep = "\t")
    names = data.columns.tolist()
    renames = [ codedict[i] if i in list(codedict.keys()) else i for i in names ]
    data.columns = renames
    data.to_csv("Count_Table.txt", sep = "\t", index=False)



def Merge(Ref, Count):

    Input = pd.read_csv(Ref,   header=0, index_col=0,  sep="\t")
    PhIP  = pd.read_csv(Count, header=0, index_col=0,  sep="\t")
    merge_list = set(Input.index.tolist() + PhIP.index.tolist())
    merge_list = list(merge_list)
    merge_list.sort()

    Input_Lack = list(set(merge_list) - set(Input.index.tolist()))
    PhIP_Lack  = list(set(merge_list) - set(PhIP.index.tolist()))

    Input_Lack_Data = pd.DataFrame(Input_Lack)
    PhIP_Lack_Data = pd.DataFrame(PhIP_Lack)

    for i in Input.columns:
        Input_Lack_Data[i] = 0
    for i in PhIP.columns:
        PhIP_Lack_Data[i]  = 0

    Input_Lack_Data.index = Input_Lack_Data.iloc[:,0]
    Input_Lack_Data = Input_Lack_Data.drop(Input_Lack_Data.columns[0], axis=1)
    PhIP_Lack_Data.index  = PhIP_Lack_Data.iloc[:,0]
    PhIP_Lack_Data = PhIP_Lack_Data.drop(PhIP_Lack_Data.columns[0], axis=1)

    Input_Full_Data = pd.concat([Input,Input_Lack_Data],axis=0)
    Input_Full_Data = Input_Full_Data.loc[merge_list,:]
    PhIP_Full_Data = pd.concat([PhIP,PhIP_Lack_Data],axis=0)
    PhIP_Full_Data = PhIP_Full_Data.loc[merge_list,:]

    Result = pd.concat([Input_Full_Data, PhIP_Full_Data], axis=1)
    Result = Result.loc[Result.Input > 0,:]
    Result.to_csv("Count.txt.gz", sep="\t", index_label="Peptide", compression="gzip")


def main():
    usage = "usage: python %prog [options] -1 [Files with #1 mates] -2 [Files with #1 mates] -b [barcode file] [OPTIONAL_FLAGS]"
    parser = OptionParser(usage = usage)
    
    #required flags
    parser.add_option("-1","--read1", dest="read1",nargs = 1, default=None, help = "Enter a .gz file with #1 mates")
    parser.add_option("-2","--read2", dest="read2",nargs = 1, default=None, help = "Enter a .gz file with #2 mates")
    parser.add_option("-b","--barcode", dest="barcode",nargs = 1, default=None, help = "Enter a .tsv file with standard barcode information")
    
    #optional flags
    #parser.add_option("-p","--threads", dest="threads",nargs = 1, default=8, help = "Enter a number of CPU cores to use")
    
    #RETRIEVING FLAGS
    (options,args) = parser.parse_args()

    if not options.read1 or not options.read2:
        print('Please input sequencing reads file in fastq format.')
        parser.print_help()
        exit()

    if not options.barcode:
        print('Please input barcode in tsv format.')
        parser.print_help()
        exit()

    

    ### Input Data Region
    R1=options.read1
    R2=options.read2
    r1="Temp_R1.fq.gz"
    r2="Temp_R2.fq.gz"
    barcode=options.barcode
    
    threads = 6
    print("Read1: " + options.read1)
    print("Read2: " + options.read2)
    print("Number of threads: " + str(threads))
    
    ### Make Temp Dirs
    if not os.path.exists("0_Raw/"):
        os.system('mkdir 0_Raw/')
    if not os.path.exists("1_Extract/"):
        os.system('mkdir 1_Extract/')
    if not os.path.exists("2_Sample"):
        os.system('mkdir 2_Sample/')
    if not os.path.exists("3_Table"):
        os.system('mkdir 3_Table')

    ### step-0 Prepare Barcode Information
    #os.system('cp ' + barcode + '   ' + BARCODE)
    Barcode_Prepare(barcode)
    os.system('mv Barcode_info.txt 0_Raw/')
    os.system('mv Barcode.txt 0_Raw/')
    os.system('mv Index.txt 0_Raw/')
    #os.system('mv ' + barcode + " 0_Raw/")

    ### Step-1 Clean and Extract
    os.system('cp ' + R1 + '  ' + r1)
    os.system('cp ' + R2 + '  ' + r2)
    os.system('fastp -w ' + str(threads) + ' -i ' + r1 + ' -I ' + r2 + ' -n 0 -e 20 -U --umi_loc=per_read --umi_len=8 --overlap_diff_limit 0 -m --merged_out merge.fq.gz --disable_trim_poly_g  --length_required 150')
    os.system('mv fastp.html 1_Extract/')
    os.system('mv fastp.json 1_Extract/')
    os.system('rm ' + r1)
    os.system('rm ' + r2)

    #os.system('cutadapt  -j ' + str(threads) + ' ' + ' -g ^ATGCTCGGGGATCCGAATTCT...TGAAAGCTTGCGGCCGCACTCGAGTAAC$ -o clean.fq.gz -m 168 -M 168 merge.fq.gz')
    os.system('cutadapt  -j ' + str(threads) + ' ' + ' -g ^ATGCTCGGGGATCCGAATTCCGCTGCG...GATTACAAGGACGACGACGACAAGTAAAAGCTTGCGGCCGCACTCGAGTAAC$ -o clean.fq.gz -m 168 -M 168 merge.fq.gz')
    os.system('seqkit translate -j ' + str(threads) + ' clean.fq.gz -o clean.aa.gz')
    os.system("zcat clean.aa.gz | awk 'NR%2==1{split($1,a," + '"' + ":" + '"' + ");b=a[length(a)];sub(" + '"' + "_" + '"' + "," + '"' + '"' + ",b);getline l2;if(b !~ /N/) {print b " + '"' + "\t" + '"'  + "0" + '"' +  "\t" + '"' +"l2}}' > clean.seq")
    os.system("echo -n Raw Read Pairs:' '  >> 1_Extract/stat.txt; zcat "+R1+" | wc -l |awk '{print $1/4}' >> 1_Extract/stat.txt")
    os.system("echo -n Clean Read Pairs:' '  >> 1_Extract/stat.txt; wc -l clean.seq|awk '{print $1}' >> 1_Extract/stat.txt")
    os.system("awk '{print $1}' clean.seq |egrep -o " + '"' + "\\b[[:alpha:]]+\\b" + '"' +   "| awk 'BEGIN{print " + '"' + "Word\tCount" + '"' +"}{ count[$0]++ } END{ for(ind in count) { print ind" + '"' + "\t" + '"' + "count[ind]} }'|sort -k 2,2nr  >> 1_Extract/Frequency.txt")
    os.system("cp 1_Extract/Frequency.txt ./")
    os.system('rm clean.aa.gz ')
    os.system('rm clean.fq.gz ')
    os.system('rm merge.fq.gz ')
    
    
    os.system('''awk '{if(NR!=1){print $2}}' 0_Raw/Barcode.txt > Raw.text''')
    os.system('grep -w -F -f 0_Raw/Index.txt clean.seq > peptide_raw.seq')
    os.system("echo -n With Stop Condons:' '     >> 1_Extract/stat.txt; cat peptide_raw.seq|wc -l|awk '{print $1}' >> 1_Extract/stat.txt")
    os.system("grep -v '*' peptide_raw.seq > Peptide.seq")
    os.system('rm peptide_raw.seq')
    os.system("echo -n Without Stop Condons:' '  >> 1_Extract/stat.txt; cat Peptide.seq|wc -l|awk '{print $1}' >> 1_Extract/stat.txt")
    os.system('rm clean.seq')

    ## Step-2 Make Table
    os.system("awk '{print $3 >> $1" + '"' + '.tmp.TXT' + '"' + " }'  Peptide.seq")
    #os.system("t=0;for i in $(cat  0_Raw/Index.txt);do egrep -o " + '"' + "\\b[[:alpha:]]+\\b" + '"' + " $i'.tmp.TXT' | awk '{ count[$0]++ } END{ for(ind in count) { print ind" + '"' + "\t" + '"' + "count[ind] } }'|sort -r -n -k2 >> 2_Sample/$i'.txt' & t=$(($t+1)) ; if [[ $t -gt " +  str(threads)  + " ]];then wait ;t=0;fi; done" )
    os.system('for i in $(cat 0_Raw/Index.txt); do egrep -o "\\b[[:alpha:]]+\\b" $i\'.tmp.TXT\' | awk \'{ count[$0]++ } END{ for(ind in count) { print ind"\\t"count[ind] } }\' | sort -r -n -k2 >> 2_Sample/$i\'.txt\'; done')
    os.system("awk '{print $3}' Peptide.seq|sort -u > peptide")
    os.system('rm Peptide.seq')
    Frequency()
    Replace("Raw_Table.txt")
    os.system('rm Raw_Table.txt')
    os.system('rm *.tmp.TXT')
    os.system('rm 2_Sample/* ')
    Merge("HP15w_Scan.txt", "Count_Table.txt")
    os.system('rm Count_Table.txt')
    
if __name__ == "__main__":
    main()
