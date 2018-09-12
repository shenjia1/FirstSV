

import sys
import os
import re
import time
import argparse
import json
from collections import OrderedDict

prog_name = "FirstSV V 1.0"

TEMP_STR = re.compile(r"\d+([SMHIDP])")
TEMP_INT = re.compile(r"(\d+)[SMHIDP]")

base_reverse = {
                 "A" : "T",
                 "C" : "G",
                 "T" : "A",
                 "G" : "C",
               }

sv = OrderedDict()
sv1 = OrderedDict()
sv2 = OrderedDict()
srregion = {}
refgene = {}
drunit = {}
config = OrderedDict()
gene_list = []
bed_gene = []

#readin config_file to the dictionary
def parser_config(config_dir):
    if not os.path.exists(config_dir):
        print "%s is not exists! EXIT!" % config_dir
        exit(0)
    with open(config_dir,"r") as config_file:
        for eachline in config_file.readlines():
            if eachline[0] == "#":
                continue
            if eachline == "" or eachline == "\n":
                continue
            line = eachline.strip().split("=")
            config[line[0]] = line[1].lstrip()

#deal with the cigar, readin pair cigar pattern, output "NC/CN" and the number of bp used to correct the position of breakpoint
#"23M45S" --> [(23,M),(45,S)] -(by the fuction of del_cigar_tmp)-> ["NC",23]
def del_cigar(cigar1,cigar2):
    #store_user = []
    
    tempstr = re.findall(TEMP_STR,cigar1)
    tempint = re.findall(TEMP_INT,cigar1)
    cigar1_temp = [(tempint[i],tempstr[i]) for i in range(0,len(tempstr))]
    tempstr = re.findall(TEMP_STR,cigar2)
    tempint = re.findall(TEMP_INT,cigar2)
    cigar2_temp = [(tempint[i],tempstr[i]) for i in range(0,len(tempstr))]
    cigar1_order = del_cigar_tmp(cigar1_temp)
    cigar2_order = del_cigar_tmp(cigar2_temp)
    if not cigar1_order:
        return False
    if not cigar2_order:
        return False
    store_user = [cigar1_order,cigar2_order]
    return store_user

def del_cigar_2(cigar1,cigar2):
    tempstr = re.findall(TEMP_STR,cigar1)
    tempint = re.findall(TEMP_INT,cigar1)
    cigar1_temp = [(tempint[i],tempstr[i]) for i in range(0,len(tempstr))]
    tempstr = re.findall(TEMP_STR,cigar2)
    tempint = re.findall(TEMP_INT,cigar2)
    cigar2_temp = [(tempint[i],tempstr[i]) for i in range(0,len(tempstr))]
    cigar1_order = del_cigar_tmp(cigar1_temp)
    cigar2_order = del_cigar_tmp(cigar2_temp)
    cigar1_dis = del_cigar_tmp2(cigar1_temp)

    if not cigar1_order:
        return False
    if not cigar2_order:
        return False
    if not cigar1_dis:
        return False

    store_user = [cigar1_order,cigar2_order,cigar1_dis]
    return store_user

def del_cigar_tmp(cigar_temp):
    if cigar_temp == []:
        return False
    if (cigar_temp[0][1] == 'S' or cigar_temp[0][1] == 'H') and (cigar_temp[-1][1] == 'S' or cigar_temp[-1][1] == 'H') :
        if int(cigar_temp[0][0]) < 10 and int(cigar_temp[-1][0]) >= 10:
            return ["NC",int(cigar_temp[1][0])]
        elif int(cigar_temp[-1][0]) < 10 and int(cigar_temp[0][0]) >= 10:
            return ["CN",0]
        else:
            return False
    else:
        if cigar_temp[0][1] == 'S' or cigar_temp[0][1] == 'H': 
            return ["CN",0]
        elif cigar_temp[-1][1] == 'S' or cigar_temp[-1][1] == 'H':
            return ["NC",int(cigar_temp[0][0])]
        else:
            return False

#calculate the breakpoint 
def del_cigar_tmp2(cigar_temp):
    if cigar_temp == []:
        return False
    if int(cigar_temp[0][0]) < 10:
        return int(cigar_temp[1][0])
    else:
        return int(cigar_temp[0][0])

#calculate the depth at the breakpoint by "samtools mpileup"
def cal_depth(chrN,pos,bam):
    samtools = config["samtools "]
    cmd = "%s mpileup -Ar %s:%s-%s %s -Q 20 -q 20" % (samtools,chrN,pos,pos,bam)
    fp = os.popen(cmd,"r")
    rdepth = 0
    for eachline in fp.readlines():
        line = eachline.strip().split()
        rdepth = int(line[3])
    fp.close()
    return rdepth

#readin simplerepeatbed file to the dictionary {srregion}
def load_bed(bedfile):
    with open(bedfile,"r") as bfile:
        for eachline in bfile.readlines():
            if eachline == "" or eachline[0] == "#":
                continue
            line = eachline.strip().split()
            chrN = line[0]
            spos = int(line[1])
            epos = int(line[2])
            if srregion.has_key(chrN):
                srregion[chrN].extend(range(spos,epos))
            else:
                srregion[chrN] = range(spos,epos)

def load(genefile):
    with open(genefile) as json_file:
        data = json.load(json_file)
        return data

def dict2list(data):

    for each in data:
        gene_list.append(each)
        chrN = data[each]["gene_range"].split(":")[0]
        start = int(data[each]["gene_range"].split(":")[1].split("-")[0])
        end = int(data[each]["gene_range"].split(":")[1].split("-")[1])
        bed_gene.append([chrN,start,end])
    return gene_list,bed_gene


#use refseq to annoate the breakpoint 

def is_inbed(c_str,start,end):
    return lambda x:((x[0]==c_str)and (x[1]<=start) and (x[2]>=end))

def anno_refgene(chrN,start,end):
    run = filter(is_inbed(chrN,start,end),bed_gene)
    if run == []:
        return "?"
    else:
        return gene_list[bed_gene.index(run[0])]



def seq_filter(samtools,chr1,pos1,ref):
    pos1_r = int(pos1) - 150 if int(pos1)-150>0 else 0
    pos1_l = int(pos1) + 150 if int(pos1)+150>0 else int(pos1) + 150
    sflag = "%s:%d-%d" % (chr1,pos1_r,pos1_l)
    files = os.popen("%s faidx %s %s" % (samtools,ref,sflag))
    segment ="".join([i.strip() for i in files.readlines()[1:]])
    files.close()
    repeat_flag = True
    for i in xrange(len(segment) - 150):
        if get_seq_complexity(segment[i:i + 150],3) < 0.3:
            repeat_flag = False
    if repeat_flag and cal_gc(segment) > 0.2 and cal_gc(segment) < 0.8:
        return False
    else:
        return True
   
def seq_cal(samtools,chr1,pos1,ref):
    pos1_r = int(pos1) - 25 if int(pos1)-25>0 else 0
    pos1_l = int(pos1) + 25 if int(pos1)+25>0 else int(pos1) + 25
    rflag = "%s:%d-%d" % (chr1,pos1_r,pos1_l)
    files = os.popen("%s faidx %s %s" % (samtools,ref,rflag))
    rsegment ="".join([i.strip() for i in files.readlines()[1:]])
    files.close()
    return cal_gc(rsegment),get_seq_complexity(rsegment)

def call_sv(sv_filename,sam,bp_threshold,sp_threshold,mp_threshold,samtools,ref,hd,model):
    if  not hd:
        import numpy as np
        from sklearn.neural_network import MLPClassifier
        from sklearn.externals import joblib
        clf = MLPClassifier()
        clf = joblib.load(model)
        #print clf.predict([[0.5,0.5,0.5,0.5]])
        ##########################################
        ############        ?        #############
        ##########################################
    sam_file = open(sam,"r")
    sv_file = open(sv_filename,"w")
    for eachline in sam_file.readlines():
        line = eachline.strip().split("\t")
        if "SA" not in line[11]: #useful screads
            continue
        if int(line[4]) < mp_threshold: #MAP Score<20
            continue
        chr1 = line[2]
        pos1 = line[3]
        pattern = line[11].split(":")[2].split(",")
        chr2 = pattern[0]
        pos2 = pattern[1]
        maps = pattern[4]
        cigar1 = line[5]
        cigar2 = pattern[3]
        store_order = del_cigar(cigar1,cigar2)
        if not store_order: 
            continue
        if int(maps) < mp_threshold:
            continue
        pos1 = str(int(pos1) + store_order[0][1])  #correct the breakpoint
        pos2 = str(int(pos2) + store_order[1][1])
        strand1 = ">" if store_order[0][1] > 0 else "<"
        strand2 = ">" if store_order[1][1] > 0 else "<"
        
        if chr1 == chr2 and abs(int(pos2)-int(pos1)) < bp_threshold: #discard the distance of two breakpoints < 2Kbp
            continue
        if sv.has_key((chr2,pos2,chr1,pos1)):
            continue
        if sv.has_key((chr1,pos1,chr2,pos2)):
            sv[(chr1,pos1,chr2,pos2)][0] = sv[(chr1,pos1,chr2,pos2)][0] + 1
        else:
            sv[(chr1,pos1,chr2,pos2)] = [0,strand1,strand2,store_order[0][0],store_order[1][0]]
    for key in sv.keys():
        chr1_c = key[0]
        pos1_c = key[1][:-3] if len(key[1]) > 3 else "-"
        chr2_c = key[2]
        pos2_c = key[3][:-3] if len(key[3]) > 3 else "-"
        if (drunit.has_key((chr1_c,pos1_c,chr2_c,pos2_c))):
            sv[key][0] = sv[key][0] + drunit[(chr1_c,pos1_c,chr2_c,pos2_c)]/2
        elif (drunit.has_key((chr2_c,pos2_c,chr1_c,pos1_c))):
            sv[key][0] = sv[key][0] + drunit[(chr2_c,pos2_c,chr1_c,pos1_c)]/2
        if sv[key][0] < sp_threshold:   #useful softclip reads supports >= 5
            continue
        if key[0] == "chrM" or key[2] == "chrM":  #exclude the chrM
            continue
        if int(key[1]) in srregion[key[0]] or int(key[3]) in srregion[key[2]]:  #only one breakpoint in simple repeat region will be discarded
            continue
        #gc_filter and repeat_segment_filter 2017/09/28

        if not hd:
            lgc,lcom = seq_cal(samtools,key[0],key[1],ref)
            rgc,rcom = seq_cal(samtools,key[2],key[3],ref)
            if clf.predict([[lgc,rgc,lcom,rcom]])[0] == 0:
                continue
        else:
            if seq_filter(samtools,key[0],key[1],ref) or seq_filter(samtools,key[2],key[3],ref):
                continue
        
        if sv1.has_key((key[0],key[1])):
            sv1[(key[0],key[1])].append((key[2],key[3],sv[key][0],sv[key][1],sv[key][2],sv[key][3],sv[key][4]))
        else:
            sv1[(key[0],key[1])] = [(key[2],key[3],sv[key][0],sv[key][1],sv[key][2],sv[key][3],sv[key][4])]
    sv_file.write("#gene1chr\tgene1pos\tgene1dir\tgene1order\tgene2chr\tgene2pos\tgene2dir\tgene2order\tbpdistance\tbpsupports\tother\n")

    #deal with the (1 to more),we consider that if one breakpoint is the same and the other is not the same (also not in the region), then we take them to the false positive.

    for each in sv1:
        if len(sv1[each]) > 2:     #1 to more case
            for one in sv1[each]:
                if one[-1] == one[2]: #merge supports equal to supports, the breakpoint will be discarded.(here maybe a bug)
                    continue
                else:
                    distance = 99999
                    if each[0] == one[0]:
                        distance = int(one[1]) - int(each[1])
                    sv_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (each[0],each[1],one[3],one[5],one[0],one[1],one[4],one[6],distance,str(one[2])))
            continue
        distance = 99999
        if each[0] == sv1[each][0][0]:
            distance = int(sv1[each][0][1]) - int(each[1])
        sv_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (each[0],each[1],sv1[each][0][3],sv1[each][0][5],sv1[each][0][0],sv1[each][0][1],sv1[each][0][4],sv1[each][0][6],distance,str(sv1[each][0][2])))
    sam_file.close()
    sv_file.close()
    sv.clear()   
    sv1.clear()


#readin the sv file to dictionary{sv}
def load_sv(sv_file):
    with open(sv_file,"r") as fp:
        for eachline in fp.readlines():
            if eachline[0] == "#" or eachline == "":
                continue
            line = eachline.strip().split("\t")
            sv2[eachline] = [(line[0],int(line[1])),(line[4],int(line[5]))]

#calculate the gc
def cal_gc(seq):
    _c = seq.count("c") + seq.count("C")
    _g = seq.count("g") + seq.count("G")
    return float((_c+_g))/float(len(seq))

#calculate the complexity of sequence
def get_seq_complexity(seq,kmer):
    info = {}
    for i in range(kmer,len(seq)):
        info[seq.upper()[i-kmer:i]] = ""
    if kmer == 3:
        return float(len(info))/48
    else:
        return float(len(info))/float(len(seq)-2)

#the result recalled will be check in or not in the novel sv results $chr2;202146659;+;NC;chr2;202149443;+;CN;2784;552$
def find_sv(chr1,pos1,chr2,pos2,seq,dis,fasta,del_c,insert_t): #modified bug repeat results by shenjia 20180103
    fasta_f = open(fasta,"a")
    for key,value in sv2.items():
        print value
        hey = key.strip().split()
        if chr1==value[0][0] and chr2==value[1][0]:
            if abs(int(pos1)-value[0][1]) < 10 and abs(int(pos2)-value[1][1]) < 10:
                if hey[6] == ">":
                    pos2 = str(int(pos2)-del_c)
                else:
                    pos2 = str(int(pos2)+del_c)
                fasta_f.write(">"+hey[0]+";"+pos1+";"+hey[2]+";"+hey[3]+";"+hey[4]+";"+pos2+";"+hey[6]+";"+hey[7]+";"+hey[8]+";"+hey[9]+";"+str(dis)+";"+insert_t+"\n"+seq+"\n")
                break
        if chr1==value[1][0] and chr2==value[0][0]:
            if abs(int(pos2)-value[0][1]) < 10 and abs(int(pos1)-value[1][1]) < 10:
                if hey[6] == ">":
                    pos2 = str(int(pos2)-del_c)
                else:
                    pos2 = str(int(pos2)+del_c)
                fasta_f.write(">"+hey[0]+";"+pos2+";"+hey[2]+";"+hey[3]+";"+hey[4]+";"+pos1+";"+hey[6]+";"+hey[7]+";"+hey[8]+";"+hey[9]+";"+str(dis)+";"+insert_t+"\n"+seq+"\n")
                break
    fasta_f.close()

def deal_with_sv(line,bam):
    [chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert] = line.strip().split("\t")
    print line
    gene_anno1 = anno_refgene(chr1,int(pos1),int(pos1))
    gene_anno2 = anno_refgene(chr2,int(pos2),int(pos2))
    depth = max(cal_depth(chr1,int(pos1),bam),cal_depth(chr2,int(pos2),bam))
    return gene_anno1,gene_anno2,depth

# moudle
 
def add_insert(flag1,cigar1,flag2,cigar2,seq):
    if flag1 == "+" and flag2 == "+":
        first = seq
        second = "".join([int(i[0])*i[1] for i in treat_cigar(cigar1)])
        three = "".join([int(i[0])*i[1] for i in treat_cigar(cigar2)])
    elif flag1 == "+" and flag2 == "-":
        first = seq
        second = "".join([int(i[0])*i[1] for i in treat_cigar(cigar1)])
        three = "".join([int(i[0])*i[1] for i in cigar_reverse(treat_cigar(cigar2))])
    elif flag1 == "-" and flag2 == "+":
        first = read_reverse(seq)
        second = "".join([int(i[0])*i[1] for i in cigar_reverse(treat_cigar(cigar1))])
        three = "".join([int(i[0])*i[1] for i in treat_cigar(cigar2)])
    elif flag1 == "-" and flag2 == "-":
        first = read_reverse(seq)
        second = "".join([int(i[0])*i[1] for i in cigar_reverse(treat_cigar(cigar1))])
        three = "".join([int(i[0])*i[1] for i in cigar_reverse(treat_cigar(cigar2))])
    insert = []
    i_range = []
    sdf = ""
    df = ""    

    for i in xrange(len(first)):
        if second[i]==three[i]:
            sdf = sdf + first[i]
            df = df + second[i]
            i_range.append(i)
        elif sdf != "" and second[i]!=three[i] and i_range[-1] == i-1:
            if (i_range[0]!=0 and i_range[-1]!=len(first)-1):
                insert.append([i_range[0],i_range[-1],sdf,df])
    return insert




def treat_cigar(cigar):
    tempstr = re.findall(TEMP_STR,cigar)
    tempint = re.findall(TEMP_INT,cigar)
    cigar_temp = [(tempint[i],tempstr[i]) for i in range(0,len(tempstr))]
    return cigar_temp


def read_reverse(seq):
    temp = list(seq)
    temp.reverse()
    return "".join([base_reverse[one] for one in temp])

def cigar_reverse(cigar):
    cigar.reverse()
    return cigar



# add_insert end

 
def extract_fasta(sam,fasta,bam,sv_file):
    load_sv(sv_file)
    sam_file = open(sam,"r")
    for eachline in sam_file.readlines():
        line = eachline.strip().split("\t")
        if len(line) < 16:
            continue
        chr1 = line[2]
        pos1 = line[3]
        flag1 = "-" if int(line[1]) & 16 else "+"
        pattern = line[15].split(":")[2].split(",")
        chr2 = pattern[0]
        pos2 = pattern[1]
        flag2 = pattern[2]
        cigar1 = line[5]
        cigar2 = pattern[3]
        seq = line[9]
        store_order = del_cigar_2(cigar1,cigar2)
        if not store_order:
            continue
        if cal_gc(seq) > 0.8 or cal_gc(seq) < 0.2:
            continue
        if len(seq) < 152:
            if get_seq_complexity(seq,3) < 0.3:
                continue
        else:
            for i in xrange(len(seq) - 150):
                if get_seq_complexity(seq[i:i + 150],3) < 0.3:
                    continue
        pos1 = str(int(pos1) + store_order[0][1])
        pos2 = str(int(pos2) + store_order[1][1])
        # deal with insert/deletion around breakpoints
        if ("H" in cigar1 or "H" in cigar2): continue
        ins_r = add_insert(flag1,cigar1,flag2,cigar2,seq)
        del_l = 0
        print ins_r
        if len(ins_r) > 1:
            insert = "many_blocks"
        elif ins_r == []:
            insert = "-"
        elif "M" in ins_r[0][-1] and "S" in ins_r[0][-1]:
            insert = "many_blocks"
        elif "M" in ins_r[0][-1]:
            insert = "-"
            del_l = ins_r[0][1] - ins_r[0][0] +1         # store the length of the repeat pattern 
        elif "S" in ins_r[0][-1]:
            insert = ins_r[0][2]
        else:
            insert = "many_blocks"

        print insert
        find_sv(chr1,pos1,chr2,pos2,seq,store_order[2],fasta,del_l,insert);
    sam_file.close() 

def extract_sr(samtools,bam,sam,region):
    cmd1 = " -L %s" % region if (region != "") and (os.path.exists(region)) else ""
    cmd = " %s view -F 1024 -q 20 %s %s| awk \'{if (($6~/S/)||($6~/H/))print}\' > %s" % (samtools,cmd1,bam,sam)
    os.system(cmd)   

def cal_dr_sp(samtools,bam,stat,region):
    cmd1 = " -L %s" % region if (region != "") and (os.path.exists(region)) else ""
    cmd = " %s view -F 2 -F 1024 %s %s | grep -v \"SA\" | awk '{a=$7;if($7 == \"=\")a=$3;print $3\"\t\"$4\"\t\"a\"\t\"$8}\' | uniq -c > %s" % (samtools,cmd1,bam,stat)
    os.system(cmd)
    with open(stat,"r") as fd:
        for eachline in fd.readlines():
            line = eachline.strip().split()
            chr1 = line[1]
            chr2 = line[3]
            support = line[0]
            pos1 = line[2][:-3] if len(line[2]) > 3 else "-"
            pos2 = line[4][:-3] if len(line[4]) > 3 else "-"
            if (drunit.has_key((chr1,pos1,chr2,pos2))):
                drunit[(chr1,pos1,chr2,pos2)] = drunit[(chr1,pos1,chr2,pos2)] + int(support)
            elif (drunit.has_key((chr2,pos2,chr1,pos1))):
                drunit[(chr2,pos2,chr1,pos1)] = drunit[(chr2,pos2,chr1,pos1)] + int(support)
            else:
                drunit[(chr1,pos1,chr2,pos2)] = int(support)

def sam2fasta(sam,fasta):
    fasta_f = open(fasta,"w")
    with open(sam,"r") as sam_f:
        for eachline in sam_f.readlines():
            line = eachline.strip().split()
            if int(line[1]) & 2048 :
                continue
            fasta_f.write(">%s%s\n%s\n"%(line[0],line[1],line[9]))
    fasta_f.close()

def fd_lasm_run(fd_lasm,bam,sv_file,fastq):
    cmd = "%s %s %s -c 1,8 > %s" % (fd_lasm,bam,sv_file,fastq)
    os.system(cmd)

def bwa_aln(bwa,fastq,ref,sam):
    cmd = "%s mem  %s %s > %s" % (bwa,ref,fastq,sam)
    os.system(cmd)

def blat_aln(blat,assembly_fa,database_fa,pslx):
    cmd = "%s %s %s -out=pslx %s" % (blat,database_fa,assembly_fa,pslx)
    os.system(cmd)

def classify(chr1,pos1,dir1,order1,chr2,pos2,dir2,order2):
    if chr1 != chr2:
        return "CTX"
    elif dir1 == dir2:
        return "INV"
    elif int(pos1) < int(pos2) and dir1 == ">":
        return "DEL"
    elif int(pos1) > int(pos2) and dir1 == "<":
        return "DEL"
    elif int(pos1) < int(pos2) and dir1 == "<":
        return "INV"
    elif int(pos1) > int(pos2) and dir1 == ">":
        return "INV"
    else:
        return "-"
    

def print_report(bam,pslx,out):
    sv_result = OrderedDict()
    with open(pslx,"r") as fp:
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        line = fp.readline()
        while(line):
            print line
            match = int(line.strip().split("\t")[0])
            info = line.strip().split("\t")[9]
            zlen = int(line.strip().split("\t")[14])
            start = int(line.strip().split("\t")[11])
            end = int(line.strip().split("\t")[12])
            [chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,dis,insert] = info.split(";")
            if match > zlen * 0.9 and int(dis) > start and int(dis) < end :
                if sv_result.has_key((chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert)):
                     sv_result[(chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert)] = sv_result[(chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert)] + 1
                else:
                     sv_result[(chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert)] = 1
            line = fp.readline()
    with open(out,"w") as fo:
        fo.write("#chr1	pos1	gene1	dir1	chr2	pos2	gene2	dir2	bpdistance	type	softclipreads_support	blat_support	breakpoint_depth\n")
        for key,value in sv_result.items():
            (chr1,pos1,dir1,order1,chr2,pos2,dir2,order2,bpdistance,bpsupports,insert) = key
            gene1,gene2,depth = deal_with_sv("\t".join(key),bam)
            sv_type = classify(chr1,pos1,dir1,order1,chr2,pos2,dir2,order2)
            fo.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr1,pos1,gene1,dir1,chr2,pos2,gene2,dir2,str(bpdistance),sv_type,bpsupports,str(value),str(depth),insert))


def fusion_anno(chrN,pos,dirN,gene,gene_data):
    if gene_data[gene]["gene_ord"] == "-":
        exonlist = gene_data[gene]["gene_exon"][::-1]
        exonrg = [gene_data[gene]["gene_exon"][-1][0],gene_data[gene]["gene_exon"][0][0]]
    else:
        exonlist = gene_data[gene]["gene_exon"]
        exonrg = [gene_data[gene]["gene_exon"][0][0],gene_data[gene]["gene_exon"][-1][0]]
    gene_start = int(gene_data[gene]["gene_range"].split(":")[1].split("-")[0])
    gene_end = int(gene_data[gene]["gene_range"].split(":")[1].split("-")[1])
    if int(pos)>gene_end or int(pos)<gene_start:
        return 0
    pos_r = "%d-%d"%(gene_start,int(pos)) if dirN == ">" else "%d-%d"%(int(pos),gene_end)
    for i in xrange(len(exonlist)-1):
        if int(pos) < int(exonlist[i][1]):
            pass
        elif int(pos) < int(exonlist[i][2]) and int(pos) >= int(exonlist[i][1]):
            exon_pos = "%d+%d"%(int(exonlist[i][1]),abs(int(exonlist[i][1])-int(pos))) if abs(int(exonlist[i][1])-int(pos)) <  abs(int(exonlist[i][2])-int(pos))  else "%d-%d"%(int(exonlist[i][2]),abs(int(exonlist[i][2])-int(pos)))
            exon_t = exonlist[i][0]
            exon_r = "E%s-E%s"%(exonrg[0],exon_t) if dirN == ">" else "%s-%s"%(exon_t,exonrg[1])
            flag = "exon"
            break
        elif int(pos) < int(exonlist[i+1][1]) and int(pos) >= int(exonlist[i][2]):
            exon_pos = "%d-%d"%(int(exonlist[i+1][1]),abs(int(exonlist[i+1][1])-int(pos))) if abs(int(exonlist[i+1][1])-int(pos)) <  abs(int(exonlist[i][2])-int(pos))  else "%d+%d"%(int(exonlist[i][2]),abs(int(exonlist[i][2])-int(pos)))
            exon_t = exonlist[i+1][0] if abs(int(exonlist[i+1][1])-int(pos)) <  abs(int(exonlist[i][2])-int(pos)) else exonlist[i][0]
            exon_r = "E%s-E%s"%(exonrg[0],exon_t) if dirN == ">" else "%s-%s"%(exon_t,exonrg[1])
            flag = "intron"
            break
        elif int(pos) < int(exonlist[i+1][2]) and int(pos) >= int(exonlist[i+1][1]):
            exon_pos = "%d+%d"%(int(exonlist[i+1][1]),abs(int(exonlist[i+1][1])-int(pos))) if abs(int(exonlist[i+1][1])-int(pos)) <  abs(int(exonlist[i+1][2])-int(pos))  else "%d-%d"%(int(exonlist[i+1][2]),abs(int(exonlist[i+1][2])-int(pos)))
            exon_t = exonlist[i+1][0]
            exon_r = "E%s-E%s"%(exonrg[0],exon_t) if dirN == ">" else "%s-%s"%(exon_t,exonrg[1])
            flag = "exon"
            break
    return pos_r,exon_r,exon_t,exon_pos,gene_data[gene]["gene_ord"],flag,gene_data[gene]["gene_id"]


def sv_filter2(out,out_filtered,anno,gene_data):
    filtered = open(out_filtered,"w")
    anno_file = open(anno,"w")
    filtered.write("#chr1 pos1    gene1   dir1  chr2    pos2    gene2   dir2    bpdistance      type    softclipreads_support   blat_support    breakpoint_depth        filter_tag\n")
    anno_file.write("#Gene_Name Last_Observed_Exon      Inferred_Breakpoint     Head_breakpoint_location/Strand Gene_Name Last_Observed_Exon      Inferred_Breakpoint     Head_breakpoint_location/Strand       Remaining_sequence      Remaining_Protein\n")
    with open(out,"r") as sv_file:
        for eachline in sv_file.readlines():
            if eachline[0] == "#":
                continue
            line = eachline.strip().split()
            gene1,gene2 = line[2],line[6]
            if int(line[10]) < 2 or int(line[11]) < 2 or int(line[12]) < 5:  #modified by shenjia 20180126
                continue
            else:
                continue
            b_anno_1 = fusion_anno(line[0],line[1],line[3],line[2],gene_data)
            b_anno_2 = fusion_anno(line[4],line[5],line[7],line[6],gene_data)
            if b_anno_1==0 or b_anno_2 ==0:
                print "Maybe something with wrong in Fusion Annotation!Attention!"
                continue
            anno_file.write("%s\t%s:%s\t%s\t%s:%s/%s\t"%(gene1,b_anno_1[5],b_anno_1[2],b_anno_1[3],line[0],line[1],b_anno_1[4]))
            anno_file.write("%s\t%s:%s\t%s\t%s:%s/%s\t"%(gene2,b_anno_2[5],b_anno_2[2],b_anno_2[3],line[4],line[5],b_anno_2[4]))
            anno_file.write("%s_%s{%s:%s}|%s_%s{%s:%s}\t%s|%s\n"%(gene1,b_anno_1[-1],line[0],b_anno_1[0],gene2,b_anno_2[-1],line[4],b_anno_2[0],b_anno_1[1],b_anno_2[1]))
    anno_file.close()
    filtered.close()


def sv_filter(out,out_filtered,file_inhouse,file_cosmic):
    cosmic = {}
    inhouse = {}
    into_database(file_cosmic,cosmic)
    into_database(file_inhouse,inhouse)
    filtered = open(out_filtered,"w")
    filtered.write("#chr1 pos1    gene1   dir1	chr2    pos2    gene2   dir2	bpdistance	type    softclipreads_support   blat_support    breakpoint_depth	insert_seq	filter_tag\n")
    with open(out,"r") as sv_file:
        for eachline in sv_file.readlines():
            if eachline[0] == "#":
                continue
            line = eachline.strip().split()
            gene1,gene2 = line[2],line[6]
            if int(line[10]) < 2 or int(line[11]) < 2 or int(line[12]) < 5:
                continue
            if inhouse.has_key((gene1,gene2)):
                filtered.write(eachline.strip()+"\tinhouse\n")
            elif cosmic.has_key((gene1,gene2)) and not inhouse.has_key((gene1,gene2)):
                filtered.write(eachline.strip()+"\tcosmic\n")
            else:
                continue
    filtered.close()

def print_time():
    return time.asctime(time.localtime(time.time()))

def rm_file(filename):
    cmd = "rm %s" % filename
    os.system(cmd)


def main():
    parser = argparse.ArgumentParser(description="This tool is used to detect structural variants from a BAM file using local assembly.", prog=prog_name)
    parser.add_argument('--bamfile', '-b',  required=True, metavar='bamfile',  type=str, help='BAM file created by BWA or other alignment software.')
    parser.add_argument('--outputdir', '-o', required=True, metavar='outdir', type=str, help='Output direction for results.')
    parser.add_argument('--sampleid', '-i', required=True, metavar='sampleid', type=str, help='Sample ID.')
    parser.add_argument('--configinfo', '-c', required=True, metavar='config_file', type=str, help='Config file direction.')
    parser.add_argument('--targetregion', '-r', metavar='targetregion', type=str, default="", help='Detect Region.')
    parser.add_argument('--bpthreshold', '-p', metavar='bpthreshold', type=int, default=2000, help='Breakpoint Mini Distance.')
    parser.add_argument('--spthreshold', '-s', metavar='spthreshold', type=int, default=5, help='Breakpoint Mini Support Counts.')
    parser.add_argument('--mpthreshold', '-m', metavar='mpthreshold', type=int, default=20, help='Mini Mapping Quality.')
    parser.add_argument('--hard_filter', '-f', metavar='hard_filter', type=bool, default=True, help='Weather to using hard filter.(EXPERIMENTAL)')
    args = parser.parse_args()
    parser_config(args.configinfo)
    samtools = config["samtools "]
    bwa = config["bwa "]
    blat = config["blat "]
    ref = config["ref_hg19 "]
    fd_lasm = config["lasm_lite "]
    simplerepeatbed = config["simplerepeatbed "]
    genefile = config["genelist4sv "]
    our_model = config["our_model "]
    bam = args.bamfile
    outdir = args.outputdir
    sampleid = args.sampleid
    bp_threshold = args.bpthreshold
    sp_threshold = args.spthreshold
    mp_threshold = args.mpthreshold
    t_region = args.targetregion
    hard_filter = args.hard_filter
    sam = "%s/%s.sam" % (outdir,sampleid)
    stat = "%s/%s.dr.stat" % (outdir,sampleid)
    fasta = "%s/%s.fasta" % (outdir,sampleid)
    sv_filename = "%s/%s.sv" % (outdir,sampleid)
    fastq_lasm = "%s/%s.assembly.fastq" % (outdir,sampleid)
    sam_lasm = "%s/%s.assembly.sam" % (outdir,sampleid)
    fasta_lasm = "%s/%s.assembly.fasta" % (outdir,sampleid)
    pslx = "%s/%s.pslx" % (outdir,sampleid)
    out = "%s/%s.result" % (outdir,sampleid)
    filtered_out = "%s/%s_filtered_fusions.xls" % (outdir,sampleid)
    anno_out = "%s/%s_anno_fusions.xls" % (outdir,sampleid)
    load_bed(simplerepeatbed)
    gene_data = load(genefile)
    dict2list(gene_data)
    print "Step1: extract softclip reads from bamfile"
    print ">>> start at %s" % print_time()
    extract_sr(samtools,bam,sam,t_region)
    cal_dr_sp(samtools,bam,stat,t_region)
    print "<<< end at %s" % print_time()
    print "Step2: transfer samfile to fastafile"
    print ">>> start at %s" % print_time()
    sam2fasta(sam,fasta)
    print "<<< end at %s" % print_time()
    print "Step3: call sv from samfile"
    print ">>> start at %s" % print_time()
    call_sv(sv_filename,sam,bp_threshold,sp_threshold,mp_threshold,samtools,ref,hard_filter,our_model)
    print "<<< end at %s" % print_time()
    print "Step4: local assembly"
    print ">>> start at %s" % print_time()
    fd_lasm_run(fd_lasm,bam,sv_filename,fastq_lasm)
    print "<<< end at %s" % print_time()
    print "Step5: realignment"
    print ">>> start at %s" % print_time()
    bwa_aln(bwa,fastq_lasm,ref,sam_lasm)
    print "<<< end at %s" % print_time()
    print "Step6: extract sv from samfile created by realignment"
    print ">>> start at %s" % print_time()
    extract_fasta(sam_lasm,fasta_lasm,bam,sv_filename)
    print "<<< end at %s" % print_time()
    print "Step7: blat sv fasta to reference to calculate the support counts"
    print ">>> start at %s" % print_time()
    blat_aln(blat,fasta_lasm,fasta,pslx) 
    print "<<< end at %s" % print_time()
    print "Step8: print report"
    print ">>> start at %s" % print_time()
    print_report(bam,pslx,out)
    rm_file(sam)
    rm_file(fasta)
    sv_filter2(out,filtered_out,anno_out,gene_data)
    print "<<< end at %s" % print_time()


if __name__ == "__main__":
    main()
