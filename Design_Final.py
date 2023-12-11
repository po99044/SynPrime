#%%
import pandas as pd
import numpy as np
from scipy import stats
import numpy as np
import math
import os
import random
import Endo_REF
import synony_dic
import Pam_disrupt
import Gen_Syn


def buffer(RTPBS, totallen):
    remn_buflen = totallen- len(RTPBS)

    nucleo_lst = ["A","T","G","C"]
    nucleo_dic = {"A":["A","G","C"], "T":["T","G","C"],"G":["T","A","G"],"C":["T","C","A"]}

    inter_lst = []

    if  remn_buflen > len(RTPBS):
        y = (remn_buflen % len(RTPBS))
        x = (remn_buflen // len(RTPBS))
        for k in range(remn_buflen):
            remn_buf = RTPBS * x + RTPBS[:y]        
            sub_lst = nucleo_dic[remn_buf[k]]
            a = sub_lst.copy()
            s = random.choice(a)
            if k<2:
                inter_lst.append(s)
            else:
                if inter_lst[k-2] != inter_lst[k-1]:
                    inter_lst.append(s)
                elif s == inter_lst[k-1]:
                    a.remove(s)
                    s = random.choice(a)
                    inter_lst.append(s)
                else: 
                    inter_lst.append(s)
        buff = ''.join(inter_lst)
        

    else:
        for k in range(remn_buflen):
            ran_num = len(RTPBS) - remn_buflen + 1 
            case = random.randint(1,ran_num)    
            remn_buf = RTPBS[case-1:case-1+remn_buflen]
            sub_lst = nucleo_dic[remn_buf[k]]
            a = sub_lst.copy()
            s = random.choice(a)
            if k<2:
                inter_lst.append(s)
            else:
                if inter_lst[k-2] != inter_lst[k-1]:
                    inter_lst.append(s)
                elif s == inter_lst[k-1]:
                    a.remove(s)
                    s = random.choice(a)
                    inter_lst.append(s)
                else: 
                    inter_lst.append(s)
        buff = ''.join(inter_lst)

    return buff

def complement(seq):
    nu_dic = {'A': 'T', 'G':'C', 'T': 'A', 'C': 'G'}
    seq_len = len(seq)
    lst = []
    for i in range(seq_len):
        lst.append(nu_dic[seq[i]])
    
    return ''.join(lst)

def reverse_seq(seq):
    seq_list=list(seq)
    seq_list.reverse()
    seq = ''.join(seq_list)
    return seq

def reverse_complement(seq):
    a = reverse_seq(seq)
    b = complement(a)
    
    return b

base_dic = {"A":["T","C","G"], "T":["A","C","G"], "C":["T","A","G"],"G":["A","C","T"] }
nu_dic = {'A': 'T', 'G':'C', 'T': 'A', 'C': 'G'}

## Variable ##
exo = "E0"
rem1 = ""  # CDS 시작부분 Codon 부족부분 있으면 표기
rem2 = "" # CDS 끝부분 Codon 부족부분 있으면 표기

OUT_DIR = exo
os.mkdir("AR_Modified_HEK/%s" % OUT_DIR)

# %%

files = os.listdir("AR_DeepPRIME_HEK/AR_%s"  % exo) 
for file in files:


    df = pd.read_parquet("AR_DeepPRIME_HEK/AR_%s/%s" % (exo, file))
    if len(df) == 0: 
        print (file)
        pass
    else:

        ref = Endo_REF.exon(exo)  # 60 + CDS + 60
        exon_len = len(Endo_REF.exon(exo)) - 120
        alpha = len(rem1)
        beta = len(rem2)
        ref2 = ref[:60] + rem1 + ref[60:60+exon_len] + rem2 + ref[60+exon_len:]
        cds = ref2[60:-60]


        # Input form1

        for i in df.index:
            a = df.loc[i,"Edited74_On"]
            pbs_len = df.loc[i,'PBSlen']
            rtpbs_len = df.loc[i,"RT-PBSlen"]
            rt_len = rtpbs_len - pbs_len
            rha = df.loc[i,"RHA_len"]
            ID = df.loc[i,"ID"]
            edit_pos = 60 + int(ID.split('_')[2][1:-1]) - 5  # 2 original, 3 modify
            df.loc[i,"EDIT_POS"] = edit_pos

            # RTPBS
            df.loc[i,"RTPBS_R"] = a.replace('x','')
            df.loc[i,"RTPBS"] = reverse_complement(df.loc[i,"RTPBS_R"])
            RTPBS = df.loc[i,"RTPBS"]


            pbs_r = df.loc[i,"RTPBS_R"][:pbs_len]
            pam_site = a.find(pbs_r) + pbs_len + 3  # PAM start site based on WT74
            sgRNA = df.loc[i,'WT74_On'][pam_site - 19:pam_site] 


            df.loc[i,"sgRNA"] = sgRNA
            df.loc[i,"LHA"] = int(df.loc[i,"RTlen"])- int(df.loc[i,"RHA_len"]) -1
            


            df.loc[i,"PBS"] = df.loc[i,"RTPBS"][rt_len:]
            df.loc[i,"RTPBS"] = RTPBS[:rha] + RTPBS[rha].lower() + RTPBS[rha+1:]
            
            # strand
            if sgRNA in ref:
                df.loc[i,"strand"] = "F"
            else:
                df.loc[i,"strand"] = "R"



            if df.loc[i,"strand"] =="F":
                    pam_loc = ref.find(sgRNA) + 19  # pam start loc
                    df.loc[i,"PAM_LOC"] = pam_loc

            elif df.loc[i,"strand"] =="R":
                    pam_loc = ref.find(reverse_complement(sgRNA)) -3  # pam start loc
                    df.loc[i,"PAM_LOC"] = pam_loc

            
            df.loc[i,"PAMdist_fromedit"] = abs(edit_pos-pam_loc)  # PAM start site. Edit pos으로부터상대적인. NGG에서 N 위치.

        df = df.sort_values(by=["PE2max-e_score"], ascending= False)

        df = df[['ID', 'WT74_On', 'Edited74_On', 'PBSlen', 'RTlen', 'RT-PBSlen',
            'Edit_pos', 'RHA_len', 'PE2max-e_score', 'EDIT_POS', 'RTPBS_R', 'RTPBS',
            'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit']]


        df = Pam_disrupt.pam_disrupt(df,alpha, ref2,ref,cds)

        if len(df.columns) != 25:
            
            if df[df["RHA_len"] >=6].shape[0] == 0:
                df = Gen_Syn.syn_gen(df, alpha, beta, cds, ref, ref2)
                result = df
                
                result = result.sort_values(by="PE2max-e_score", ascending=False)
                try:
                    result = result.dropna(subset = ["SynonyRTPBS"], axis =0)
                    result = result.reset_index()
                    result["synony pam disrupt"] = "F"
                    result["Intended pam disrupt"] = "F"
                    result["Intend_to_pamdisrupt_dist"] = "F"
                    result["Pam disrupt codon"] = "F"
                    

                    result =result[['ID', 'WT74_On', 'Edited74_On', 'PBSlen', 'RTlen', 'RT-PBSlen',
            'Edit_pos', 'RHA_len', 'PE2max-e_score', 'EDIT_POS', 'RTPBS_R', 'RTPBS',
            'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit',
            'SynonyRTPBS', 'synony pam disrupt', 'Intended pam disrupt',
            'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]
                    result.to_csv("AR_Modified_HEK/%s/%s.csv" % (OUT_DIR, file) )
                except KeyError : print(file)         
            else:
                
                df = df[df["RHA_len"] >=6]
                df = Gen_Syn.syn_gen(df, alpha, beta, cds, ref, ref2)
                result = df
                
                result = result.sort_values(by="PE2max-e_score", ascending=False)

                try:
                    result = result.dropna(subset = ["SynonyRTPBS"], axis =0)
                    result = result.reset_index()
                    result["synony pam disrupt"] = "F"
                    result["Intended pam disrupt"] = "F"
                    result["Intend_to_pamdisrupt_dist"] = "F"
                    result["Pam disrupt codon"] = "F"


                    result =result[['ID', 'WT74_On', 'Edited74_On', 'PBSlen', 'RTlen', 'RT-PBSlen',
            'Edit_pos', 'RHA_len', 'PE2max-e_score', 'EDIT_POS', 'RTPBS_R', 'RTPBS',
            'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit',
            'SynonyRTPBS', 'synony pam disrupt', 'Intended pam disrupt',
            'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]
                    result.to_csv("AR_Modified_HEK/%s/%s.csv" % (OUT_DIR, file) )
                except KeyError : print(file)

        else:
            
            df1 = df.loc[(df["synony pam disrupt"] == "T"),] ## PAM disrupting completed rows
            df2 = df.loc[(df["synony pam disrupt"] != "T"),] ## need to further evoke synonymous rows


            df2 = df2[df2["RHA_len"] >= 6]

            ## Not-Pam disrupting guides
            df2 = Gen_Syn.syn_gen(df2, alpha, beta, cds, ref, ref2)


            result = pd.concat([df1, df2])

            result = result.sort_values(by="PE2max-e_score", ascending=False)
            try:

                result = result.dropna(subset = ["SynonyRTPBS"], axis =0)

                result = result.reset_index()

                result =result[['ID', 'WT74_On', 'Edited74_On', 'PBSlen', 'RTlen', 'RT-PBSlen',
        'Edit_pos', 'RHA_len', 'PE2max-e_score', 'EDIT_POS', 'RTPBS_R', 'RTPBS',
        'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit',
        'SynonyRTPBS', 'synony pam disrupt', 'Intended pam disrupt',
        'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]



                result.to_csv("AR_Modified_HEK/%s/%s.csv" % (OUT_DIR, file) )
            
            except KeyError:
                print(file)
            
        