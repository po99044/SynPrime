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
exo = "E5"
rem1 = ""  # CDS 시작부분 Codon 부족부분 있으면 표기
rem2 = "" # CDS 끝부분 Codon 부족부분 있으면 표기

OUT_DIR = exo
os.mkdir("Output/%s" % OUT_DIR)

# %%

files = os.listdir("Input/SNCA_%s"  % exo) 
for file in files:


    df = pd.read_parquet("Input/SNCA_%s/%s" % (exo, file))
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
            pbs_len = df.loc[i,'PBS_len']
            rtpbs_len = df.loc[i,"RT-PBS_len"]
            rt_len = rtpbs_len - pbs_len
            rha = df.loc[i,"RHA_len"]
            ID = df.loc[i,"ID"]
            edit_pos = 60 + int(ID.split('_')[2][1:-1]) - 5
            df.loc[i,"Edit_pos"] = edit_pos

            RTPBS = df.loc[i,"RT-PBS"]
            df.loc[i,"RTPBS_R"] = reverse_complement(RTPBS)

            pbs_r = df.loc[i,"RTPBS_R"][:pbs_len]
            sgRNA = df.loc[i,'Spacer'][1:]
            df.loc[i,"sgRNA"] = sgRNA
            df.loc[i,"LHA"] = int(df.loc[i,"RTT_len"])- int(df.loc[i,"RHA_len"]) -1
            


            df.loc[i,"PBS"] = df.loc[i,"RT-PBS"][rt_len:]
            df.loc[i,"RT-PBS"] = RTPBS[:rha] + RTPBS[rha].lower() + RTPBS[rha+1:]
            
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

        df = df.sort_values(by=["PE2max_score"], ascending= False)
        df.columns = ['ID', 'PE2max_score', 'Spacer', 'RTPBS', 'PBSlen', 'RTlen',
            'RT-PBSlen', 'EDIT_POS', 'Edit_len', 'RHA_len', 'Target', 'RTPBS_R',
            'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit']

        df = Pam_disrupt.pam_disrupt(df,alpha, ref2,ref,cds)



        if len(df.columns) != 25:
            
            if df[df["RHA_len"] >=6].shape[0] == 0:
                df = Gen_Syn.syn_gen(df, alpha, beta, cds, ref, ref2)
                result = df
                
                result = result.sort_values(by="PE2max_score", ascending=False)
                try:
                    result = result.dropna(subset = ["SynonyRTPBS"], axis =0)
                    result = result.reset_index()
                    result = result[result["RHA_len"] >= 1]
                    result["synony pam disrupt"] = "F"
                    result["Intended pam disrupt"] = "F"
                    result["Intend_to_pamdisrupt_dist"] = "F"
                    result["Pam disrupt codon"] = "F"
                    

                    result =result[['ID', 'PE2max_score', 'Spacer', 'RTPBS', 'PBSlen', 'RTlen', 'RT-PBSlen',
       'EDIT_POS', 'Edit_len', 'RHA_len', 'Target', 'RTPBS_R', 'sgRNA', 'LHA',
       'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit', 'SynonyRTPBS',
       'synony pam disrupt', 'Intended pam disrupt',
       'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]
                    result.to_csv("Output/%s/%s.csv" % (OUT_DIR, file) )
                except KeyError : print(file)         
            else:
                
                df = df[df["RHA_len"] >=6]
                df = Gen_Syn.syn_gen(df, alpha, beta, cds, ref, ref2)
                result = df
                
                result = result.sort_values(by="PE2max_score", ascending=False)

                try:
                    result = result.dropna(subset = ["SynonyRTPBS"], axis =0)
                    result = result.reset_index()
                    result = result[result["RHA_len"] >= 1]
                    result["synony pam disrupt"] = "F"
                    result["Intended pam disrupt"] = "F"
                    result["Intend_to_pamdisrupt_dist"] = "F"
                    result["Pam disrupt codon"] = "F"


                    result =result[['ID', 'PE2max_score', 'Spacer', 'RTPBS', 'PBSlen', 'RTlen', 'RT-PBSlen',
       'EDIT_POS', 'Edit_len', 'RHA_len', 'Target', 'RTPBS_R', 'sgRNA', 'LHA',
       'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit', 'SynonyRTPBS',
       'synony pam disrupt', 'Intended pam disrupt',
       'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]
                    result.to_csv("Output/%s/%s.csv" % (OUT_DIR, file) )
                except KeyError : print(file)

        else:
            
            df1 = df.loc[(df["synony pam disrupt"] == "T"),] ## PAM disrupting completed rows
            df2 = df.loc[(df["synony pam disrupt"] != "T"),] ## need to further evoke synonymous rows


            df2 = df2[df2["RHA_len"] >= 6]

            ## Not-Pam disrupting guides
            df2 = Gen_Syn.syn_gen(df2, alpha, beta, cds, ref, ref2)


            result = pd.concat([df1, df2])

            result = result.sort_values(by="PE2max_score", ascending=False)
            try:

                result = result.dropna(subset = ["SynonyRTPBS"], axis =0)

                result = result.reset_index()
                result = result[result["RHA_len"] >= 1]
                result =result[['ID', 'PE2max_score', 'Spacer', 'RTPBS', 'PBSlen', 'RTlen', 'RT-PBSlen',
       'EDIT_POS', 'Edit_len', 'RHA_len', 'Target', 'RTPBS_R', 'sgRNA', 'LHA',
       'PBS', 'strand', 'PAM_LOC', 'PAMdist_fromedit', 'SynonyRTPBS',
       'synony pam disrupt', 'Intended pam disrupt',
       'Intend_to_pamdisrupt_dist', 'Pam disrupt codon', 'SynonyCodon', 'seq']]



                result.to_csv("Output/%s/%s.csv" % (OUT_DIR, file) )
            
            except KeyError:
                print(file)
            
