from matplotlib.lines import lineStyles
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from moepy import lowess
import seaborn as sns 
import matplotlib.font_manager as fm
from datetime import datetime
from sklearn.metrics import roc_curve
from sklearn import metrics
from sklearn.metrics import PrecisionRecallDisplay, precision_recall_curve, average_precision_score
from scipy import stats
import numpy as np
import math
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
import os
from statannot import add_stat_annotation
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
font_path = 'C:\\WINDOWS\\Fonts\\arial.ttf'
font_name = fm.FontProperties(fname=font_path).get_name()
plt.rcParams['font.family'] = font_name
plt.rcParams['font.size'] = 20
import random
import Endo_REF
import synony_dic

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


def pam_disrupt(df,alpha,ref2,ref, cds):

    for i in df.index:
        strand = df.loc[i,"strand"]
        edit_pos = int(df.loc[i,"EDIT_POS"])
        rtpbs = df.loc[i,"RTPBS"]
        rha = df.loc[i,"RHA_len"]
        ID = df.loc[i,"ID"]
        nu = ID.split('_')[2][-1]  ## 82의 경우만 3, original 2
        
        if strand == "F":
            pam_loc=  int(df.loc[i,"PAM_LOC"] + 2)
            # PAM, Edit codon number 설정
            # codon 내에서의 위치관계설정
            pam_codon_loc = math.ceil((pam_loc+1  + alpha - 60)/3) 
            pam_rem = (pam_loc+1  +alpha)%3
            edit_codon_loc = math.ceil((edit_pos+1  + alpha - 60)/3) 
            
            if pam_codon_loc <=1:
                df.loc[i,"synony pam disrupt"] ="F" 
                df.loc[i,"Intended pam disrupt"] ="F"            
            elif pam_codon_loc > len(cds)/3:
                df.loc[i,"synony pam disrupt"] ="F" 
                df.loc[i,"Intended pam disrupt"] ="F"
            else:
                
                if pam_codon_loc == edit_codon_loc:  # 굳이 Intended edit codon과
                                                    # PAM disrupting codon을 나눈 이유는
                                                    # 다른 코돈에 Synony 유발하는 컨셉.

                    if pam_rem != 1 :   # 2 or 0이면 PAM disrupt하면서 Synony 만들면서 다른 코돈이 불가..
                        df.loc[i,"synony pam disrupt"] ="F" 
                        df.loc[i,"Intended pam disrupt"] ="T"

                    
                    elif pam_rem == 1:  # Intend codon 앞 코돈에 유발가능
                        pam_codon = ref2[(pam_loc+alpha-3):(pam_loc+alpha)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)
                        lst = synony_dic.synony_code(pam_codon)
                        llst = []
                        com = ref2[(pam_loc+alpha-3):(pam_loc+alpha-1)]
                        
                        if len(lst) != 0:
                            for k in ["A","C","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if int(rha+ edit_pos-pam_loc +1) > len(r_lst) -1 : pass
                                else:
                                    r_lst[int(rha+ edit_pos-pam_loc +1)] = nu_dic[modi_codon[-1]].lower()

                                    
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc - 1 # RHA Ok
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc-1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    

                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="T"  


                        elif len(lst) == 0:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="T"  


                elif pam_codon_loc < edit_codon_loc :

                    if pam_rem ==2 :
                        df.loc[i,"synony pam disrupt"] ="F" 
                        df.loc[i,"Intended pam disrupt"] ="F"

                    elif pam_rem ==1 : # XNG G

                        pam_codon = ref2[(pam_loc+alpha-3):(pam_loc+alpha)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)
                        lst = synony_dic.synony_code(pam_codon)
                        llst = []
                        com = ref2[(pam_loc+alpha-3):(pam_loc+alpha-1)]
                        
                        if len(lst) != 0:
                            for k in ["A","C","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if int(rha+ edit_pos-pam_loc+1)> len(r_lst) -1 : pass
                                else:    
                                    r_lst[int(rha+ edit_pos-pam_loc+1)] = nu_dic[modi_codon[-1]].lower()
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc -1 # RHA ok
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter
                                    

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc -1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                       
                                    
                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"  
                                

                        elif len(lst) == 0:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F"  
                            

                    elif pam_rem ==0 : # NGG

                        pam_codon = ref2[(pam_loc+alpha-2):(pam_loc+alpha+1)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)
                        lst = synony_dic.synony_code(pam_codon)
                        llst = []
                        com = ref2[(pam_loc+alpha-2):(pam_loc+alpha)]
                        
                        if len(lst) != 0:
                            for k in ["A","C","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if int(rha+ edit_pos-pam_loc) > len(r_lst) -1: pass
                                else:    
                                    r_lst[int(rha+ edit_pos-pam_loc)] = nu_dic[modi_codon[-1]].lower()
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc # RHA ok
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter
                                    

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                       
                                    

                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"  
                                
                        elif len(lst) == 0:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F" 
                            

                    print("2nd")
                    
                    
                elif pam_codon_loc > edit_codon_loc :


                    if pam_rem ==2 :
                        df.loc[i,"synony pam disrupt"]="F" 
                        df.loc[i,"Intended pam disrupt"] ="F"
                        
                    
                    elif pam_rem ==1 : # XNG G
                        if pam_codon_loc - edit_codon_loc ==1 :  # 앞코돈 = Intended
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="T" 
                            
                        else:
                            pam_codon = ref2[(pam_loc+alpha-3):(pam_loc+alpha)]
                            pam_codon_letter = Endo_REF.gene_code(pam_codon)
                            lst = synony_dic.synony_code(pam_codon)
                            llst = []
                            com = ref2[(pam_loc+alpha-3):(pam_loc+alpha-1)]
                            
                            if len(lst) != 0:
                                for k in ["A","C","T"]:
                                    if com+k in lst: llst.append(com+k)
                                    else: pass

                                if len(llst) != 0:   
                                    modi_codon = random.choice(llst) ## 
                                    r_lst = list(rtpbs)
                                    if int(rha- (pam_loc-edit_pos) +1) > len(r_lst) -1: pass
                                    else:    
                                        r_lst[int(rha- (pam_loc-edit_pos) +1)] = nu_dic[modi_codon[-1]].lower()
                                        dist = edit_pos-pam_loc
                                        synonyrtpbs = ''.join(r_lst)
                                        df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                        df.loc[i,"synony pam disrupt"] ="T" 
                                        df.loc[i,"Intended pam disrupt"] ="F"  
                                        df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                        df.loc[i,"Pam disrupt codon"] = pam_codon_loc -1
                                        df.loc[i,"RHA_len"] = rha- (pam_loc-edit_pos) +1
                                        df.loc[i,"SynonyCodon"] = pam_codon_letter

                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[pam_loc-1] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq    
                                        

                                else:
                                    df.loc[i,"synony pam disrupt"] ="F" 
                                    df.loc[i,"Intended pam disrupt"] ="F"
                                    

                            elif len(lst) == 0:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"  
                                

                    elif pam_rem ==0 : # NGG

                        pam_codon = ref2[(pam_loc+alpha-2):(pam_loc+alpha+1)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)
                        lst = synony_dic.synony_code(pam_codon)
                        llst = []
                        com = ref2[(pam_loc+alpha-2):(pam_loc+alpha)]
                        
                        if len(lst) != 0:
                            for k in ["A","C","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if  int(rha- (pam_loc-edit_pos))  > len(r_lst) -1: pass
                                else:    
                                    r_lst[int(rha- (pam_loc-edit_pos))] = nu_dic[modi_codon[-1]].lower()
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc
                                    df.loc[i,"RHA_len"] = rha- (pam_loc-edit_pos)
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                      


                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"
                                            
                            
                        elif len(lst) == 0:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F"  
                            

                    print("3rd")

    ###### Reverse strand ###########
        elif strand == "R":
            # pam_loc CCN에서 '첫 C' (NGG의 G 위치임)
            # PAM, Edit codon number 설정
            # codon 내에서의 위치관계설정
            pam_loc=  int(df.loc[i,"PAM_LOC"])
            pam_codon_loc = math.ceil((pam_loc+1  + alpha - 60)/3) 
            pam_rem = (pam_loc+1  +alpha)%3
            edit_codon_loc = math.ceil((edit_pos+1  + alpha - 60)/3) 
            if pam_codon_loc <=1: 
                df.loc[i,"synony pam disrupt"] ="F" 
                df.loc[i,"Intended pam disrupt"] ="F"                  
            elif pam_codon_loc > len(cds)/3:
                df.loc[i,"synony pam disrupt"] ="F" 
                df.loc[i,"Intended pam disrupt"] ="F"                  
            else:
                    
                
                if pam_codon_loc == edit_codon_loc:  # 굳이 Intended edit codon과
                                                    # PAM disrupting codon을 나눈 이유는
                                                    # 다른 코돈에 Synony 유발하는 컨셉.

                    if pam_rem != 0 :  ## CCN condon = intended edit codon --> pass
                        df.loc[i,"synony pam disrupt"] ="F" 
                        df.loc[i,"Intended pam disrupt"] ="T"
                        
                    
                    elif pam_rem == 0:  # Intend codon 뒤 코돈에 유발가능 
                        ## CUA -> UUA // CUG -> UUG (LEU) 만 가능

                        pam_codon = ref2[(pam_loc+alpha+1):(pam_loc+alpha+4)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)

                        if pam_codon == "CTA":
                                r_lst = list(rtpbs)
                                if int(rha+ 1 + (pam_loc-edit_pos)) > len(r_lst) -1: pass
                                else:
                                    r_lst[int(rha+ 1 + (pam_loc-edit_pos))] = "t"
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist # RHA ok
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter


                                    ref_lst = list(ref)
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc+1] = "t" # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                      
                                    

                                    # RHA 길이는 변함없음 ok
                        elif pam_codon == "CTG":
                                r_lst = list(rtpbs)
                                if int(rha+ 1 + (pam_loc-edit_pos)) > len(r_lst) -1: pass
                                else:
                                    r_lst[int(rha+ 1 + (pam_loc-edit_pos) )] = "t"
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter

                                    ref_lst = list(ref)
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc+1] = "t" # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                       
                        else:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F"  
                            

                if pam_codon_loc > edit_codon_loc:

                    if pam_rem == 0:  # PAM codon과 PAM +1 코돈에 유발가능 
                        ## CUA -> UUA // CUG -> UUG (LEU) 만 가능
                        pam_codon1 = ref2[(pam_loc+alpha-2):(pam_loc+alpha+1)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon1)

                        pam_codon2 = ref2[(pam_loc+alpha+1):(pam_loc+alpha+4)]
                        pam_codon_letter2 = Endo_REF.gene_code(pam_codon2)

                                
                        lst = synony_dic.synony_code(pam_codon1)
                        llst = []
                        com = ref2[(pam_loc+alpha-2):(pam_loc+alpha)]
                        
                        if len(lst) != 0:
                            for k in ["A","G","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if int(rha + (pam_loc-edit_pos))> len(r_lst) -1  : pass
                                else:

                                    r_lst[int(rha + (pam_loc-edit_pos))] = modi_codon[-1].lower()
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc # RHA ok
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                                     
                                    
                        elif pam_codon2 == "CTA":
                                r_lst = list(rtpbs)
                                if int(rha+ 1 + (pam_loc-edit_pos)) > len(r_lst) -1 : pass
                                else:
                                    r_lst[int(rha+ 1 + (pam_loc-edit_pos) )] = "t"
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist # RHA ok
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter
                                    
                            
                                    ref_lst = list(ref)
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc+1] = "t" # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                      
                                    
                                    
                                    # RHA 길이는 변함없음 ok
                        elif pam_codon2 == "CTG":
                                r_lst = list(rtpbs)
                                if int(rha+ 1 + (pam_loc-edit_pos)) > len(r_lst) -1 : pass
                                else:
                                    r_lst[int(rha+ 1 + (pam_loc-edit_pos) )] = "t"
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter
                                
                                    ref_lst = list(ref)
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc+1] = "t" # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                      
                                    
                                    
                        else:
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F"  
                            

                    elif pam_rem == 2: 
                        # XC'C' NXX

                        pam_codon = ref2[(pam_loc+alpha-1):(pam_loc+alpha+2)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon)
                        lst = synony_dic.synony_code(pam_codon)
                        llst = []
                        com = ref2[(pam_loc+alpha-1):(pam_loc+alpha+1)]
                        
                        if len(lst) != 0:
                            for k in ["A","G","T"]:
                                if com+k in lst: llst.append(com+k)
                                else: pass

                            if len(llst) != 0:   
                                modi_codon = random.choice(llst) ## 
                                r_lst = list(rtpbs)
                                if int(rha + (pam_loc-edit_pos)) +1 > len(r_lst) -1:pass
                                else:
                                    r_lst[int(rha + (pam_loc-edit_pos)) +1 ] = modi_codon[-1].lower()
                                    dist = edit_pos-pam_loc
                                    synonyrtpbs = ''.join(r_lst)
                                    df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                    df.loc[i,"synony pam disrupt"] ="T" 
                                    df.loc[i,"Intended pam disrupt"] ="F"  
                                    df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                    df.loc[i,"Pam disrupt codon"] = pam_codon_loc # RHA ok
                                    df.loc[i,"SynonyCodon"] = pam_codon_letter
                                    

                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[pam_loc +1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                     
                                    
                                    
                        else : 
                            df.loc[i,"synony pam disrupt"] ="F" 
                            df.loc[i,"Intended pam disrupt"] ="F"     
                            

                    elif pam_rem == 1:            
                        df.loc[i,"synony pam disrupt"] ="F" 
                        df.loc[i,"Intended pam disrupt"] ="F"  
                        

                if pam_codon_loc < edit_codon_loc:
                    if pam_rem == 1:            
                        df.loc[i,"synony pam disrupt"] ="F" 
                        df.loc[i,"Intended pam disrupt"] ="F"  
                        
                    
                    elif pam_rem == 0:

                        # PAM codon과 PAM +1 코돈에 유발가능 
                        ## CUA -> UUA // CUG -> UUG (LEU) 만 가능
                        pam_codon1 = ref2[(pam_loc+alpha-2):(pam_loc+alpha+1)]
                        pam_codon_letter = Endo_REF.gene_code(pam_codon1)

                        pam_codon2 = ref2[(pam_loc+alpha+1):(pam_loc+alpha+4)]
                        pam_codon_letter2 = Endo_REF.gene_code(pam_codon2)

                        lst = synony_dic.synony_code(pam_codon1)
                        llst = []
                        com = ref2[(pam_loc+alpha-2):(pam_loc+alpha)]

                        if edit_codon_loc - pam_codon_loc > 1:
                            
                            if len(lst) != 0:
                                for k in ["A","G","T"]:
                                    if com+k in lst: llst.append(com+k)
                                    else: pass

                                if len(llst) != 0:   
                                    modi_codon = random.choice(llst) ## 
                                    r_lst = list(rtpbs)
                                    if int(rha + (pam_loc-edit_pos)) > len(r_lst) -1: pass
                                    else:
                                        r_lst[int(rha + (pam_loc-edit_pos))] = modi_codon[-1].lower()
                                        dist = edit_pos-pam_loc
                                        synonyrtpbs = ''.join(r_lst)
                                        df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                        df.loc[i,"synony pam disrupt"] ="T" 
                                        df.loc[i,"Intended pam disrupt"] ="F"  
                                        df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                        df.loc[i,"Pam disrupt codon"] = pam_codon_loc
                                        df.loc[i,"RHA_len"] = rha - (edit_pos-pam_loc)
                                        df.loc[i,"SynonyCodon"] = pam_codon_letter    
                                        

                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[pam_loc] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                          
                                            
                

                            elif pam_codon2 == "CTA":
                                    r_lst = list(rtpbs)
                                    if int(rha+ 1 + (pam_loc-edit_pos) + 1)> len(r_lst) -1: pass
                                    else:
                                        r_lst[int(rha+ 1 + (pam_loc-edit_pos))] = "t"
                                        dist = edit_pos-pam_loc
                                        synonyrtpbs = ''.join(r_lst)
                                        df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                        df.loc[i,"synony pam disrupt"] ="T" 
                                        df.loc[i,"Intended pam disrupt"] ="F"  
                                        df.loc[i,"Intend_to_pamdisrupt_dist"] = dist 
                                        df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                        df.loc[i,"RHA_len"] = rha - (edit_pos-pam_loc) +1
                                        df.loc[i,"SynonyCodon"] = pam_codon_letter

                                        ref_lst = list(ref) 
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[pam_loc+1] = "t" # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                          
                                        

                                        # RHA 길이는 변함없음 ok
                            elif pam_codon2 == "CTG":
                                    r_lst = list(rtpbs)
                                    if int(rha+ 1 + (pam_loc-edit_pos) + 1) > len(r_lst) -1: pass
                                    else:
                                        r_lst[int(rha+ 1 + (pam_loc-edit_pos))] = "t"
                                        dist = edit_pos-pam_loc
                                        synonyrtpbs = ''.join(r_lst)
                                        df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                        df.loc[i,"synony pam disrupt"] ="T" 
                                        df.loc[i,"Intended pam disrupt"] ="F"  
                                        df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                        df.loc[i,"Pam disrupt codon"] = pam_codon_loc +1
                                        df.loc[i,"RHA_len"] = rha - (edit_pos-pam_loc) +1
                                        df.loc[i,"SynonyCodon"] = pam_codon_letter
        
                                        ref_lst = list(ref) 
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[pam_loc+1] = "t" # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                             
                            
                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"  
                                
                        if edit_codon_loc - pam_codon_loc == 1:
                            
                            if len(lst) != 0:
                                for k in ["A","G","T"]:
                                    if com+k in lst: llst.append(com+k)
                                    else: pass

                                if len(llst) != 0:   
                                    modi_codon = random.choice(llst) ## 
                                    r_lst = list(rtpbs)
                                    if int(rha + (pam_loc-edit_pos)) > len(r_lst) -1: pass
                                    else:
                                        r_lst[int(rha + (pam_loc-edit_pos))] = modi_codon[-1].lower()
                                        dist = edit_pos-pam_loc
                                        synonyrtpbs = ''.join(r_lst)
                                        df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                        df.loc[i,"synony pam disrupt"] ="T" 
                                        df.loc[i,"Intended pam disrupt"] ="F"  
                                        df.loc[i,"Intend_to_pamdisrupt_dist"] = dist
                                        df.loc[i,"Pam disrupt codon"] = pam_codon_loc
                                        df.loc[i,"RHA_len"] = rha - (edit_pos-pam_loc)
                                        df.loc[i,"SynonyCodon"] = pam_codon_letter

                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[pam_loc] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                          
                                        
                            else:
                                df.loc[i,"synony pam disrupt"] ="F" 
                                df.loc[i,"Intended pam disrupt"] ="F"  

        
        
    return df

