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


def syn_gen(df, alpha, beta, cds, ref, ref2):

    for i in df.index:
        strand = df.loc[i,"strand"]
        pam_loc=  int(df.loc[i,"PAM_LOC"] + 2)  # NGG에서 마지막 G의 위치
        edit_pos = int(df.loc[i,"EDIT_POS"])
        rtpbs = df.loc[i,"RTPBS"]
        rha = df.loc[i,"RHA_len"]
        lha = df.loc[i,"LHA"]
        pbs_len = len(df.loc[i,"PBS"])
        ID = df.loc[i,"ID"]
        nu = ID.split('_')[2][-1]
        
    ## 인트론부분 alpha or beta 있을 때 상황 다시 체크하기. RTPBS랑 Seq
    ## 3' 쪽은 +beta,  5' 쪽은 -alpha가 답
        
        try:
            
            if strand == "F":  # Forward strand여서 5' 쪽 위주로 생각
                rtpbs_terminal = int(edit_pos - lha - pbs_len)
                edit_codon_loc = math.ceil((edit_pos+1  + alpha - 60)/3) 
                edit_rem = (edit_pos+1 - 60 + alpha) %3   # OOO  1 2 0 순서
                nick_loc = pam_loc -6

                if edit_pos-55 in [0,1,2,3,4]:  #5'은 codon3,4 only
                    k = edit_pos - 55
                    if k == 0:
                        Syn_codon3 = ref2[edit_pos+5:edit_pos+8]
                        Syn_codon4 = ref2[edit_pos+8:edit_pos+11]
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                    
                        if (len(lst) != 0) & (edit_pos+5> nick_loc):  #CODON3
                            
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2 + alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 7-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+ 5-alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        if (len(llst) == 0) & (edit_pos+8> nick_loc):
                            
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2 + alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 10-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 8-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif k ==1:
                        Syn_codon3 = ref2[edit_pos+4:edit_pos+7]
                        Syn_codon4 = ref2[edit_pos+7:edit_pos+10]

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                                
                        if (len(lst) != 0) & (edit_pos+4> nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1 + alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 6-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+ 4-alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a

                        if (len(llst) == 0) & (edit_pos+7> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1 + alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 9-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 7-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif k ==2:
                        Syn_codon3 = ref2[edit_pos+3:edit_pos+6]
                        Syn_codon4 = ref2[edit_pos+6:edit_pos+9] 

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                    
                    

                        if (len(lst) != 0) & (edit_pos+3> nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 + alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 5-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+ 3-alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        if (len(llst) == 0) & (edit_pos+6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 + alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 8-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 6-alpha ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif k ==3:
                        Syn_codon3 = ref2[edit_pos+2:edit_pos+5]
                        Syn_codon4 = ref2[edit_pos+5:edit_pos+8]


                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                        
                    

                        if (len(lst) != 0) & (edit_pos+2> nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1 + alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 4-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+ 2-alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        if (len(llst) == 0) & (edit_pos+5> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1 + alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 7-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 5-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif k ==4:
                        Syn_codon3 = ref2[edit_pos+1:edit_pos+4]
                        Syn_codon4 = ref2[edit_pos+4:edit_pos+7]

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                        
                        
                        if (len(lst) != 0) & (edit_pos+1> nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2 + alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+3-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+ 1-alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        if (len(llst) == 0) & (edit_pos+4> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2 + alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 6-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 4-alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                elif edit_pos-55 in [len(ref)-120 + i for i in [5,6,7,8,9]]: # 3'은 codon 1,2로만..
                    k = edit_pos - 55

                    if k == len(ref)-120 + 5:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-3:edit_pos+alpha+beta]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-6:edit_pos+alpha+beta-3]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                        

                        if (len(lst) != 0) & (edit_pos+alpha+beta-3> nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 - beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                elif modi_codon =="TAA": 
                                    a += 1
                                    r_lst[int(rha +1 + a)] = "t"
                                    ref_lst = list(ref)
                                    syn_nu = "t"
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -2+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                    llst.append(1)

                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass   

                                else:
                                    a += 0
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -1+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass                            
                                    # r_lst.append(ref2[:]) ## 길이연장
                                    
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3+beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass             
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        if (len(llst) == 0) & (edit_pos+alpha+beta-6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0 - beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                         
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += 2
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                    
                    elif k == len(ref)-120 + 6:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-4:edit_pos+alpha+beta-1]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-7:edit_pos+alpha+beta-4]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                        

                        if (len(lst) != 0) & (edit_pos+alpha+beta-4> nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1 - beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -2+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                                                
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass                            
                                    
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4+beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        if (len(llst) == 0) & (edit_pos+alpha+beta-7> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1 - beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += 2
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -7+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                    elif k == len(ref)-120 + 7:      
                        Syn_codon1 = ref2[edit_pos+alpha+beta-5:edit_pos+alpha+beta-2]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-8:edit_pos+alpha+beta-5]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+alpha+beta-5> nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2 - beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass                            
                                    
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5+beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            

                        if (len(llst) == 0) & (edit_pos+alpha+beta-8> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2 - beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += 2
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -8+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                    elif k == len(ref)-120 + 8:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-6:edit_pos+alpha+beta-3]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-9:edit_pos+alpha+beta-6]
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = [] ##
                        
                        if (len(lst) != 0) & (edit_pos+alpha+beta-6> nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 3 - beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass

                                else:
                                    a += 0
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -6+beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                
                                llst.append(1)
            
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        if (len(llst) == 0) & (edit_pos+alpha+beta-9> nick_loc): 
                        
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 3 - beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -7+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += 2
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -9+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass    
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                    elif k == len(ref)-120 + 9:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-7:edit_pos+alpha+beta-4]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-10:edit_pos+alpha+beta-7]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst =[]
                        

                        if (len(lst) != 0) & (edit_pos+alpha+beta-7> nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 4 - beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass                            
                                    
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -7+beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        if (len(llst) == 0) & (edit_pos+alpha+beta-10> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 4 - beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -8+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += 2
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -10+beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                elif edit_rem == 1:
                    Syn_codon_loc1 = edit_codon_loc - 1
                    edit_pos = edit_pos + alpha  ##
                    Syn_codon1 = ref2[edit_pos-3:edit_pos]
                    Syn_codon2 = ref2[edit_pos-6:edit_pos-3]
                    Syn_codon3 = ref2[edit_pos+3:edit_pos+6]
                    Syn_codon4 = ref2[edit_pos+6:edit_pos+9]
                    edit_pos = edit_pos-alpha


                    if 2< edit_codon_loc < len(cds)/3 -1 : 
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-3 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                        
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                            
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst = []
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                  
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:

                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                        
                                        llst.append(1)
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +3> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)
                            llst = []

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                        
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos + 3 ] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq       
                                    llst.append(1) 
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                        elif (len(llst) == 0) & (edit_pos +6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon4[0]:
                                if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                else:

                                    a += -5
                                    r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos + 8] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                            
                                    
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq        
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3
                        
                    elif edit_codon_loc == 2:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-3 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                          
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                if alpha>=1 : pass
                                else:

                                    a += 2
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                              
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용


                        elif (len(llst) == 0) & (edit_pos +3> nick_loc):

                            lst = synony_dic.synony_code(Syn_codon3)
                            llst = []

                            if len(lst) != 0:  #CODON3
        
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                              
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                              
                                    llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a


                        elif (len(llst) == 0) & (edit_pos +6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)

                                #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+8] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                              
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+6 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3
                                                
                    elif edit_codon_loc == 1:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst= []

                        if (len(lst) != 0) & (edit_pos+3 > nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 5 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                          
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += -3
                                r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq     
                                llst.append(1)
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        elif (len(llst) == 0) & (edit_pos +6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+8 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos+ 6] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif edit_codon_loc == len(cds)/3-1:
                    
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-3 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst = []
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                    llst.append(1)                           
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:

                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos - 6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq         
                                        llst.append(1)                          
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +3> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += -5
                                        r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos+ 5 ] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                       
                                        
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos+3 ] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a
                        
                    elif edit_codon_loc == len(cds)/3:
                    
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-3 > nick_loc): #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3 ] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -6> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos-4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                               
                                    
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:

                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos-6 ] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                   
                                        
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                elif edit_rem == 2:
                    
                    Syn_codon_loc1 = edit_codon_loc - 1
                    edit_pos = edit_pos + alpha
                    Syn_codon1 = ref2[edit_pos-4:edit_pos-1]
                    Syn_codon2 = ref2[edit_pos-7:edit_pos-4]
                    Syn_codon3 = ref2[edit_pos+2:edit_pos+5]
                    Syn_codon4 = ref2[edit_pos+5:edit_pos+8]
                    edit_pos = edit_pos-alpha

                    if 2< edit_codon_loc < len(cds)/3 -1: 
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-4 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)

                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)

                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -7> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst= []
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)                           
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq      
                                        llst.append(1)                            
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +2> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)
                            llst= []

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq    
                                    llst.append(1)        

                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                        elif (len(llst) == 0) & (edit_pos +5> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon4[0]:
                                if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                else:
                                    a += -5
                                    r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +7] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                      
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif edit_codon_loc ==2: 
                        
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-4 > nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)

                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 2
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)                           
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos +2> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)
                            llst = []

                            if len(lst) != 0:  #CODON3

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  

                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a


                        elif (len(llst) == 0) & (edit_pos +5> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +7] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                              
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3      

                    elif edit_codon_loc ==1: 
                    
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+2 > nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq       
                                llst.append(1)                   
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += -3
                                r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        elif (len(llst) == 0) & (edit_pos +5> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +7] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   

                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a -3

                    elif edit_codon_loc == len(cds)/3-1:
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-4 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq       
                                llst.append(1)                  
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq     
                                llst.append(1)                   
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -7> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq      
                                    llst.append(1)                         
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq      
                                        llst.append(1)                          
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +2> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                            
                                        a += -5
                                        r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +4] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                       
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                    elif edit_codon_loc == len(cds)/3:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-4 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq     
                                llst.append(1)                      
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    
                                else: pass                        
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq       
                                llst.append(1)                 
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -7> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                             
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                              
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass  

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                elif edit_rem == 0:
                    edit_pos = edit_pos + alpha    ##
                    Syn_codon_loc1 = edit_codon_loc - 1

                    Syn_codon1 = ref2[edit_pos-5:edit_pos-2]
                    Syn_codon2 = ref2[edit_pos-8:edit_pos-5]
                    Syn_codon3 = ref2[edit_pos+1:edit_pos+4]
                    Syn_codon4 = ref2[edit_pos+4:edit_pos+7]
                    edit_pos = edit_pos-alpha

                    if 1<edit_codon_loc <len(cds)/3 -1: 
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-5 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq    
                                llst.append(1)                       
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq     
                                llst.append(1)                   
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -8> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst =[]
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq       
                                    llst.append(1)                        
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq        
                                        llst.append(1)                        
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +1> nick_loc):
                            llst = []
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq     
                                    llst.append(1)                              
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                        elif (len(llst) == 0) & (edit_pos +4> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            
                            #CODON4
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon4[0]:
                                if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                else:
                                    a += -5
                                    r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                       
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif edit_codon_loc ==2: 
                        
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-5 > nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq      
                                llst.append(1)                     
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                if alpha>=1:pass
                                else:
                                    a += 2
                                    r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq     
                                    llst.append(1)                       

                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                    else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos +1> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)
                            llst =[]

                            if len(lst) != 0:  #CODON3

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon3[0]:
                                    a += -5
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq      
                                    llst.append(1)                         
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                        elif (len(llst) == 0) & (edit_pos+4> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                               
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3
                        
                    elif edit_codon_loc ==1: 
                    
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+1 > nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon3[0]:
                                a += -5
                                r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq 
                                llst.append(1)                          
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += -3
                                r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a


                        elif (len(llst) == 0) & (edit_pos +4> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon4[0]:
                                a += -5
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq                           
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += -3
                                r_lst[int(rha-3 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                            df.loc[i,"RHA_len"] = rha + a -3

                    elif edit_codon_loc == len(cds)/3 -1:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-5 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -8> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst = []
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                              
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        llst.append(1)
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass  
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)

                        elif (len(llst) == 0) & (edit_pos +1> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += -5
                                        r_lst[int(rha + a)] = nu_dic[modi_codon[-1]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +3] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq                                       
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += -3
                                    r_lst[int(rha + a)] = nu_dic[modi_codon[0]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha + a

                    elif edit_codon_loc == len(cds)/3:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-5 > nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[-1]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)

                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += 2
                                r_lst[int(rha +1 + a)] = nu_dic[modi_codon[0]].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(reverse_complement(ref[rtpbs_terminal-3:rtpbs_terminal]))
                                else: pass
                                
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용

                        elif (len(llst) == 0) & (edit_pos -8> nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[-1]].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    if pbs_len+lha<11 :
                                        r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                    else: pass                                   
                                    
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (Syn_codon_loc1 ==2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha +1 + 3 + a)] = nu_dic[modi_codon[0]].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        if pbs_len+lha<11 :
                                            r_lst.append(reverse_complement(ref[rtpbs_terminal-6:rtpbs_terminal]))
                                        else: pass       
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)



            elif strand == "R":  
                rtpbs_terminal = int(edit_pos + lha + pbs_len)
                
                edit_codon_loc = math.ceil((edit_pos+1  + alpha - 60)/3) 
                edit_rem = (edit_pos+1 - 60 + alpha) %3   # OOO  1 2 0 순서
                nick_loc = pam_loc + 6
                ## Splice
                if edit_pos-55 in [0,1,2,3,4]:  #5'은 codon3,4 only
                    k = edit_pos - 55
                    if k == 0:
                        Syn_codon3 = ref2[edit_pos+5:edit_pos+8]
                        Syn_codon4 = ref2[edit_pos+8:edit_pos+11]
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+7< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 - alpha   

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 7
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +7 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 5
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5 - alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass                            

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        if (len(llst) == 0) & (edit_pos+10< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 - alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += 10
                                r_lst[int(rha-3 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +10 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                               
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += 8
                                r_lst[int(rha-3 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +8 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                               
                            

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif k ==1:
                        Syn_codon3 = ref2[edit_pos+4:edit_pos+7]
                        Syn_codon4 = ref2[edit_pos+7:edit_pos+10]

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst=[]
                    
                        if (len(lst) != 0) & (edit_pos+6< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1 - alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 7
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 5
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4 - alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass                               
                                    

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        if (len(llst) == 0) & (edit_pos+9< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1 - alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += 10
                                r_lst[int(rha-3 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +9 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += 8
                                r_lst[int(rha-3 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +7 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif k ==2:
                        Syn_codon3 = ref2[edit_pos+3:edit_pos+6]
                        Syn_codon4 = ref2[edit_pos+6:edit_pos+9] 

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+5< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2 - alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 7
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 5
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +3 - alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        if (len(llst) == 0) & (edit_pos+8< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2 - alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += 10
                                r_lst[int(rha-3 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +8 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += 8
                                r_lst[int(rha-3 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                if pbs_len+lha<11 :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif k ==3:
                        Syn_codon3 = ref2[edit_pos+2:edit_pos+5]
                        Syn_codon4 = ref2[edit_pos+5:edit_pos+8]


                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+4< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -3 - alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 7
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4 - alpha] = syn_nu # synony site #### 여기 alpha 빼줘야지..
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 5
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +2 - alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        if (len(llst) == 0) & (edit_pos+7< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -3 - alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += 10
                                r_lst[int(rha-3 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +7  - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                           
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += 8
                                r_lst[int(rha-3 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+ 7])
                                else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif k ==4:
                        Syn_codon3 = ref2[edit_pos+1:edit_pos+4]
                        Syn_codon4 = ref2[edit_pos+4:edit_pos+7]

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []
                        

                        if (len(lst) != 0) & (edit_pos+3< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -4 - alpha

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 7
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq 
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                if alpha >=1: pass
                                else:
                                    a += 5
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +1 - alpha] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        if (len(llst) == 0) & (edit_pos+6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            #CODON4

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -4 - alpha

                            if modi_codon[0] == Syn_codon4[0]:
                                a += 10
                                r_lst[int(rha-3 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +6 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass                           
                                
                            elif modi_codon[-1] == Syn_codon4[-1]:
                                a += 8
                                r_lst[int(rha-3 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4 - alpha] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                else: pass   
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 


                elif edit_pos-55 in [len(ref)-120 + i for i in [5,6,7,8,9]]: # 3'은 codon 1,2로만..
                    k = edit_pos - 55

                    if k == len(ref)-120 + 5:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-3:edit_pos+alpha+beta]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-6:edit_pos+alpha+beta-3]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+alpha+beta-1< nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0 + beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -1+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3+ beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha-1+a
                            
                        if (len(llst) == 0) & (edit_pos+alpha+beta-4< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0 + beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq 
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += -2
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha-1 -3 +a
                        
                    elif k == len(ref)-120 + 6:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-4:edit_pos+alpha+beta-1]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-7:edit_pos+alpha+beta-4]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+alpha+beta-2< nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1 + beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -2 + beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4 + beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq

                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha-1+a
                        if (len(llst) == 0) & (edit_pos+alpha+beta-5< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1 + beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += -2
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -7 + beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq


                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha-1 -3 +a

                    elif k == len(ref)-120 + 7:      
                        Syn_codon1 = ref2[edit_pos+alpha+beta-5:edit_pos+alpha+beta-2]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-8:edit_pos+alpha+beta-5]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                        

                        if (len(lst) != 0) & (edit_pos+alpha+beta-3< nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2 + beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5+ beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha-1+a
                        if (len(llst) == 0) & (edit_pos+alpha+beta-6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2 + beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += -2
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -8+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha-1 -3 +a

                    elif k == len(ref)-120 + 8:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-6:edit_pos+alpha+beta-3]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-9:edit_pos+alpha+beta-6]
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+alpha+beta-4< nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -3 + beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq 
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -6+ beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha-1+a
                        if (len(llst) == 0) & (edit_pos+alpha+beta-7< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -3 + beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -7+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += -2
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -9+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha-1 -3 +a

                    elif k == len(ref)-120 + 9:
                        Syn_codon1 = ref2[edit_pos+alpha+beta-7:edit_pos+alpha+beta-4]
                        Syn_codon2 = ref2[edit_pos+alpha+beta-10:edit_pos+alpha+beta-7]

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []
                    

                        if (len(lst) != 0) & (edit_pos+alpha+beta-5< nick_loc):  #CODON1

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -4 + beta

                            if modi_codon[0] == Syn_codon1[0]:
                                if beta >=1: pass
                                else:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -7+ beta] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq

                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha-1+a
                        if (len(llst) == 0) & (edit_pos+alpha+beta-8< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2

                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -4 + beta

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -8+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    a += -2
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -10+ beta] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha-1 -3 +a


                elif edit_rem ==1:
                    Syn_codon_loc1 = edit_codon_loc - 1
                    edit_pos = edit_pos + alpha
                    Syn_codon1 = ref2[edit_pos-3:edit_pos]
                    Syn_codon2 = ref2[edit_pos-6:edit_pos-3]
                    Syn_codon3 = ref2[edit_pos+3:edit_pos+6]
                    Syn_codon4 = ref2[edit_pos+6:edit_pos+9]
                    edit_pos = edit_pos-alpha

                    if 1< edit_codon_loc < len(cds)/3 -1: 

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst= []

                        if (len(lst) != 0) & (edit_pos+5< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 3

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+8< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst=  []

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 6

                                if modi_codon[0] == Syn_codon4[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq   
                                        llst.append(1)
                                        if pbs_len+lha<11  :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                        else: pass                                
                                        
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 


                        elif (len(llst) == 0) & (edit_pos-1< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon1)
                            llst= []

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    a += -2
                                    r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                                df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-4< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq


                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a 

                    elif edit_codon_loc ==2:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst =[]

                        if (len(lst) != 0) & (edit_pos+5< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 3

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+8< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst = []

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 6

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +8] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos-1< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon1)

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq   
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    if alpha>=1: pass
                                    else:
                                        a += -2
                                        r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -3] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha -1 + a

                    elif edit_codon_loc ==1:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+5< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 3

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq   
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+8< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 6

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +8] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+ 7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3 -1:
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-1< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq 
                                llst.append(1) 
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-4< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst = []
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a


                        elif (len(llst) == 0) & (edit_pos+5< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 3

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +5] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq  
                                        if lha<3 :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                        else: pass                                    
                                        
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if lha<3 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst= []

                        if (len(lst) != 0) & (edit_pos-1< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 0

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-4< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 0

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a


                elif edit_rem ==2:
                    Syn_codon_loc1 = edit_codon_loc - 1
                    edit_pos = edit_pos + alpha
                    Syn_codon1 = ref2[edit_pos-4:edit_pos-1]
                    Syn_codon2 = ref2[edit_pos-7:edit_pos-4]
                    Syn_codon3 = ref2[edit_pos+2:edit_pos+5]
                    Syn_codon4 = ref2[edit_pos+5:edit_pos+8]
                    edit_pos = edit_pos-alpha
                    
                    if 1< edit_codon_loc < len(cds)/3 -1: 
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst= []

                        if (len(lst) != 0) & (edit_pos+4< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+7< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst = []

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 5

                                if modi_codon[0] == Syn_codon4[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        llst.append(1)
                                        if pbs_len+lha<11  :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                        else: pass                                
                                        
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 


                        elif (len(llst) == 0) & (edit_pos-2< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon1)
                            llst = []

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    a += -2
                                    r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                                df.loc[i,"RHA_len"] = rha -1 + a 

                        elif (len(llst) == 0) & (edit_pos-5< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a

                    elif edit_codon_loc ==2:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst=[]

                        if (len(lst) != 0) & (edit_pos+4< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+7< nick_loc): 
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst = []

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 5

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +7] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos-2< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon1)

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    if alpha>=1: pass
                                    else:
                                        a += -2
                                        r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -4] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha -1 + a

                    elif edit_codon_loc ==1:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst =[]

                        if (len(lst) != 0) & (edit_pos+4< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 2

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+7< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 5

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +7] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3 -1:
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst =[]

                        if (len(lst) != 0) & (edit_pos-2< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-5< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst = []
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a


                        elif (len(llst) == 0) & (edit_pos+4< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 2

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +4] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        if lha<3 :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                        else: pass                                    
                                                    
                                        
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +2] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq                                
                                    if lha<3 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3:
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst = []

                        if (len(lst) != 0) & (edit_pos-2< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -1

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -2] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -4] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-5< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -1

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -7] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a


                elif edit_rem ==0:
                    Syn_codon_loc1 = edit_codon_loc - 1
                    edit_pos = edit_pos + alpha
                    Syn_codon1 = ref2[edit_pos-5:edit_pos-2]
                    Syn_codon2 = ref2[edit_pos-8:edit_pos-5]
                    Syn_codon3 = ref2[edit_pos+1:edit_pos+4]
                    Syn_codon4 = ref2[edit_pos+4:edit_pos+7]
                    edit_pos = edit_pos-alpha


                    if 1< edit_codon_loc < len(cds)/3 -1:
                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+3< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst = []
                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 4

                                if modi_codon[0] == Syn_codon4[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +6] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq  
                                        llst.append(1)
                                        if pbs_len+lha<11 :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                        else: pass                                
                                        
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 


                        elif (len(llst) == 0) & (edit_pos-3< nick_loc):
                            llst = []
                            lst = synony_dic.synony_code(Syn_codon1)

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    a += -2
                                    r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -5] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                                df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a

                    elif edit_codon_loc ==2:

                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+3< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)

                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 


                        elif (len(llst) == 0) & (edit_pos+6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)
                            llst = []

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 4

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    llst.append(1)
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                        

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos-3< nick_loc):

                            lst = synony_dic.synony_code(Syn_codon1)

                            if len(lst) != 0:  #CODON1
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2

                                if modi_codon[0] == Syn_codon1[0]:
                                    a += 0
                                    r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -3] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                elif modi_codon[-1] == Syn_codon1[-1]:
                                    if alpha>=1: pass
                                    else:
                                        a += -2
                                        r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -5] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq


                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 
                                df.loc[i,"RHA_len"] = rha -1 + a

                    elif edit_codon_loc ==1:


                        lst = synony_dic.synony_code(Syn_codon3)
                        llst = []

                        if (len(lst) != 0) & (edit_pos+3< nick_loc):  #CODON3

                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = 1

                            if modi_codon[0] == Syn_codon3[0]:
                                a += 2
                                r_lst[int(rha + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass                        
                                
                            elif modi_codon[-1] == Syn_codon3[-1]:
                                a += 0
                                r_lst[int(rha + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos +1] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)
                                if pbs_len+lha<11  :
                                    r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                else: pass
                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                        elif (len(llst) == 0) & (edit_pos+6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon4)

                            if len(lst) != 0:  #CODON4
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 4

                                if modi_codon[0] == Syn_codon4[0]:
                                    a += 2
                                    r_lst[int(rha + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            
                                    
                                elif modi_codon[-1] == Syn_codon4[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +4] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if pbs_len+lha<11  :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+7])
                                    else: pass                            

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3 -1:
                        lst = synony_dic.synony_code(Syn_codon1)
                        llst =[]

                        if (len(lst) != 0) & (edit_pos-3< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            llst =[]
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                    llst.append(1)
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq
                                        llst.append(1)

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a


                        elif (len(llst) == 0) & (edit_pos+3< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon3)

                            if len(lst) != 0:  #CODON3
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = 1

                                if modi_codon[0] == Syn_codon3[0]:
                                    if (beta >=1)&(edit_codon_loc == len(cds)/3-2): pass
                                    else:
                                        a += 2
                                        r_lst[int(rha + a)] = modi_codon[-1].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[-1]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos +3] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq  
                                        if lha<3 :
                                            r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                        else: pass                                    
                                        
                                elif modi_codon[-1] == Syn_codon3[-1]:
                                    a += 0
                                    r_lst[int(rha + a)] = modi_codon[0].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[0]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos +1] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq
                                    if lha<3 :
                                        r_lst.append(ref[rtpbs_terminal+1:rtpbs_terminal+4])
                                    else: pass
                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) 

                    elif edit_codon_loc == len(cds)/3:

                        lst = synony_dic.synony_code(Syn_codon1)
                        llst =[]

                        if (len(lst) != 0) & (edit_pos-3< nick_loc):  #CODON1
        
                            modi_codon = random.choice(lst) 
                            r_lst = list(rtpbs)
                            a = -2

                            if modi_codon[0] == Syn_codon1[0]:
                                a += 0
                                r_lst[int(rha -1 + a)] = modi_codon[-1].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[-1]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -3] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq  
                                llst.append(1)
                            elif modi_codon[-1] == Syn_codon1[-1]:
                                a += -2
                                r_lst[int(rha -1 + a)] = modi_codon[0].lower()
                                ref_lst = list(ref)
                                syn_nu = modi_codon[0]
                                ref_lst[edit_pos] = nu  # original
                                ref_lst[edit_pos -5] = syn_nu # synony site
                                seq = ''.join(ref_lst)
                                df.loc[i,"seq"] = seq
                                llst.append(1)

                            synonyrtpbs = ''.join(r_lst)
                            df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                            df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon) # 나중 검증용
                            df.loc[i,"RHA_len"] = rha -1 + a

                        elif (len(llst) == 0) & (edit_pos-6< nick_loc):
                            lst = synony_dic.synony_code(Syn_codon2)
                            
                            if len(lst) != 0:  #CODON2
            
                                modi_codon = random.choice(lst) 
                                r_lst = list(rtpbs)
                                a = -2

                                if modi_codon[0] == Syn_codon2[0]:
                                    a += 0
                                    r_lst[int(rha -1 -3 + a)] = modi_codon[-1].lower()
                                    ref_lst = list(ref)
                                    syn_nu = modi_codon[-1]
                                    ref_lst[edit_pos] = nu  # original
                                    ref_lst[edit_pos -6] = syn_nu # synony site
                                    seq = ''.join(ref_lst)
                                    df.loc[i,"seq"] = seq  
                                elif modi_codon[-1] == Syn_codon2[-1]:
                                    if (alpha >=1) & (edit_codon_loc ==2): pass
                                    else:

                                        a += -2
                                        r_lst[int(rha -1 -3 + a)] = modi_codon[0].lower()
                                        ref_lst = list(ref)
                                        syn_nu = modi_codon[0]
                                        ref_lst[edit_pos] = nu  # original
                                        ref_lst[edit_pos -8] = syn_nu # synony site
                                        seq = ''.join(ref_lst)
                                        df.loc[i,"seq"] = seq

                                synonyrtpbs = ''.join(r_lst)
                                df.loc[i,"SynonyRTPBS"] = synonyrtpbs
                                df.loc[i,"SynonyCodon"] = Endo_REF.gene_code(modi_codon)
                                df.loc[i,"RHA_len"] = rha -1 -3 + a

        except IndexError: 
            print("list index out of range")
        
   
    return df

