# %%
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
import Pam_disrupt
import Gen_Syn

## Example 
exo = "E5"

result = pd.DataFrame(columns = ['ID', 'PE2max_score', 'Spacer', 'RTPBS', 'PBSlen',
       'RTlen', 'RT-PBSlen', 'EDIT_POS', 'Edit_len', 'RHA_len', 'Target',
       'RTPBS_R', 'sgRNA', 'LHA', 'PBS', 'strand', 'PAM_LOC',
       'PAMdist_fromedit', 'SynonyRTPBS', 'synony pam disrupt',
       'Intended pam disrupt', 'Intend_to_pamdisrupt_dist',
       'Pam disrupt codon', 'SynonyCodon', 'seq'])



files = os.listdir("Output/%s"  % exo) 
for file in files:


    df = pd.read_csv("Output/%s/%s" % (exo, file), index_col=0)
    df = df.dropna(subset = ["synony pam disrupt", "Intended pam disrupt"], axis =0)
    sgrna_dic = {}  # spacer 다른 것도 하나 고를거고, #Top 1 + pam disrupting 1 + other spacer , # Top1이 pam disrupting이면 + Pam not disrupting + other spacer, # other spacer가 없다면 Top2 + other.
    for i in df.index: 
        sgrna = df.loc[i,"sgRNA"]
        if sgrna not in sgrna_dic:
            sgrna_dic[sgrna] = 0 
        else: pass
    # result = pd.concat([result, df.head(3)])
    
    if len(sgrna_dic) == 1 :  # spacer 1개
        
        result = pd.concat([result, df.head(3)])


    else:

        d2f = df[df["sgRNA"] == list(sgrna_dic.keys())[1]]
        
        result = pd.concat([result, df.head(2), d2f.head(1)])
    


result.to_csv("SCNA_Library/%s.csv"  %exo)
