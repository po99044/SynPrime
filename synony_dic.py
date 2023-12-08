
def synony_code(c):
    
    dic = {

        "TTT":["TTC"],
        "TTC":["TTT"], 

        "TTA": ["TTG","CTA"], 
        "TTG": ["TTA","CTG"], 
        "CTT": ["CTC", "CTA","CTG"], 
        "CTC": ["CTT", "CTA","CTG"], 
        "CTA": ["TTA","CTT","CTC","CTG"], 
        "CTG": ["TTG","CTT","CTC","CTA"], 

        "TCT": ["TCC","TCA","TCG"], 
        "TCC": ["TCT","TCA","TCG"],
        "TCA": ["TCT","TCC","TCG"],
        "TCG": ["TCT","TCC","TCA"],
        "AGT": ["AGC"],
        "AGC": ["AGT"],
        
        "CCT":["CCC","CCA","CCG"], 
        "CCC":["CCT","CCA","CCG"], 
        "CCA":["CCT","CCC","CCG"], 
        "CCG":["CCT","CCC","CCA"],

        "TAT" : ["TAC"], 
        "TAC" : ["TAC"], 

        "TAA":["TAG","TGA"],
        "TAG":["TAA"],
        "TGA":["TAA"],

        "CAT":["CAC"],
        "CAC":["CAT"], 
 
        "CAA": ["CAG"], 
        "CAG": ["CAA"], 

        "TGT":["TGC"],
        "TGC":["TGT"],

        "TGG":[], # W

        "CGT":["CGC","CGA","CGG"],
        "CGC":["CGT","CGA","CGG"], 
        "CGA":["CGT","CGC","CGG","AGA"], 
        "CGG":["CGT","CGC","CGA","AGG"], 
        "AGA":["CGA","AGG"], 
        "AGG":["CGG","AGA"],  

        "ATT" : ["ATC","ATA"],
        "ATC" : ["ATT","ATA"],
        "ATA" : ["ATT","ATC"],

        "ATG":[],  # M
        "GTT":["GTC","GTA","GTG"], 
        "GTC":["GTT","GTA","GTG"], 
        "GTA":["GTT","GTC","GTG"], 
        "GTG":["GTT","GTC","GTA"], 


        "ACT": ["ACC","ACA","ACG"],
        "ACC": ["ACT","ACA","ACG"],
        "ACA": ["ACT","ACC","ACG"],
        "ACG": ["ACT","ACC","ACA"],


        "GCT":["GCC","GCA","GCG"], 
        "GCC":["GCT","GCA","GCG"], 
        "GCA":["GCT","GCC","GCG"], 
        "GCG":["GCT","GCC","GCA"], 
        
        "AAT":["AAC"], 
        "AAC":["AAT"], 

        "AAA" : ["AAG"],
        "AAG" : ["AAA"],


        "GAT":["GAC"], 
        "GAC":["GAT"], 
        
        "GAA":["GAG"], 
        "GAG":["GAA"], 
        
        "GGT":["GGC","GGA","GGG"],
        "GGC":["GGT","GGA","GGG"],
        "GGA":["GGT","GGC","GGG"],
        "GGG":["GGT","GGC","GGA"]

        

    }

    return dic[str(c)]