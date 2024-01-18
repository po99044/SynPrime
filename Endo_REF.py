

## Genetic code

def gene_code(b):
    dic = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S","TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CTT":"L","CTC":"L","CTA":"L",
        "CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T","AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
        "GTT":"V","GTC":"V","GTA":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A","GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
    }

    return dic[b]

def synony_code(c):
    
    dic = {
        "F":["TTT","TTC",], "L": ["TTA","TTG","CTT","CTC", "CTA","CTG"], "S": ["TCT","TCC","TCA","TCG","AGT","AGC"], 
        "P":["CCT","CCC","CCA","CCG"], "Y" : ["TAT","TAC"], "*":["TAA","TAG","TGA"],
        "H":["CAT","CAC"], "Q": ["CAA","CAG"], "C":["TGT","TGC"],"W":["TGG"],
        "R":["CGT","CGC","CGA","CGG","AGA","AGG"], "I" : ["ATT","ATC","ATA"],
        "M":["ATG"], "V":["GTT","GTC","GTA","GTG"], "T": ["ACT","ACC","ACA","ACG"],
        "A":["GCT","GCC","GCA","GCG"], "N":["AAT","AAC"], "K" : ["AAA","AAG"],
        "D":["GAT","GAC"], "E":["GAA","GAG"], "G":["GGT","GGC","GGA","GGG"]

        

    }

    return dic[c]

