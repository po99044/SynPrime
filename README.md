# SynPrime
Provide SynPrime library design


The article 'Saturation resistance profiling of EGFR variants against tyrosine kinase inhibitors using prime editing' is based on this design concept. (https://www.biorxiv.org/content/10.1101/2023.12.03.569825v1)


**Principle**

1.Add additional synonymous snv adjacent to edited snv

2.additional synonymous snv is in the other codon different with edited snv

3.synonymous snv is randomly selected in the synonymous codon (i.e., Leucine : TTA, TTG, CTT, CTC, CTA, CTG -> in case of CTX, X can be T, C, A or G; that could be randomly selected)

4.Intron SNV 5bp near the CDS : synonymous snv is located in the CDS codon

5. Synonymous SNV could be inserted into the LHA, RHA, or PAM-disrupting regions

6. Suggest pegRNAs according to the DeepPrime Scores
7. 
