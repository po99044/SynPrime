# PEER-seq
Provide PEER-seq library design (Old Synprime in biorxiv)


The article 'Saturation profiling of drug-resistant genetic variants using prime editing' is based on this design concept. (https://www.nature.com/articles/s41587-024-02465-z)


**Principle**

1.Add additional synonymous snv adjacent to edited snv

2.additional synonymous snv is in the other codon different with edited snv

3.synonymous snv is randomly selected in the synonymous codon (i.e., Leucine : TTA, TTG, CTT, CTC, CTA, CTG -> in case of CTX, X can be T, C, A or G; that could be randomly selected)

4.Intron SNV 5bp near the CDS : synonymous snv is located in the CDS codon

5.Synonymous SNV could be inserted into the LHA, RHA, or PAM-disrupting regions

6.Suggest pegRNAs according to the DeepPrime Scores 
