# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 15:26:32 2015

@author: Glaucia
"""
#####Sequence Manipulation Exercise#######
#1.
dnaSeq="aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"
print(dnaSeq)
#2.
#Number of nucleotids in the sequence
print(len(dnaSeq))#I already took of the last two bases
#3.
#RNA equivalent
rna = dnaSeq.replace("t","u")#replacing timine with uracile to get RNA
print(rna)
#4.
#Reverse complement of my Sequence
####Using uppercase to avoid unintended substitutions.
revComplement1 = dnaSeq.replace("t","A")
revComplement2 = revComplement1.replace("a","T")
revComplement3 = revComplement2.replace("c","G")
revComplement4 = revComplement3.replace("g","C")
####Printing the reverse complement of my Sequence (Final result):
revComplement_FINAL=revComplement4[::-1]#To reverse the complement
print(revComplement_FINAL)
#5.
#Extracting the bases correponding to the 13rd and 14th codons:
print(dnaSeq[36:42])#indexing the bases position
####13rd
print(dnaSeq[36:39])
####14th
print(dnaSeq[39:42])
#6.
#Creating a function to translate the nucleotide sequence to amino acids 
def DNAtranslator(dnaSeq):
    dnaSeq_codons=[dnaSeq[x:x+3] for x in range(0,len(dnaSeq),3)]#transforming the sequence in a list (splicing at each 3)
    for codon in dnaSeq_codons:#defining the loop
        if codon in ["ttt","ttc"]:#defining if statements for each codon and aminoacid
            print("F"),
        elif codon in ["tta","ttg","ctt","ctc","cta","ctg"]:
            print("L"),
        elif codon in ["att","atc"]:
            print("I"),
        elif codon in ["ata","atg"]:
            print("M"),
        elif codon in ["gtt","gtc","gta","gtg"]:
            print("V"),
        elif codon in ["tct","tcc","tca","tcg"]:
            print("S"),
        elif codon in ["cct","ccc","cca","ccg"]:
            print("P"),
        elif codon in ["act","acc","aca","acg"]:
            print("T"),
        elif codon in ["gct","gcc","gca","gcg"]:
            print("A"),
        elif codon in ["tat","tac"]:
            print("Y"),
        elif codon in ["taa","tag","aga","agg"]:
            print("*"),
        elif codon in ["cat","cac"]:
            print("H"),
        elif codon in ["caa","cag"]:
            print("Q"),
        elif codon in ["aat","aac"]:
            print("N"),
        elif codon in ["aaa","aag"]:
            print("K"),
        elif codon in ["gat","gac"]:
            print("D"),
        elif codon in ["gaa","gag"]:
            print("E"),
        elif codon in ["tgt","tgc"]:
            print("C"),
        elif codon in ["tga","tgg"]:
            print("W"),
        elif codon in ["cgt","cgc","cga","cgg"]:
            print("R"),
        elif codon in ["agt","agc"]:
            print("S"),
        elif codon in ["ggt","ggc","gga","ggg"]:
            print("G"),
        else:
            print("error")#to detect possible errors in my function
#7.
#Using the function to translate my DNA sequence            
amino_acids=DNAtranslator(dnaSeq)



