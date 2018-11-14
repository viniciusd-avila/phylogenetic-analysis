import pandas as pd
import numpy as np

#Nossos datapoints serão a contagem das 'palavras' (codons) de 1, 2, 3 ou 4 letras, em cada um dos fragmentos com FragSize letras do DNA.

#Data points. Create a new matrix with same size fragments of a DNAseq
def DNA_frag(DNAseq, FragSize):
    pos = 0
    frags = []
    for i in range(1, int(len(DNAseq)/Frag_Size)+2): #1+1=2, 1 para arredondar para cima e +1 para o range alcançar toda sequência
        new_frag = DNAseq[pos:i*Frag_Size]
        frags.append(new_frag)
        i += 1
        pos = pos + Frag_Size
    np.transpose(frags)
    return frags

#Usaremos a função para quebrar nossa sequencia de DNA em fragmentos de tamanho 300 e contaremos a frequência de codons usando a função Codon_Freq

def Codon_Freq(Frags, Codon_Size):
    Amino_acids = 'ACTG' 
    Codons = [] #criar lista de códons daquele tamanho . (Ex: size 2 : aa, ac, at, ag, ca, cc, ct, cg,ta, tc, tg, tt...)
    
    #criar matriz com frequencia dos codons de tamanho Codon_Size
    
    
    
    