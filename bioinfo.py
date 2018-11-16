import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import itertools
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

lstamn = ['A','C','G','T']

amn_list = [['Adenine','A'],
 ['Cytosine','C'],
 ['Guanine','G'],
 ['Thymine','T'],
 ['Uracil','U'],
 ['Purine','A','G'],
 ['Pyrimidines','T','U'],
 ['Non-identified','N']]

plot_list = ['Adenine','Cytosine','Guanine','Thymine','Non-identified']

def defcodons(n):
    return list(map((lambda element: ''.join(list(element))), itertools.product(lstamn,repeat=n)))

def amino_acids_count(df, amn_list):
    for amn in amn_list:
        df[amn[0]] = 0.0
        for letter in amn[1:]:
            df[amn[0]] = df[amn[0]] + df['Sequence'].apply(lambda seq:seq.count(letter))
        df[amn[0]] = df[amn[0]]/df['Length']
    return df

def split_DNA(seq,n):
    return [seq[i:i+n] for i in range(0, len(seq), n)]

def gen_df2(seq,n):
    df2 = pd.DataFrame(split_DNA(seq,24))
    df2.columns = ['Fragment']
    for codon in defcodons(n):
        df2[codon] = df2['Fragment'].apply(lambda frag: float(frag.count(codon)))
    return df2

def gen_df(infile):
    with open(infile, 'r') as myfile:
        data=myfile.read()
    
    data = data.split('>')
    dic = {}

    for i in range(len(data)):
        data[i] = data[i].split('\n',maxsplit=1)
        if len(data[i]) > 1:
            dic[data[i][0]] = data[i][1].replace('\n','')
    
    df = pd.DataFrame.from_dict(dic,orient='index')
    df.columns = ['Sequence']
    df = df.reset_index()
    df['Accession number'] = df['index'].apply(lambda name: name.split(' ')[0])
    df['Species'] = df['index'].apply(lambda name: name.split(' ')[1] + ' ' + name.split(' ')[2])
    df = df[['Accession number','Species','Sequence']]
    df = df.set_index('Species')
    df.reset_index(level=0, inplace=True)
    df['Species'] = df['Species'].apply(lambda s: s.replace(' ','-'))
    df['Length'] = df['Sequence'].apply(lambda seq: len(seq))
    df = amino_acids_count(df, amn_list)
    return df 

def amino_plot(plot_list):
    return df[plot_list].plot(kind='bar', title ="V comp", figsize=(15, 10), legend=True, stacked=True,fontsize=12)

def tree_node_names(dnd_file, df):
    with open(dnd_file, 'r') as myfile:
        text=myfile.read()
    rep = df.set_index('Accession number').to_dict()['Species']
    for i, j in rep.items():
        text = text.replace(i, j)
    return text

def pca(df):
    x = df.loc[:, df.columns != 'Fragment'].values
    y = df.loc[:, df.columns == 'Fragment'].values
    x = StandardScaler().fit_transform(x)
    
    pca = PCA(n_components=2)
    
    principalComponents = pca.fit_transform(x)
    
    principalDf = pd.DataFrame(data = principalComponents, 
            columns = ['principal component 1', 'principal component 2'])
    
    return pd.concat([principalDf, df[['Fragment']]], axis = 1)

def plot_pca(finalDf):
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)
    
    ax.scatter(finalDf.loc[:, 'principal component 1'], finalDf.loc[:, 'principal component 2'])
    
    ax.grid()
