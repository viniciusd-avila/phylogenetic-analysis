import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo
from sklearn.cluster import DBSCAN
import itertools
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans

aminoacids = ['A','C','G','T']

#['AAA','AAC','AGT',...]
codons = list(map((lambda element: ''.join(list(element))), itertools.product(aminoacids,repeat=3)))

def clustalw_alignment(filename):
    cline = ClustalwCommandline("clustalw2", infile = filename)
    cline()
    return AlignIO.read(filename + ".aln","clustal")

def pca_2d(df):
    #input: codon frequency dataframe
    #output: 2d numpay array of principal components
    X = df.loc[:, codons].values
    X = StandardScaler().fit_transform(X)
    y = df.loc[:, codons].values
    pca = PCA(n_components=2)
    return pca.fit_transform(X)

def t_SNE(df):
    #input: codon frequency dataframe
    #output: 2d numpay array of t-SNE distribution
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300) 
    return tsne.fit_transform(df.loc[:,codons].values)

def kmeans_plot(k,X):
    #input: k and a 2d-array
    #output: K-Means cluster plot 
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(X)
    y_means = kmeans.predict(X)

    plt.scatter(X[:,0],X[:,1],c=y_means,s=50,cmap='viridis')

    centers = kmeans.cluster_centers_
    plt.scatter(centers[:,0],centers[:,1],c='black',s=200,alpha=0.5)

def dbscan_plot(X):
    #input: 2d-array
    #output: DBScan cluster plot, number of clusters

    db = DBSCAN(eps=0.3, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool) 
    core_samples_mask[db.core_sample_indices_] = True 
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    unique_labels = set(labels) 
    colors = [plt.cm.Spectral(each)           
        for each in np.linspace(0, 1, len(unique_labels))] 
    for k, col in zip(unique_labels, colors):     
        if k == -1:         
        # Black used for noise.         
            col = [0, 0, 0, 1]      
            
    class_member_mask = (labels == k)      
    xy = X[class_member_mask & core_samples_mask]     
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                        markeredgecolor='k', markersize=14)      
    xy = X[class_member_mask & ~core_samples_mask]     
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col),
                    markeredgecolor='k', markersize=6)  
    
    plt.title('Estimated number of clusters: %d' % n_clusters_) 
    plt.show()

def clean_data(infile):
    #input: fasta file of DNA fragments from multiple species
    #output: dictionary of fragments for each species
    with open(infile,'r') as myfile:
        data=myfile.read()

    data = data.split('>')
    dic = {}

    for i in range(len(data)):
        data[i] = data[i].split('\n',maxsplit=1)
        if len(data[i])>1:
            dic[data[i][0]] = data[i][1].replace('\n','')
    return dic

def gen_species(dic,dim_red=False):
    #input: dictionary of fragments of multiple species
    #output: list of Species object instances
    species_list = []  
    for key in dic:     
        x = Species(key.split(' ')[1] + ' ' + key.split(' ')[2], key.split(' ')[0],dic[key],dim_red)     
        species_list.append(x)
    return species_list

def species_name_list(species_list):
    #input: list of Species object instances
    #output: list of their names
    res = []
    for species in species_list:
        res.append(species.name)
    return res

def tree_nodes_names(dnd_file,species_list):
    #input: tree file in dnd format, species list
    #output: swaps accession number in the tree with the species scientific name
    dic = {}
    for species in species_list:
        dic[species.accession_number] = species.name.replace(' ','-')
    
    with open(dnd_file,'r') as myfile:
        text=myfile.read()

    for i,j in dic.items():
        text = text.replace(i,j)
    return text


def attr_alignment(species_list,alignment):
    #input: list of species object instances; clustal results
    #output: None, changes the object instances, adds attribute of aligned sequence
    for species in species_list:
        for recs in alignment:
            if species.accession_number == recs.id:
                species.aligned_seq = str(recs.seq)

def lavenshtein_dist(stringA,stringB):
    #input: two strings
    #output: lavenshtein (edit) distance between them
    m = len(stringA)
    n = len(stringB)
    d = np.zeros(shape=(m,n))
    d[0,:] = np.arange(1,n+1)
    d[:,0] = np.arange(1,m+1)
    for j in range(1,n):
        for i in range(1,m):
            if stringA[i] == stringB[j]:
                d[i, j] = d[i-1, j-1]
            else:
                d[i,j] = min(d[i-1,j] + 1,d[i,j-1],d[i-1,j-1] + 1)
    return d[m-1,n-1] - 1

def laven_dist_species(species_list,alignment):
    #input: list of species object instance, alignment sequences
    #output: matrix (numpy 2d array) with the lavenshtein distances between each species
    attr_alignment(species_list,alignment)
    length = len(species_list)

    dm = np.empty(shape=(length,length))
    dm[:] = np.nan
    for i in range(length):
        for j in range(length):
            if i == j:
                dm[i,j] = 0.0
            elif np.isnan(dm[i,j]):
                dm[i,j] = dm[j,i] = lavenshtein_dist(species_list[i].aligned_seq,species_list[j].aligned_seq) 
    pd.DataFrame(dm).to_csv("lavenshtein_distance_species.csv")
    return dm

def codon_freq_multi_species(species_list):
    length = len(species_list)
    
    df = pd.DataFrame()
    
    label = 0
    for species in species_list:
        df_tmp = species.codon_freq_df.copy()
        y = np.empty(len(df_tmp.index))
        y[:] = label
        df_tmp['label'] = y
        df_tmp.reset_index(drop=True,inplace=True)
        df.reset_index(drop=True,inplace=True)
        df = pd.concat([df,df_tmp],axis=0)
        label = label + 1
    df.reset_index(drop=True,inplace=True)
    return df
        

class Species():
    def __init__(self,name,accession_number,sequence,red_dimension=False):
        self.name = name
        self.accession_number = accession_number
        self.sequence = sequence
        self.split_seq = self.split_sequence(24)
        self.codon_freq_df = self.codon_freq()
        if red_dimension:
            self.principal_components = pca_2d(self.codon_freq_df)
            self.tsne_results = t_SNE(self.codon_freq_df)

    def split_sequence(self,n):
        res = [self.sequence[i:i+n] for i in range(0,len(self.sequence),n)]
        #res = list(filter(lambda x: set(x) != set('N'),res))
        return res

    def codon_freq(self):
        df = pd.DataFrame()
        df['Fragment'] = self.split_seq
        for cdn in codons:
            df[cdn] = df['Fragment'].apply(lambda frgmnt: frgmnt.count(cdn))
        return df