import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import itertools
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans

aminoacids = ['A','C','G','T']

codons = list(map((lambda element: ''.join(list(element))), itertools.product(aminoacids,repeat=3)))

def pca_2d(df):
    X = df.loc[:, codons].values
    X = StandardScaler().fit_transform(X)
    y = df.loc[:, codons].values
    pca = PCA(n_components=2)
    return pca.fit_transform(X)

def t_SNE(df):
    tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300) 
    return tsne.fit_transform(df.loc[:,codons].values)

def kmeans_plot(self,k,X):
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(X)
    y_means = kmeans.predict(X)

    plt.scatter(X[:,0],X[:,1],c=y_means,s=50,cmap='viridis')

    centers = kmeans.cluster_centers_
    plt.scatter(centers[:,0],centers[:,1],c='black',s=200,alpha=0.5)

def dbscan_plot(self,X):
    db = DBSCAN(eps=0.3, min_samples=10).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool) 
    core_samples_mask[db.core_sample_indices_] = True 
    labels = db.labels_
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
    with open(infile,'r') as myfile:
        data=myfile.read()

    data = data.split('>')
    dic = {}

    for i in range(len(data)):
        data[i] = data[i].split('\n',maxsplit=1)
        if len(data[i])>1:
            dic[data[i][0]] = data[i][1].replace('\n','')
    return dic

def gen_species(dic):
    species_list = []  
    for key in dic:     
        x = draft.Species(key.split(' ')[1] + ' ' + key.split(' ')[2], key.split(' ')[0],dic[key])     
        species_list.append(x)
    return species_list

class Species():
    def __init__(self,name,accession_number,sequence):
        self.name = name
        self.accession_number = accession_number
        self.sequence = sequence
        self.split_seq = self.split_sequence(24)
        self.codon_freq_df = self.codon_freq()
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
