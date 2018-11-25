# BioinfoAnalysis

## 1. Introdução

A bioinformática é uma área multdisciplinar que surgiu recentemente motivada pela utilização de ferramentas computacionais para a análise de dados genéticos,
bioquímicos e de biologia molecular. Seu principal objetivo é estudar e desvendar a grande quantidade de dados obtidos a partir de seqüências de DNA usando e desenvolvento novas técnicas computacionais.

### 1.1. Motivação
Durante a realização do Curso de Verão de Bioinformática 2018, na Universidade de São Paulo fui apresentada à área em seus diversos aspectos. Com a popularização, desenvolvimento e refinamento dos métodos de sequenciamento do DNA houve um enorme aumento da quantidade de dados biológicos disponíveis. Inspirados pelo o impacto nas descobertas causado pela utilização de métodos de Data Science nas mais diversas areas, decidimos assim explorar sequencias de DNA de diferentes espécies. 

### 1.2. Dados
Foram utilizados amostras de DNA TRIM5alpha genômico de 17 espécies. O TRIM5alpha é uma proteína fator de restrição de retrovírus que medeia a infecção precoce por bloqueio. 

Nos chamados Macacos do Velho Mundo, Cercopithecidae, o  TRIM5α foi isolado como uma proteína de macaco rhesus responsável por bloquear a infecção pelo HIV-1. Assim estes macacos não podem ser infectados com o HIV-1, o vírus que causa a AIDS em humanos.
A versão humana do TRIM5α não tem como alvo o HIV-1, mas pode inibir as cepas do vírus da leucemia murina (MLV), bem como o vírus da anemia infecciosa eqüina (EIAV).

Segue a lista das 17 espécies estudadas:

* Chlorocebus 
* Papio anubis 
* Pan troglodytes 
* Colobus guereza 
* Pygathrix nemaeus 
* Symphalangus syndactylus 
* Gorilla gorilla 
* Alouatta sara 
* Callithrix pygmaea 
* Pongo pygmaeus 
* Erythrocebus patas 
* Pithecia pithecia 
* Ateles geoffroyi 
* Saimiri sciureus 
* Saguinus labiatus 
* Callicebus donacophilus 
* Lagothrix lagotricha 


Os dados usados foram baixados do NCBI (National Center for Biotechnology Information Search database) e podem ser encontrados neste link: goo.gl/HYqnxD. Foi escolhido o PopSet usado no trabalho "Positive selection of primate TRIM5alpha identifies a critical species-specific retroviral restriction domain." [1] que inseriu no GenBank os 17 arquivos sob os números de acesso AY843504, AY843505, ... até  AY843520.

### 1.3. Pacotes e versões

Para execução do código encontrado neste repositório foram utilizados os seguintes pacotes:
* Biopython 1.7
* Matplotlib 3.0.1
* Numpy 1.15.2
* Pandas 0.23.4
* Pygraphviz 1.5
* CLUSTAL 2.1 Multiple Sequence Alignments



## 2. Relação entre espécies

As espécies são classificadas com base nas suas caracteristicas compartilhadas. Quanto maior o número de características compartilhadas, maior a probabilidade de que as espécies estejam fortemente relacionadas entre si.
As caracteristicas anatômicas compartilhadas, como a presença de uma coluna vertebral, é um exemplo do que pode-se ser utilizado para estudar relações evolutivas. Todavia nem todas as características podem ser fácilmente interpretadas.

### 2.1. Relação baseada em sequências moleculares

As sequências moleculares,(DNA, RNA, proteinas) também são usadas para estudar as relações evolutivas. Assim como na comparação anatômica, a comparação de sequêcias também busca semelhanças e diferenças para deduzir relações. 
se realizada apropriadamente, esta comparação pode ser mais objetiva e menos ambigua que a comparação anatômica.

### 2.2. Diferenças nas sequências e nas relações evolutivas
Os organismos relacionados evolutivamente possuem um ancestral comum, com uma sequencia ancestral de DNA. A medida que os organismos evoluem e divergem, suas sequencias de DNA vão acumulando diferenças, chamadas de mutações.
Os tipos mais comuns de mutações são conehcidos como SNP (single nucleotide polymorphism) e Indels

### 2.3. Árvore Filogenética
A construção da árvore filogenética das 17 espécies estudadas encontra-se no notebook 'PhylogeneticTree.ipynb'. O alinhamento e a árvore foram feitas utilizando a biblioteca Biopython e uma das visualizações exploradas com o Matplotlib encontra-se a seguir:

![image](tree.png)

## 3. Usando Machine Learning para decifrar o genoma
Machine Learning (em português, Aprendizado de Máquina) explora o estudo e construção de algoritmos para coletar dados, aprender com eles, e então fazer uma determinação ou predição sobre os dados. Dentre os vários tipos de aprendizado, iremos explorar o aprendizado por representação PCA (principal components analysis).

### 3.1. Breve resumo do método PCA 
O método de redução da dimensionalidade do PCA é um método de redução de dimensionalidade linear. Ele funciona projetando um número de variáveis correlacionadas em um número (menor) de variáveis não correlacionadas, chamadas componentes principais.
O primeiro componente principal é responsável pela maior parte da variabilidade dados possíveis, e cada componente sucessor é responsável por tanto a variabilidade restante quanto possível. O algoritmo resolve os autovalores e autovetores de uma matriz simétrica quadrada com somas de quadrados e cruz produtos. O autovetor associado ao maior autovalor tem a mesma direção do primeiro componente principal. O autovetor associado ao segundo maior autovalor determina a direção do segundo componente principal. A soma dos autovalores é igual ao traço da matriz quadrada e o número máximo de autovetores é igual ao número de linhas (ou colunas) desta matriz.

### 3.2. 
codons/matriz/pca

## Referências
* [1] Sawyer SL, Wu LI, Emerman M, Malik HS. Positive selection of primate
TRIM5alpha identifies a critical species-specific retroviral restriction domain. 
Proc Natl Acad Sci U S A. 2005 Feb 22;102(8):2832-7. Epub 2005 Feb 2. PubMed
PMID: 15689398; PubMed Central PMCID: PMC549489.

## Resources

* [Ensembl](https://www.ensembl.org/)
* _The Phylogenetic Handbook: a Practical Approach to Phylogenetic Analysis and Hypothesis Testing_,
Philippe Lemey, Marco Salemi, and Anne-Mieke Vandamme (eds.). 2009.
* [Positive selection of primate _TRIM5α_ identifies a critical species-specific retroviral restriction domain.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC549489/)

