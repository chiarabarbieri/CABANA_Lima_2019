---
title: "Analisis y filogenia de poblaciones humanas"
author: "chiara barbieri"
date: "28/11/2019"
output: rmarkdown::github_document
---

## A dataset of Human Diversity

En este ejercicio vamos a trabajar con un dataset de SNP chip que incluye diferentes poblaciones. El objetivo es comprender la variación genética humana a través de la historia poblacional. Elegimos un SNP chip array diseñado para los estudios de historia humana, para maximizar la información sobre diversidad humana y reconstruir eventos demográficos. El nombre de lo SNP chip es Human Origins (Affymetrix). El array incluye SNPs que se encuentran en poblaciones de diferentes continentes, para minimizar el ascertainment bias effect. 


Trabajamos con PLINK para mirar el dataset y correr simples análisis, y ADMIXTURE para reconstruir componentes de ancestralidad entre los individuos.   

El dataset incluye 100 individuos y 14 poblaciones de Africa y Medio Oriente. El dato esta publicado en [Patterson et al. 2012](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2012_Patterson_AncientAdmixture_Genetics.pdf) 
and [Lazaridis et al. 2014](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2014_Nature_Lazaridis_EuropeThreeAncestries.pdf) .



```{r echo=FALSE ,message=FALSE, warning = FALSE}
infoo<-read.table("infopopExercis.txt", as.is = T, header=T,quote = "", sep="\t")
library("maps")
library("ggrepel")
map.world<-map_data(map="world")
gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-30,40), xlim=c(0,50))+
  geom_label(data=infoo,aes(x=lon, y=lat,label =PopName)) +
geom_label_repel(force = 1)
gg
```

Trabajamos local, desde el terminal. Puedes descargar la carpeta con todo los datos desde GitHub

```
git clone https://github.com/chiarabarbieri/Bio373_Blockcourse.git
```


___________________________

# PLINK

Disponible en   <https://www.cog-genomics.org/plink/1.9/>.
Es una herramienta para genome-wide association studies (GWAS) y investigación de genética poblacional.
PLINK trabaja con comandos y “flag”. Cada flag empieza con dos guiónes. El flag puede ser seguido por parámetros. PLINK fue desarrollado para estudios médicos, y los formatos utilizan informaciones como pedigree o fenotipo que no nos interesan para nuestras análisis. 


### Formatos diferentes para los input files

PLINK puede tomar varios input files, también .vcf. El formato nativo de PLINK consiste en tablas de individuos y variant calls. 

El formato se encuentra en dos variantes: binario (bed + bim + fam) y texto (ped + map).

El .ped incluye ID, pedigree (optional) + tabla de genotipo. El .map es básicamente la lista de SNPs con posición cromosómica y alleles. Los dos o tres archivos deben tener el mismo nombre (solo varia el sufijo). 

*.bed* for binary and *.ped* for text files:
Contiene la información del genotipo. Una línea para cada individuo.

El *.ped* file contiene también información para cada individuo. En la forma binaria, esta información esta contenida en el file .fam. 
Las primeras seis columnas del *.ped*  y del  *.fam* son las mismas:
     
     Family ID
     Individual ID
     Paternal ID
     Maternal ID
     Sex (1=male; 2=female; other=unknown)
     Phenotype (for association studies)



*.map* para archivos de textos:
Lista de marcadores geneticos.
Cada línea del .map describe un marcador y debe contener 4 columnas:

     chromosome (1-22, X, Y or 0 if unplaced)
     rs# or snp identifier
     Genetic distance (Centimorgans)
     Base-pair position (bp units)
     

*.bim* para archivos binarios:
Cada línea del  .bim describe un marcador y contiene seis columnas. Es una forma extendida del .map file con dos columnas mas con los alelos. 
     


### desde un formato hasta otro

La línea de comando de PLINK es

```
plink --file TuArchivo --flag modifiers que hacen algo con tu fichero
```

donde * TuArchivo * es la raíz del nombre compartida entre los dos ficheros .ped y .map. Si tu utilizas el flag --bfile, en cambio, vas a llamar los tres ficheros binarios .bed, .bin y .fam.

Ejemplo: pasar entre un formato y otro y mirar a las diferencias en el terminal.

```
plink --bfile HumanDataHO --recode
```

La documentación online describe otros flags para manipular ficheros. 

Otras herramientas útiles son: hacer un subset de SNPs, un subset de individuos, merge entre datasets. 

## Herramientas básicas de genética de poblacion

Generar simple summary statistics. Proporción de datos faltante. Diversidad dentro población y entre poblaciones. 

### datos faltantes

Con el flag --missing. Exploramos los outputs. Como esta la proporción de datos faltante para marcador genético y para individuos? 


```
plink --bfile HumanDataHO --missing
```
Trazamos un grafico en R. 


```{r echo=FALSE}
library("ggplot2")
aa<-read.table("missing.imiss", header=T)
pdf("missing.pdf")

ggplot(aa, aes(FID,F_MISS))+
         geom_boxplot()+
     theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
   labs(title = "percentage of missing data per individual")
   dev.off()

```


### Heterozygosity - inbreeding - consanguinity 





Con otro comando en PLINK miramos a F, el grado de consanguineidad. Si es el caso, podemos eliminar individuos con un F muy alto. 


*--het* computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates (i.e. ([observed hom. count] - [expected count]) / ([total observations] - [expected count])) to plink.het. 


```
plink --bfile HumanDataHO --het
```

Visualizamos las diferencias en heterocigosidad entre poblaciones.


```{r message=FALSE, warning = FALSE}
het<-read.table("plink.het", header=T)
pdf("het.pdf")
ggplot(het, aes(FID,F))+
         geom_boxplot()+
     theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
   labs(title = "homozygosity per individual")
   dev.off()
```


```{r  message=FALSE, warning = FALSE}

infoo<-read.table("infopopExercis.txt", as.is = T, header=T,quote = "", sep="\t")
infoo$het<-NA
for (i in 1:nrow(infoo)){
  temp<-het[which(het$FID==infoo$PopName[i]),]
  infoo$het[i]<-mean(temp$F)
}
library("maps")
map.world<-map_data(map="world")
pdf("het_map.pdf")
gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)
gg<- gg+coord_quickmap(ylim=c(-30,40), xlim=c(0,50))
gg + geom_point(data=infoo, aes(x=lon, y=lat, color=het), size=5 )+
  scale_color_gradient(low = "blue", high = "red") +
  ggtitle("Intensity of homozygosity in each population")
   dev.off()

```

La diversidad sigue el modelo de migraciones humanas Out of Africa?


___________________________

## PCA

With the modifier *--pca*. PLINK extracts the top 20 principal components of the variance-standardized relationship matrix. 
The results consist in a *.eigenvec* file with the coordinates for each individual in rows and eigenvectors in columns, and a *.eigenval* which explains how much variance there is in the data for that given vector. 

```{r message=FALSE, warning = FALSE}
eigenvec<-read.table("plink.eigenvec")
eigenval<-read.table("plink.eigenval")

library("ggplot2")
library("RColorBrewer")
library("colorRamps")
pdf("pca.pdf")
gg<-ggplot(eigenvec,aes(V3,V4, color=V1))+
  geom_point()+
 labs(x = eigenval[1,], y=eigenval[2,],title = "PCA analysis dimension 1 vs. 2")+
scale_color_manual(values = colorRampPalette(brewer.pal(12,  "Accent"))(12))
gg
   dev.off()

```

One group is clearly an outlier. Repeat the analysis after excluding this population.

___________________________

# ADMIXTURE analysis


Identifying ancestry components shared between individuals of a set of populations.

ADMIXTURE is a software that works similarly to Structure, but with faster computation time. Download and manual from https://www.genetics.ucla.edu/software/admixture/. It takes plink format input files.



### Pruning

First, we prune the dataset for excluding the SNPs in Linkage, with Plink. The resulting file will have less SNPs, and the computation will be faster. 

The settings define window size, step and the r2 threshold.

```
plink --bfile HumanDataHO --indep-pairwise 200 25 0.4 --out x.tmp
plink --bfile HumanDataHO --extract x.tmp.prune.in --recode12 --out HumanDataHO_pruned
```

How many SNPs are left after pruning?


### ADMIXTURE run

Now the proper ADMIXTURE run. The following commands will run ADMIXTURE for each *K* (number of ancestry blocks) desired. One value of *K* will be more supported by the analysis: the one with the lowest associated cross-validation error. This *K* will be considered as the best representation of the actual data variation.

```
typeset -i run=0
for K in 2 3 4 5; do  # select a meaningful series of K - the more Ks, the longer the run obviously
admixture -s time --cv HumanDataHO_pruned.ped $K | tee log.K${K}.RUN1.out;
mv HumanDataHO_pruned.$K.P K$K.Run1.P;
mv HumanDataHO_pruned.$K.Q K$K.Run1.Q;
done
```
For each run there are three output: .out, .P, and .Q

Now we will elaborate the outputs with a mix of bash commands and *R*.

### Cross-validation error to identify the best value of *K*

```
grep -h CV log*out > CV.txt
```
plot the distribution of values associated to each K in R.

```{r message=FALSE, warning = FALSE}

read.table("CV.txt")->cv  

# set which one was the smallest and the largest K that you ran
minK<-2
maxK<-10

ordine<-c()
for (k in minK:maxK){
  ordine[k]<-paste("(K=",k,"):",sep="", collapse = "")
}

ordine <- ordine[!is.na(ordine)]

library(ggplot2)

p <- ggplot(cv, aes(x=V3, y=V4)) + 
  geom_boxplot()+ 
  ggtitle("values associated to each K")+
  scale_x_discrete(limits=ordine)
  
p
```

I previously run 5 iterations for each K and determined which of the 5 runs has the highest likelihood. This is to exclude the chance that some run did not perform correctly.

_______________


### Plotting the ADMIXTURE results for each K

prepare to plot: information to plot on the admixture bars, from your external info population file


```

MYINFO<-read.table("infoAdmixtureExercis.txt", header=T, as.is=T)


table(MYINFO$population)->pops
namespop<-unique(MYINFO$population)

my.labels <- vector()   ## plotting pop labels instead of sample ids
for (k in 1:length(namespop)){
  paste("^",namespop[k],"$",sep="")->a
  length(grep(a, MYINFO$population)) -> my.labels[k]
}

labels.coords<-vector()  ### where to plot pop labels
labels.coords[1]<-my.labels[1]/2
for (i in 1:(length(my.labels)-1)) {
  labels.coords[i]+(my.labels[i]/2+my.labels[i+1]/2)->labels.coords[i+1]
}
z<-vector()
z[1]<-my.labels[1]
for (i in 1:(length(my.labels)-1)){
  z[i]+my.labels[i+1]->z[i+1]
}

# select a color palette
# you can use colorbrewer. put together a number of colours equal Kmax.

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
colorchoice<-sample(col_vector, maxK)
#pie(rep(1,maxK), col=coso)

# now plot for each K
K<-2 # chose the K to plot. Start with the lowest.

valuesToplot<-read.table(paste("K",K,".Run1.Q", sep="", collapse = ""))

valuesToplotSort<-valuesToplot[MYINFO$oderAdmix,]

#pdf(paste("AdmixtureForK",K,".pdf", sep="", collapse = ""),pointsize=8, height=3.5)

barplot(t(as.matrix(valuesToplotSort)), col=colorchoice[1:K], axisnames=F, axes=F, space=0,border=NA)
axis(3, labels = FALSE, tick=F)
for (i in c(0,z)){
  lines(x=c(i,i),y=c(0,3), lwd=0.7, col="white")
}
text(labels.coords, par("usr")[3] + 1.03 , srt = 45, adj = 0, labels = namespop, cex=0.7, xpd = TRUE)
#dev.off()

```

![](admixtureFigure.png)





Look at patterns across populations. Do they follow a geographic structure? Is there a sign of Admixture?

... 
________________________________


# SECOND PART: PHYLOGENIES

### Handling sequence data and a matrix of distances

Handling genetic alignments, genetic distance, trees, genetic distance between populations.
Sequence data from Barbieri et al. 2013 "Unraveling the Complex Maternal History of Southern African Khoisan Populations"" (mitochondrial sequences mtDNA)
Original sequences are deposited at <http://www.ncbi.nlm.nih.gov/popset?DbFrom=nuccore&Cmd=Link&LinkName=nuccore_popset&IdsFromResult=444211048>

This part of the tutorial will be performed in R.

## 1. Sequence alignment

Load packages in R with the analysis that you will need.

```{r, message=FALSE, warning=FALSE}
library(ape, adegenet)
library(phangorn)
library (lattice)
library(MASS)

```

Import the sequence data
```{r}
seq<-read.dna("reducedSet100sequences.fasta", format="fasta")
```
check the alignment!
``` {r seq}
print(seq)
```


# 2. Tree based on individual distances  - NJ

Now we calculate a simple genetic distance between sequences.

```{r}
dist.matrix <- dist.dna(seq) # Calculate distance matrix
m5<- as.dist(dist.matrix, diag=F, upper=F)
```
We create a Neighbour Joining tree from the matrix of pairwise distance
```{r}
treeRED<-nj(m5)   #create the object neighbour joining tree 
treeRED$edge.length[treeRED$edge.length < 0] = 0.002 #little trick to avoid negative branches
```
Plot the tree and add a legend
```
plot.phylo(treeRED, type="u", tip.col="blue", cex=0.3 ) #plot the tree as unrooted

```
here the tree without colors

```{r, echo=FALSE}
plot(treeRED, type="u", cex=0.3 ) #plot the tree as unrooted
```

## 3. Maximum Parsimony analysis on sequences
```
fasta.phyDat  <- as.phyDat(seq) # convert to phangorn data type
njtree <- nj(dist.matrix) # Calculate NJ tree
parsimony(njtree, fasta.phyDat) # Determine the tree length of NJ tree
tree.optim <- optim.parsimony(njtree, fasta.phyDat) # more settings: weight matrix, algorithm
```

Parsimony Ratchet algorithm
```
tree.pratchet <- pratchet(fasta.phyDat, maxit=100, k=10, rearrangements="SPR")
```
 Root tree: I have an outlier, a sequence from a Neandertal bone remain.
```
tree.pratchet <- root(tree.pratchet, match("Neandertal", tree.pratchet$tip.label), resolve.root=T)
```
Annotate with the number of mutations per branch
```
tree.pratchet <- acctran(tree.pratchet, fasta.phyDat)
```
 Plot tree
```
plot(tree.pratchet, tip.col="blue", cex=0.3)
edgelabels(tree.pratchet$edge.length)
```
 Write the tree to a file
```
write.tree(tree.pratchet, file="output.tree")
```
now you can open it with FigTree or SplitsTree if you want!

download the software FigTree <http://tree.bio.ed.ac.uk/software/figtree/>


Plus: check the number of mutations. How old is the root of the tree? Soares et al. (2009) calculated the mutation rate of the whole mtDNA genome molecule to be one mutation every 3624 years. <http://www.sciencedirect.com/science/article/pii/S0002929709001633/>


# 4. Trees based on population distance
### Neighbour Joining tree populations
Visualize population distances with a tree (unrooted). Find outliers.
I will use an external file to add a color coding to the populations, corresponding to their Language Family.

```{r}
popOrder <- read.table("popOrderColor", as.is=T)
mat = read.table("matrix_mtDNA_popdist.txt",sep="\t", header=T) #matrix of distances between populations
colnames(mat)<-rownames(mat) 
 	replace(mat,mat <0, 0) -> m3## does not accept negative values

m4<-m3+0.0001      ## does not accept values==0
diag(m4)=0
m5<- as.dist(mat, diag=F, upper=F)
treeTest<-nj(m5)
treeTest$edge.length[treeTest$edge.length < 0] = 0.002
rownames(popOrder)<-popOrder[,1]
popOrder2<-popOrder[treeTest$tip.label,]
plot.phylo(treeTest, type="u", tip.col=popOrder2[,2] )
nomi<-c("Bantu","Khoe","Kx'a","Tuu")
legend("topleft",nomi,text.col=c("darkmagenta","blue","green","red"))
```

### UPGMA tree populations

```{r}
upgmatree<-upgma(m5)
plot.phylo(upgmatree,  tip.col=popOrder2[,2] )
legend("topleft",nomi,text.col=c("darkmagenta","blue","green","red"))
```

What are the differences between the two trees?
