##### Clear all variables ########
rm(list = ls())
gc()

##### Libraries #####
source("http://bioconductor.org/biocLite.R")
require(igraph)
require(plotrix)
library(STRINGdb)

##### Load data ######

### Insert the path to the folder (it will be use for loading and saving all the files)
folder= getwd()

### Load diferential expression genes from proteomics analysis
a = read.table("Proteomics", header = T)
head(a)

### Load diferential expression genes from transcripts data
b = read.table("Transcripts", header = T)
head(b)

colnames(a) = c("Fold Change", "Protein")
colnames(b) = c("Fold Change", "Protein")

dim(a)
dim(b)

##### String Database #####

# Load String Database
string_db = STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="")

# Get the GO ID  
annotations = string_db$get_annotations()
head(annotations)


#Select the maximum absolute values of fold change, it will be used for normalize the both distributions at the same scale
#maximo = max(abs(min(as.numeric(a$`Fold Change`))), max(as.numeric(a$`Fold Change`)),
#             abs(min(as.numeric(b$`Fold Change`))), max(as.numeric(b$`Fold Change`)))
#maximo = ceiling(maximo)

#You can define maximo as well
maximo = 4

#hist(as.numeric(a$`Fold Change`))
#hist(as.numeric(b$`Fold Change`))

# For every value of the fold change, we add the ceiling of the maximum (called maximo) and then dive by 2 times of the maximo, 
# so the values will be distributed between 0 and 1
a$`Fold Change` = (a$`Fold Change`+maximo)/(2*maximo) 
b$`Fold Change` = (b$`Fold Change`+maximo)/(2*maximo) 

a$`Fold Change`[a$`Fold Change`>1] = 1
b$`Fold Change`[b$`Fold Change`>1] = 1

a$`Fold Change`[a$`Fold Change`<0] = 0.0001
b$`Fold Change`[b$`Fold Change`<0] = 0.0001

#hist(a$`Fold Change`)
#hist(b$`Fold Change`)


##### Select the DEs #####

allGenes = c(as.character(a$Protein),as.character(b$Protein))
allGenes = unique(allGenes)
allGenes = as.data.frame(allGenes)

# Get the String Database IDs 
allGenes = string_db$map(allGenes, "allGenes", removeUnmappedRows = TRUE)
links = string_db$get_interactions(allGenes[,2])

### Choose here the types of interactions, in this case, only interactions from experiments are choosen.
interactions = unique(c(which(links$experiments!=0)))
links = links[interactions, c(1,2,10)]
dim(links)
head(links)
summary(links$experiments)
head(links)

### Select the minimum confidence level, in this case, the minimum is set to 0.400.
links = links[links$experiments>399,]
dim(links)

links = links[,c(1,2)]

# Change the String ID for gene names
for(i in 1:nrow(links)){
  links[i,1] = allGenes[which(allGenes$STRING_id==links[i,1])[1],1]
  links[i,2] = allGenes[which(allGenes$STRING_id==links[i,2])[1],1]
}

head(links)

# Create graph 
graphObj = graph_from_data_frame(links, directed=F)

# The first group will be "vermelho" (red in portuguese)
head(a)
a$Protein = as.character(a$Protein)

a = a[a$Protein %in% V(graphObj)$name,]

vermelho = a
head(vermelho)

# The second group will be "azul" (blue in portuguese)
head(b)
b$Protein = as.character(b$Protein)

b = b[b$Protein %in% V(graphObj)$name,]
dim(b)

azul = b[,c(1, 2)]
head(azul)

# Get the names of the vertices in the graph object (genes that have at least one conection) 
nomes = V(graphObj)$name
nomes = as.data.frame(nomes)

#Get the fold change information of the genes/proteins 
nomes = merge(nomes, vermelho,  by.x = "nomes", by.y = "Protein", all.x = T)
nomes = merge(nomes, azul,  by.x = "nomes", by.y = "Protein", all.x = T)

head(nomes)

nomes[is.na(nomes[,2]),2] = 0
nomes[is.na(nomes[,3]),3] = 0

head(nomes)

# For each gene, we classify it by "trascripts", "proteins"  and "both". The variable is 

classGene = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,]
  
  if(d[2]==0){
    classGene = c(classGene, "transcripts")
  }else if(d[3]==0){
    classGene = c(classGene, "proteins")
  }else{
    classGene = c(classGene, "both")
  }
}

head(nomes)

nomes = cbind(nomes, classGene)

nomes[1:20,]

# Order "nomes" based on the graph object
order = c()

for(i in 1:length(V(graphObj)$name)){
  order = c(order, which(as.character(nomes[,1])==V(graphObj)$name[i]))  
}

order

head(nomes)
nomes = nomes[order,]
head(nomes)
nomes$classGene = as.character(nomes$classGene)

rm(azul, vermelho, classGene, i, order)

##### GO information #####

# Get the GO information
go =  string_db$get_enrichment(allGenes[,2])
head(go)

head(annotations)
head(go)

# Merge all informations
allGenes = merge(allGenes, annotations, by.x = "STRING_id", by.y = "STRING_id")
head(allGenes)

allGenes = merge(allGenes, go, by.x = "term_id", by.y = "term_id")
head(allGenes)

dim(allGenes)

#If you desire to save the all the GO information about the genes. 
#allGenes$pvalue = format(allGenes$pvalue, scientific = FALSE)
#allGenes$pvalue_fdr = format(allGenes$pvalue_fdr, scientific = FALSE)
#write.table(allGenes, "info.csv", row.names = F)

#If you want to have a list of all functions with at least one gene highly associated
#lista = (unique(allGenes[allGenes$pvalue_fdr<0.01,10]))
#lista

# Change the colors classification by format
formats = rep("pie", length(as.character(nomes$classGene)))

# Define colors of the vertices, if the fold change is positive, it should be in a scale red, if it's  negative, 
# it should be in a scale of blue. Vertices with two fold change values will be set as white for now.

head(nomes)
classGene = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,3]
  
  if(d > 0.5){
    classGene = c(classGene, rgb(1, (1-d)*2, (1-d)*2))
  } else if(d < 0.5 && d != 0){
    classGene = c(classGene, rgb((d*2), (d*2), 1))
  } else if(d == 0.5){
    classGene = c(classGene, rgb(1,1,1))
  } else if(d == 0){
    classGene = c(classGene, rgb(0,0,0))
  }
}


colorFinal = classGene
classGene = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,2]
  
  if(d > 0.5){
    classGene = c(classGene, rgb(1, (1-d)*2, (1-d)*2))
  } else if(d < 0.5 && d != 0){
    classGene = c(classGene, rgb((d*2), (d*2), 1))
  } else if(d == 0.5){
    classGene = c(classGene, rgb(1,1,1))
  } else if(d == 0){
    classGene = c(classGene, rgb(0,0,0))
  }
}

colorFinal = cbind(colorFinal, classGene)

classGene = list()

for(i in 1:nrow(colorFinal)){
  classGene[[i]] = colorFinal[i,]
}

classGene

head(colorFinal)
head(nomes)

for(i in 1:nrow(nomes)){
  values <- lapply(1:i, function(x) c(5,5))  
}

# Choose biological function here
bio = "actin filament organization"


# Identify which genes has the selected function
head(nomes)

aux = allGenes[allGenes$term_description==bio,3]
aux = as.character(aux)
aux = as.character(nomes[nomes$nomes %in% aux,1])

for(i in 1:length(aux)){
  aux[i] = which(V(graphObj)$name==aux[i])
}

aux = as.numeric(aux)
aux2 = rep("azure4", length(V(graphObj)$name))
aux2[aux] = "red"


### Plot ###
#Plot will be saved in your current directory with the name of the biological function
png(paste(bio, ".png", sep=""), width = 1224, height = 1000)

par(mar=c(5.1, 5, 6, 8.2), xpd=TRUE)

par(fig=c(0, 1, 0, 1), xpd=TRUE)
#, layout = newlayout
if (interactive()) {
  plot(graphObj, vertex.shape=formats, vertex.pie=values,
       vertex.pie.color= classGene, vertex.label.dist=0.25, vertex.label.degree=pi/2, vertex.label.cex=0.8, 
       vertex.label.color = aux2, vertex.size=7, main = "Gene Network", edge.color = "darkgrey")
}

testcol<-color.gradient(c(1,1,0),c(0,1,0),c(0,1,1),nslices=100)
color.legend(1.3,1,1.5,0.3, c("Up[4+]", "","Zero", "", "Down[4-]"), testcol, gradient="y", cex = 0.55)
text(1.35, 1.1,"Fold Change", cex=0.7, xpd=T)

par(fig=c(0.7, 1, 0, 0.3), new=TRUE)

slices <- c(50,50)
lbls <- c("Proteins", "Transcripts")
pie(slices, labels = lbls, clockwise = T, col = c("white","lightgray"), radius = 0.5)

dev.off()

