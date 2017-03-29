##### Clear all variables ########
rm(list = ls())
gc()

##### Libraries #####
source("http://bioconductor.org/biocLite.R")
require(igraph)
require(plotrix)
library(STRINGdb)

##### Load expression information ######

### Insert the path to the folder (it will be use for loading and saving all the files)
folder= "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/metaanalysis/new/" 

### Load diferential expression genes from proteomics analysis
expP = read.table(paste(folder, "General_parameters_1D25.txt", sep = ""), header = T, sep = ",")
head(expP)
dim(expP)

### Load diferential expression genes from transcripts data
expT = read.table(paste(folder ,"General_parameters_20160302_LundData_transcriptomics_freeze_exons.txt", sep = ""), header = T, sep = ",")
head(expT)
dim(expT)

head(expP)
head(expT)

a = expP[,c(2,7)]
b = expT[,c(2,7)]

head(a)
head(b)

colnames(a) = c("Fold Change", "Protein")
colnames(b) = c("Fold Change", "Protein")


##### String Database #####

# Load String Database
string_db = STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="")

# Get the GO ID  
annotations = string_db$get_annotations()
head(annotations)


# Select the maximum absolute values of fold change, it will be used for normalize the both distributions at the same scale
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

head(expP)
head(expT)
dim(expT)

#which(expP$adj.p<0.05)


# Get the genes from the both groups
allGenes = c(as.character(expP[expP$adj.p<0.05,7]),as.character(expP[expT$adj.p<0.01,7]))
allGenes = unique(allGenes)
allGenes = as.data.frame(allGenes)

write.table(allGenes$allGenes, "/home/cicconella/Desktop/all01", quote = F, row.names = F, col.names = F)

# Get the String Database IDs 
allGenes = string_db$map(allGenes, "allGenes", removeUnmappedRows = TRUE)
links = string_db$get_interactions(allGenes[,2])

### Choose here the types of interactions
interactions = unique(c(which(links$experiments!=0)))
links = links[interactions, c(1,2,10)]
dim(links)
head(links)
summary(links$experiments)
head(links)

links = links[links$experiments>399,]
dim(links)

links = links[,c(1,2)]

# Change the String ID for gene names
for(i in 1:nrow(links)){
  links[i,1] = allGenes[which(allGenes$STRING_id==links[i,1])[1],1]
  links[i,2] = allGenes[which(allGenes$STRING_id==links[i,2])[1],1]
}

head(links)

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

# For each gene, we classify it by "trascripts", "proteins"  and "both" 

cores = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,]
  
  if(d[2]==0){
    cores = c(cores, "transcripts")
  }else if(d[3]==0){
    cores = c(cores, "proteins")
  }else{
    cores = c(cores, "both")
  }
}

head(nomes)

nomes = cbind(nomes, cores)

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
nomes$cores = as.character(nomes$cores)



rm(azul, vermelho, cores, i, order)

##### GO information #####

# Get the GO information
go =  string_db$get_enrichment(allGenes[,2])
head(go)

head(annotations)
head(go)

# Merge all informations
allGenes = merge(allGenes, annotations, by.x = "STRING_id", by.y = "STRING_id")
#head(allGenes)

allGenes = merge(allGenes, go, by.x = "term_id", by.y = "term_id")
head(allGenes)

dim(allGenes)

allGenes$pvalue <- format(allGenes$pvalue, scientific = FALSE)
allGenes$pvalue_fdr <- format(allGenes$pvalue_fdr, scientific = FALSE)

#write.table(allGenes, "/home/cicconella/Desktop/gogenesProt05.csv", row.names = F)


dim(allGenes[allGenes$pvalue_fdr<0.05,])

lista = (unique(allGenes[allGenes$pvalue_fdr<0.01,10]))

lista



# Change the colors classification by format
formats = rep("pie", length(as.character(nomes$cores)))


# Define colors of the vertices, if the fold change is positive, it should be in a scale red, if it's  negative, 
# it should be in a scale of blue. Vertices with two fold change values will be set as white for now.

head(nomes)
cores = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,3]
  
  if(d > 0.5){
    cores = c(cores, rgb(1, (1-d)*2, (1-d)*2))
  } else if(d < 0.5 && d != 0){
    cores = c(cores, rgb((d*2), (d*2), 1))
  } else if(d == 0.5){
    cores = c(cores, rgb(1,1,1))
  } else if(d == 0){
    cores = c(cores, rgb(0,0,0))
  }
}


colorFinal = cores
cores = c()

for(i in 1:nrow(nomes)){
  d = nomes[i,2]
  
  if(d > 0.5){
    cores = c(cores, rgb(1, (1-d)*2, (1-d)*2))
  } else if(d < 0.5 && d != 0){
    cores = c(cores, rgb((d*2), (d*2), 1))
  } else if(d == 0.5){
    cores = c(cores, rgb(1,1,1))
  } else if(d == 0){
    cores = c(cores, rgb(0,0,0))
  }
}

colorFinal = cbind(colorFinal, cores)

cores = list()

for(i in 1:nrow(colorFinal)){
  cores[[i]] = colorFinal[i,]
}

cores

head(colorFinal)
head(nomes)

for(i in 1:nrow(nomes)){
  values <- lapply(1:i, function(x) c(5,5))  
}

threshold = 0.1

if(threshold==0.1){
  layoutGraph = layout_with_fr(graphObj)
  set = as.character(nomes$nomes) 
  newlayout = layoutGraph
}else{
  dim((layoutGraph))
  length(set)
  
  nomes$nomes = as.character(nomes$nomes)
  
  newlayout = cbind(set, layoutGraph)
  
  newlayout = merge(newlayout, nomes, by.x = "set", by.y = "nomes") 
  
  head(newlayout)
  head(nomes)
  
  order = c()
  
  for(i in 1:length(V(graphObj)$name)){
    order = c(order, which(as.character(newlayout$set)==V(graphObj)$name[i]))  
  }
  
  order
  
  head(newlayout)
  
  newlayout = newlayout[order,]
  
  head(newlayout)
  
  newlayout = newlayout[,c(2,3)]
  
  
  newlayout = cbind(as.numeric(as.character(newlayout[,1])),as.numeric(as.character(newlayout[,2])))
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


png("/home/cicconella/Desktop/actin.png", 
    width = 1224, height = 1000)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

par(fig=c(0, 1, 0, 1), xpd=TRUE)
#, layout = newlayout
if (interactive()) {
  plot(graphObj, vertex.shape=formats, vertex.pie=values,
       vertex.pie.color= cores, vertex.label.dist=0.25, vertex.label.degree=pi/2, vertex.label.cex=0.8, 
       vertex.label.color = aux2, vertex.size=7, main = "Gene Network", edge.color = "darkgrey")
}

testcol<-color.gradient(c(1,1,0),c(0,1,0),c(0,1,1),nslices=100)
color.legend(1.3,1,1.5,0.3, c("Up[4+]", "","Zero", "", "Down[4-]"), testcol, gradient="y", cex = 0.55)
text(1.35, 1.1,"Fold Change", cex=0.7, xpd=T)

par(fig=c(0.6, 1, 0, 0.5), new=TRUE)

slices <- c(50,50)
lbls <- c("Proteins", "Transcripts")
pie(slices, labels = lbls, clockwise = T, col = c("white","lightgray"), radius = 0.7)

dev.off()


bio = "cytoskeleton organization"

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

png("/home/cicconella/Desktop/cyto.png", 
    width = 1224, height = 1000)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

par(fig=c(0, 1, 0, 1), xpd=TRUE)
#, layout = newlayout
if (interactive()) {
  plot(graphObj, vertex.shape=formats, vertex.pie=values,
       vertex.pie.color= cores, vertex.label.dist=0.25, vertex.label.degree=pi/2, vertex.label.cex=0.65, 
       vertex.label.color = aux2, vertex.size=7, main = "Gene Network", edge.color = "darkgrey")
}

testcol<-color.gradient(c(1,1,0),c(0,1,0),c(0,1,1),nslices=100)
color.legend(1.3,1,1.5,0.3, c("Up", "","Zero", "", "Down"), testcol, gradient="y", cex = 0.55)
text(1.35, 1.1,"Fold Change", cex=0.7, xpd=T)

par(fig=c(0.6, 1, 0, 0.5), new=TRUE)

slices <- c(50,50)
lbls <- c("Proteins", "Transcripts")
pie(slices, labels = lbls, clockwise = T, col = c("white","lightgray"), radius = 0.2)

dev.off()


lista = (unique(allGenes[allGenes$pvalue_fdr<0.01,10]))


# Choose biological function here



# Identify which genes has the selected function
head(nomes)

aux = allGenes[allGenes$term_description==bio,3]
aux = as.character(aux)
aux = as.character(nomes[nomes$nomes %in% aux,1])

for(i in 1:length(aux)){
  aux[i] = which(V(graphObj)$name==aux[i])
}

aux = as.numeric(aux)
aux2 = rep("gray", length(V(graphObj)$name))
aux2[aux] = "black"

png(paste(folder, "bio/", bio,".png", sep=""), width = 960, height = 760)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

par(fig=c(0, 1, 0, 1), xpd=TRUE)

if (interactive()) {
  plot(graphObj, vertex.shape=formats, vertex.pie=values,
       vertex.pie.color= cores, vertex.label.dist=0.25, vertex.label.degree=pi/2, vertex.label.cex=0.65, 
       vertex.label.color = aux2, vertex.size=7, main = "Gene Network", edge.color = "darkgrey", 
       layout = newlayout)
}

testcol<-color.gradient(c(1,1,0),c(0,1,0),c(0,1,1),nslices=100)
color.legend(1.3,1,1.5,0.3, c("Up", "","Zero", "", "Down"), testcol, gradient="y", cex = 0.55)
text(1.35, 1.1,"Fold Change", cex=0.7, xpd=T)

par(fig=c(0.6, 1, 0, 0.5), new=TRUE)

slices <- c(50,50)
lbls <- c("Proteins", "Transcripts")
pie(slices, labels = lbls, clockwise = T, col = c("white","lightgray"), radius = 0.5)

dev.off()
}

