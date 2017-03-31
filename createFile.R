
##### Load expression information ######

### Insert the path to the folder (it will be use for loading and saving all the files)
folder= getwd() 

### Preparing files for analises
expP = read.table(paste(folder, "/genParameterProt", sep = ""), header = T)
head(expP)
dim(expP)
summary(expP)


### Load diferential expression genes from transcripts data
expT = read.table(paste(folder ,"/genParameterTrans", sep = ""), header = T)
head(expT)
dim(expT)
summary(expT)
###Check the header of both before doing that
colnames(expT) = colnames(expP)


# Get the genes from the both groups
allGenes = c(as.character(expP[expP$adj.p<0.05,7]),as.character(expT[expT$adj.p<0.05,7]))
allGenes = unique(allGenes)
allGenes = as.data.frame(allGenes)
dim(allGenes)

prot = expP[which(expP$gene %in% allGenes[,1]),]
trans = expT[which(expT$gene %in% allGenes[,1]),] 

a = prot[,c(2,7)]
b = trans[,c(2,7)]

colnames(a) = c("Fold Change", "Gene")
colnames(b) = c("Fold Change", "Gene")

write.table(a, "Proteomics-output")
write.table(b, "Transcripts-output")

