
##### Load expression information ######

### Insert the path to the folder (it will be use for loading and saving all the files)
folder= "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/metaanalysis/new/" 

folder= "/media/cicconella/8AA6013CA6012A71/Documents and Settings/Nina/Dropbox/metaanalysis/new/" 

### Preparing files for analises
expP = read.table(paste(folder, "General_parameters_1D25.txt", sep = ""), header = T, sep = ",")
head(expP)
dim(expP)
summary(expP)


a = sample(c(1:dim(expP)[1]))[1:100]


a = expP[a,]

summary(a)

write.table(a, "genParameterProt", row.names = F)

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