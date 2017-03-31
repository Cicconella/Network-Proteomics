# Network-Proteomics

This directory contains two R scripts:

### network.R 

It will create a network of genes based on StringDB (http://string-db.org/) information. The input is two files containing  two columns (Fold Change and Gene). As the example:
```
"Fold Change" "Protein"
-0.18 "CFH"
-0.34 "LAP3"
0.00 "CD99"
0.18 "LASP1"
-0.05 "AK2"
```

### createFile.R

It will create the files to be used on network.R. As input it receives the files with all the parameters of the dataset and it will select only the genes with the desired adjust p-value.

```
"X" "logFC" "logCPM" "LR" "PValue" "adj.p" "Associated.Gene.Name" "Description"
"ENSG00000105514" -0.134335254064683 4.42349961902604 0.205467365607031 0.650343473489778 0.906432279929866 "RAB3D" "RAB3D, member RAS oncogene family [Source:HGNC Symbol;Acc:9779]"
"ENSG00000091136" 0.222660203856437 7.84300864342524 0.767899182485522 0.380867742652381 0.791050284305124 "LAMB1" "laminin, beta 1 [Source:HGNC Symbol;Acc:6486]"
"ENSG00000255197" -0.709764358179691 2.43547840441288 3.2398881068154 0.0718655461747584 0.410048673924121 "RP11-750H9.5" ""
"ENSG00000241158" 0.162111727693081 2.86341830996078 0.235728157108815 0.627308634359797 0.898896453298309 "ADAMTS9-AS1" "ADAMTS9 antisense RNA 1 [Source:HGNC Symbol;Acc:40625]"
```


 
