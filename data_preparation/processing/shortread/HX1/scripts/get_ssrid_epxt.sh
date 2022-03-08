esearch -db sra -query "SRX3207163 [Project]" | efetch -format xml | 
xtract -pattern EXPERIMENT_PACKAGE -element PRIMARY_ID