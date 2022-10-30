data_file=read.table("GSE156922_read_counts.txt",row.names  = 1,header = TRUE )
head(data_file)

cpmatrix=data_file
for(i in 1:ncol(data_file)){
  cpmatrix[,i]=(data_file[,i]/sum(data_file[,i]))*1000000
}

#Calculate a log of cpm
logcpm=log2(cpmatrix+1)
saveRDS(logcpm,file="logCPM.rds")
summary(logcpm)

#Calculate a z score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
z_score

#Calculate variance using log 
variance = apply(logcpm, 1, var)
variance = sort(variance,decreasing = T)
top50 = variance[1:50]
pmat = z_score[names(top50),]

#Create a heatmap
library(ComplexHeatmap)
Heatmap(pmat)

#To identify genes which are differential in tumor vs control samples

data_file1=matrix(NA,ncol=4,nrow = nrow(logcpm))
rownames(data_file1)=rownames(logcpm)
colnames(data_file1)=c('meanTumor','meanControl','pvalue','log2FC')
data_file1
for(i in 1:nrow(logcpm)){
  vector1 = as.numeric(logcpm[i, 1:3])
  
  vector2 = as.numeric(logcpm[i, 4:6])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  data_file1[i,1]=res$estimate[[1]]
  data_file1[i,2]=res$estimate[[2]]
  data_file1[i,3]=res$p.value
  data_file1[i,4]=data1[i,1]-data_file1[i,2]
  
}

data_file1=as.data.frame(data_file1)
num=which(is.nan(data_file1$pvalue))
data_file1[num,'pvalue']=1

library(EnhancedVolcano) #importing library
pdf('volcano.pdf',width=10,height=10)
EnhancedVolcano(data_file1,lab = rownames(data_file1),x = 'log2FC' ,y ='pvalue',pCutoff = 3.5)
dev.off()
saveRDS(data_file1,file = "C:/Users/Admin/OneDrive/Desktop/CB/DEG_data.rds")
S=readRDS('DEG_data.rds')
View(S)

