cancer=read.csv("C:/Users/Admin/OneDrive/Documents/cancer genomics/GSE164064__seq.csv",sep=",",header = T,row.names = 1)
cancer #printing the file

#Create a count per matrix
cpmatrix=data_file
for(i in 1:ncol(data)){
  cpmatrix[,i]=(data[,i]/sum(data[,i]))*1000000
}

#Calculate a log of cpm
logcpm=log2(cpmatrix+1)
saveRDS(logcpm,file="logC.rds")
summary(logcpm)

#Calculate a z score
library(matrixStats)
z_score = (logcpm - rowMeans(logcpm))/rowSds(as.matrix(logcpm))[row(logcpm)]
z_score

#Calculate variance using log 
var_c1 = apply(logcpm, 1, var)
var_c2 = sort(var_c1,decreasing = T)
top50 = var_c2[1:50]
pmat = z_score[names(top50),]
data_mat=data.matrix(pmat)
data_mat


#install.packages("ComplexHeatmap")
#Create a heatmap
library(ComplexHeatmap)    #importing library
pdf('data.pdf',width=10,height=10)
Heatmap(data_mat)
saveRDS(logCPM)
dev.off()


#To identify genes which are differential in tumor vs control samples

data1=matrix(NA,ncol=4,nrow = nrow(logcpm))
rownames(data1)=rownames(logcpm)
colnames(data1)=c('meanTumor','meanControl','pvalue','log2FC')
data1
for(i in 1:nrow(logcpm)){
  vector1 = as.numeric(logcpm[i, 1:3])
  
  vector2 = as.numeric(logcpm[i, 4:6])
  
  res=t.test(vector1, vector2, paired = F, alternative = "two.sided")
  data1[i,1]=res$estimate[[1]]
  data1[i,2]=res$estimate[[2]]
  data1[i,3]=res$p.value
  data1[i,4]=data1[i,1]-data1[i,2]
  
}

data1=as.data.frame(data1)
num=which(is.nan(data1$pvalue))
data1[num,'pvalue']=1

library(EnhancedVolcano) #importing library
#pdf('volcano.pdf',width=10,height=10)
EnhancedVolcano(data1,lab = rownames(data1),x = 'log2FC' ,y ='pvalue',pCutoff = 3.5)

