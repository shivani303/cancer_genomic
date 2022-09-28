#DIFFERENTIAL GENE EXPRESSION
data1= readRDS("log.RDS")

mat= matrix(NA,ncol=4,nrow=nrow(data1))
rownames(mat)=rownames(data1)
colnames(mat)= c("control","test","pval","logfc2")

data2 = round(data1, 4)

for (i in 1:nrow(data2)){
  cdata = as.numeric(data2[i, 1:3])
  tdata = as.numeric(data2[i, 4:19])
  
  if(length(unique(cdata))==1){
    next
  }
  
  if(length(unique(tdata))==1){
    next
  }
  t= t.test(cdata,tdata,paired = F,alternative = 'two.sided')
  mat[i,1]= t$estimate[[1]]
  mat[i,2]= t$estimate[[2]]
  mat[i,3]= t$p.value
  mat[i,4]= mat[i,1]- mat[i,2]
}
summary(data1[,1])
mat
#volcano plot

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

library(EnhancedVolcano)
library(ggplot2)
library(ggrepel)
mat= as.data.frame(mat)
num= which(is.na(mat$pval))
num
mat[num,'pval']= 1

EnhancedVolcano(mat,lab=rownames(mat),x='logfc2',y='pval')

