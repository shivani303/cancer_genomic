# LEUKEMIA CANCER DATA
data= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/leukemia.csv")
data
dim(data) **checking dimensions**
countdata <- data[,-1]
countdata[,20] = NULL
head(countdata)

**Store Genename as rownames**

rownames(countdata) <- data[,1]
head(countdata)

# Convert counts to DGEList object
**calling a library**

library(edgeR) 
library(limma)
y <- DGEList(countdata)
print(head(y))

**See what slots are stored in y**

names(y)

# Create count matrix 

cpm_data=cpm(y)
cpm_data

log_data<- cpm(countdata, log=TRUE)
log_data


**saving log file for further use**

saveRDS(log_data,file="log.RDS")

View(log_data)

**coverting to matrix**

mat= as.matrix(log_data)
mat
 # check wheather data is normalised or not
boxplot(log_data)
![image](https://user-images.githubusercontent.com/66779651/193062901-4bb7f64a-ef27-40ca-ae4b-daf94d14f219.png)


# Z-value caculations

for(i in 1:nrow(mat)){
  vec= as.numeric(mat[i,])
  mat[i, 1:ncol(mat)] = (vec-mean(vec))/sd(vec)
}
View(mat)

# Heatmap
library(ComplexHeatmap)
library(circlize)

**ploting only 200 genes**

hm= Heatmap(mat[1:200,],col= colorRamp2(c(-2,0,2),c("green","white","yellow")))
![image](https://user-images.githubusercontent.com/66779651/193064306-fa782403-4a21-4874-ae85-27d0d75c6dfb.png)

# variance

**apply variance formula on log data**

var_data=apply(log_data,1,var)

sort_var=sort(var_data,decreasing = TRUE)

head(sort_var)

**extracted top 50 genes from sorted matrix**

top_gene= sort_var[1:50]
top_gene

mat1=mat[names(top_gene),]
mat1

names(top_gene)

library(ComplexHeatmap)
library(circlize)

h1= Heatmap(mat1,col= colorRamp2(c(-2,0,2),c("green","white","yellow")))
![image](https://user-images.githubusercontent.com/66779651/193064700-09b56a4f-a418-44b7-a78d-4b31cbc8c1ca.png)

# DIFFERENTIAL GENE EXPRESSION 
# T_TEST

**read previous log file**

data1= readRDS("log.RDS")

**creating empty matrix with 4 columns**

mat= matrix(NA,ncol=4,nrow=nrow(data1))
rownames(mat)=rownames(data1)
colnames(mat)= c("control","test","pval","logfc2")

**round off data to 4th decimal postion**

data2 = round(data1, 4)

**divide data itno control and test** 

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

# volcano plot

**installing package**

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

**calling library**

library(EnhancedVolcano)

library(ggplot2)

library(ggrepel)

mat= as.data.frame(mat)

**changing NA values in p_values to 1**

num= which(is.na(mat$pval))
num

mat[num,'pval']= 1

**ploting**

EnhancedVolcano(mat,lab=rownames(mat),x='logfc2',y='pval')

![image](https://user-images.githubusercontent.com/66779651/193067461-24dafa70-8872-4a8b-b2bd-bf1dbe5f2568.png)
