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
