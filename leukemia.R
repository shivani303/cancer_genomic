
data= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/leukemia.csv")
data
dim(data)
countdata <- data[,-1]
countdata[,20] = NULL
head(countdata)

# Store Genename as rownames
rownames(countdata) <- data[,1]
head(countdata)

#Convert counts to DGEList object
library(edgeR)
library(limma)
y <- DGEList(countdata)
print(head(y))

# See what slots are stored in y
names(y)
cpm_data=cpm(y)
cpm_data
log_data<- cpm(countdata, log=TRUE)
log_data
saveRDS(log_data,file="log.RDS")
View(log_data)
mat= as.matrix(log_data)
mat
boxplot(log_data) # data is normalised

#Z= (value - mean)/ (Standard Deviation)
for(i in 1:nrow(mat)){
  vec= as.numeric(mat[i,])
  mat[i, 1:ncol(mat)] = (vec-mean(vec))/sd(vec)
}
View(mat)
hm= Heatmap(mat,col= colorRamp2(c(-2,0,2),c("green","white","yellow")))
#z_value=apply(mat,1, function(x) (x - mean(x)) / sd(x))
#z_value
#summary(z_value)

#variance
var_data=apply(log_data,1,var)
sort_var=sort(var_data,decreasing = TRUE)
head(sort_var)
top_gene= sort_var[1:50]
top_gene
mat1=mat[names(top_gene),]
mat1
names(top_gene)
library(ComplexHeatmap)
library(circlize)
h1= Heatmap(mat1,col= colorRamp2(c(-2,0,2),c("green","white","yellow")))

#HEATMAP




