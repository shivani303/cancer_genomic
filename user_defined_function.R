data= read.csv("C:/Users/shiva/Desktop/shivani college/rgitbt/cancer genomics/leukemia.csv")
data
dim(data)
countdata <- data[,-1]
countdata[,20] = NULL
head(countdata)

# Store Genename as rownames
rownames(countdata) <- data[,1]
head(countdata)

#user-defined function to return log_cpm value
log_cpm= function(data){
  mat=data
  rownames(mat)=rownames(data)
  for(i in 1:ncol(data)){
    mat[,i]= (data[,i]/sum(data[,i]))*1000000
  }
  l2cpm= log2(mat+1)
  return(l2cpm)
  
}

log_cpm(countdata)
