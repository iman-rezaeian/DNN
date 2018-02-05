NetUnion2<-function(matrix1, matrix2){
  if((row.names(matrix1) != colnames(matrix1)) || (nrow(matrix1) != ncol(matrix1))){
    print("Matrix1 must be square and have same row/column names")
    return(0)
  }
  if((row.names(matrix2) != colnames(matrix2)) || (nrow(matrix2) != ncol(matrix2))){
    print("Matrix2 must be square and have same row/column names")
    return(0)
  }
  #ordering the graph
  if(!(nrow(matrix1)==1 & ncol(matrix1)==1))
    matrix1<-matrix1[sort(row.names(matrix1)), sort(colnames(matrix1))]
  if(!(nrow(matrix2)==1 & ncol(matrix2)==1))
    matrix2<-matrix2[sort(row.names(matrix2)), sort(colnames(matrix2))]
  #this part transforms mat1 and mat2 into a 1/0 graph
  #need to put in error for in matrix1 or 2 is not square or is null
  shared.names<-intersect(colnames(matrix1), colnames(matrix2))
  #this is the total names, and the length of the final product
  total.names<-sort(union(colnames(matrix1), colnames(matrix2)))
  #need to get the missing from each
  missing2.names<-setdiff(colnames(matrix1), colnames(matrix2))
  missing1.names<-setdiff(colnames(matrix2), colnames(matrix1))
  #get the number of samples
  l1<-length(colnames(matrix1))
  l2<-length(colnames(matrix2))
  #this goes through and addes rows/col of 0 for missing in matrix2
  if(length(missing2.names) >0){
    temp <- matrix(0, nrow = l2, ncol = length(missing2.names))
    temp<-as.data.frame(temp)
    colnames(temp)<-missing2.names
    m2.temp<-cbind(matrix2, temp)
    temp2 <- matrix(0, nrow = length(missing2.names), ncol = length(total.names))
    temp2<-as.data.frame(temp2)
    colnames(temp2)<-colnames(m2.temp)
    rownames(temp2)<-missing2.names
    m2<-rbind(m2.temp, temp2)
  } else{
    m2<-matrix2
  }
  #this goes through and addes rows/col of 0 for missing in matrix1
  if(length(missing1.names) >0){
    temp <- matrix(0, nrow = l1, ncol = length(missing1.names))
    temp<-as.data.frame(temp)
    colnames(temp)<-missing1.names
    m1.temp<-cbind(matrix1, temp)
    temp2 <- matrix(0, nrow = length(missing1.names), ncol = length(total.names))
    temp2 <- as.data.frame(temp2)
    colnames(temp2) <- colnames(m1.temp)
    rownames(temp2)<-missing1.names
    m1<-rbind(m1.temp, temp2)
  } else{
    m1<-matrix1
  }
  matrix1b<-as(m1[total.names, total.names],"matrix")
  rownames(matrix1b)<-colnames(matrix1b)<-row.names(m1)
  matrix2b<-as(m2[total.names, total.names],"matrix")
  rownames(matrix2b)<-colnames(matrix2b)<-row.names(m2)
  
  matrix3=matrix(0, nrow(matrix1b),nrow(matrix1b))
  for (i in 1:nrow(matrix1b)) {
    for(j in 1:nrow(matrix2b)) {
      if((matrix1b[i,j]==-1 & matrix2b[i,j]==0)|(matrix1b[i,j]==0 & matrix2b[i,j]==-1))
        matrix3[i,j]=-1
      else if ((matrix1b[i,j]==0 & matrix2b[i,j]==1)|(matrix1b[i,j]==1 & matrix2b[i,j]==0))
        matrix3[i,j]=1
      else if(matrix1b[i,j]==1 & matrix2b[i,j]==1)
        matrix3[i,j]=1
      else if(matrix1b[i,j]==-1 & matrix2b[i,j]==-1)
        matrix3[i,j]=-1
      else
        matrix3[i,j]=0
      
      
    }
  }
  rownames(matrix3)<-colnames(matrix3)<-row.names(matrix1b)
  matrix3
}