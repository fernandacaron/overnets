#degree -> indegree (i.e., number of connections that point inward at a vertex),
#outdegree (i.e., number of connections that originate at a vertex and point 
#outward to other vertices, or freeman (i.e., total).
#         digraph = relationship between two vertices are mutual
metRangeOver <- function(matrix){
  library(sna)
  
  deg<-degree(matrix, gmode = "graph")
  bet<-betweenness(matrix, gmode = "graph")
  clo<-closeness(matrix, gmode = "graph")
  
  res<-list()
  res$degree<-deg
  res$betweeness<-bet
  res$closeness<-clo
  
  return(res)
}

