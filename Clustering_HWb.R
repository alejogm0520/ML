library(readxl)
library(ggplot2)
library(StatMatch)
library(reshape2)
library(dendextend)
library(GMD)
#library(NbClust)
library(cluster)
library(ICGE)
library(mclust)
crimen <- read_excel("./Analisis Multivariado/crimen.xlsx")
# crimen <- read_excel("~/Descargas/crimen.xlsx")
View(crimen)

hcd <- function(X, labs, metric, meth, ps){
  Distances <- dist(X, method=metric, diag = TRUE, upper = TRUE, p = ps)
  hc <- hclust(Distances,method = meth)
  hc$labels <- labs
  my_list <- list("hcluster" = hc, "dist" = Distances)
  return(my_list) 
}

dendo <- function(hc, nc, title='Dendogram Cluster'){
  dend <- as.dendrogram(hc)
  dend <- color_branches(dend, k=nc, col=rainbow(nc))
  par(mar=c(4,4,4,5)+0.1)
  plot(dend, horiz = TRUE, main=title)
}


WSS <- function(distance, labs){
  if (var(labs)== 0) {
    wss <- 0
  }
  else{
    a <- deltas(distance, labs)
    a <- a/nrow(a)
    
    wss <- sum(a)/2  
  }
  return(wss)
}

WDD <- function(distance, clust, k){
  labs <- cutree(clust,k)
  
  if (var(labs)== 0) {
    wss <- 0
  }
  else{
    a <- deltas(distance, labs)
    a <- a/nrow(a)
    
    wss <- sum(a)/2
  }
  wss <- abs(wss)*k^2
  return(wss)
}
WTT <- function(distance, clust, k){
  labs <- cutree(clust,k)
  
  if (var(labs)== 0) {
    wss <- 0
  }
  else{
    a <- deltas(distance, labs)
    a <- a/nrow(a)
    wss <- sum(a)/2
  }
  wss <- wss
  return(wss)
}

DIFK <- function(distance, clust, k){
  
  d <- 7
  lab <- cutree(clust,k-1)
  W.k0 <- WSS(distance, lab)
  lab <- cutree(clust,k)
  W.k1 <- WSS(distance, lab)
  difK <- W.k0*((k-1)^(2/d)) -  W.k1*((k)^(2/d))
  return(difK)
}

KL <- function(distance, clust, k){
  kl <- abs( DIFK(distance, clust, k)/DIFK(distance, clust, k+1))
  return(kl)
}


nClust <- function(distance, hcluster, nMax, title='Cluster Number'){
  klindex <- sapply(1:nMax, function(x) WTT(hc$dist, hc$hcluster, x))
  #clust <- cutree(hcluster, k=1:nMax)
  #Wss <- sapply(1:nMax, function(x)  css(distance, array(clust[,x]))$totwss)
  plot(2:nMax, klindex[2:nMax], ylab = 'Indice de WSS', xlab = 'Número de Cluster', main=title)
  #abline(v=which.min(diff(Wss))+1, lty=1, col = 'red')
  v2=which.max(diff(append(0,diff(append(0,klindex))))[1:(nMax-2)])
  abline(v=v2, lty=1, col = 'red')
  return(v2)
  }

nClust3 <- function(distance, hcluster, nMax, title='Cluster Number'){
  klindex <- sapply(2:nMax, function(x) KL(distance, hcluster, x))
  plot(2:nMax, log(klindex), ylab = 'Indice KL', xlab = 'Número de Cluster', main=title)
  v=which.min(klindex[3:nMax])+3
  abline(v=v, lty=1, col = 'red')
  return(v)
  
}

nClust4 <- function(distance, hcluster, nMax, title='Cluster Number'){
  klindex <- sapply(1:nMax, function(x) WDD(distance, hcluster, x))
  plot(2:nMax, log1p(klindex)[2:nMax], ylab = 'Indice M', xlab = 'Número de Cluster', main=title)
#  abline(v=which.max(diff((diff((klindex))))), lty=1, col = 'red')
#  abline(v=which.max(diff(append(diff((klindex)),0))), lty=1, col = 'blue')
  v=which.max(diff(append(0,diff(append(0,klindex))))[1:(nMax-3)])
  abline(v=v, lty=1, col = 'red')
  return(v)
}


nClust2 <- function(X, distanc, hcm,  nMax, ind,  title='Cluster Number'){
  res<-NbClust(X, distance = distanc, min.nc=1, max.nc=nMax, 
               method = hcm, index = ind)
  plot(1:7, res$All.index, ylab = 'Índice de KL', xlab = 'Número de Cluster', main=title)
  abline(v=res$Best.nc[1], lty=1, col = 'red')
}



par(mfrow = c(1,3))

## Single Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'single', 2)

### Número de Cluster

cairo_ps(paste('Single','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Single','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Single','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
#nClust2(crimen[,2:8],"euclidean", "single",7, 'marriot', "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Single','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Vinculo Unico')
dev.off()


## Complete Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'complete', 2)

### Número de Cluster

cairo_ps(paste('Complete','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Complete','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Complete','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Complete','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Vinculo Completo')
dev.off()



## Average Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'average', ps=2)

### Número de Cluster

cairo_ps(paste('Average','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Average','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Average','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Average','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Vinculo Promedio')
dev.off()


## Ward Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'ward.D2', ps=2)

### Número de Cluster

cairo_ps(paste('Ward','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Ward','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Ward','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Ward','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Ward')
dev.off()



## Centroid Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'centroid', ps=2)

### Número de Cluster

cairo_ps(paste('Centroid','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Centroid','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Centroid','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Centroid','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster,res$mids[which.max(res$counts)]+0.75, title='Dendograma Centroide')
dev.off()



## Median Link - Euclideana
hc <- hcd(crimen[,2:8], crimen$Ciudad, 'euclidean', 'median', ps=2)

### Número de Cluster

cairo_ps(paste('Median','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Median','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Median','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Median','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Mediana')
dev.off()



## Beta Flexible B=3 - Euclideana
hc$dist <- dist(crimen[,2:8], method = "euclidean", diag = TRUE, upper = TRUE, p = 2)
hc$hcluster <- agnes(crimen[,2:8], diss = inherits(crimen[,2:8], "dist"), metric = "euclidean", method = "flexible", par.method = 3)
hc$hcluster <- as.hclust(hc$hcluster)
### Número de Cluster

cairo_ps(paste('Beta','_','NC','_','WS','.eps',sep=""), width = 5, height = 5)
r1<-nClust(hc$dist, hc$hcluster, 8, 'Método de WSS')
dev.off()
cairo_ps(paste('Beta','_','NC','_','KL','.eps',sep=""), width = 5, height = 5)
r2<-nClust3(hc$dist, hc$hcluster, 8, "Método de Krzanowski-Lai")
dev.off()
cairo_ps(paste('Beta','_','NC','_','MA','.eps',sep=""), width = 5, height = 5)
r3<-nClust4(hc$dist, hc$hcluster, 8, "Método de Marriot")
dev.off()

### Dendograma
hc$hcluster$height <- log(hc$hcluster$height)
res <- hist(c(r1,r2,r3))

cairo_ps(paste('Beta','_','DN','.eps',sep=""), width = 19, height = 5)
dendo(hc$hcluster, res$mids[which.max(res$counts)], title='Dendograma Beta Flexible')
dev.off()


## No Jerarquicos

dioxido <- read_excel("./Analisis Multivariado/dioxido.xlsx")
