library(igraph) # Greedy, greedy2, dsatur, rlf
library(tmaptools) # map_coloring (as_adj_list)
library(BiocManager) # sequential_vertex_coloring (graph::graphAM)

vertex.color.greedy <- function(g, seq=V(g)){
  grps <- list()
  for(i in seq){
    if(length(grps)==0){grps[[1]]<-i}else{
      for(j in 1:length(grps)){
        if(sum(grps[[j]]%in%neighbors(g,i))==0){grps[[j]]<- c(grps[[j]],i);break}
      }
      if(sum(i%in%unlist(grps))==0){grps[[(length(grps)+1)]] <- i}
    }
  }
  grps2 <- list(no._de_grupos=length(grps),grupos=grps)
  return(grps2)
}

vertex.color.greedy.2 <- function(g, seq=V(g),times=1){
  n=0
  while(n<times){
    grps <- list()
    for(i in seq){
      if(length(grps)==0){grps[[1]]<-i}else{
        for(j in 1:length(grps)){
          if(sum(grps[[j]]%in%neighbors(g,i))==0){grps[[j]]<- c(grps[[j]],i);break}
        }
        if(sum(i%in%unlist(grps))==0){grps[[(length(grps)+1)]] <- i}
      }
    }
    for(i in 1:length(grps)){
      if(length(grps[[i]])>1){grps[[i]] <- sample(grps[[i]])}
    }
    seq <- unlist(sample(grps),use.names = F)
    n=n+1
  }
  
  grps2 <- list(no._de_grupos=length(grps),grupos=grps)
  return(grps2)
}

vertex.color.dsatur <- function(g){
  grps <- list()
  verts <- V(g)
  d.satur <- rep(0,times=length(verts))
  while(length(verts)!=0){
    d <- data.frame(as.vector(verts),d.satur,deg=degree(g,verts))
    seq <- d[order(d$d.satur,d$deg,decreasing=T),1]
    for(i in seq){
      if(length(grps)==0){grps[[1]]<-i;verts<-verts[verts!=i]}else{
        for(j in 1:length(grps)){
          if(sum(grps[[j]]%in%neighbors(g,i))==0){grps[[j]]<- c(grps[[j]],i);verts<-verts[verts!=i];break}
        }
        if(sum(i%in%unlist(grps))==0){grps[[(length(grps)+1)]] <- i;verts<-verts[verts!=i]}
      }
    }
    d.satur <- c()
    for(i in verts){
      s <- 0
      for(j in grps){
        if(sum(j%in%neighbors(g,i))>0){s <- s+1}
      }
      d.satur <- c(d.satur,s)
    }
  }
  grps2 <- list(no._de_grupos=length(grps),grupos=grps)
  return(grps2)
}

vertex.color.RLF <- function(g){
  grps <- list()
  while(sum(V(g)%in%unlist(grps))!=length(V(g))){
    X <- V(g)[!V(g)%in%unlist(grps)]
    grp.new <- c()
    while(length(X)>0){
      if(length(X)==1){v<-X}else{v <- sample(X,1,replace=F)}
      grp.new <- c(grp.new,v)
      X <- X[X!=v]
      X <- X[!X%in%neighbors(g,v)]
    }
    grps[[length(grps)+1]] <- grp.new
  }
  grps2 <- list(no._de_grupos=length(grps),grupos=grps)
  return(grps2)
}

############## SIMULATION ################
##### GREEDY #####
p=0.20
set.seed(6969)
tempos.greedy <- c()
tamanhos.greedy <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.greedy(g,seq=sample(V(g)))
  fin <- Sys.time()
  tempos.greedy <- c(tempos.greedy, fin-ini) ; tamanhos.greedy <- c(tamanhos.greedy, sol[[1]])
}

greedy.2 <- cbind(tempos.greedy, tamanhos.greedy)

p=0.50
set.seed(6969)
tempos.greedy <- c()
tamanhos.greedy <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.greedy(g,seq=sample(V(g)))
  fin <- Sys.time()
  tempos.greedy <- c(tempos.greedy, fin-ini) ; tamanhos.greedy <- c(tamanhos.greedy, sol[[1]])
}

greedy.5 <- cbind(tempos.greedy, tamanhos.greedy)

p=0.90
set.seed(6969)
tempos.greedy <- c()
tamanhos.greedy <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.greedy(g,seq=sample(V(g)))
  fin <- Sys.time()
  tempos.greedy <- c(tempos.greedy, fin-ini) ; tamanhos.greedy <- c(tamanhos.greedy, sol[[1]])
}

greedy.9 <- cbind(tempos.greedy, tamanhos.greedy)

comp.greedy <- data.frame(vertex=seq(10,500,10),tem.2=greedy.2[,1],sol.2=greedy.2[,2],
                          tem.5=greedy.5[,1],sol.5=greedy.5[,2],
                          tem.9=greedy.9[,1],sol.9=greedy.9[,2])
write.table(comp.greedy, "comparison_greedy.txt", col.names = TRUE, row.names = FALSE)

##### DSATUR #####
p=0.20
set.seed(6969)
tempos.dsatur <- c()
tamanhos.dsatur <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.dsatur(g)
  fin <- Sys.time()
  tempos.dsatur <- c(tempos.dsatur, fin-ini) ; tamanhos.dsatur <- c(tamanhos.dsatur, sol[[1]])
}

dsatur.2 <- cbind(tempos.dsatur, tamanhos.dsatur)

p=0.50
set.seed(6969)
tempos.dsatur <- c()
tamanhos.dsatur <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.dsatur(g)
  fin <- Sys.time()
  tempos.dsatur <- c(tempos.dsatur, fin-ini) ; tamanhos.dsatur <- c(tamanhos.dsatur, sol[[1]])
}

dsatur.5 <- cbind(tempos.dsatur, tamanhos.dsatur)

p=0.90
set.seed(6969)
tempos.dsatur <- c()
tamanhos.dsatur <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.dsatur(g)
  fin <- Sys.time()
  tempos.dsatur <- c(tempos.dsatur, fin-ini) ; tamanhos.dsatur <- c(tamanhos.dsatur, sol[[1]])
}

dsatur.9 <- cbind(tempos.dsatur, tamanhos.dsatur)

comp.dsatur <- data.frame(vertex=seq(10,500,10),tem.2=dsatur.2[,1],sol.2=dsatur.2[,2],
                          tem.5=dsatur.5[,1],sol.5=dsatur.5[,2],
                          tem.9=dsatur.9[,1],sol.9=dsatur.9[,2])
write.table(comp.dsatur, "comparison_dsatur.txt", col.names = TRUE, row.names = FALSE)

##### RLF #####
p=0.20
set.seed(6969)
tempos.rlf <- c()
tamanhos.rlf <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.RLF(g)
  fin <- Sys.time()
  tempos.rlf <- c(tempos.rlf, fin-ini) ; tamanhos.rlf <- c(tamanhos.rlf, sol[[1]])
}

rlf.2 <- cbind(tempos.rlf, tamanhos.rlf)

p=0.50
set.seed(6969)
tempos.rlf <- c()
tamanhos.rlf <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.RLF(g)
  fin <- Sys.time()
  tempos.rlf <- c(tempos.rlf, fin-ini) ; tamanhos.rlf <- c(tamanhos.rlf, sol[[1]])
}

rlf.5 <- cbind(tempos.rlf, tamanhos.rlf)

p=0.90
set.seed(6969)
tempos.rlf <- c()
tamanhos.rlf <- c()
for(n in seq(10,500,10)){
  g <- sample_gnp(n,p)
  ini <- Sys.time()
  sol <- vertex.color.RLF(g)
  fin <- Sys.time()
  tempos.rlf <- c(tempos.rlf, fin-ini) ; tamanhos.rlf <- c(tamanhos.rlf, sol[[1]])
}

rlf.9 <- cbind(tempos.rlf, tamanhos.rlf)

comp.rlf <- data.frame(vertex=seq(10,500,10),tem.2=rlf.2[,1],sol.2=rlf.2[,2],
                          tem.5=rlf.5[,1],sol.5=rlf.5[,2],
                          tem.9=rlf.9[,1],sol.9=rlf.9[,2])
write.table(comp.rlf, "comparison_rlf.txt", col.names = TRUE, row.names = FALSE)

########### COMPARISON ############
greedy <- read.table("comparison_greedy.txt", header = T)
dsatur <- read.table("comparison_dsatur.txt", header = T)
rlf <- read.table("comparison_rlf.txt", header = T)

#### Tempos para p=0.2 ####
plot(dsatur$vertex, dsatur$tem.2, type='l', lwd=3, col="red", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação do Tempo Computacional para p=0.2")
lines(greedy$vertex, greedy$tem.2, type='l', lwd=3, col='blue')
lines(rlf$vertex, rlf$tem.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Tempos para p=0.5 ####
plot(dsatur$vertex, dsatur$tem.5, type='l', lwd=3, col="red", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação do Tempo Computacional para p=0.5")
lines(greedy$vertex, greedy$tem.5, type='l', lwd=3, col='blue')
lines(rlf$vertex, rlf$tem.5, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Tempos para p=0.9 ####
plot(dsatur$vertex, dsatur$tem.9, type='l', lwd=3, col="red", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação do Tempo Computacional para p=0.9")
lines(greedy$vertex, greedy$tem.9, type='l', lwd=3, col='blue')
lines(rlf$vertex, rlf$tem.9, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução para p=0.2 ####
plot(greedy$vertex, greedy$sol.2, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de cores da solução", main="Comparação da Solução para p=0.2")
lines(dsatur$vertex, dsatur$sol.2, type='l', lwd=3, col='red')
lines(rlf$vertex, rlf$sol.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução para p=0.5 ####
plot(greedy$vertex, greedy$sol.5, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de cores da solução", main="Comparação da Solução para p=0.5")
lines(dsatur$vertex, dsatur$sol.5, type='l', lwd=3, col='red')
lines(rlf$vertex, rlf$sol.5, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução para p=0.9 ####
plot(greedy$vertex, greedy$sol.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de cores da solução", main="Comparação da Solução para p=0.9")
lines(dsatur$vertex, dsatur$sol.9, type='l', lwd=3, col='red')
lines(rlf$vertex, rlf$sol.9, type='l', lwd=3, col='black')
legend("topleft", legend=c('GREEDY',"DSatur","RLF"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Tempos Greedy ####
plot(greedy$vertex, greedy$tem.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação dos Tempos - GREEDY")
lines(greedy$vertex, greedy$tem.5, type='l', lwd=3, col='red')
lines(greedy$vertex, greedy$tem.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Tempos DSatur ####
plot(dsatur$vertex, dsatur$tem.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação dos Tempos - DSatur")
lines(dsatur$vertex, dsatur$tem.5, type='l', lwd=3, col='red')
lines(dsatur$vertex, dsatur$tem.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Tempos RLF ####
plot(rlf$vertex, rlf$tem.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Tempo Computacional (segundos)", main="Comparação dos Tempos - RLF")
lines(rlf$vertex, rlf$tem.5, type='l', lwd=3, col='red')
lines(rlf$vertex, rlf$tem.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução Greedy ####
plot(greedy$vertex, greedy$sol.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de Cores", main="Comparação das Soluções - GREEDY")
lines(greedy$vertex, greedy$sol.5, type='l', lwd=3, col='red')
lines(greedy$vertex, greedy$sol.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução DSatur ####
plot(dsatur$vertex, dsatur$sol.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de Cores", main="Comparação das Soluções - DSatur")
lines(dsatur$vertex, dsatur$sol.5, type='l', lwd=3, col='red')
lines(dsatur$vertex, dsatur$sol.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)

#### Solução RLF ####
plot(rlf$vertex, rlf$sol.9, type='l', lwd=3, col="blue", xlab="Número de Vértices",
     ylab="Número de Cores", main="Comparação das Soluções - RLF")
lines(rlf$vertex, rlf$sol.5, type='l', lwd=3, col='red')
lines(rlf$vertex, rlf$sol.2, type='l', lwd=3, col='black')
legend("topleft", legend=c('p=0.9',"p=0.5","p=0.2"), col=c('blue','red','black'), lty=1, lwd=3, cex=0.7)
