library(tidyverse)
library(cluster)
library(clusterSim)
library(ade4)

# Linking long-term dietary patterns with gut microbial enterotypes
# Dirichlet multinomial mixtures: generative models for microbial metagenomics
# https://enterotype.embl.de/enterotypes.html
# https://github.com/abikoushi/enigma

m3 <- read_delim('/Users/congliu/prs/mv34/HZKS01_1148_mv34_otu_qiime_silva99v2_m3_std.txt', delim = '\t') %>% 
  dplyr::select(-taxonomy) %>% filter(str_detect(`OTU ID`, 'g__'))

OTU_ID <- dplyr::select(m3, `OTU ID`)
foo <- m3 %>% dplyr::select(-`OTU ID`) %>% mutate_each(function(x){x/sum(x)})
bar <- bind_cols(OTU_ID, foo) %>% 
  extract(`OTU ID`, into = 'OTU_ID', regex = 'g__([A-Z].*)')
write_delim(bar, '/Users/congliu/prs/hbrm_Alveoli/genus.txt', delim = '\t')

# 11年文献 -------------------------------------------------------------------

#df <- read.table("/Users/congliu/prs/enterotype/MetaHIT_SangerSamples.genus.txt", header=T, row.names=1, dec=".", sep="\t")
#df <- df[-1, ]
df <- read.table('/Users/congliu/prs/hbrm_Alveoli/genus.txt', header = T, row.names = 1, sep = '\t', dec = ".")

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
#df <- noise.removal(df, percent=0.01)

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data.dist=dist.JSD(df)

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

#data.cluster=pam.clustering(data.dist, k=3)
#确定最优的k，Calinski-Harabasz index
#nclusters = index.G1(t(df), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL
for (k in 1:20) {
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(df),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
data.cluster=pam.clustering(data.dist, k=3)
obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
obs.silhouette

obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

###BCA 验证聚类情况 对于短序列测序情况 推荐去除
noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}
#data.denoized=noise.removal(df, percent=0.01)
obs.pca=dudi.pca(data.frame(t(df)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
#s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))


# 对比 ----------------------------------------------------------------------

foo <- pam(as.dist(data.dist), k=3, diss=TRUE)
summary(foo)
sample_type <- as.data.frame(foo$clustering)

get_top10 <- function(tax){
  sample_type_1 <- subset(sample_type, foo$clustering==tax)
  df_1 <- df[rownames(sample_type_1)]
  df_1['total'] <- rowSums(df_1)
  bar <- df_1[order(df_1$total, decreasing = T), ]
  df_top10 <- subset(bar[1:10, ], select = total)
  return(df_top10)
}
get_top10(1)

get_top10_m3 <- function(tax){
  sample_type_1 <- subset(sample_type, foo$clustering==tax)
  df_1 <- m3[rownames(sample_type_1)]
  df_1['total'] <- rowSums(df_1)
  df_1 <- bind_cols(dplyr::select(m3, `OTU ID`), df_1)
  bar <- df_1[order(df_1$total, decreasing = T), ]
  df_top10 <- subset(bar[1:10, ], select = c(`OTU ID`,total))
  df_last <- subset(bar[11:nrow(bar),], select = c(`OTU ID`,total))
  df.final <- rbind(df_top10, c("Others", colSums(df_last[,2:length(df_last)])))
  df.final <- df.final %>% mutate_at(vars(total), as.double) %>% 
    mutate_at(vars(total), funs(./sum(.)))
  return(df.final)
}
get_top10_m3(1)

define_enterotype <- function(m3){
  m3_rate <- m3 %>% arrange(desc(PRS003170019)) %>% mutate_if(is.double, ~(./sum(.))) %>% 
    rename(OTU_ID=`OTU ID`)
  for (i in colnames(m3_rate)[2:length(m3_rate)]) {
    one_sample <- m3_rate[c('OTU_ID', i)]
    B <- filter(one_sample, OTU_ID == 'g__Bacteroides')[[1,2]]
    P <- filter(one_sample, OTU_ID == 'g__Prevotella')[[1,2]]
    if((B > 0.4) && (B > P)){
      print(paste(i, 'TypeB', sep = '\t'))
      cat(paste(i, 'TypeB', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype.txt", sep='\n',append = TRUE)
    } else if ((P >= 0.3) && (P >= B)){
      print(paste(i, 'TypeP', sep = '\t'))
      cat(paste(i, 'TypeP', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype.txt", sep='\n',append = TRUE)
    } else {
      print(paste(i, 'TypeIII', sep = '\t'))
      cat(paste(i, 'TypeIII', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype.txt", sep='\n',append = TRUE)
    }
  }
}

define_enterotype(m3)

define_enterotype_2 <- function(m3){
  m3_rate <- m3 %>% arrange(desc(PRS003170019)) %>% mutate_if(is.double, ~(./sum(.))) %>% 
    rename(OTU_ID=`OTU ID`)
  for (i in colnames(m3_rate)[2:length(m3_rate)]) {
    one_sample <- m3_rate[c('OTU_ID', i)]
    B <- filter(one_sample, OTU_ID == 'g__Bacteroides')[[1,2]]
    P <- filter(one_sample, OTU_ID == 'g__Prevotella')[[1,2]]
    if((B > 0.2) && (B > P)){
      print(paste(i, 'TypeB', sep = '\t'))
      cat(paste(i, 'TypeB', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype_2.txt", sep='\n',append = TRUE)
    } else if ((P >= 0.3) && (P >= B)){
      print(paste(i, 'TypeP', sep = '\t'))
      cat(paste(i, 'TypeP', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype_2.txt", sep='\n',append = TRUE)
    } else {
      print(paste(i, 'TypeIII', sep = '\t'))
      cat(paste(i, 'TypeIII', sep = '\t'), file = "/Users/congliu/prs/hbrm_Alveoli/enterotype_2.txt", sep='\n',append = TRUE)
    }
  }
}
define_enterotype_2(m3)
# 18年文献 -------------------------------------------------------------------















