##Figure
#library(corrplot)
##feature selsction
library(tidyverse)
library(Seurat)
library(ggplot2)
library(umap)
library(M3Drop)
library(ggpubr)

klein.dat <- read.csv("klein.csv",header = T,row.names = 1)
klein.label <- read.csv("klein_label.csv",row.names = NULL)
klein.label <- as.character(as.matrix(klein.label))
l<- factor(klein.label)


set.seed(1234567)
####
raw <- function(dat = klein.dat){
  set.seed(1234567)
  X <- t(klein.dat)
  Y.umap_raw <- umap(X)
  a<- as.data.frame(Y.umap_raw$layout)
  a$l <- l
  return(a)
}


###Random
RAN<- function(dat = klein.dat,n_feature = 2000){
  set.seed(1234567)
  klein_ran <- klein.dat[sample(rownames(klein.dat),n_feature),]
  
  X <- t(klein_ran)
  a <- as.data.frame(umap(X)$layout)
  a$l <- l
  return(a)
}


####HVG
HVG<- function(dat = klein.dat,n_feature = 2000){
  set.seed(1234567)
  seurat_klein <- CreateSeuratObject(counts = klein.dat)
  seurat_klein <- FindVariableFeatures(seurat_klein,nfeatures = n_feature)
  seurat_hvg<- klein.dat[which(VariableFeatures(seurat_klein)%in% rownames(klein.dat)),]
  
  X <- t(seurat_hvg)
  print(dim(X))
  a <- as.data.frame(umap(X)$layout)
  a$l <- l
  return(a)
}


###HEG
HEG<- function(dat = klein.dat,n_feature = 2000){
  set.seed(1234567)
  ep <- apply(klein.dat,1,mean)
  top_heg <- as.matrix(head(ep,n_feature))
  klein_heg <- klein.dat[rownames(top_heg),]
  
  X <- t(klein_heg)
  a <- as.data.frame(umap(X)$layout)
  a$l <- l
  return(a)
}


###giniclust

calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
  set.seed(1234567)
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind = is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x = x[!na.ind]
  n = length(x)
  mu = mean(x)
  N = if (unbiased) n * (n - 1) else n * n
  ox = x[order(x)]
  dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
  dsum / (mu * N)
}

GINI<- function(dat=klein.dat,n_feature=2000){
  set.seed(1234567)
  
  gini = apply(as.data.frame(klein.dat), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
  GiniIndex = as.data.frame(cbind(1:dim(klein.dat)[1], gini))
  top_gini <- rownames(head(GiniIndex[order(GiniIndex[,2]),],n_feature))
  klein_gini <- klein.dat[top_gini,]
  
  X <- t(klein_gini)
  a <- as.data.frame(umap(X)$layout)
  a$l <- l
  return(a)


}


#####m3DROP
M3 <- function(dat = klein.dat, mt_threshold=0.05){
  set.seed(1234567)
  keep <- rownames(M3DropFeatureSelection(klein.dat,mt_method="none",mt_threshold=mt_threshold,suppress.plot=T))
  klein_m3 <-klein.dat[keep,]
  
  X<- t(klein_m3)
  a <- as.data.frame(umap(X)$layout)
  a$l <- l
  return(a)
}

##########ggplot2

n_feature = 2000
mt_threshold=0.05
r_1 <- raw()
r_2 <- RAN(n_feature = n_feature)
r_3 <- HVG(n_feature = n_feature)
r_4 <- HEG(n_feature = n_feature)
r_5 <- GINI(n_feature = n_feature)
r_6 <- M3(mt_threshold= mt_threshold)

p_1 <- ggplot(data=r_1, mapping=aes(r_1[,1],y=r_1[,2],colour=factor(r_1[,3])))+
  geom_point()+  
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Without feature selection")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_2 <- ggplot(data=r_2, mapping=aes(r_2[,1],y=r_2[,2],colour=factor(r_2[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 2,000 genes randomly")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_3 <- ggplot(data=r_3, mapping=aes(r_3[,1],y=r_3[,2],colour=factor(r_3[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 2,000 genes using high variable genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_4 <- ggplot(data=r_4, mapping=aes(r_4[,1],y=r_4[,2],colour=factor(r_4[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 2,000 genes using high expression genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_5 <- ggplot(data=r_5, mapping=aes(r_5[,1],y=r_5[,2],colour=factor(r_5[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 2,000 genes using Gini index")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_6 <- ggplot(data=r_6, mapping=aes(r_6[,1],y=r_6[,2],colour=factor(r_6[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 2,000 genes using M3drop")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))




n_feature = 1000
mt_threshold=0.01

r_8 <- RAN(n_feature = n_feature)
r_9 <- HVG(n_feature = n_feature)
r_10 <- HEG(n_feature = n_feature)
r_11 <- GINI(n_feature = n_feature)
r_12 <- M3(mt_threshold=  mt_threshold)

#p_7 <- ggplot(data=r_1, mapping=aes(r_1[,1],y=r_1[,2],colour=factor(r_1[,3])))+
#  geom_point()
p_8 <- ggplot(data=r_8, mapping=aes(r_8[,1],y=r_8[,2],colour=factor(r_2[,3])))+
  geom_point()  +
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 1,000 genes randomly")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_9 <- ggplot(data=r_9, mapping=aes(r_9[,1],y=r_9[,2],colour=factor(r_3[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 1,000 genes using high variable genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_10 <- ggplot(data=r_10, mapping=aes(r_10[,1],y=r_10[,2],colour=factor(r_4[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 1,000 genes using high expression genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_11 <- ggplot(data=r_11, mapping=aes(r_11[,1],y=r_11[,2],colour=factor(r_5[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Selecte 1,000 genes using Gini index")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_12 <- ggplot(data=r_12, mapping=aes(r_12[,1],y=r_12[,2],colour=factor(r_6[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 1,000 genes using M3drop")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))



n_feature = 500
mt_threshold=0.005

r_14 <- RAN(n_feature = n_feature)
r_15 <- HVG(n_feature = n_feature)
r_16 <- HEG(n_feature = n_feature)
r_17 <- GINI(n_feature = n_feature)
r_18 <- M3(mt_threshold=  mt_threshold)

#p_13 <- ggplot(data=r_1, mapping=aes(r_1[,1],y=r_1[,2],colour=factor(r_1[,3])))+
#  geom_point()
p_14 <- ggplot(data=r_14, mapping=aes(r_14[,1],y=r_14[,2],colour=factor(r_2[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 500 genes randomly")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_15<- ggplot(data=r_15, mapping=aes(r_15[,1],y=r_15[,2],colour=factor(r_3[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 500 genes using high variable genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_16 <- ggplot(data=r_16, mapping=aes(r_16[,1],y=r_16[,2],colour=factor(r_4[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 500 genes using high expression genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_17 <- ggplot(data=r_17, mapping=aes(r_17[,1],y=r_17[,2],colour=factor(r_5[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 500 genes using Gini index")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_18 <- ggplot(data=r_18, mapping=aes(r_18[,1],y=r_18[,2],colour=factor(r_6[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 500 genes using M3drop")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))


n_feature = 200
mt_threshold=0.001

r_20 <- RAN(n_feature = n_feature)
r_21 <- HVG(n_feature = n_feature)
r_22 <- HEG(n_feature = n_feature)
r_23 <- GINI(n_feature = n_feature)
r_24 <- M3(mt_threshold=  mt_threshold)

#p_19 <- ggplot(data=r_1, mapping=aes(r_1[,1],y=r_1[,2],colour=factor(r_1[,3])))+
#  geom_point()
p_20 <- ggplot(data=r_20, mapping=aes(r_20[,1],y=r_20[,2],colour=factor(r_2[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 200 genes randomly")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_21 <- ggplot(data=r_21, mapping=aes(r_21[,1],y=r_21[,2],colour=factor(r_3[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 200 genes using high variable genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_22 <- ggplot(data=r_22, mapping=aes(r_22[,1],y=r_22[,2],colour=factor(r_4[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 200 genes using high expression genes")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_23 <- ggplot(data=r_23, mapping=aes(r_23[,1],y=r_23[,2],colour=factor(r_5[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 200 genes using Gini index")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))
p_24 <- ggplot(data=r_24, mapping=aes(r_24[,1],y=r_24[,2],colour=factor(r_6[,3])))+
  geom_point()+
  theme_classic()+
  labs(y="UMAP2")+
  labs(x="UMAP1")+
  labs(color = "  Cell type")+
  ggtitle("Select 200 genes using M3drop")+
  theme(plot.title = element_text(hjust = 0.5,size=14,face="bold"))

img <- ggarrange(p_1,p_2,p_3,p_4,p_5,p_6,
                 p_1,p_8,p_9,p_10,p_11,p_12,
                 p_1,p_14,p_15,p_16,p_17,p_18,
                 p_1,p_20,p_21,p_22,p_23,p_24,
                 nrow=4,ncol=6)

img

ggsave("D:/BIB_figure/S2.eps", width = 80, height = 60, units = "cm")

dev.off()






