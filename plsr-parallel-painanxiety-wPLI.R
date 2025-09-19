rm(list = ls())

# Load the 'pls' package
library(pls)
library(npreg)
library(caret)
library(parallel)

#setwd("/Users/majid/Projects/PTSD/Connectivity/GroupLevel/")

#importing data
load("DataFormatR.RData")

rois <- dimnames(AEC_Conn)[[2]]

eff.comp.numb <- function(data){
  dm <- dim(data)[1]
  ncompo <- vector()
  error <- vector()
  for(ncomp_ in 1:5){
    err_ <- vector()
    k <- 0
    for(repea in 1:5){
      folds <- createFolds(seq(dm),k = 10)
      for(fold in folds){
        k <- k + 1
        traindata <- data[-fold,]
        testdata <- data[fold,]
        pls_model <- plsr(Y~., data=traindata, ncomp = ncomp_ , scale=T , center=T)
        pred <- predict(pls_model,testdata , ncomp=ncomp_)
        act_ <- testdata$Y
        pre_ <- as.numeric(pred)
        err_[k] <- mean(abs(act_ - pre_))
      }
    }
    error[ncomp_] <- mean(err_)
  }
  return(which.min(error))
}

resample.data <- function(data){
  Y_ <- data$Y
  Y_ <- Y_[sample(length(Y_))]
  redata <- data
  redata$Y <- Y_
  return(redata)
}

plsr.loading <- function(data,resnum){
  # Fit a PLSR model with optimized number of components
  ncomp_ <- eff.comp.numb(data)
  pls_model <- plsr(Y ~ ., data = data, ncomp = ncomp_)
  loadings <- coef(pls_model, type = "loading", comp = ncomp_)
  actu_loads <- as.vector(loadings)
  
  #resampling
  resa_loads <- matrix(ncol = 630,nrow = resnum)
  for(re in 1:resnum){
    data_ <- resample.data(data)
    ncomp_ <- eff.comp.numb(data_)
    pls_model_ <- plsr(Y ~ ., data = data_, ncomp = ncomp_)
    loadings_ <- coef(pls_model_, type = "loading", comp = ncomp_)
    resa_loads[re,] <- as.vector(loadings_)
    print(paste("resampleing" , re) )
  }
  
  pv <- vector()
  for(i in 1:630){
    pv[i] <- sum(abs(resa_loads[,i]) >= abs(actu_loads[i]) , na.rm=T)/resnum
  }
  
  coef <- pcoef <- matrix(ncol = 36,nrow = 36)
  colnames(coef) <- rownames(coef) <- 
    colnames(pcoef) <- rownames(pcoef) <- rois
  
  coef[upper.tri(coef)][] <- actu_loads[]
  pcoef[upper.tri(pcoef)][] <- pv[]
  
  loadingvec <- list(coef,pcoef)
  
  return(loadingvec)
}

pls.modeling <- function(fqb){
  
  X <- matrix(ncol = 630, nrow = 99)
  for(i in 1:99){
    X[i,] <- Connec[,,fqb,i][upper.tri(Connec[,,fqb,i])]
  }
  # Example data (replace this with your own data)
  set.seed(123)
  data <- data.frame(Y = YY, X)
  dim(data)
  
  output <- plsr.loading(data,resnu)
  return(output)
}

#####
##wpli
#pain conn 

Connec =wPLI_Conn
YY=Pheno$Pain * Pheno$Anxiety
Z1=Pheno$Age
Z2=Pheno$Sex

#regressing out age and sex
data_ <- data.frame(Z1,Z2,YY)
mdl <- lm(YY~Z1+Z2 , data_)
YY <- mdl$residuals

for(i in 1:36){
  for(j in 1:36){
    for(k in 1:7){
      data_ <- data.frame(Z1,Z2,Connec[i,j,k,])
      colnames(data_)[3] <- "Conn"
      mdl <- lm(Conn~Z1+Z2 , data_)
      Conn_ <- mdl$residuals
      Connec[i,j,k,] <- Conn_
    }
  }
}

##

resnu = 1000

ncores <- 7

cl <- makeCluster(ncores)
clusterExport(cl, varlist = c("pls.modeling","plsr.loading","eff.comp.numb",
                              "createFolds","plsr","resample.data",
                              "Connec","YY", "resnu","rois" ))

PLSR_painanxiety_wPLI_noagesex <- parLapply(cl = cl,1:7 , pls.modeling )

stopCluster(cl)

rm(list = ls()[ls() != "PLSR_painanxiety_wPLI_noagesex"])
save.image("PLSR_painanxiety_wPLI_noagesex.RData")

