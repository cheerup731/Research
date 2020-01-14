
### 3 indicators, c=9
start_time <- Sys.time()
## True value (x) 
# number of categories (3, 6, 9)
C <- c(9, 6, 3)

## indicators (y)
# number of indicators (2, 3)
ni <- 3

# missing proportion (0, 0.05, 0.1)
MP <- c(0.1, 0.05, 0)
  
# error proportion (0.05, 0.15, 0.25)
EP <- c(0.25, 0.15, 0.05)

## senrio: 
# sample size (1000)
n <- 1000

# numbers of imputation/ bootstrap
M <- 10

# number of simulation 
S <- 1000
#s <- 1
#ii <- 4

MPEP <- expand.grid(MP = MP, EP = EP, C = C)

Result <- list()

# packages
library(Amelia)
library(lavaan)

###############################################################

for (ii in 1:nrow(MPEP)) {
  cat("ii=",ii,"\n\n")
   c <- MPEP[ii, 3]
  mp <- MPEP[ii, 1]
  ep <- MPEP[ii, 2]
  
  truth <- 1/c
  truth_UC<- mean(1:c)
  Bias    <- matrix(data = NA, nrow = S, ncol = c)
  Bias_MO <- matrix(data = NA, nrow = S, ncol = c)
  
  Ori_Bias_UC <- matrix(data = NA, nrow = S, ncol = 1)
  Bias_UC <- matrix(data = NA, nrow = S, ncol = 1)
  
  RMSE <- matrix(data = NA, nrow = S, ncol = 1)
  Ori_RMSE_UC <- matrix(data = NA, nrow = S, ncol = 1)
  RMSE_MO <- NA
  RMSE_UC <- NA
  
   
# begin of simulation
for (s in 1:S) {
cat("s=",s,"\n\n")
set.seed(s*ii)

#create truth
x <- sample(1:c, size = n, replace = TRUE) 
data <- as.data.frame(x)

##### create indicators y #####

for (i in 1:ni) {
  
  data[,i+1] <- x
  
  data[,i+1][sample(1:n, size = (n*(mp+ep)), replace = FALSE)] <- NA  # take out missing and error members
  
}


# fill in errors

errormember <- matrix(data = NA, nrow = ni, ncol = n*ep)

for (i in 1:ni) {
  
  errormember[i,] <- sample( which(is.na(data[ ,i+1])) , size = n*ep, replace = FALSE)
  
  for (j in 1:round(n*ep)) {
    
    data[ ,i+1][ errormember[i,j] ] <- sample( c(1:c)[ -x[ errormember[i,j] ] ], size = 1, replace = TRUE)
    
  }
}

name <- NA

for (i in 1:ni) {
  
  name[i] <- paste0("y", i)
  
}

names(data)<- c("x", name)

##### original bias #####

Bias[s,]<- table(data[ ,2]) / (n*(1-mp)) - truth
  
RMSE[s,]<- sum(Bias[s,]^2)

Ori_Bias_UC[s]<- mean(data[ ,2],na.rm = T)-truth_UC

Ori_RMSE_UC[s]<- Ori_Bias_UC[s]^2
 

##### MO #####

da <- data[,-1]

# calculate error variance by CFA
fit <- cfa(model = 'x =~ y1 + y2 + y3 ', data = da) #**
errorsd <- sqrt(fit@ParTable$est[4])

#MO
mo <- moPrep(da, y1~y1, error.sd = errorsd)
mo.out <- amelia(x = mo, m = M, p2s=0, incheck=FALSE)

# save the imputated values 
da$w1 <- mo.out[[1]]$imp1$y1  
da$w2 <- mo.out[[1]]$imp2$y1
da$w3 <- mo.out[[1]]$imp3$y1
da$w4 <- mo.out[[1]]$imp4$y1
da$w5 <- mo.out[[1]]$imp5$y1
da$w6 <- mo.out[[1]]$imp6$y1 
da$w7 <- mo.out[[1]]$imp7$y1
da$w8 <- mo.out[[1]]$imp8$y1
da$w9 <- mo.out[[1]]$imp9$y1
da$w10 <- mo.out[[1]]$imp10$y1

# average before cutoff
mean_MO <-  apply( da[,-c(1,2)], MARGIN = 2, FUN = mean, na.rm = TRUE)
Bias_UC[s] <- mean(mean_MO) - truth_UC
RMSE_UC[s] <- Bias_UC[s]^2 

theta_MO <- matrix(data = NA, nrow = M, ncol = c) # storing each imputation's parameter

for (i in 1:M) {
  
  # cut continous imputed values into groups
  da[ ,ni+i] <- cut(da[ ,ni+i], breaks = c(-Inf, 1.5:(c - 0.5), +Inf), labels = c(1:c))   
  
  # calculate parameter
  theta_MO[i, ] <- table(da[ ,(ni+i)])/n                                          

}


#pooled mean
p_mean <- apply(theta_MO, MARGIN = 2, FUN = mean)
Bias_MO[s,] <- p_mean - truth
RMSE_MO[s] <- sum( Bias_MO[s,]^2 )


} #end of simulation


Result[[ii]] <- list(
  S = S, c= c, ep = ep, mp = mp,
Ori_Bias_cut = apply( Bias, MARGIN = 2, FUN = mean ),
Bias_cut = apply( Bias_MO, MARGIN = 2, FUN = mean ),
Ori_Bias_uncut = apply( Ori_Bias_UC, MARGIN = 2, FUN = mean ),
Bias_uncut = apply( Bias_UC, MARGIN = 2, FUN = mean ),
Ori_RMSE_cut = sqrt(mean(RMSE)),
RMSE_cut = sqrt(mean(RMSE_MO)),
Ori_RMSE_uncut = sqrt(mean(Ori_RMSE_UC)),
RMSE_uncut = sqrt(mean(RMSE_UC)) )

}


capture.output(Result, file = "amelia_test_result.txt")
Result
end_time <- Sys.time()
end_time - start_time
