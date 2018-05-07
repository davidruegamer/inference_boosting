library(mboost)
library(devtools)
# install_github("davidruegamer/iboost")
library(iboost)

# define settings for boosting and
# inference as well as simulations

# step-length and nr of iterations for boosting (fixed)
nu = 0.1
mstop = 40
# alpha level (fixed)
alpha = 0.05

# nr of true effects
pAst <- 4
# nr of null effects
# nr of observations
# signal-to-noise ratio
settings <- expand.grid(p0 = c(22, 4),
                        n = c(25, 100),
                        SNR = c(1, 4)
)
# p0 = 8 settings for CV-based simulation
# so save some time by excluding some combinations
settings <- settings[-c(4,8),]
settings <- rbind(settings, data.frame(p0=4, n=25, SNR=1))

# simulation iterations
simReps <- 1000
# nr of samples used in sampling approach
B <- 1000

# function to generate list of baselearners
# depending on number observations and effects
makeBLlist <- function(this_n, this_p, seed = 201089)
{
  
  set.seed(seed)
  
  xList <- lapply(1:this_p, function(i){
    scale(rnorm(this_n), scale = F)
  })
  blList <- lapply(xList, function(x){
    bols(x, intercept=FALSE)
  })
  return(#list(xList = xList,
    blList = blList)#)
  
}


for(this_set in 1:nrow(settings)){
  
  set.seed(42)
  
  n <- settings$n[this_set]
  p0 <- settings$p0[this_set]
  p <- p0 + pAst
  
  # simulate covariates
  x1 <- scale(rnorm(n), scale = F)
  x2 <- scale(rnorm(n), scale = F)
  x3 <- scale(rnorm(n), scale = F)
  x4 <- scale(rnorm(n), scale = F)
  
  bx1 <- bols(x1, intercept=FALSE)
  bx2 <- bols(x2, intercept=FALSE)
  bx3 <- bols(x3, intercept=FALSE)
  bx4 <- bols(x4, intercept=FALSE)
  
  ### simulate true expectation
  mu = as.numeric(4 * x1 - 3 * x2 + 2 * x3 - 1 * x4)
  mu <- mu - mean(mu)
  
  sdmu <- sd(mu)
  
  this_var <- (sdmu / settings$SNR[this_set])^2

  blList <- makeBLlist(this_n = n, this_p = p0)  
  blList <- append(list(bx1, bx2, bx3, bx4), blList)
  names(blList) <- paste0("bx", 1:p)  
  
  if(settings$SNR[this_set]==4) mstop <- 150
  
  res <- mclapply(1:simReps, function(nr){
    
    set.seed(nr)
    
    y <- mu + rnorm(n, 0, sqrt(this_var))
    modOrg <- mboost_fit(blList, y, control = boost_control(nu = nu, mstop = mstop), offset = 0)
    # X <- do.call("cbind", lapply(modOrg$baselearner, extract, "design"))

    # do cross-validation for selected settings
    if(this_set==7){ 
    
      fixFolds <- cv(weights = model.weights(modOrg),
                     type = "kfold", B = 10)
      cvr <- cvrisk(modOrg, folds = fixFolds, papply = lapply)
      modf <- modOrg[mstop(cvr)]
      
      # define corresponding refit function
      modFun <- function(y){
        
        mod <- mboost_fit(response = y,
                          blg = blList,
                          offset = 0, 
                          control = boost_control(mstop = mstop,
                                                  nu = nu))
        cvr <- cvrisk(mod, folds = fixFolds, papply = lapply)
        return(mod[mstop(cvr)])
        
      }
      
    }else{
      
      modf <- modOrg
      modFun <- function(y) mboost_fit(blList, y, control = boost_control(nu = nu, mstop = mstop), 
                                       offset = 0)
      
    }
    
    # save selection
    (selC <- unique(selected(modf)))
    
    vars <- c(var(resid(modf)),
              var(resid(lm(y ~ -1+getDesignmat(modf)))), 
              this_var,
              var(y)*(n-1)/n)
              
    

    
    # do the polyhedron for selected settings
    if(this_set==4){ 
      
      nrMethods <- 2 
      Ups <- getUpsilons(modf)
      
    }else{
      
      nrMethods <- 1
      Ups <- NULL
      
    }
    
    
    sols <- lapply(1:nrMethods, function(mnr) 
      iboost(obj = modf, 
             method = c("unifsamp", "analytic")[mnr],
             var = vars, 
             varForSampling = this_var,
             alpha = alpha, 
             B = B, 
             refit.mboost = modFun, 
             ncore = 1,
             Ups = Ups, 
             checkBL = FALSE))
    
    len <- length(selC)
    sols <- do.call("rbind", lapply(sols, function(ss) 
      do.call("rbind", lapply(ss, "[[", "resDF"))))
    sols$vartype <- rep(rep(c("residvar", 
                              "lmresidvar",
                              "truevar",
                              "empyvar"), each = len), nrMethods)
    sols$methods <- rep(c("unifsamp", "analytic")[1:nrMethods], each = len*4)
    sols$covariates <- rep(sort(selC), nrMethods*4)
    sols$sel <- paste(sort(selC), collapse="_")
    sols$simnr <- nr
    sols$mstop <- mstop(modf)
    
    return(sols)
    
  }, mc.cores=25)
  
  saveRDS(res, file=paste0("sim_results/linear/simPCnew_",this_set,".RDS"))
  
}

