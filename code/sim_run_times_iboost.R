library(mboost)
library(devtools)
# install_github("davidruegamer/iboost")
library(iboost)

if(!dir.exists("sim_results/linear_run_time/"))
  dir.create("sim_results/linear_run_time/")

# define settings for boosting and
# inference as well as simulations
pAst = 4
# step-length and nr of iterations for boosting (fixed)
nu = 0.1
# alpha level (fixed)
alpha = 0.05

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

# nr of true effects
for(p0 in rev(c(5, 50, 108))){
  # nr of null effects
  # nr of observations
  for(n in c(100, 1000, 10000)){
    # signal-to-noise ratio
    
    # define number of samples and mstop adaptively
    
    B = 500
    mstop = pmin(p0*100, 10000)
    
    set.seed(42)
    
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
    
    this_var <- (sdmu / 1)^2 # SNR = 1
    
    blList <- makeBLlist(this_n = n, this_p = p0)  
    blList <- append(list(bx1, bx2, bx3, bx4), blList)
    names(blList) <- paste0("bx", 1:p)  
    
    res <- mclapply(1:5, function(nr){
      
      set.seed(nr)
      
      y <- mu + rnorm(n, 0, sqrt(this_var))
      modOrg <- mboost_fit(blList, y, control = boost_control(nu = nu, mstop = mstop), offset = 0)
      # X <- do.call("cbind", lapply(modOrg$baselearner, extract, "design"))
      
      modf <- modOrg
      
      fixFolds <- cv(weights = model.weights(modOrg),
                     type = "kfold", B = 5)
      cvr <- cvrisk(modOrg, folds = fixFolds, papply = mclapply)
      modf <- modOrg[mstop(cvr)]
      
      # define corresponding refit function
      modFun <- function(y){
        
        mod <- mboost_fit(response = y,
                          blg = blList,
                          offset = 0, 
                          control = boost_control(mstop = mstop,
                                                  nu = nu))
        cvr <- cvrisk(mod, folds = fixFolds, papply = mclapply)
        return(mod[mstop(cvr)])
        
      }
      
      
      # save selection
      (selC <- unique(selected(modf)))
      
      nrMethods <- 1
      Ups <- NULL
      
      st <- system.time(
        sols <- iboost(obj = modf, 
                       method = "impsamp",
                       var = this_var, 
                       varForSampling = this_var,
                       alpha = alpha, 
                       B = B, 
                       refit.mboost = modFun, 
                       which = 1,
                       ncore = 1,
                       Ups = Ups, 
                       checkBL = FALSE)
      )
              
      res <- data.frame(ss = sum(sols$dist[[1]]$logvals),
                        time_whole = st["elapsed"],
                        time = sols$dur,
                        simnr = nr,
                        mstop = mstop(modf),
                        sel = paste(sort(selC), collapse="_"))

      return(res)
      
    }, mc.cores = 5)
    
    saveRDS(res, file=paste0("sim_results/linear_run_time/sim_p_",p,"_n_",n,".RDS"))
    
  }
  
}