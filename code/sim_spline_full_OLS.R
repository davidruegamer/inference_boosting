library(devtools)
library(ggplot2)
library(data.table)
# install_github("selective-inference/R-software/selectiveInference")
library(selectiveInference)
# install_github("davidruegamer/iboost")
library(iboost)

# fixed settings
alpha <- 0.05
n <- 300
nu <- 0.1
mstop <- 50
commonDf <- 6
commonKnots <- 14
simReps <- 500
B <- 1000
nrCores <- 25

set.seed(201089)

# covariates
x1 <- scale(rnorm(n), scale = F)
x2 <- scale(rnorm(n), scale = F)
x3 <- scale(rnorm(n), scale = F)
x4 <- scale(rnorm(n), scale = F)
x5 <- scale(rnorm(n), scale = F)
x6 <- scale(rnorm(n), scale = F)
x7 <- scale(rnorm(n), scale = F)
x8 <- scale(rnorm(n), scale = F)
x9 <- scale(rnorm(n), scale = F)
x10 <- scale(rnorm(n), scale = F)
x11 <- scale(rnorm(n), scale = F)
x12 <- scale(rnorm(n), scale = F)
x13 <- scale(rnorm(n), scale = F)
x14 <- scale(rnorm(n), scale = F)
x15 <- scale(rnorm(n), scale = F)

# baselearner
bx1 <- bbs(x1, knots = commonKnots, df = commonDf)
bx2 <- bbs(x2, knots = commonKnots, df = commonDf)
bx3 <- bbs(x3, knots = commonKnots, df = commonDf)
bx4 <- bbs(x4, knots = commonKnots, df = commonDf)
bx5 <- bbs(x5, knots = commonKnots, df = commonDf)
bx6 <- bbs(x6, knots = commonKnots, df = commonDf)
bx7 <- bbs(x7, knots = commonKnots, df = commonDf)
bx8 <- bbs(x8, knots = commonKnots, df = commonDf)
bx9 <- bbs(x9, knots = commonKnots, df = commonDf)
bx10 <- bbs(x10, knots = commonKnots, df = commonDf)
bx11 <- bbs(x11, knots = commonKnots, df = commonDf)
bx12 <- bbs(x12, knots = commonKnots, df = commonDf)
bx13 <- bbs(x13, knots = commonKnots, df = commonDf)
bx14 <- bbs(x14, knots = commonKnots, df = commonDf)
bx15 <- bbs(x15, knots = commonKnots, df = commonDf)

# list of covariates / baselearner
xList <- list(x1,x2,x3,x4,x5,
              x6,x7,x8,x9,x10,
              x11,x12,x13,x14,x15)
blListOrg <- list(bx1, bx2, bx3, bx4, bx5,
                  bx6, bx7, bx8, bx9, bx10,
                  bx11, bx12, bx13, bx14, bx15)

# true linear predictor
mu <- as.numeric(sin(2*x1) + 0.5*x2^2)
mu <- mu - mean(mu)
sdmu <- sd(mu)

settings <- expand.grid(list(SNR = c(0.5, 1)))#, 
                             #ind = c(5, 15)))

for(this_set in 1:nrow(settings)){

  sigma = sdmu / settings$SNR[this_set]
  # ind = settings$ind[this_set]
  ind = 15
  
  blList <- blListOrg[1:ind]
  
  res <- mclapply(1:simReps, function(nr){
    
    set.seed(nr)
    
    y <- as.numeric(scale(mu + rnorm(n, 0, sigma), scale=F))
    modOrg <- mboost_fit(blList, offset = 0, response = y, 
                         control = boost_control(nu = nu, mstop = mstop))
    
    # fixFolds <- cv(weights = model.weights(modOrg),
    #                type = "bootstrap", B = 5)
    # cvr <- cvrisk(modOrg, folds = fixFolds, papply = lapply)
    # (stopiter <- mstop(cvr))
    # # stopiter = 8
    # modf <- modOrg[stopiter]
    (selC <- unique(selected(modOrg)))
    if(!2%in%selC | !1%in%selC) return(NULL)
    
    
    # define corresponding refit function
    modFun <- function(y){
      
      mod <- mboost_fit(response = y,
                        blg = blList,
                        offset = 0, 
                        control = boost_control(mstop = mstop,
                                                nu = nu))
      # cvr <- cvrisk(mod, folds = fixFolds, papply = lapply)
      # stopiter <- mstop(cvr)
      # # stopiter = 8
      # mod <- mod[stopiter]
      return(mod)
      
    }

    
    
    ss <- sort(selC)
    
    vars <- c(var(resid(modOrg)),
              sigma^2,
              var(y)*(n-1)/n)
    
    s1 <- Sys.time()
    
    sols <- try(iboost(obj = modOrg, 
             method = "impsamp",
             var = vars, 
             varForSampling = sigma^2,
             alpha = alpha, 
             B = B, 
             refit.mboost = modFun, 
             ncore = 1,
             computeCI = FALSE,
             checkBL = FALSE))
    
    s2 <- Sys.time() - s1
    
    if(class(sols)=="try-error"){
      
      reti <- data.frame(lower = rep(NA, length(ss)),
                         mean = rep(NA, length(ss)),
                         upper = rep(NA, length(ss)),
                         pval = rep(NA, length(ss)),
                         lowtrunc = rep(NA, length(ss)),
                         uptrunc = rep(NA, length(ss)),
                         vartype = rep(NA, length(ss)))
      reti$covariate <- ss
      reti$selection <- paste(ss, collapse = "_")
      reti$simnr <- nr
      
    }else{
      
      reti <- lapply(sols, "[[", "resDF")
      for(i in 1:3){
        
        reti[[i]]$vartype <- c("Variance of boosting residuals", 
                               "True variance",
                               "Empirical response variance")[i]
        reti[[i]]$covariate <- ss
        reti[[i]]$selection <- paste(ss, collapse = "_")
        reti[[i]]$simnr <- nr
        reti[[i]]$nrSamples <- sapply(sols[[i]]$dist, function(x) sum(x$logvals))
        
      }
      
      reti <- do.call("rbind", reti)
  
    }

    return(reti)
    
  }, mc.cores=nrCores)
  
  res <- res[!sapply(res, is.null)]
  # res <- res[sapply(res, class)!="try-error"]
  res <- rbindlist(res)
  
  saveRDS(res, file=paste0("sim_results/splines/OLS_setting_",this_set,".RDS"))
  
}

res <- lapply(1:nrow(settings), function(set) readRDS(paste0("sim_results/splines/OLS_setting_", set, ".RDS")))

for(i in 1:nrow(settings)){

  # res[[i]] <- rbindlist(res[[i]][sapply(res[[i]], class)!="try-error"])
  res[[i]]$SNR <- settings$SNR[i]
  # res[[i]]$nrCov <- settings$ind[i]

}

res <- do.call("rbind", res)

# res$cov[!res$cov%in%1:2] <- "noise"

res$covariate[!res$covariate%in%c(1:2)] <- "noise"

table(grepl(res$selection, pattern = "1_2"))
# => use all observations

# subset <- res$selection%in%names(table(res$selection)[table(res$selection)>100])

g <- ggplot(res, aes(sample = pval, colour = factor(SNR))) +
  geom_abline(slope = 1, intercept = 0, linetype=2) +
  geom_qq(distribution = stats::qunif, size = 1, geom="line") +
  theme_bw() + facet_grid(covariate ~ vartype) + # coord_flip() +
  xlab("Expected Quantiles") + ylab("Observed p-values") # +
  
g

saveRDS(g, file="plot_objects/splines/refit_w_OLScov_15covs.RDS")
