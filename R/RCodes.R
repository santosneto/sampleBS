#'@name rbs
#'
#'@aliases rbs
#'
#'@title Random generation for the Birnbaum-Saunders distribution.
#'
#'@description A function for generating values from the Birnbaum-Saunders distribution.
#'
#'
#'@param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#'@param alpha vector of shape parameter values.
#'@param beta vector of scale parameter value.
#'
#'@details The density function of the Birnbaum-Saunders distribution used in the function \code{rbs()} is
#'
#'\deqn{f_{X}(x|\alpha, \beta) = \frac{1}{\sqrt{2\,\pi}}\,\exp\left[-\frac{1}{2\alpha^2} \left(\frac{x}{\beta}+ \frac{\beta}{x}-2\right)\right]\frac{(x+\beta)}{2\alpha \sqrt{\beta x^{3}}}}
#'
#'@return A sample of size n from the Birnbaum-Saunders distribution.
#' 
#'@note If X is Birnbaum-Saunders distributed then  
#' 
#' \eqn{X = (\beta/4)(\alpha Z + \sqrt{(\alpha Z)^2 + 4})^2,} 
#' 
#' where Z follows a standard normal distribution.
#' 
#'@author Eliardo G. Costa \email{eliardocosta@ccet.ufrn.br} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#' 
#'@examples 
#' x <- rbs(n=10, a = 10, b = 2.0) 
#' x   
#'      
#' @export
#' 
#' @importFrom stats rnorm
 
rbs <- function(n=1.0, alpha=0.5, beta=1.0) {
  if (n == 1) {
    x <- numeric()
    for (i in 1:length(alpha)) {
      z <- rnorm(1)
      x[i] <- (beta[i]/4)*(alpha[i]*z + sqrt((alpha[i]*z)^2 + 4))^2
    }
  } else if (n > 1 && length(alpha) == 1 && length(beta) == 1) {
    z <- rnorm(n)
    x <- (beta/4)*(alpha*z + sqrt((alpha*z)^2 + 4))^2
  }
  return(x)
}

#'@name logp.beta
#'@aliases logp.beta
#'@title x
#'@description x
#'
#'
#' @param beta scale parameter.
#' @param x vector of observed values of the model.
#' @param a1 hyperparameter of the prior distribution for beta. 
#' @param b1 hyperparameter of the prior distribution for beta.
#' @param a2 hyperparameter of the prior distribution for \code{shape^2}.
#' @param b2 hyperparameter of the prior distribution for \code{shape^2}.
#'
#'@return The log da marginal posterior of beta
#'
#'@export
logp.beta <- function(beta, x, a1, b1, a2, b2) {
  n <- length(x)
  if (beta > 0) {
    out <- (-(n + a1 + 1) * log(beta) - b1/beta + 
              sum(log((beta/x)^(1/2) + (beta/x)^(3/2))) - 
              ((n + 1)/2 + a2) * log(sum(1/2 * (x/beta + beta/x - 2)) + b2))
  } else {
    out <- -Inf
  }
  return(out)
}


#'@name rbeta.post
#'
#'@aliases rbeta.post
#'
#'@title Random generation for the joint posterior distribution of the Birnbaum-Saunders/inverse-gamma model.
#'
#'@description A function for generating values from the joint posterior distribution of the Birnbaum-Saunders/inverse-gamma model.
#'
#'
#' @param N number of observations.
#' @param x vector of observed values of the model.
#' @param a1 hyperparameter of the prior distribution for beta. 
#' @param b1 hyperparameter of the prior distribution for beta.
#' @param a2 hyperparameter of the prior distribution for \code{shape^2}.
#' @param b2 hyperparameter of the prior distribution for \code{shape^2}.
#' @param burnin a positive constant for the sampling method. 
#' @param thin a positive constant for the sampling method.
#' @param start a positive constant for the sampling method. 
#' @param varcov a positive constant for the sampling method. 
#' @param scale a positive constant for the sampling method.  
#'
#'
#' @return A random sample of the joint posterior distribution of the model Birnbaum-Saunders/inverse-gamma model. 
#' 
#' @references 
#' 
#' Wang, M., Sun, X. and Park, C. (2016) Bayesian analysis of Birnbaum-Saunders distribution via the generalized ratio-of-uniforms methods. Comput. Stat. 31: 207--225.
#' 
#'@author Eliardo G. Costa \email{eliardocosta@ccet.ufrn.br} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#'@export
#' 
#'@importFrom LearnBayes  rwmetrop 
 
rbeta.post <- function(N, x, a1, b1, a2, b2, burnin, thin, start,
                       varcov, scale = 1) {
  prop <- list(var = varcov, scale = scale) # parametros da dist. proposta
  sam.post <- rwmetrop(logpost = logp.beta, proposal = prop, 
                       start = start, m = burnin + N*thin, 
                       x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
  beta <- as.vector(sam.post$par[burnin + (1:N)*thin,])
  out <- list(sam = beta, accept = sam.post$accept)
  return(out)
}

#'@name bss.dt.bs
#'
#'@aliases bss.dt.bs
#'
#'@title Bayesian sample size in a decision-theoretic approach under the Birbaum-Saunders/inverse-gamma model.
#'
#'@description A function to obtain the optimal Bayesian sample size via a decision-theoretic approach for estimating the mean of the Birbaum-Saunders distribution.
#'
#' @param loss L1 (Absolute loss), L2 (Quadratic loss), L3 (Weighted loss) and L4 (Half loss) representing the loss function used. The default is absolute loss function.
#' @param a1 hyperparameter of the prior distribution for beta. The default is 3.
#' @param b1 hyperparameter of the prior distribution for beta. The default is 2.
#' @param a2 hyperparameter of the prior distribution for \code{shape^2}. The default is 3.
#' @param b2 hyperparameter of the prior distribution for \code{shape^2}. The default is 2.
#' @param cost a positive real number representing the cost of colect one observation. The default is 0.010.
#' @param rho a number in (0, 1). The probability of the credible interval is \eqn{1-rho}. Only
#' for loss function L3. The default is 0.95. 
#' @param gam a positive real number connected with the credible interval when using loss
#' function L4. The default is 0.5.
#' @param nmax a positive integer representing the maximum number for compute the Bayes risk.
#' Default is 100.
#' @param nlag a positive integer representing the lag in the n's used to compute the Bayes risk. Default is 10.
#' @param nrep a positive integer representing the number of samples taken for each \eqn{n}.
#' @param lrep a positive integer representing the number of samples taken for \eqn{S_n}. Default is 100.
#' @param npost a positive integer representing the number of values to draw from the posterior distribution of the mean. Default is 100.
#' @param plots Boolean. If TRUE (default) it plot the estimated Bayes risks and the fitted curve.
#' @param prints Boolean. If FALSE (default) the output is a list. 
#' @param nburn a positive constant for the sampling method.
#' @param thin a positive constant for the sampling method.
#' @param scale a positive constant for the sampling method. 
#' @param diag.name x
#' @param path.diag x
#' @param save.plot Boolean. If TRUE, the plot is saved to an external file. The default is FALSE.    
#' @param ... Currently ignored.
#' 
#'
#' @return An integer representing the optimal sample size.
#' 
#' @references 
#'Costa, E.G., Paulino, C.D., and Singer, J. M. (2019). Sample size determination to evaluate ballast water standards: a decision-theoretic approach. Tech. rept. University of Sao Paulo. 
#'
#'@author Eliardo G. Costa \email{eliardocosta@ccet.ufrn.br} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#'@examples  
#' #bss.dt.bs(loss="L1", plot=TRUE, lrep=10, npost=10)
#'
#' @export
#' @importFrom LearnBayes rigamma laplace
#' @import ggplot2
#' @importFrom graphics par plot
#' @importFrom stats lm acf plot.ts median quantile rnorm var
#' @importFrom grDevices cairo_pdf pdf dev.off
#' @importFrom pbmcapply pbmcmapply
#' @importFrom parallel detectCores
#' @import dplyr 
#' @import magrittr
#' @importFrom graphics legend points

bss.dt.bs <- function(loss = 'L1', a1 = 8, b1 = 50, a2 = 8, b2 = 50, 
                      cost = 0.01, rho = 0.05, gam = 1, nmax = 2E3, 
                      nlag = 2E2, nrep = 6L, lrep = 1E2, npost = 5E2, 
                      nburn = 5E2, thin = 20L, scale = 1L,
                      plots = TRUE, prints = TRUE, save.plot = FALSE,
                      path.diag = getwd(), diag.name='plot', ...) 
{
  cl <- match.call()
  ns <- rep(seq(2, nmax, by = nlag), each = nrep)
  nprint <- ns[nrep+1]
  a <- b <- NULL
  
  if (loss == 'L1') { # absolute loss
    
    loops <- pbmcmapply(function(k){ 
      lapply(X=1:lrep,function(i,n){  
        alpha2 <- LearnBayes::rigamma(n = 1, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- LearnBayes::rigamma(n = 1, a = a1, b = b1)
        x <- rbs(n = n, alpha = alpha, beta = beta)
        lapla <- laplace(logpost = logp.beta, mode = median(x), x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2)
        beta.pos <- rbeta.post(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2,burnin = nburn, thin = thin, start = median(x), varcov = lapla$var, scale = scale)
        
        if (n == nprint & i == 1) {
          graph_name <- paste(path.diag, "/diag_", diag.name, ".pdf", sep = "")
          pdf(graph_name)
          par(mfrow = c(2, 1))
          plot.ts(beta.pos$sam, xlab = "iteration", ylab = "")
          acf(beta.pos$sam, main = "")
          dev.off()
          par(mfrow = c(1, 1))
        }
        
        t_k <- lapply(seq_len(length(beta.pos$sam)), function(i) sum(0.5*(x/beta.pos$sam[i] + beta.pos$sam[i]/x - 2))) %>% unlist() #ok!
        lam.pos <- mapply(LearnBayes::rigamma,n = 1, a = (n + 1)/2 + a2, b = b2 + t_k) 
        theta.pos <- beta.pos$sam*(1 + lam.pos/2) 
        
        medi.pos <- median(theta.pos)
        loss <-  mean(abs(theta.pos - medi.pos)) + cost*n
        accept <- beta.pos$accept
        
        
        c(loss,accept)
      },n=k) %>% unlist() %>% matrix(ncol=2,nrow=lrep,byrow = TRUE) %>% apply(MARGIN = 2,mean)
    }, k=ns,mc.cores = detectCores()-2) %>% unlist() %>% matrix(ncol=2,nrow=length(ns),byrow = TRUE)  
    
  }else if (loss == 'L2') { # quadratic loss ATUALIZANDO
    loops <- pbmcmapply(function(k){ 
      lapply(X=1:lrep,function(i,n){
        alpha2 <- LearnBayes::rigamma(n = 1, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- LearnBayes::rigamma(n = 1, a = a1, b = b1)
        x <- rbs(n = n, alpha = alpha, beta = beta)
        lapla <- laplace(logpost = logp.beta, mode = median(x), x = x, a1 = a1, b1 = b1,
                         a2 = a2, b2 = b2)
        beta.pos <- rbeta.post(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                               burnin = nburn, thin = thin, start = median(x), 
                               varcov = lapla$var, scale = scale)
        if (n == nprint & i == 1) {
          graph_name <- paste(path.diag, "/diag_", diag.name, ".pdf", sep = "")
          pdf(graph_name)
          par(mfrow = c(2, 1))
          plot.ts(beta.pos$sam, xlab = "iteration", ylab = "")
          acf(beta.pos$sam, main = "")
          dev.off()
          par(mfrow = c(1, 1))
        }
        
        t_k <- lapply(seq_len(length(beta.pos$sam)), function(i) sum(0.5*(x/beta.pos$sam[i] + beta.pos$sam[i]/x - 2))) %>% unlist() #ok!
        lam.pos <- mapply(LearnBayes::rigamma,n = 1, a = (n + 1)/2 + a2, b = b2 + t_k) #ok!
        theta.pos <- beta.pos$sam*(1 + lam.pos/2) 
        loss <- var(theta.pos) + cost*n
        accept <-  beta.pos$accept
        c(loss,accept)
      },n=k) %>% unlist() %>% matrix(ncol=2,nrow=lrep,byrow = TRUE) %>% apply(MARGIN = 2,mean)
    }, k=ns,mc.cores = detectCores()-2) %>% unlist() %>% matrix(ncol=2,nrow=length(ns),byrow = TRUE)  
  } else if (loss == 'L3') { # loss function for interval inference depending on rho
    loops <- pbmcmapply(function(k){ 
      lapply(X=1:lrep,function(i,n){
        alpha2 <- LearnBayes::rigamma(n = 1, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- LearnBayes::rigamma(n = 1, a = a1, b = b1)
        x <- rbs(n = n, alpha = alpha, beta = beta)
        lapla <- laplace(logpost = logp.beta, mode = median(x), x = x, a1 = a1, b1 = b1,
                         a2 = a2, b2 = b2)
        beta.pos <- rbeta.post(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                               burnin = nburn, thin = thin, start = median(x), 
                               varcov = lapla$var, scale = scale)
        if (n == nprint & i == 1) {
          graph_name <- paste(path.diag, "/diag_", diag.name, ".pdf", sep = "")
          pdf(graph_name)
          par(mfrow = c(2, 1))
          plot.ts(beta.pos$sam, xlab = "iteration", ylab = "")
          acf(beta.pos$sam, main = "")
          dev.off()
          par(mfrow = c(1, 1))
        }
        
        t_k <- lapply(seq_len(length(beta.pos$sam)), function(i) sum(0.5*(x/beta.pos$sam[i] + beta.pos$sam[i]/x - 2))) %>% unlist() #ok!
        lam.pos <- mapply(LearnBayes::rigamma,n = 1, a = (n + 1)/2 + a2, b = b2 + t_k) #ok!
        theta.pos <- beta.pos$sam*(1 + lam.pos/2) 
        qs <- quantile(theta.pos, probs = c(rho/2, 1 - rho/2))
        loss <-  sum(theta.pos[which(theta.pos > qs[2])])/npost - sum(theta.pos[which(theta.pos < qs[1])])/npost + cost*n
        accept <-  beta.pos$accept
        c(loss,accept)
      },n=k) %>% unlist() %>% matrix(ncol=2,nrow=lrep,byrow = TRUE) %>% apply(MARGIN = 2,mean)
    }, k=ns,mc.cores = detectCores()-2) %>% unlist() %>% matrix(ncol=2,nrow=length(ns),byrow = TRUE)  
  } 
  else if(loss == 'L4'){ # loss function for interval inference depending on gamma
    loops <- pbmcmapply(function(k){ 
      lapply(X=1:lrep,function(i,n){
        alpha2 <- LearnBayes::rigamma(n = 1, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- LearnBayes::rigamma(n = 1, a = a1, b = b1)
        x <- rbs(n = n, alpha = alpha, beta = beta)
        lapla <- laplace(logpost = logp.beta, mode = median(x), x = x, a1 = a1, b1 = b1,
                         a2 = a2, b2 = b2)
        beta.pos <- rbeta.post(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2,
                               burnin = nburn, thin = thin, start = median(x), 
                               varcov = lapla$var, scale = scale)
        if (n == nprint & i == 1) {
          graph_name <- paste(path.diag, "/diag_", diag.name, ".pdf", sep = "")
          pdf(graph_name)
          par(mfrow = c(2, 1))
          plot.ts(beta.pos$sam, xlab = "iteration", ylab = "")
          acf(beta.pos$sam, main = "")
          dev.off()
          par(mfrow = c(1, 1))
        }
        
        t_k <- lapply(seq_len(length(beta.pos$sam)), function(i) sum(0.5*(x/beta.pos$sam[i] + beta.pos$sam[i]/x - 2))) %>% unlist() #ok!
        lam.pos <- mapply(LearnBayes::rigamma,n = 1, a = (n + 1)/2 + a2, b = b2 + t_k) #ok!
        theta.pos <- beta.pos$sam*(1 + lam.pos/2)
        
        loss <-  2*sqrt(gam*var(theta.pos)) + cost*n
        accept <- beta.pos$accept
        c(loss,accept)
      },n=k) %>% unlist() %>% matrix(ncol=2,nrow=lrep,byrow = TRUE) %>% apply(MARGIN = 2,mean)
    }, k=ns,mc.cores = detectCores()-2) %>% unlist() %>% matrix(ncol=2,nrow=length(ns),byrow = TRUE)
  }
  
  accept <- mean(loops[,2])
  risk <- loops[,1]
  Z <- log(risk - cost*ns)
  fit <- lm(Z ~ I(log(ns + 1)))
  E <- as.numeric(exp(fit$coef[1]))
  G <- as.numeric(-fit$coef[2])
  nmin <- ceiling((E*G/cost)^(1/(G + 1)) - 1)
  
  if (plots == TRUE) {
    
    if(save.plot == FALSE){ 
      par(mar=c(4.1,4.1,0.2,0.2))
      plot(ns, risk, xlim = c(0, max(ns) + 1), ylim = c(min(risk) - 0.5, max(risk) + 0.5), xlab = "n", ylab = "TC(n)",pch=19)
      curve <- function(x) {cost*x + E/(1 + x)^G}
      plot(function(x)curve(x), 0, max(ns) + 1, col = "blue", add = TRUE)
      points(nmin,curve(nmin),pch=19,col=2)
      legend("top",pch=19,legend = paste("Optimal Sample Size = ",nmin,sep = ''), col=2, bty ='n')
    }else{
      
      if(loss == 'L1'|| loss == 'L2'){
        file.name <- paste('case','_',loss,'_',a,'_',b,'_',cost,'.pdf', sep='')
      } else if(loss == 'L3'){
        file.name <- paste('case','_',loss,'_',a,'_',b,'_',cost,'.pdf', sep='')
      } else{
        file.name <- paste('case','_',loss,'_',a,'_',b,'_',cost,'.pdf', sep='')
      }
      pdf(file.name)
      par(mar=c(4.1,4.1,0.2,0.2))
      plot(ns, risk, xlim = c(0, max(ns) + 1), ylim = c(min(risk) - 0.5, max(risk) + 0.5), xlab = "n", ylab = "TC(n)",pch=19)
      curve <- function(x) {cost*x + E/(1 + x)^G}
      plot(function(x)curve(x), 0, max(ns) + 1, col = "blue", add = TRUE)
      points(nmin,curve(nmin),pch=19,col=2)
      legend("top",pch=19,legend = paste("Optimal Sample Size = ",nmin,sep = ''), col=2, bty ='n')
      dev.off()
      
    }
  }
  
  if(prints == TRUE)
  {  
    # Output
    cat("\nCall:\n")
    print(cl)
    cat("\nSample size:\n")
    cat("n  = ", nmin, "\n")
    cat("---------------\n")
    cat("accept = ", accept , "\n")
  }else{ 
    out <- list(n = nmin, risk=risk, cost=cost, loss = loss, E=E, G=G,accept=accept )
    
    return(out)
  }
}

