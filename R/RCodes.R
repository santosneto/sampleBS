#' Random generation for the Birnbaum-Saunders distribution.
#'
#'@description The random generation for the Birnbaum-Saunders distribution.
#'
#'@usage rbs(n=1.0, alpha = 0.5, beta = 1.0)
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
#'@author Eliardo Costa \email{eliardo@ccet.ufrn} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#' 
#'@examples 
#' x <- rbs(n=10, alpha = 10, beta = 2.0) 
#' x   
#'      
#' @export


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


#' Random generation for the posterior distribution of the Birnbaum-Saunders/inverse-gamma model.
#'
#'@description The function for random generation of the posterior distribution of the Birnbaum-Saunders/inverse-gamma model.
#'
#'
#'@usage rpost.bs(N, x, a1, b1, a2, b2, r)
#'
#' @param N number of observations.
#' @param x vector of observed values of the model.
#' @param a1 hyperparameter of the prior distribution for beta. 
#' @param b1 hyperparameter of the prior distribution for beta.
#' @param a2 hyperparameter of the prior distribution for \code{alpha^2}.
#' @param b2 hyperparameter of the prior distribution for alpha^2.
#' @param r a positive constant for the sampling method. 
#'
#' @note We adapted the available R script (Supplementary Material) of Wang et. al. (2016). Wang et. al. (2016) shows that
#' a good choice of \code{r} should be an integer that is greater than or equal 1 and that the Bayesian estimation of the BS distribution
#' is insensitive to the choice of \code{r}. Wang et. al. (2016) also shows that the best combinations of the hyperparameters should be taken \eqn{a1=b1=a2=b2 \neq 10^{-3}} 
#' in order to use diffuse proper priors and obtain more stable results.
#'
#' @return A random sample of the posterior distribution of the model Birnbaum-Saunders/inverse-gamma. 
#' 
#' @references 
#' 
#' Wang, M., Sun, X. and Park, C. (2016) Bayesian analysis of Birnbaum-Saunders distribution via the generalized ratio-of-uniforms methods. Comput. Stat. 31: 207--225.
#' 
#'@author Eliardo Costa \email{eliardo@ccet.ufrn} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#' @export
#' 
#' @importFrom stats median optimize quantile rnorm runif var
 

rpost.bs <- function(N, x, a1, b1, a2, b2, r) 
  {
  
  n <- length(x)
  a1 <- a1
  b1 <- b1
  a2 <- a2
  b2 <- b2
  
  betaLog <- function(b) {
   (1/(r + 1))*( -(n + a1 + 1)*log(b) - b1/b + (n/2)*log(b) + sum(x+b) - (3/2)*sum(log(x)) - ((n + 1)/2 + a2)*((n/2)*(mean(x)/b + mean(1/x)*b -2) + b2   ) )   
  }
  betafLog <- function(b) {
  log(b) + (r/(r + 1)) * (-(n + a1 + 1)*log(b) - b1/b + (n/2)*log(b) + sum(x+b) - (3/2)*sum(log(x)) - ((n + 1)/2 + a2)*((n/2)*(mean(x)/b + mean(1/x)*b -2) + b2   ) )
  }
  
  a.max <- optimize(betaLog, lower = 0, upper = 1E2, maximum = TRUE)$objective
  b.max <- optimize(betafLog, lower = 0, upper = 1E2, maximum = TRUE)$objective
  
  print(a.max)
  print(b.max)
  
  a.val <- b.val <- rep(0, N)
  for (j in 1:N) {
    U <- runif(1, 0, exp(a.max))
    V <- runif(1, 0, exp(b.max))
    rho <- V/(U^r)
    while (log(U) > betaLog(rho)) {
      U <- runif(1, 0, exp(a.max))
      V <- runif(1, 0, exp(b.max))
      rho <- V/(U^r)
    }
    b.val[j] <- rho
    a.val[j] <- rigamma(1, n/2 + a2, sum(x/rho + rho/x - 2)/2 + b2)
  }
  #cred.a <- emp.hpd(sqrt(a.val), conf = cred.level)
  #cred.b <- emp.hpd(b.val, conf = cred.level)
  alpha <- sqrt(a.val) #c(median(sqrt(a.val)), sd(sqrt(a.val)), cred.a)
  beta <-  b.val #c(median(b.val), sd(b.val), cred.b)
  output <- cbind(alpha, beta)
  #colnames(output) <- c("Median", "SD", "Lower", "Upper")
  return(output)
}

#' Bayesian sample size in a decision-theoretic approach for the Birbaum-Saunders/inverse-gamma model.
#'
#'@description Function for Bayesian sample size determination via decision-theoretic approach for the Birbaum-Saunders/inverse-gamma model.
#'
#'@usage bss.dt.bs(loss = "L1", a1 = 1E-3, b1 = 1E-3, a2 = 1E-3, 
#'                 b2 = 1E-3, c=0.005, rho = NULL, gam = NULL,
#'                 nmax = 1E2, nlag = 1E1, nrep = 1E2, lrep = 1E2,
#'                 npost = 1E2, plot = FALSE, ...)
#'
#' @param loss L1 (Absolute loss), L2 (Quadratic loss), L3 (Weighted loss) and L4 (Half loss) representing the loss function used. The default is absolute loss function.
#' @param a1 hyperparameter of the prior distribution for beta. The default is 1E-3.
#' @param b1 hyperparameter of the prior distribution for beta. The default is 1E-3.
#' @param a2 hyperparameter of the prior distribution for alpha^2. The default is 1E-3.
#' @param b2 hyperparameter of the prior distribution for alpha^2. The default is 1E-3.
#' @param c a positive real number representing the cost of colect one aliquot. The default is 0.005.
#' @param rho a number in (0, 1). The probability of the credible interval is \eqn{1-rho}. Only
#' for lost function 1.
#' @param gam a positive real number connected with the credible interval when using lost
#' function 2.
#' @param nmax a positive integer representing the maximum number for compute the Bayes risk.
#' Default is 100.
#' @param nlag a positive integer representing the lag in the n's used to compute the Bayes risk. Default is 10.
#' @param nrep a positive integer representing the number of samples taken for each \eqn{n}.
#' @param lrep a positive integer representing the number of samples taken for \eqn{S_n}. Default is 100.
#' @param npost a positive integer representing the number of values to draw from the posterior distribution of the mean. Default is 100.
#' @param plot Boolean. If TRUE (default) it plot the estimated Bayes risks and the fitted curve.
#' @param ... Currently ignored.
#'
#' @return An integer representing the sample size.
#' 
#' @references 
#'Costa, E.G., Paulino, C.D., and Singer, J. M. (2019). Sample size determination to evaluate ballast water standards: a decision-theoretic approach. Tech. rept. University of Sao Paulo. 
#'
#'
#'@author Eliardo Costa \email{eliardo@ccet.ufrn} and Manoel Santos-Neto \email{manoel.ferreira@ufcg.edu.br}
#'
#'@examples  
#'#n_bs <- bss.dt.bs(plot=TRUE)
#'#print(n_bs)
#'  
#' @export
#' @importFrom LearnBayes rigamma 


bss.dt.bs <- function(loss = 'L1', a1 = 1E-3, b1 = 1E-3, a2 = 1E-3, b2 = 1E-3, c = 0.005, rho = NULL, gam = NULL,
                      nmax = 1E2, nlag = 1E1, nrep = 1E2, lrep = 1E2, npost = 1E2, plot = FALSE, ...) {

  cl <- match.call()
  ns <- rep(seq(3, nmax, by = nlag), each = nrep)
  if (loss == 'L2') { # quadratic loss
    risk <- sapply(ns, function(n) {
      loss <- sapply(seq_len(lrep), function(j) {
        alpha2 <- rigamma(n = n, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- rigamma(n = n, a = a1, b = b1)
        x <- rbs(n = 1, alpha = alpha, beta = beta)
        post.xn <- rpost.bs(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2, r = max(1/(a1+a2+(1/2)) + 1E-3,1.0))
        mu.post <- post.xn[, 2]*(1 + post.xn[, 1]^2/2)
        out.loss <- var(mu.post) + c*n
        return(out.loss)
      })
      out.risk <- mean(loss)
      return(out.risk)
    })
  }
  if (loss == 'L1') { # absolute loss
    risk <- sapply(ns, function(n) {
      loss <- sapply(seq_len(lrep), function(j) {
        alpha2 <- rigamma(n = n, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- rigamma(n = n, a = a1, b = b1)
        print(alpha)
        print(beta)
        x <- rbs(n = 1, alpha = alpha, beta = beta)
        post.xn <- rpost.bs(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2, b2 = b2, r = max(1/(a1+a2+(1/2)) + 1E-3,1.0)) 
        mu.post <- post.xn[, 2]*(1 + post.xn[, 1]^2/2)
        med.post <- median(mu.post)
        out.loss <- mean(abs(mu.post - med.post)) + c*n
        return(out.loss)
      })
      out.risk <- mean(loss)
      return(out.risk)
    })
  }
  if (loss == 'L3') { # loss function for interval inference depending on rho
    risk <- sapply(ns, function(n) {
      loss <- sapply(seq_len(lrep), function(j) {
        alpha2 <- rigamma(n = n, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- rigamma(n = n, a = a1, b = b1)
        x <- rbs(n = 1, alpha = alpha, beta = beta)
        post.xn <- rpost.bs(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2,b2 = b2, r = max(1/(a1+a2+(1/2)) + 1E-3,1.0)) 
        mu.post <- post.xn[, 2]*(1 + post.xn[, 1]^2/2)
        qs <- quantile(mu.post, probs = c(rho/2, 1 - rho/2))
        out.loss <- sum(mu.post[which(mu.post > qs[2])])/npost - sum(mu.post[which(mu.post < qs[1])])/npost + c*n
        return(out.loss)
      })
      out.risk <- mean(loss)
      return(out.risk)
    })
  }
  if (loss == 'L4') { # loss function for interval inference depending on gamma
    risk <- sapply(ns, function(n) {
      loss <- sapply(seq_len(lrep), function(j) {
        alpha2 <- rigamma(n = n, a = a2, b = b2)
        alpha <- sqrt(alpha2)
        beta <- rigamma(n = n, a = a1, b = b1)
        x <- rbs(n = 1, alpha = alpha, beta = beta)
        post.xn <- rpost.bs(N = npost, x = x, a1 = a1, b1 = b1, a2 = a2,b2 = b2, r = max(1/(a1+a2+(1/2)) + 1E-3,1.0))
        mu.post <- post.xn[, 2]*(1 + post.xn[, 1]^2/2)
        out.loss <- 2*sqrt(gam*stats::var(mu.post)) + c*n
        return(out.loss)
      })
      out.risk <- mean(loss)
      return(out.risk)
    })
  }
  
  Y <- log(risk - c*ns)
  fit <- stats::lm(Y ~ I(log(ns + 1)))
  E <- as.numeric(exp(fit$coef[1]))
  G <- as.numeric(-fit$coef[2])
  nmin <- ceiling((E*G/c)^(1/(G + 1))-1)
  
  if (plot == TRUE) {
    plot(ns, risk, xlim = c(0, nmax), xlab = "n", ylab = "TC(n)")
    curve <- function(x) {c*x + E/(1 + x)^G}
    plot(function(x) curve(x), 0, nmax, col = "blue", add = TRUE)
    graphics::abline(v = nmin, col = "red")
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nSample size:\n")
  cat("n  = ", nmin, "\n")
}
