#' Interface for the BayesMultMeta class
#'
#' The BayesMultMeta package implements two methods of constructing Markov
#' chains to assess the posterior distribution of the model parameters, namely
#' the overall mean vector \eqn{\mathbf{\mu}} and the between-study covariance matrix
#' \eqn{\mathbf{\Psi}}, of the generalized marginal multivariate random effects models.
#' The Bayesian inference procedures are performed when the model parameters are
#' endowed with the Berger and Bernardo reference prior
#' \insertCite{berger1992development}{BayesMultMeta} and the Jeffreys prior
#' \insertCite{1946RSPSA.186..453J}{BayesMultMeta}. This is achieved by
#' constructing Markov chains using the Metropolis-Hastings algorithms developed
#' in \insertCite{bodnar2021objective}{BayesMultMeta}. The convergence
#' properties of the generated Markov chains are investigated by the rank plots
#' and the split-\eqn{\hat{R}} estimate based on the rank normalization, which are
#' proposed in \insertCite{vehtari2021rank}{BayesMultMeta}.
#'
#' @param X A \eqn{p \times n} matrix which contains \eqn{n} observation vectors
#' of dimension \eqn{p}
#' @param U A \eqn{pn \times pn} block-diagonal matrix which contains the
#' covariance matrices of observation vectors.
#' @param N Length of the generated Markov chain.
#' @param burn_in Number of burn-in samples
#' @param likelihood Likelihood to use. It currently supports "normal" and
#' "t".
#' @param prior Prior to use. It currently supports "reference" and
#' "jeffrey".
#' @param algorithm_version One of "mu" or "Psi". Both algorithms samples the
#' same quantities.
#' @param d Degrees of freedom for the t-distribution when the "t" option is
#' used for the likelihood.
#'
#' @return a BayesMultMeta class which contains simulations from the MCMC
#' inference procedure as well as many of the input parameters. The elements
#' 'psi' and 'mu' in the list contains simulations from the posterior
#' distribution. All other elements are input parameters to the class.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' dataREM<-mvmeta::hyp
#' # Observation matrix X
#' X<-t(cbind(dataREM$sbp,dataREM$dbp))
#' p<-nrow(X)  # model dimension
#' n<-ncol(X)  # sample size
#' # Matrix U
#' U<-matrix(0,n*p,n*p)
#' for (i_n in 1:n) {
#'   Use<-diag(c(dataREM$sbp_se[i_n],dataREM$dbp_se[i_n]))
#'   Corr_mat<-matrix(c(1,dataREM$rho[i_n],dataREM$rho[i_n],1),p,p)
#'   U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)]<- Use%*%Corr_mat%*%Use
#' }
#'
#' bmgmr_run <- BayesMultMeta(X, U, 1e2, burn_in = 100,
#'                    likelihood = "normal", prior="jeffrey",
#'                    algorithm_version = "mu")
#' summary(bmgmr_run)
#' plot(bmgmr_run)
#'
#' @export
BayesMultMeta <- function(X, U, N, burn_in, likelihood, prior, algorithm_version, d=NULL) {
  assertthat::assert_that(likelihood %in% c("normal", "t"))
  assertthat::assert_that(prior %in% c("jeffrey", "reference"))
  assertthat::assert_that(algorithm_version %in% c("mu", "Psi"))

  if (algorithm_version == "mu") {
    if (likelihood == "normal" & prior == "jeffrey") {
      simulations <- sample_post_nor_jef_marg_mu(X, U, N + burn_in)
    }else if (likelihood == "normal" & prior == "reference") {
      simulations <- sample_post_nor_ref_marg_mu(X, U, N + burn_in)
    }else if (likelihood == "t" & prior == "jeffrey") {
      simulations <- sample_post_t_jef_marg_mu(X, U*(d-2)/d, d, N + burn_in)
    }else if (likelihood == "t" & prior == "reference") {
      simulations <- sample_post_t_ref_marg_mu(X, U*(d-2)/d, d, N + burn_in)
    }
  }else{
    if (likelihood == "normal" & prior == "jeffrey") {
      simulations <- sample_post_nor_jef_marg_Psi(X, U, N + burn_in)
    }else if (likelihood == "normal" & prior == "reference") {
      simulations <- sample_post_nor_ref_marg_Psi(X, U, N + burn_in)
    }else if (likelihood == "t" & prior == "jeffrey") {
      simulations <- sample_post_t_jef_marg_Psi(X, U*(d-2)/d, d, N + burn_in)
    }else if (likelihood == "t" & prior == "reference") {
      simulations <- sample_post_t_ref_marg_Psi(X, U*(d-2)/d, d, N + burn_in)
    }
  }

  structure(list(mu=simulations[[1]],
       psi=simulations[[2]],
       X=X,
       U=U,
       N=N,
       p=nrow(X),
       burn_in=burn_in,
       likelihood=likelihood,
       prior=prior,
       algorithm_version=algorithm_version,
       d=d
       ), class="BayesMultMeta")
}

#' Summary statistics from the posterior of a BayesMultMeta class
#'
#' @param object BayesMultMeta class
#' @param alpha Significance level used in the computation of the credible interval.
#' @param ... not used
#'
#' @returns a list with summary statistics
#' @export
summary.BayesMultMeta <- function(object, alpha=0.95, ...) {
  Gp<-duplication_matrix(object$p)
  Lp<-Gp%*%solve(t(Gp)%*%Gp)

  list("mu"=bayes_inference(object$mu[,(object$burn_in+1):(object$burn_in+object$N)], alpha),
       "psi"=bayes_inference(object$psi[,(object$burn_in+1):(object$burn_in+object$N)], alpha)%*%Lp)
}

#' Summary statistics from a posterior distribution
#'
#' Given a univariate sample drawn from the posterior distribution, this
#' function computes the posterior mean, the posterior median, the posterior
#' standard deviation, and the limits of the \eqn{(1-\alpha)} probability-symmetric
#' credible interval.
#'
#' @param x Univariate sample from the posterior distribution of a parameter.
#' @param alp Significance level used in the computation of the credible interval
#'
#' @return a matrix with summary statistics
bayes_inference <- function(x,alp){
  x_mean<-apply(x,1,mean)
  x_med<-apply(x,1,quantile,probs=0.5)
  x_ql<-apply(x,1,quantile,probs=alp/2)
  x_qu<-apply(x,1,quantile,probs=1-alp/2)
  x_sd<-sqrt(apply(x,1,var))
  rbind(x_mean,x_med,x_sd,x_ql,x_qu)
}


#' Plot a BayesMultMeta object
#'
#' This function produces the trace plots of the constructed Markov chains.
#'
#' @param x a BayesMultMeta object
#' @param ... additional arguments
#'
#' @return No return value, produces trace plots
#'
#' @export
plot.BayesMultMeta <- function(x, ...) {
  for (var in c("mu", "psi")) {
    df <- x[[var]]
    ylab <- function(p) ifelse(var == "mu", parse(text=paste0("mu[", p, "]")),
                               parse(text=paste0("psi[", p, "]")))
    for (idx in 1:nrow(df)) {
      plot(y=df[idx,], x=1:ncol(df), type="l", ylab=ylab(idx),
           xlab="index")
      abline(v=x$burn_in, col = "lightgray", lty = 3)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}

#' Computes the ranks within the pooled draws of Markov chains
#'
#' The function computes the ranks within the pooled draws of Markov
#' chains. Average ranks are used for ties.
#'
#' @param MC An \eqn{N \times M} matrix with N draws in each of M constructed
#' Markov chains.
#'
#' @return a matrix with the ranks from the MCMC procedure
#' @export
#'
#' @examples
#' dataREM<-mvmeta::hyp
#' # Observation matrix X
#' X<-t(cbind(dataREM$sbp,dataREM$dbp))
#' p<-nrow(X) # model dimension
#' n<-ncol(X) # sample size
#' # Matrix U
#' U<-matrix(0,n*p,n*p)
#' for (i_n in 1:n) {
#'   Use<-diag(c(dataREM$sbp_se[i_n],dataREM$dbp_se[i_n]))
#'   Corr_mat<-matrix(c(1,dataREM$rho[i_n],dataREM$rho[i_n],1),p,p)
#'   U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)]<- Use%*%Corr_mat%*%Use
#' }
#' # Generating M Markov chains for mu_1
#' M<-4 # number of chains
#' MC <-NULL
#' for (i in 1:M) {
#' chain <-  BayesMultMeta(X, U, 1e2, burn_in = 1e2,
#'                           likelihood = "t", prior="jeffrey",
#'                           algorithm_version = "mu",d=3)
#'   MC<- cbind(MC,chain$mu[1,])
#' }
#' ranks<-MC_ranks(MC)
#' id_chain <- 1
#' hist(ranks[,id_chain],breaks=25,prob=TRUE, labels = FALSE, border = "dark blue",
#'   col = "light blue", main = expression("Chain 1,"~mu[1]), xlab = expression(),
#'   ylab = expression(),cex.axis=1.2,cex.main=1.7,font=2)
#'
MC_ranks <- function(MC) {
  Np<-nrow(MC)
  M<-ncol(MC)
  matrix(rank(c(MC),ties.method = "average"),Np,M)
}

#' Computes the split-\eqn{\hat{R}} estimate based on the rank normalization
#'
#' The function computes the split-\eqn{\hat{R}} estimate based on the rank
#' normalization.
#'
#' @inheritParams MC_ranks
#'
#' @return a value with the the split-\eqn{\hat{R}} estimate based on the rank
#' normalization
#'
#' @examples
#' dataREM<-mvmeta::hyp
#' # Observation matrix X
#' X<-t(cbind(dataREM$sbp,dataREM$dbp))
#' p<-nrow(X) # model dimension
#' n<-ncol(X) # sample size
#' # Matrix U
#' U<-matrix(0,n*p,n*p)
#' for (i_n in 1:n) {
#'   Use<-diag(c(dataREM$sbp_se[i_n],dataREM$dbp_se[i_n]))
#'   Corr_mat<-matrix(c(1,dataREM$rho[i_n],dataREM$rho[i_n],1),p,p)
#'   U[(p*(i_n-1)+1):(p*i_n),(p*(i_n-1)+1):(p*i_n)]<- Use%*%Corr_mat%*%Use
#' }
#' # Generating M Markov chains for mu_1
#' M<-4 # number of chains
#' MC <-NULL
#' for (i in 1:M) {
#'   chain <-  BayesMultMeta(X, U, 1e2, burn_in = 1e2,
#'                           likelihood = "t", prior="jeffrey",
#'                           algorithm_version = "mu",d=3)
#'   MC<- cbind(MC,chain$mu[1,])
#' }
#' split_rank_hatR(MC)
#'
split_rank_hatR<-function(MC) {
  Np<-nrow(MC)
  M<-ncol(MC)

  if (Np%%2) {
    simpleError('The length of Markov chain should be an even number')
  } else {
    ranks<-MC_ranks(MC)
    x_ranks_full<-qnorm((ranks-3/8)/(M*Np+1/4))
    x_ranks<-cbind(x_ranks_full[1:(Np/2),],x_ranks_full[(Np/2+1):Np,])
    means<-apply(x_ranks,2,mean)
    BR<-Np*var(means)/2
    vars<-apply(x_ranks,2,var)
    WR<-sum(vars)/(2*M)
    sqrt(1-2/Np+(2/Np)*BR/WR)
  }
}
