#' Computing multivariate Gaussian density
#'
#' `dmnorm` computes the multivariate joint density Gaussian distribution
#'
#' @param y a vector (1 x d) or matrix (n x d) of multivariate observations
#' @param mu a vector (1 x d) or matrix (n x d) of multivariate observations
#' @param Sigma variance-covariance matrix (n x d)
#' @param log1 TRUE if logaritmith transformation is used. Default is FALSE.
#'
#' @return Multivariate joint density Gaussian distribution
#' @export
#' @import dplyr
#' @import Matrix
#' @import stringr
#' @author Joaquín Martínez-Minaya <\email{jmarminaya@@gmail.com}>
dmnorm <- function(y, mu, Sigma, log1 = FALSE){
  k <- ncol(Sigma)
  y <- t(y)
  mu <- t(mu)
  dmn <- exp((-1/2)*diag(t(y-mu)%*%solve(Sigma)%*%(y-mu)))/sqrt(((2*pi)^k)*det(Sigma))
  if(log1 == TRUE){
    log(dmn)
  }else{
    dmn
  }
}



#' Extracting posterior distribution of the linear predictor a of the covariance matrix from an object obtained from inla.posterior.sample
#'
#' `extract_lp_sigma` Extracts posterior distribution of the linear predictor a of the covariance matrix from an object obtained from inla.posterior.sample
#'
#' @param xx a resulting object of applying the function `inla.posterior.sample`
#'
#' @return A list with two objects: 1. Linear predictor substracting id.z. 2. Posterior Covariance matrix
#' @export
#' @import stringr
#' @import dplyr
#' @author Joaquín Martínez-Minaya <\email{jmarminaya@@gmail.com}>
extract_lp_sigma <- function(xx){
  hyper_like <- xx$hyperpar %>% names(.) %>%
    stringr::str_starts("Precision for the Gaussian observations") %>%
    xx$hyperpar[.]
  len <- length(hyper_like)

  hyper_cov <- xx$hyperpar %>% names(.) %>%
    stringr::str_starts("Precision for id.z") %>%
    xx$hyperpar[.]

  Sigma_post <- diag(1/hyper_like) +
    matrix(1/hyper_cov, ncol = len, nrow = len)

  #Predictor
  predictor <- xx$latent %>% rownames(.) %>%
    stringr::str_starts("APredictor") %>% xx$latent[.]
  reffect <- xx$latent %>% rownames(.) %>%
    stringr::str_starts("id.z") %>% xx$latent[.] %>% rep(., len)

  pred1 <- (predictor - reffect) %>% matrix(., ncol = len, byrow = FALSE)
  list(predictor = pred1,
       Sigma = Sigma_post)
}


#' Computing dic for a Multivariate Gaussian likelihood
#'
#' `dic.mult` computes DIC for a Multivariate Gaussian likelihood
#'
#' @param inf a list with nsim elements. Each element of the list is a result of applying extract_lp_sigma.
#' @param y is a matrix in N x D with the response variable.
#' @return e.dev, dev.e, p.eff, dic
#' @export
#' @import dplyr
#' @author Joaquín Martínez-Minaya <\email{jmarminaya@@gmail.com}>
dic.mult <- function(inf, y){
  nsim <- length(inf)
  #Expected deviance
  e.dev <- parallel::mclapply(inf,
                    function(x){
                      a <- -2*dmnorm(mu = x$predictor,
                                     Sigma = x$Sigma,
                                     y = y,
                                     log1 = TRUE) %>% sum(.)}) %>%
    unlist(.) %>%
    mean(.)

  #Deviance of the expected value
  pred_mean <- inf %>%
    lapply(., `[[`, "predictor") %>%
    Reduce("+", .)/nsim

  Sigma_post <- inf %>%
    lapply(., `[[`, "Sigma") %>%
    Reduce("+", .)/nsim

  dev.e <- sum(-2*dmnorm(y = as.matrix(y),
                         mu = pred_mean,
                         Sigma = Sigma_post,
                         log1 = TRUE))

  p.eff <- e.dev - dev.e
  dic1 <- p.eff + e.dev
  data.frame(e.dev = e.dev,
             dev.e = dev.e,
             p.eff = p.eff,
             dic   = dic1)
}


#' Computing waic for a Multivariate Gaussian likelihood
#'
#' `waic.mult` computes WAIC for a Multivariate Gaussian likelihood
#'
#' @param inf a list with nsim elements. Each element of the list is a result of applying extract_lp_sigma.
#' @param y is a matrix in N x D with the response variable.
#' @return p.eff.waic, waic
#' @export
#' @import dplyr
#' @author Joaquín Martínez-Minaya <\email{jmarminaya@@gmail.com}>
waic.mult <- function(inf, y){
  eval.p <- sapply(inf, function(x){
    dmnorm(mu = x$predictor,
           Sigma =x$Sigma,
           y = y,
           log1 = TRUE)})

  p.eff.waic <- apply(eval.p, 1, var) %>%
    sum(.)


  eval.p2 <- sapply(inf, function(x){
    dmnorm(mu = x$predictor,
           Sigma =x$Sigma,
           y = y,
           log1 = FALSE)})

  lppd <- eval.p2 %>%
    apply(., 1, mean) %>%
    log(.) %>%
    sum(.)

  WAIC.1 <- -2*(lppd - p.eff.waic)

  data.frame(p.eff.waic = p.eff.waic,
             waic       = WAIC.1)
}
