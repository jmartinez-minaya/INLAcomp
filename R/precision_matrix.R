#' Computing the structured matrix R, such as Q = tau*R
#'
#' `precision_matrix` Main function to compute R
#'
#' @param term one of besag, rw1, rw2 or iid
#' @param n number of points
#' @param H inla graph, just required if term = "besag"
#'
#' @return Structured matrix R, such as Q = tau*R
#'
#' @examples
#' precision_matrix(term = "rw1", n = 10)
#' precision_matrix(term = "rw2", n = 15)
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
precision_matrix <- function(term = "besag", n, H = NULL){
  if(term == "besag")
  {
    if(is.null(H)){
      stop("Be careful! a graph is required! H \n")
    }else{
      #Besag term
      S <- H$n
      Q.xi <- matrix(0, H$n, H$n)
      for(i in 1:H$n)
      {
        Q.xi[i,i] <- H$nnbs[[i]]
        Q.xi[i, H$nbs[[i]]] <- -1
      }
      
      Q <- Q.xi
      
    }
  }else if(term == "rw1"){
    D1 <- diff(diag(n), differences = 1)
    Q <- t(D1)%*%D1
  }else if(term == "rw2"){
    D2 <- diff(diag(n), differences = 2)
    Q <- t(D2)%*%D2
  }else{
    Q  <- diag(n)
  }
  Q
}  


#' Computing structured matrix for space-time interactions in areal data
#' such as Q = tau*R
#'
#' `precision_matrix` Main function to compute R
#'
#' @param term a vector of two elements. In each one should include one of besag, rw1, rw2 or iid
#' @param n a vector with two elements with the number of points of each effect
#' @param H a list of two elements indicating inla graph, just required if term = "besag"
#'
#' @return list with the structured matrix of the term 1, term 2 and
#' the interaction between effect term1 and term2 
#' 
#' @examples
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
precision_int <- function(term = c("besag", "rw2"), n, H){
  Q1 <- precision_matrix(term = term[1],
                         n    = n[1],
                         H    = H[[1]])
  Q2 <- precision_matrix(term = term[2],
                         n    = n[2],
                         H    = H[[2]])
  a <- list(Q1, Q2, kronecker(Q1, Q2))
  names(a) <- c(term[1], term[2], paste0(term[1], "x", term[2]))
  a
}


