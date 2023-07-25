#' Selecting models based on DIC, WAIC, LCPO
#'
#' `Bdiclcpomodel_stack` is a function for selecting models based on DIC, WAIC, LCPO
#'
#'
#'
#' @param variables: names of the terms to be included in the model selection process
#' @param spatial
#' @param spatial.short
#' @param friends parameter for computing the CPO when we assume that some data points are friends. 
#' @param n: number of models that we want to fit
#' @param ...: more arguments for the inla function
#'
#' @return A list of three elements:
#' 1: best modelS ordered by DIC
#' 2: best models ordered by WAIC
#' 3. best models ordered by LCPO
#' @export
#' @author Joaquín Martínez-Minaya <\email{jomarminaya@@gmail.com}>
Bdiclcpomodel_comp <- function(variables, 
                              n = 256,
                              shared.var = FALSE,
                              spatial,
                              spatial.short,
                              friends = NULL,
                              ...)
{
  if(shared.var == TRUE){
    variables2 <- variables
    variables <- paste0(variables, "_s")
  }else{
    #To include in the compositional model
    variables2 <- lapply(variables, function(x){
      paste0("f(", "id.", x, ",", x, ",",
             "model   = 'iid',",
             "initial = log(1/1000),",
             "fixed   = TRUE)")
    }) %>% unlist(.)
  }

  # #Intercept
  # var_intercept <- variables2[1]
  # variables <- variables[-1]
  # variables2 <- variables2[-1]
  
  #T?rminos que usaremos
  sel.terms <- switch('terms', terms = variables)
  sel.terms2 <- switch('terms', terms = variables2)
  
  if(!is.na(spatial)){
    sel.terms <- c(sel.terms, spatial.short, "dir")
    sel.terms2 <- c(sel.terms2, 
                    spatial,
                    "  f(id.z, 
    model = 'iid', 
    hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(0.1, 0.01))))")
  }else{
    sel.terms <- c(sel.terms, "dir")
    sel.terms2 <- c(sel.terms2,
                    "  f(id.z, 
    model = 'iid', 
    hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(0.1, 0.01))))")
    
  }
  
  
  # todas las combinaciones de los m elementos de v
  comb.terms <- function(m, v = sel.terms) {
    intercept1 <- paste0('resp ~ -1')
    v2 <- v[-length(v)]
    
    if(m==0){
      return(paste(intercept1, v[length(v)],
                          sep =' + '))
    }else {
      combis <- apply(combn(v2, m), 2, paste, collapse =' + ')
      return(paste(intercept1, combis, v[length(v)], 
                   sep =' + '))
    }
  }

  #Lista con todos los modelos posibles
  f.list <- unlist(sapply(0:(length(sel.terms)-1), comb.terms, v = sel.terms))
  f.list2 <- unlist(sapply(0:(length(sel.terms2)-1), comb.terms, v = sel.terms2))
  
  f.list.tot <- cbind(f.list, f.list2)
  # lanzamos cada uno de los modelos guardados en el objeto 'f.list' y nos quedamos con el DIC
  dic <- numeric()
  waic <- numeric()
  LCPO <- numeric()
  
  dic.o <- numeric()
  waic.o <- numeric()
  LCPO.o <- numeric()
  res <- list()
  for(i in 1:length(f.list.tot[,1])){
  #for(i in 1:5){
    res[[i]] =  tryCatch(inla(eval(parse(text= f.list.tot[,2][i])),...),
                         error=function(e){})

    
    if(is.null(res[[i]])){
      dic[i] <- NA
      LCPO[i] <- NA
      waic[i] <- NA
  
      dic.o[i] <- NA
      LCPO.o[i] <- NA
      waic.o[i] <- NA
    }else{
      #cat("--- Sampling for computing DIC and WAIC --- \n")
      xx <- inla.posterior.sample(1000, res[[i]])
      inf <- mclapply(xx, INLAcomp::extract_lp_sigma)
      
      #cat("--- Computing new DIC and WAIC --- \n")
      #New DIC
      dic.mod1 <- INLAcomp::dic.mult(inf, y = data[, c(paste0("alr.gc", 1:3))])
      
      #New WAIC
      waic.mod1 <- INLAcomp::waic.mult(inf, y = data[, c(paste0("alr.gc", 1:3))])
      
      #cat("--- Computing new LCPO --- \n")
      # New LCPO
      LCPO[i] <- INLA::inla.group.cv(result = res[[i]],
                                num.level.sets = -1,
                                strategy = "posterior",
                                friends = friends) %>%
        .[["cv"]] %>% log(.) %>% mean(.) %>% -.
      
      
      dic[i] <- dic.mod1$dic
      waic[i]<- waic.mod1$waic
      
      #Obsolete DIC, LCPO and WAIC
      dic.o[i] <- res[[i]]$dic$dic
      LCPO.o[i] = -mean(log(res[[i]]$cpo$cpo))
      waic.o[i]<- res[[i]]$waic$waic
      
      
      
      # LGOCV[i] <-    INLA::inla.group.cv(result = res[[i]],
      #                                    num.level.sets = 1,
      #                                    strategy = "posterior") %>%
      #   .[["cv"]] %>% log(.) %>% mean(.) %>% -.
      
    }
        print(c(f.list.tot[i,1], dic[i], waic[i],
            LCPO[i], dic.o[i], waic.o[i], LCPO.o[i]))
  }


  models_total <-data.frame(f.list.tot[,1][1:n], 
                          dic[1:n], 
                          waic[1:n], 
                          LCPO[1:n],
                          dic.o[1:n], 
                          waic.o[1:n], 
                          LCPO.o[1:n])
  colnames(models_total)<-c("Models", "Dic", "Waic", "LCPO", "Dic.o", "Waic.o", "LCPO.o")

  # #mostramos los modelos ordenados SEGÚN EL waic
  # modelos_waic<-data.frame(f.list.tot[,1][order(waic)[1:n]], 
  #                          dic[order(waic)[1:n]], 
  #                          waic[order(waic)[1:n]], 
  #                          LCPO[order(waic)[1:n]],
  #                          LGOCV[order(waic)[1:n]])
  # colnames(modelos_waic)<-c("Models", "Dic", "Waic", "LCPO", "LGOCV")
  # 
  # 
  # #mostramos los modelos ordenados SEGÚN EL LCPO
  # modelos_lcpo<-data.frame(f.list.tot[,1][order(LCPO)[1:n]], dic[order(LCPO)[1:n]], 
  #                          waic[order(LCPO)[1:n]], 
  #                          LCPO[order(LCPO)[1:n]],
  #                          LGOCV[order(LCPO)[1:n]])
  # colnames(modelos_lcpo)<-c("Models", "Dic", "Waic", "LCPO", "LGOCV")
  # 
  # 
  # #mostramos los modelos ordenados SEGÚN EL LCPO
  # modelos_lgocv<-data.frame(f.list.tot[,1][order(LGOCV)[1:n]], 
  #                           dic[order(LGOCV)[1:n]], 
  #                          waic[order(LGOCV)[1:n]], 
  #                          LCPO[order(LGOCV)[1:n]],
  #                          LGOCV[order(LGOCV)[1:n]])
  # colnames(modelos_lcpo)<-c("Models", "Dic", "Waic", "LCPO", "LGOCV")
  # 
  # 
  # modelos<-list(modelos_dic, modelos_waic, modelos_lcpo,
  #               modelos_lgocv)
  # names(modelos)<-c("Models using DIC", "Models using WAIC", 
  #                   "Models using LCPO",
  #                   "Models using LGOCV")
  models_total

}
