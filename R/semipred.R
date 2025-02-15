semipred <- function(semireg, modelterm=NULL, covariate, sm_term=NULL, contr=NULL,
                     covariateV=NULL, boundary=NULL, level=0.05, trans=NULL, trellis=TRUE, 
                     scales=c("fixed", "free", "free_x", "free_y"), 
                     plotord=NULL, ci=TRUE, point=TRUE, jitterv=0, threeD=FALSE, prt=TRUE) {
  
  stopifnot(inherits(semireg, "semireg"), !is.null(semireg$CovMat))  
  scales <- as.character(scales)
  scales <- match.arg(scales)
  
  if (is.null(modelterm) || all(modelterm%in%c("NULL", ""))) {  
    modelterm <- covariate
  }
  
  mod_df <- semireg$data
  if (!is.null(contr)) {
    trans <- NULL
    stopifnot("'modelterm' must be a single term"=length(unlist(strsplit(modelterm, "\\:")))==1) 
    
    if(length(contr) != 2){
      stop("contr must be a character vector of length 2")
    }
    
    if(is.character(contr)){
      grp_fct <- mod_df[[modelterm]]
      grp_fct_levl <- levels(grp_fct)
      
      i = match(contr, grp_fct_levl)
      if(any(is.na(i))){
        stop(paste("The entries", paste0(contr[is.na(i)], sep = ", "), "do not match any level in ", modelterm))
      }else{
        contr = i
      }
    }
  }
  
  if (!is.null(covariateV)){
    if (length(covariate) > 1) {
      stopifnot(is.matrix(covariateV), dim(covariateV)[2] == length(covariate))
      colnames(covariateV) <- covariate
    }
    
    for (i in covariate) {
      if (!is.null(attr(mod_df, paste(i, "mean", sep="_")))) {
        if (!is.matrix(covariateV)) covariateV <- covariateV-attr(mod_df, paste(i, "mean", sep="_"))
        else covariateV[,i] <- covariateV[,i]-attr(mod_df, paste(i, "mean", sep="_"))
      }
      if (!is.null(attr(mod_df, paste(i, "sd", sep="_")))) {
        if (length(covariate) ==1) covariateV <- covariateV/attr(mod_df, paste(i, "sd", sep="_"))
        else covariateV[,i] <- covariateV[,i]/attr(mod_df, paste(i, "sd", sep="_"))
      } 
    }
  }
  
  semer <- semireg$semer 
  vcov_indN <- semireg$Cov_indN
  cov_lst <- semireg$cov_lst
  fullCovMat <- semireg$CovMat
  
  if (inherits(semer, "glmmTMB")) betaHat <- fixef(semer)$cond else betaHat <- fixef(semer)
  sig.epsHat <- sigma(semer) 
  
  sm_varsN <- setdiff(semireg$sm_vars, semireg$fomul_vars)
  
  if (length(sm_varsN) > 0){
    lm_fl <- format(formula(terms(semer)))
    for (i in sm_varsN)  lm_fl <- paste(lm_fl, "+ ",i)
    mod_lm <- lm(as.formula(lm_fl), data=mod_df)
    KK <- Kmatrix(mod_lm, modelterm, covariate, covariateV)
  }else{  
    KK <- Kmatrix(semer, modelterm, covariate, covariateV)
  }
  
  response <- response0 <- KK$response 
  K <- KK$K[, colnames(model.matrix(semer))]
  pred_df <- KK$fctnames
  if (!(response %in% names(mod_df))) {   # all.vars(str2expression(response0))]
    mod_df[["response_y"]] <- eval(parse(text=response), mod_df)
    names(pred_df)[names(pred_df)==response] <- "response_y" 
    response <- "response_y"
  }
  
  if (is.null(sm_term)) {
    sm_term <- names(as.list(semireg$smoothZ_call))[-1]
    sm_term <- intersect(names(vcov_indN), sm_term)
    sm_Cov <- cov_lst$sm_Cov
  }else{ 
    if (length(sm_term)==1) {
      sm_Cov <- cov_lst[[sm_term]] 
    }else{
      sm_term <- intersect(names(vcov_indN), sm_term)
      sm_termInds <- 1:ncol(getME(semer, "X")) 
      for (i in sm_term) {
        stmch_id <- match(i, names(vcov_indN))
        sm_termInds <- c(sm_termInds, (vcov_indN[stmch_id-1]+1):vcov_indN[stmch_id])
      }  
      sm_Cov <- fullCovMat[sm_termInds, sm_termInds]
    }   
  }
  
  u_lst <- list()
  for (i in sm_term) {
    u_indx <- grep(i, names(semireg$u_lst))
    u_lst <- c(u_lst, semireg$u_lst[u_indx])
  }
  
  knots_lst <- semireg$knots_lst[sm_term]
  range_lst <- semireg$range_lst[sm_term] 
  type_lst <- semireg$type_lst[sm_term]  
  call_lst <- as.list(semireg$smoothZ_call)[sm_term]
  
  call_lstN <- Z_lstN <- vector("list", length(call_lst))
  names(call_lstN) <- names(Z_lstN) <- sm_term 
  
  for (i in sm_term) {
    xcall_lst  <- as.list(call_lst[[i]])
    if (!is.null(xcall_lst$type) && xcall_lst$type=="smspline") xcall_lst$pred <- TRUE
    if (!is.null(xcall_lst$k)) xcall_lst$k <- NULL
    if (!is.null(range_lst[[i]]) && type_lst[[i]]!="Ztps") xcall_lst$range.x <- range_lst[[i]]
    xcall_lst$intKnots <- knots_lst[[i]]
    call_lstN[[i]] <- as.call(xcall_lst) 
  }  
  
  for (i in sm_term) {
    Z_lstN[[i]] <- with(pred_df, eval(call_lstN[[i]]))
  }
  
  Z_lstNN <- list()
  for (i in sm_term){
    if (is.list(Z_lstN[[i]])) Z_lstNN <- c(Z_lstNN, unlist(Z_lstN[i], recursive=FALSE)) else Z_lstNN <- c(Z_lstNN, Z_lstN[i])
  }
  
  Z_lstNN <- Z_lstNN[intersect(names(semireg$u_lst), names(Z_lstNN))]
  ZN <- do.call("cbind", Z_lstNN)
  D_mx <- cbind(K,ZN)
  u_hat <- do.call("c", u_lst)
  beta_u <- c(betaHat, u_hat)
  
  stderr <- sig.epsHat*sqrt(diag(D_mx%*%sm_Cov%*%t(D_mx)))
  y_hat <- as.numeric(D_mx%*%beta_u)
  LL <- y_hat - qnorm(1-level)*stderr
  UL <- y_hat + qnorm(1-level)*stderr  
  
  if (is.null(trans)) {
    if (!isLMM(semer)) {
      if (inherits(semer, "glmmTMB")) {
        family_n <- semer$modelInfo$family$family
        invfun <- semer$modelInfo$family$linkinv
        weight_fun <- function(eta, mod){
          family_fun <- mod$modelInfo$family
          family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
        } 	  
      }else{
        family_n <- slot(semer, "resp")$family$family
        invfun <- slot(semer, "resp")$family$linkinv	
        weight_fun <- function(eta, mod){
          family_fun <- slot(mod, "resp")$family
          family_fun$mu.eta(eta)^2/family_fun$variance(family_fun$linkinv(eta))
        }
      }		
      stderr_mu <- sig.epsHat*sqrt(diag(D_mx%*%sm_Cov%*%t(D_mx)))*weight_fun(y_hat, semer)
      y_hat <- invfun(y_hat)
      LL <- y_hat - qnorm(1-level)*stderr_mu
      UL <- y_hat + qnorm(1-level)*stderr_mu
      if (family_n == "binomial"){
        UL <- ifelse(UL > 1, 1, UL)	
        LL <- ifelse(LL < 0, 0, LL)	
        if (is.factor(mod_df[[response]])) {
          mod_df[[response]] <- as.numeric(mod_df[[response]]) - 1
        }else if (!is.null(dim(mod_df[, response])) && ncol(mod_df[, response])==2) {
          mod_df[["Proportion"]] <- mod_df[, response][,1]/rowSums(mod_df[, response])
          response <- "Proportion"
        }		
      } 
      pred_df$SE_mu <- stderr_mu
    }else pred_df$SE <- stderr
  }else{
    y_hat <- trans(y_hat)
    LL <- trans(LL)
    UL <- trans(UL)
    mod_df[["Trans_y"]] <- trans(mod_df[[response]])
    response <- "Trans_y"
  }
  
  pred_df[[response]] <- y_hat
  pred_df$LL <- LL
  pred_df$UL <- UL
  
  if (!is.null(attr(mod_df, "numeric_var"))) {
    for (i in attr(mod_df, "numeric_var")) {
      if (!is.null(attr(mod_df, paste(i, "sd", sep="_")))) {
        mod_df[[i]] <- mod_df[[i]]*attr(mod_df, paste(i, "sd", sep="_"))
        pred_df[[i]] <- pred_df[[i]]*attr(mod_df, paste(i, "sd", sep="_"))  
      }  
      
      if (!is.null(attr(mod_df, paste(i, "mean", sep="_")))) {
        mod_df[[i]] <- mod_df[[i]]+attr(mod_df, paste(i, "mean", sep="_"))
        pred_df[[i]] <- pred_df[[i]]+attr(mod_df, paste(i, "mean", sep="_"))  
      }  
    }
  }
  
  if (!is.null(attr(mod_df, paste(response, "sd", sep="_")))) {
    pred_df$LL <- pred_df$LL*attr(mod_df, paste(response, "sd", sep="_")) 
    pred_df$UL <- pred_df$UL*attr(mod_df, paste(response, "sd", sep="_")) 
  }
  
  if (!is.null(attr(mod_df, paste(response, "mean", sep="_")))) {
    pred_df$LL <- pred_df$LL+attr(mod_df, paste(response, "mean", sep="_")) 
    pred_df$UL <- pred_df$UL+attr(mod_df, paste(response, "mean", sep="_")) 
  }
  
  for (i in names(pred_df)){
    if (is(pred_df[,i], "AsIs")) pred_df[,i] <- as.numeric(pred_df[,i])
  }
  
  if (setequal(modelterm, covariate)) {
    if (length(covariate)==1) {
      plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response))))+  
        geom_line(linewidth=1, col="red") +
        geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.6, stat="identity", fill="#00FFFF", col="#00FFFF")+
        labs(x=covariate, y=response)+
        theme_bw()
      if (point) plt <- plt + geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response))), col="DeepPink", position = position_jitter(width = jitterv, height = jitterv))
      if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))
    }else if(length(covariate)==2) {
      if (!is.null(boundary)){
        inBdry <- HRW::pointsInPoly(pred_df[, covariate], boundary)
        pred_df <- droplevels(pred_df[inBdry,])
      }
      
      if (threeD) {
        plt <- plotly::plot_ly(pred_df, x = formula(paste("~", covariate[1])), 
                               y = formula(paste("~", covariate[2])),
                               z = formula(paste("~", response)), 
                               intensity = formula(paste("~", response)), 
                               type = 'mesh3d',
                               opacity = 0.7) 
        
        if (point) plt <- plotly::add_trace(plt, x = mod_df[[covariate[1]]], y = mod_df[[covariate[2]]], 
                                            z = mod_df[[response]], mode = "markers", type = "scatter3d", 
                                            marker = list(size = 2.5, color = "red", symbol = 104))						
      }else{		
        plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) +
          geom_contour_filled() +
          labs(x=covariate[1], y=covariate[2], z=response)
        if (point) plt <- plt + geom_point(data=mod_df, col="red", size=1.5)   
        if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))		
      }
      
    }
  }else{
    vars <- unlist(strsplit(modelterm, "\\:"))
    if (length(vars)==1) {
      if (length(covariate)==1) {
        if (is.null(contr)) {
          plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))))+
            geom_line(aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))), linewidth=1) +
            labs(x=covariate, y=response, col=modelterm) +
            theme_bw()	
          if (ci) plt <- plt + 
              geom_ribbon(aes(ymin = LL, ymax = UL, fill=eval(parse(text = modelterm))), alpha = 0.3, stat="identity")+
              labs(fill=modelterm)
          if (point) plt <- plt + 
              geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = modelterm))), position = position_jitter(width = jitterv, height = jitterv))+
              labs(col=modelterm)
          if (trellis)  plt <- plt + 
              facet_wrap(formula(paste("~", modelterm)), scales=scales) 
          if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))  
        }else{
          grp_fct <- pred_df[[modelterm]]
          grp_fct_levl <- levels(grp_fct)
          plt_df_lst <- lapply(split(pred_df, grp_fct), function(x) x[, !names(x) %in% modelterm])
          plt_df_contr <- plt_df_lst[[grp_fct_levl[contr[1]]]] - plt_df_lst[[grp_fct_levl[contr[2]]]]		  
          diff_covarn <- setdiff(names(plt_df_contr), c("SE", response, "LL", "UL"))
          plt_df_contr[,diff_covarn] <- plt_df_lst[[grp_fct_levl[contr[1]]]][, diff_covarn]
          
          D_mx_lst <- lapply(split(data.frame(as.matrix(D_mx)), grp_fct), as.matrix)
          D_contr <- D_mx_lst[[grp_fct_levl[contr[1]]]] - D_mx_lst[[grp_fct_levl[contr[2]]]]
          
          if (!isLMM(semer)) {
            stderr_contr <- sig.epsHat*sqrt(diag(D_contr%*%sm_Cov%*%t(D_contr)))*weight_fun(plt_df_contr[[response]], semer)
          }else{
            stderr_contr <- sig.epsHat*sqrt(diag(D_contr%*%sm_Cov%*%t(D_contr)))
          }
          if (!is.null(attr(mod_df, paste(response, "sd", sep="_")))) {
            plt_df_contr$LL <- plt_df_contr[[response]] - qnorm(1-level)*stderr_contr*attr(mod_df, paste(response, "sd", sep="_")) 
            plt_df_contr$UL <- plt_df_contr[[response]] + qnorm(1-level)*stderr_contr*attr(mod_df, paste(response, "sd", sep="_")) 
          }else{
            plt_df_contr$LL <- plt_df_contr[[response]] - qnorm(1-level)*stderr_contr
            plt_df_contr$UL <- plt_df_contr[[response]] + qnorm(1-level)*stderr_contr
          }          
          
          plt <- ggplot(plt_df_contr, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response))))+
            geom_line(linewidth=1, col="blue") +
            geom_hline(yintercept=0, linetype="dashed", color = "red") +
            geom_ribbon(aes(ymin = LL, ymax = UL), alpha = 0.6, stat="identity", fill="chartreuse3", col="chartreuse3") +
            ggtitle(paste("Factor", modelterm, "level", grp_fct_levl[contr[1]], "vs level", grp_fct_levl[contr[2]], "with", (1-level) * 100, "% CI")) +
            labs(x=covariate, y=paste("The difference of ", response, "\n", sep="")) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5)) 
          if (response == "response_y") plt <- plt + labs(y=paste("The difference of ", response0, "\n", sep=""))
          pred_df <- plt_df_contr
          pred_df$SE_contr <- stderr_contr
          pred_df$SE <- NULL
          pred_df$SE_mu <- NULL
        }
      }else if(length(covariate)==2) {
        if (is.null(contr)) {
          if (!is.null(boundary)){
            inBdry <- HRW::pointsInPoly(pred_df[, covariate], boundary)
            pred_df <- droplevels(pred_df[inBdry,])
          }
          
          if (threeD) {
            grp_fct <- pred_df[[modelterm]]
            plt <- lapply(split(pred_df, grp_fct), function(x) {
              plotly::plot_ly(droplevels(x), x = formula(paste("~", covariate[1])), 
                              y = formula(paste("~", covariate[2])),
                              z = formula(paste("~", response)), 
                              intensity = formula(paste("~", response)), 
                              type = 'mesh3d',
                              opacity = 0.7
              ) 
            })
            
            if (point) {
              mod_df_lst <- split(mod_df, mod_df[[modelterm]])
              for (i in names(mod_df_lst)) {
                plt[[i]] <- plotly::add_trace(plt[[i]], x = mod_df_lst[[i]][[covariate[1]]], y = mod_df_lst[[i]][[covariate[2]]], 
                                              z = mod_df_lst[[i]][[response]], mode = "markers", type = "scatter3d", 
                                              marker = list(size = 2.5, color = "red", symbol = 104)) 
                
              }
            }
          }else{		
            plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) +
              geom_contour_filled()+
              facet_wrap(formula(paste("~", modelterm)), scales=scales)+
              labs(x=covariate[1], y=covariate[2], z=response) 
            if (point) plt <- plt + geom_point(data=mod_df, col="red", size=1.5)
            if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))			
          }          
        }else{
          grp_fct <- pred_df[[modelterm]]
          grp_fct_levl <- levels(grp_fct)
          plt_df_lst <- lapply(split(pred_df, grp_fct), function(x) x[, !names(x) %in% modelterm])
          plt_df_contr <- plt_df_lst[[grp_fct_levl[contr[1]]]] - plt_df_lst[[grp_fct_levl[contr[2]]]]
          diff_covarn <- setdiff(names(plt_df_contr), c("SE", response, "LL", "UL"))
          plt_df_contr[,diff_covarn] <- plt_df_lst[[grp_fct_levl[contr[1]]]][, diff_covarn]
          
          D_mx_lst <- lapply(split(data.frame(as.matrix(D_mx)), grp_fct), as.matrix)
          D_contr <- D_mx_lst[[grp_fct_levl[contr[1]]]] - D_mx_lst[[grp_fct_levl[contr[2]]]]
          
          if (!isLMM(semer)) {
            stderr_contr <- sig.epsHat*sqrt(diag(D_contr%*%sm_Cov%*%t(D_contr)))*weight_fun(plt_df_contr[[response]], semer)
          }else{
            stderr_contr <- sig.epsHat*sqrt(diag(D_contr%*%sm_Cov%*%t(D_contr)))
          }
          if (!is.null(attr(mod_df, paste(response, "sd", sep="_")))) {
            plt_df_contr$LL <- plt_df_contr[[response]] - qnorm(1-level)*stderr_contr*attr(mod_df, paste(response, "sd", sep="_")) 
            plt_df_contr$UL <- plt_df_contr[[response]] + qnorm(1-level)*stderr_contr*attr(mod_df, paste(response, "sd", sep="_")) 
          }else{
            plt_df_contr$LL <- plt_df_contr[[response]] - qnorm(1-level)*stderr_contr
            plt_df_contr$UL <- plt_df_contr[[response]] + qnorm(1-level)*stderr_contr
          }
          
          if (!is.null(boundary)){
            inBdry <- HRW::pointsInPoly(plt_df_contr[, covariate], boundary)
            plt_df_contrN <- droplevels(plt_df_contr[inBdry,])
          }else plt_df_contrN <- plt_df_contr		
          plt_df_contrSub <- subset(plt_df_contrN, LL > 0 | UL < 0)
          
          plt <- ggplot(plt_df_contrN, aes(x=eval(parse(text = covariate[1])), y=eval(parse(text = covariate[2])), z=eval(parse(text = response)))) + 
            geom_contour_filled() +
            geom_point(data=plt_df_contrSub, col="red", size=2.5, shape=22) +
            ggtitle(paste("Factor", modelterm, "level", grp_fct_levl[contr[1]], "vs level", grp_fct_levl[contr[2]], "with Points at the Area out of", (1-level) * 100, "% CI")) +
            labs(x=covariate[1], y=paste("The difference of ", response, "\n", sep=""), z=response) +
            theme_bw() +
            theme(plot.title = element_text(hjust = 0.5)) 			
          if (response == "response_y") plt <- plt + labs(y=paste("The difference of ", response0, "\n", sep=""))
          pred_df <- plt_df_contr
          pred_df$SE_contr <- stderr_contr
          pred_df$SE <- NULL
          pred_df$SE_mu <- NULL 
        }
      }
    }else if (length(vars)==2) {	 
      if (is.null(plotord) || all(plotord%in%c("NULL", ""))) plotord <- 1:2
      fact1 <- (vars[plotord])[1]
      fact2 <- (vars[plotord])[2]
      plt <- ggplot(pred_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))))+         	  
        geom_line(aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))), linewidth=1) +
        labs(x=covariate, y=response, col=fact1) + 
        facet_wrap(formula(paste("~", fact2)), scales=scales)+		
        theme_bw()
      if (response == "response_y") plt <- plt + labs(y=paste(response0, "\n", sep=""))
      if (ci) plt <- plt +
        geom_ribbon(data=pred_df, aes(ymin = LL, ymax = UL, fill=eval(parse(text = fact1))), alpha = 0.3, stat="identity") +
        labs(fill=fact1)
      if (point) plt <- plt + 
        geom_point(data=mod_df, aes(x=eval(parse(text = covariate)), y=eval(parse(text = response)), col=eval(parse(text = fact1))), position = position_jitter(width = jitterv, height = jitterv))
    }
  }
  #     facet_grid(formula(paste(Variable1, "~", Variable2)))
  
  if (prt) {
    options(warn = -1)	  
    if (length(covariate)==2 && threeD && !setequal(modelterm, covariate)) {
      for (i in names(plt)) print(plt[[i]])
    }else print(plt)
    options(warn = 0) 
  }
  rownames(pred_df) <- NULL
  return(invisible(list(plt=plt, pred_df=pred_df)))
}


