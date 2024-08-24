library(survival)
library(MASS)
library(simsurv)
library(dplyr)
library(extRemes)
library(survival)
library(ggplot2)
library(numDeriv)
library(flexsurv)
#' simulate error from Weibull, log-logistic, and lognormal distributions for AFT model
#'
#' @param N number of observations
#' @param scale
#' @param location
#' @param type  The type of the error distribution
#' @param beta vector of regression coefficients
#' @param rateC rate of the censoring distribution
#' @param censoring_type either exponential or uniform distribution
#'
#' @return simulated dataset
#' @export
#' @import survival
#' @import MASS
#' @import simsurv
#' @examples
simaft<-function(N,scale,location,beta,rateC,ad_C = 30,type="weibull",censoring_type = "exponential"){

  Sigma = matrix(c(1,0.3,0.3,1),2,2)
  x1x2 = as.data.frame(mvrnorm(n=N,mu=c(0,0),Sigma=Sigma) )
  x1 = x1x2[,1]
  x2 = x1x2[,2]
  x3 = (rbinom(N,1,0.8))
  x4 = rep(0,N)
  x4[x3==0] = sample(c(1:4),table(x3)[1], replace = TRUE, prob = c(0.2,0.2,0.3,0.3))
  x4[x3==1] = sample(c(1:4),table(x3)[2], replace = TRUE, prob = c(0.1,0.2,0.4,0.5))

  x_data = as.data.frame(cbind(x1,x2,x3,x4))
  x_data$x3 = factor(x_data$x3)
  x_data$x4 = as.factor(x_data$x4)
  covs = model.matrix(~x1+x2+x3+x4,x_data)
  covs <- data.frame(id = 1:N, covs[,2:dim(covs)[2]])
  if (type == "weibull"){

    error<-revd(N, loc = 0, scale = scale)   #, loc = location, scale = scale
    time_ev<-exp(location+as.matrix(covs[c(2:dim(covs)[2])])%*%(beta)-1*error)

  }else if (type =="loglogistic"){
    error<-rlogis(N)
    time_ev<-exp(location+as.matrix(covs[c(2:dim(covs)[2])])%*%(beta)+scale*error)
  }else if(type =="lognormal"){
    error<-rnorm(N)
    time_ev<-exp(location+as.matrix(covs[c(2:dim(covs)[2])])%*%(beta)+scale*error)

  }

  if(censoring_type == "exponential"){
    C <- rexp(n=N, rate=rateC)

  }else if(censoring_type == "uniform"){
    C<- runif(N,0,max=rateC)
  }
  C <- ifelse(C>ad_C,ad_C,C)
  # follow-up times and event indicators
  time <- pmin(time_ev, C)
  status <- as.numeric(time_ev <= C)

  # data set
  data.frame(id=1:N,
             time=time,
             status=status,
             covs [,2:dim(covs)[2]])
}


#' pdf for generalized gamma distribution
#'
#' @param scaled_xbeta
#' @param q
#' @param type  The type of the error distribution
#' @param t
#' @param xbeta the product of coefficients and covariates
#' @param sigma scale parameter
#'
#' @return pdf
#' @export
#' @examples
pdf_eval<-function(scaled_xbeta,type,t=NA,xbeta=NA,sigma=NA,q=NA){
  if (type == "weibull"){
    res = exp(scaled_xbeta-exp(scaled_xbeta))

  }else if(type =="loglogistic"){
    res = exp(scaled_xbeta)*(1+exp(scaled_xbeta))^(-2)
  }else if(type == 'lognormal'){
    res =  dnorm(scaled_xbeta)
  }else if(type =="gengamma"){
    sigma = exp(sigma)
    res = flexsurv::dgengamma(t,xbeta,sigma,q)*sigma*t
  }
  return(res)
}

#' cdf for generalized gamma distribution
#'
#' @param scaled_xbeta
#' @param q
#' @param type  The type of the error distribution
#' @param t
#' @param xbeta the product of coefficients and covariates
#' @param sigma scale parameter
#'
#' @return pdf
#' @export
#' @examples
cdf_eval<-function(scaled_xbeta,type,t=NA,xbeta=NA,sigma=NA,q=NA){
  if (type == "weibull"){
    res = pdf_eval(scaled_xbeta,type)/exp(scaled_xbeta)
  }else if(type =="loglogistic"){
    res = (1+exp(scaled_xbeta))^(-1)
  }else if(type == "lognormal"){
    res = pnorm(scaled_xbeta,lower.tail = FALSE)
  }else if(type == "gengamma"){
    sigma = exp(sigma)
    res =  flexsurv::pgengamma(t,xbeta,sigma,q,lower.tail = FALSE)

  }

  return(res)
}

#' log-likelihood when the error term follows generalized gamma distribution
#'
#' @param par parameter
#' @param t time
#' @param d event status
#' @param X vector of covariates
#'
#' @return res log-likelihood
#' @export
#' @examples
loglik_for_gengamma<-function(par,t,d,X){ #t,d,X


  sigma = par[length(par)-1]
  q = par[length(par)]
  beta = par[1:(length(par)-2)]
  xbeta =  as.numeric(as.matrix(X) %*% beta)
  scaled_xbeta = (log(t)-xbeta)/sigma
   res = -sum(d*log(pdf_eval(scaled_xbeta,"gengamma",t,xbeta,sigma,q))-d*(sigma)+(1-d)*log(cdf_eval(scaled_xbeta,"gengamma",t,xbeta,sigma,q)))

  return(res)


}
#' individual log-likelihood when the error term follows generalized gamma distribution
#'
#' @param par parameter
#' @param t1 time
#' @param d1 event status
#' @param X1 vector of covariates
#'
#' @return res log-likelihood
#' @export
#' @examples
loglik_for_gengamma_ind <-function(par,t1,d1,X1){ #t,d,X


  sigma = par[length(par)-1]
  q = par[length(par)]
  beta = par[1:(length(par)-2)]
  xbeta =  as.numeric(t(as.matrix(X1)) %*% beta)
  scaled_xbeta = (log(t1)-xbeta)/sigma
  res = -(d1*log(pdf_eval(scaled_xbeta,"gengamma",t1,xbeta,sigma,q))-d1*(sigma)+(1-d1)*log(cdf_eval(scaled_xbeta,"gengamma",t1,xbeta,sigma,q)))

  return(res)


}

#' log-likelihood for the AFT model
#'
#' @param par parameter
#' @param t time
#' @param d event status
#' @param X vector of covariates
#' @param type  The type of the error distribution
#'
#' @return res log-likelihood
#' @export
#' @examples

loglik_eval<-function(par,t,d,X,type){

  if (type == "weibull" |type =="loglogistic"|type == "lognormal"){
    sigma = par[length(par)]
    beta = par[1:(length(par)-1)]
    xbeta = as.numeric(as.matrix(X) %*% beta)
    scaled_xbeta = (log(t)-xbeta)/sigma#log
    res = -sum(d*log(pdf_eval(scaled_xbeta,type))-d*log(sigma)+(1-d)*log(cdf_eval(scaled_xbeta,type)))
    return(res)
  }else if(type == "gengamma"){
    sigma = par[length(par)-1]
    q = par[length(par)]
    beta = par[1:(length(par)-2)]
    xbeta =  as.numeric(as.matrix(X) %*% beta)
    scaled_xbeta = (log(t)-xbeta)/sigma
    res =-sum(d*log(pdf_eval(scaled_xbeta,"gengamma",t,xbeta,sigma,q))-d*(sigma)+(1-d)*log(cdf_eval(scaled_xbeta,"gengamma",t,xbeta,sigma,q)))

    #-sum(d*log(pdf_eval(scaled_xbeta,type,q))-d*log(sigma)+(1-d)*log(cdf_eval(scaled_xbeta,type,q)))
    return(res)
  }

}

#' derivative of log-likelihood for the AFT model
#'
#' @param par parameter
#' @param t time
#' @param d event status
#' @param X vector of covariates
#' @param type The type of the error distribution
#'
#' @return res log-likelihood
#' @export
#' @examples

dloglik_aft_eval<-function(par,t,d,X,type,individual=FALSE){#par,n_gamma,degree,t,d,X,knots
  if (type == "weibull"|type == "loglogistic"|type == "lognormal"){
    sigma = par[length(par)]
    beta = par[1:(length(par)-1)]
    xbeta = as.numeric(as.matrix(X) %*% beta)
    scaled_xbeta = (log(t)-xbeta)/sigma
    if (type == "weibull"){
      ai = exp(scaled_xbeta)-d

    }else if (type == "loglogistic"){
      ai = -d+(1+d)*exp(scaled_xbeta)*(1+exp(scaled_xbeta))^(-1)
    }else if (type == "lognormal"){
      fwi =  (2*pi)^(-1/2)*exp(-scaled_xbeta^2/2)
      #ffwi<-function(par){

      # return((2*pi)^(-1/2)*exp(-scaled_xbeta^2/2))
      #}
      int_f_wi = cdf_eval(scaled_xbeta,type)
      ai = d*scaled_xbeta+(1-d)*fwi/int_f_wi
    }
    U_beta = 1/sigma*as.vector(colSums(X * (ai)))
    U_sigma =1/sigma* sum(scaled_xbeta*(ai)-d)
    res = c(U_beta,U_sigma)
    if (individual){
      res = 1/sigma* cbind(X * (ai),scaled_xbeta*(ai)-d)
    }}else if(type == "gengamma"){
      if(!individual){


        res = -grad(loglik_for_gengamma,x=par,t=t,d=d,X=X)
        #sqrt(diag(solve(hessian(loglik_for_gengamma,x=par,t=t,d=d,X=X))))
      }
      else if(individual){
        grad_ind_Res_list <- lapply(c(1:length(t)), function(i) {
          -grad(loglik_for_gengamma_ind, x=par, t=t[i], d=d[i], X=X[i, ])
        })

        res <- do.call(rbind, grad_ind_Res_list)

      }
    }

  return(res)
}
G = function(est,t,d,X,type){
  return(apply(dloglik_aft_eval(est,t,d,X,type),2,sum))

}

Godambe = function(est,t,d,X,type){

  meat.half=dloglik_aft_eval(est,t,d,X,type,TRUE)
  sum_outer_products = as.matrix(t(meat.half)) %*% as.matrix(meat.half)


  bread=ddloglik_aft_eval(est,t,d,X,type)

  t(solve(bread))%*%sum_outer_products%*%(solve(bread))

  return(  t(IF)%*%(IF) )
}
#' Get the Hessian matrix of the log-likelihood for the AFT model
#'
#' @param par parameter
#' @param t time
#' @param d event status
#' @param X vector of covariates
#' @param type The type of the error distribution
#'
#' @return res log-likelihood
#' @export
#' @examples
ddloglik_aft_eval<-function(par,t,d,X,type){#par,n_gamma,degree,t,d,X,knots
  if(type =="weibull"|type=="loglogistic"|type=="lognormal"){
    sigma = par[length(par)]
    beta = par[1:(length(par)-1)]
    xbeta = as.numeric(as.matrix(X) %*% beta)
    scaled_xbeta = (log(t)-xbeta)/sigma
    if (type == "weibull"){

      Ai = exp(scaled_xbeta)
    }else if (type =="loglogistic"){
      Ai = (1+d)* exp(scaled_xbeta)*(1+ exp(scaled_xbeta))^(-2)
    }else if (type =="lognormal"){
      fwi =  (2*pi)^(-1/2)*exp(-scaled_xbeta^2/2)

      int_f_wi = cdf_eval(scaled_xbeta,type)

      lambda_wi  = fwi/int_f_wi
      Ai = d+(1-d)*lambda_wi*(lambda_wi-scaled_xbeta)
    }
    ret_beta_all <- matrix(0,length(beta),length(beta))
    for(j in 1:length(beta)){
      for (i in j:length(beta)){
        ret_beta_all[i,j] <- 1/sigma^2*sum(X[,i]*X[,j] * Ai)
      }
    }
    for(i in 1:(length(beta)-1)){
      for(j in (1+i):length(beta) ){
        ret_beta_all[i,j]  = ret_beta_all[j,i]
      }
    }
    ret_beta_sigma <-matrix(0,1,length(beta))
    Us = dloglik_aft_eval(par,t,d,X,type)
    for(j in 1:length(beta)){
      ret_beta_sigma[,j] = 1/sigma^2*sum(X[,j]*scaled_xbeta * Ai)+1/sigma*Us[j]
    }
    ret_sigma = 1/sigma^2*sum(scaled_xbeta*scaled_xbeta*Ai+d)+2*1/sigma*Us[length(par)]
    res = cbind(rbind(ret_beta_all,ret_beta_sigma),c(ret_beta_sigma,ret_sigma))


  }else if(type =="gengamma"){
    res = hessian(loglik_for_gengamma,x=par,t=t,d=d,X=X)
  }
  return(res)
}
#' Compare the generated survival curve vs the KM survival curve
#'
#' @param data the complete data
#' @param res the results from cafta
#' @param type  The type of the error distribution
#' @param j the specific site you want to plot
#'
#' @return the plots
#' @export
#' @import ggplot2
#' @import simsurv
#' @examples
plot_survival_curve_for_renew<-function(data,res,j,type = "weibull"){
  kp<-survfit(Surv(time, status)~1,data )
  time_points <- kp$time
    covariate_names = c("x1", "x2", "x31","x42","x43","x44")
  avg_data = as.data.frame(data[1,covariate_names])
  avg_data[1,]=(colMeans(data[,covariate_names]))

  scale <- exp(res[c(1:(length(covariate_names)+1)),1]%*% as.numeric(cbind(1,avg_data[1,])))[1,1]
  shape <- 1/res[dim(res)[1],1]
  if(type == "weibull"){
    predicted_survival <- exp(-(time_points / scale) ^ shape)

  }else if (type == "loglogistic"){
    predicted_survival <- 1/(1+(time_points/ (scale))^(shape))

  }else if (type == "lognormal"){
    predicted_survival <-  plnorm(time_points, meanlog = log(scale), sdlog = 1/(shape), lower.tail = FALSE)

  }
  #exp(-(time_points / scale) ^ shape)
  km_df <- data.frame(time = kp$time, survival = kp$surv)

  # Create a data frame for the predicted survival
  pred_df <- data.frame(time = time_points, survival = predicted_survival)

  # Plot both on the same graph
  a = ggplot() +
    geom_step(data = km_df, aes(x = time, y = survival, col = "KM")) +
    geom_line(data = pred_df, aes(x = time, y = survival, col = "Predicted")) +
    labs(title =paste0("KM vs. Predicted Survival"," ", type,"Site ", j ), x = "Time", y = "Survival Probability") +
    theme_minimal() +
    scale_colour_manual("",
                        breaks = c("KM", "Predicted"),
                        values = c("blue", "red"))
  km_df <- km_df %>% arrange(time)
  pred_df <- pred_df %>% arrange(time)

  combined_df <- merge(km_df, pred_df, by = "time", suffixes = c("_km", "_pred"))

  combined_df$difference <- with(combined_df, survival_km - survival_pred)

  b <- ggplot(combined_df, aes(x = time, y = difference)) +
    geom_line(colour = "green") +  # Use green or any color of preference
    labs(title = "Difference in Survival Probabilities",
         x = "Time",
         y = paste0("Difference in Survival Probability (KM - Predicted)"," Site ", j)) +
    theme_minimal()

  b = b
  list(a = a,b=b)
}
