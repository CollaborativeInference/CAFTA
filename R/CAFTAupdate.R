#' CAFTA update function for distributed datasets
#'
#' @param B number of sites/subdatas
#' @param init initial beta estimates
#' @param p the total number of covaraites
#' @param robust whether to use robust standard error or not
#'
#' @return
#' @export
#'
#' @examples
CAFTA_update_optimx <-
  function(B, data,type, init=NA, p,robust = TRUE, intercept=TRUE){
    sum2<-diag(0,p,p);
    V_matrix<-diag(0,p,p)
    s<-0
    phi<-0
    tol=1e-8;
    max_iter=10000;
    negl<-0
    if(sum(is.na(init))>1){init<-rep(0,p)}
    betahat<-init;

    for(b in 1:B){


#      t<-d<-X<-NULL;
#      load(paste(tempdatadir,"/Simdata",b,".RData",sep=""))
      subdata= subset(data,group==b)
      t <- subdata$time
      d <- subdata$status
      X = as.matrix(subdata[,4:(4+npar-1)])

      if(intercept==TRUE){X<-as.matrix(cbind(1,X))}



      betahat_old=betahat;
      #record W in a list form rather than a matrix to simplify calculation
      #W<-invlinkdiv(X,betahat_old,type=type)
      dloglik_aft_eval_1<-function(par,t,d,X,type,betahat_old,sum2,individual=FALSE){#par,n_gamma,degree,t,d,X,knots
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
          add_term  = ((betahat_old-par)%*%sum2)
          hsigma =add_term[length(par)]
          hbeta= add_term[1:(length(par)-1)]
          res = c(U_beta+hbeta,U_sigma+hsigma)
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

        return(-res)
      }
      ddloglik_aft_eval_1<-function(par,t,d,X,type,betahat_old,sum2){#par,n_gamma,degree,t,d,X,knots
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
      loglik_restricted<-function(betahat,t,d,X,type,betahat_old,sum2){
        loglik_eval(betahat,t,d,X,type)+1/2*((betahat_old-betahat)%*%sum2%*%(betahat_old-betahat))[1,1]
      }
      if(type =="weibull"|type=="loglogistic"|type=="lognormal"){
        optimxres = optimx::optimx(betahat,loglik_restricted,gr=dloglik_aft_eval_1, t=t,d=d,X=X,type = type,betahat_old = betahat_old, sum2=sum2,method="BFGS",hessian   = TRUE,control = list(maxit = 10000,reltol=1e-12 ))

        betahat = as.numeric(optimxres[,1:p])

      }else{

        optimxres = optimx::optimx(betahat,loglik_restricted,t=t,d=d,X=X,type = type,betahat_old = betahat_old, sum2=sum2,method="BFGS",hessian   = TRUE,control = list(maxit = 10000,reltol=1e-12 ))

        betahat = as.numeric(optimxres[,1:p])
      }
      if(optimxres$convcode!=0){
        stop("optimx does not converge")
      }
      H_new<-ddloglik_aft_eval(betahat,t,d,X,type) #J2
      if(matrixcalc::is.positive.semi.definite(H_new)){
        H_new = H_new
      }else{
        H_new = as.matrix(Matrix::nearPD(Matrix::Matrix(H_new), ensureSymmetry = TRUE, )$mat)

      }


      sum2<-sum2+H_new
      meat.half=dloglik_aft_eval(betahat,t,d,X,type,TRUE)
      sum_outer_products = as.matrix(t(meat.half)) %*% as.matrix(meat.half)
      V_matrix<-V_matrix+sum_outer_products
      negl<-negl+loglik_eval(betahat,t,d,X,type)
    }


    sd <-sqrt(diag(t(solve(sum2))%*%V_matrix%*%(solve(sum2))))
    if (!robust){
      sd<-sqrt(diag(solve(sum2)))
    }
    pvalue<-2*pnorm(-abs(betahat)/sd)
    result<-cbind(betahat=betahat,sd=sd,pvalue=pvalue,negll=negl)#
    colnames(result)<-c("Estimates","Std.Errors","p-values","neg-logll")#

    return(result)

  }
