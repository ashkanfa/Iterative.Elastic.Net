library("Rcpp")
library("glmnet")
library("ggplot2")
library("dplyr")
library("coefplot")
library("rje")
library("rtkore")



#non-iterative adaptive penalized
Adapen<-function(data,lambda.delta=0.01,ols.adaptive= T, alpha1=0.05,alpha2=0.5){
    # first run
    data<-as.data.frame(data)
    x<-as.matrix(data[,-dim(data)[2]])
    
    if (ols.adaptive){
        
        ols_0<-lm(y~.-1,data=data)
        ols_0.coef<-abs(ols_0[["coefficients"]])
        pf_ols_0  <- 1/ols_0.coef
    }
    
    else{
        
        ols_0 <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha1, nlambda = 10000,
                           intercept=FALSE)
        ols_0 <- cv.glmnet(x=x, y=as.matrix(data[,"y"]), family="gaussian", offset=NULL, alpha = alpha1, lambda=rev(seq(0,max(ols_0[["lambda"]]),by=lambda.delta)),
                           intercept=FALSE)
        selected.coef<-abs(selected.feature(model=ols_0,lambda.delta=lambda.delta,keep.zero = T))
        selected.coef[selected.coef==0]=0.01
        pf_ols_0 <-1/selected.coef
    }
    #########################################################################################################################################################################
    
    cv.elas.ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha2, nlambda = 10000,
                             intercept=FALSE, penalty.factor = pf_ols_0)
    
    cv.elas.ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha2, 
                             lambda=rev(seq(0,max(cv.elas.ols[["lambda"]]),by=lambda.delta)),
                             intercept=FALSE, penalty.factor = pf_ols_0)
    
    m<- selected.feature(model=cv.elas.ols,lambda.delta=lambda.delta)
    ml<-length(m)
    n<- selected.feature(model=cv.elas.ols,lambda.delta=lambda.delta, model.selection = "lambda.1se")
    nl<-length(n)
    
    return(list("model.summary"=cv.elas.ols,"features.lambda.min"=m,"features.lambda.1se"=n ))
}


