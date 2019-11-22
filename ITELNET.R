#rm(list = ls())
#setwd("c:/users/RandKID/OneDrive/Research/Proposal/Generating Simulated Data/Generated_data")

###################################################################################################
## install and load the GGLASSO with weights from my repository
## !! the available package WLS is not working
#library("devtools")
#library("githubinstall")
#install_github("ashkanfa/gglasso",ref="weight",force = TRUE)
#library("gglasso")
library("Rcpp")
library("glmnet")
library("ggplot2")
library("dplyr")
library("coefplot")
library("rje")
library("rtkore")

###################################################################################################
###################################################################################################

Itelnet<-function(data,lambda.delta=0.01,ols.adaptive= T, alpha1=0.05,alpha2=0.5,stepp="min",l=1){
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
    #m=cv.elas.ols[["lambda.1se"]]
    #n=cv.elas.ols[["lambda.min"]]
    #########################################################################################################################################################################
    # Check if the number of feature selected by lambda.min equals to lambda.1se, if not update with the sparser model and do over until the condition satisfies
    while (abs(ml-nl)>l){
        #beta.cv.elas.ols<-extractPath(cv.elas.ols)
        #beta.cv.elas.ols<-as.matrix(beta.cv.elas.ols)
        #beta.cv.elas.ols<-t(beta.cv.elas.ols)
        #beta.cv.elas.ols<-beta.cv.elas.ols[-1,]
        #colnames(beta.cv.elas.ols)<-seq(0,max(cv.elas.ols[["lambda"]]),by=lambda.delta)
        #CV.selected.beta<- beta.cv.elas.ols[,as.character(cv.elas.ols[["lambda.min"]])]
        #CV.selected.beta<- beta.cv.elas.ols[,as.character(cv.elas.ols[["lambda.1se"]])] #big change I made from lambda.min to lambda.1se
        #CV.selected.beta<-CV.selected.beta[CV.selected.beta!=0]
        #CV.selected.beta<-selected.feature(model=cv.elas.ols,lambda.delta=lambda.delta, model.selection = "lambda.1se")
        
        if(stepp =="1se"){
            x<-as.matrix(data[,names(n)])
        } else{
            x<-as.matrix(data[,names(m)])
        }
        
        data<-as.data.frame(cbind(x,data["y"]))
        
        
        #########################################################################################################################################################################
        ## Updating the adaptive weights with reduced model
        if (ols.adaptive){
            
            ols<-lm(y~.-1,data = data)
            ols.coef<-abs(ols[["coefficients"]])
            pf_ols  <- 1/ols.coef
            #ols.coef<-extractPath(ols)
            #ols.coef<-as.matrix(ols.coef)
            #ols.coef<-t(ols.coef)
            #ols.coef<-ols.coef[-1,]
            #colnames(ols.coef)<-seq(0,max(ols[["lambda"]]),by=lambda.delta)
            #selected.coef<- abs(ols.coef[,as.character(ols[["lambda.min"]])])
        }
        else{
            
            ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha1, nlambda = 10000,
                             intercept=FALSE)
            ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha1, 
                             lambda=rev(seq(0,max(ols[["lambda"]]),by=lambda.delta)),intercept=FALSE)
            
            selected.coef<-abs(selected.feature(model=ols, lambda.delta = lambda.delta,keep.zero = T))
            selected.coef[selected.coef==0]=0.01
            pf_ols <-1/selected.coef
            
        }
        
        
        cv.elas.ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha2, nlambda = 10000,
                                 intercept=FALSE, penalty.factor = pf_ols)
        cv.elas.ols <- cv.glmnet(x=x, y=as.matrix(data["y"]), family="gaussian", offset=NULL, alpha = alpha2, lambda=rev(seq(0,max(cv.elas.ols[["lambda"]]),by=lambda.delta)),
                                 intercept=FALSE, penalty.factor = pf_ols)
        #m=cv.elas.ols[["lambda.1se"]]
        #n=cv.elas.ols[["lambda.min"]]
        m<- selected.feature(model=cv.elas.ols,lambda.delta = lambda.delta)
        ml<-length(m)
        n<- selected.feature(model=cv.elas.ols, lambda.delta = lambda.delta, model.selection = "lambda.1se")
        nl<-length(n)
     }
    #z<-selected.feature(model=cv.elas.ols,lambda.delta = 0.01)
    return(list("model.summary"=cv.elas.ols,"features.lambda.min"=m,"features.lambda.1se"=n ))
}

###################################################################################################
###################################################################################################
