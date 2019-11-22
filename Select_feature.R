###################################################################################################
###################################################################################################

# define a function to extract names of the feature selected
selected.feature<-function(model,lambda.delta=0.01,model.selection="lambda.min",keep.zero=F){
    beta.model<-extractPath(model)
    beta.model<-as.matrix(beta.model)
    beta.model<-t(beta.model)
    beta.model<-beta.model[-1,]
    colnames(beta.model)<-seq(0,max(model[["lambda"]]),by=lambda.delta)
    selected.beta<- beta.model[,as.character(model[[model.selection]])]
    if (keep.zero==F){
        selected.beta<- selected.beta[selected.beta!=0]
    }
    else 
        return(selected.beta)
}

