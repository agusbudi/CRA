rm(list=ls(all=TRUE))
#expert+spammers = 100 annotators
filename= "diabetes"

#from 0 to 5 relative experts
for(expnum in 0:5){
  
  ptm2=NULL
  
  info = NULL
  ptm2 <- proc.time()
  
  folder= paste(filename,expnum, sep="/")
  
  #matrice of truth label : Y
  rawdata=read.table(paste(filename,"/",filename,".txt",sep=""),sep=",",header=FALSE)
  rawdata=as.matrix(rawdata)
  Y = rawdata[,ncol(rawdata)]
  
  write.table(rawdata[,1:(ncol(rawdata)-1)],paste(folder,"/X-vali.txt", sep=""),row.names=FALSE,col.names=FALSE)
  write.table(Y,paste(folder,"/Yb-vali.txt", sep=""),row.names=FALSE,col.names=FALSE)
  
  X.vali=rawdata[,1:(ncol(rawdata)-1)]
  X.vali=as.matrix(X.vali)
  Yb.vali=Y
  
  #define the matrix of expert and spammer
  expt_num = expnum  #correct annotators >> sensitivity and specificity >= 0.6 (these wolley, 130)
  spam_num = 100-expt_num  #spam annotators >> randomly
  
  expt.A =matrix(0,nrow=length(Y),ncol=expt_num)
  spam.A =matrix(0,nrow=length(Y),ncol=spam_num)
  
  spam.A[,1] <- sample(x = c(0, 1), size = length(Y), replace = TRUE)
  
  
  for (t in 1:expt_num){
    #fill with the truth label as the basic prediction for expert annotators
    expt.A[,t]=Y
    
    #define the error, between 0-40%
    err <- sample(1:40, 1) #in %
    err.idx = sample(which(Y==1), round(err*length(Y[which(Y==1)])/100,0))
    expt.A[err.idx,t]=0
    
    err <- sample(1:40, 1) #in %
    err.idx = sample(which(Y==0), round(err*length(Y[which(Y==0)])/100,0))
    expt.A[err.idx,t]=1    
    
  }  
  
  #fill with the random prediction spammer annotators
  t=1
  for (t in 1:spam_num){
    #fill with the truth label as the basic prediction for spammer
    spam.A[,t]=Y
    #define random error for sensitivity
    err <- sample(20:80, 1) #in %
    err.idx = sample(which(Y==1), round(err*length(Y[which(Y==1)])/100,0))
    spam.A[err.idx,t]=0

    #define random error for specificity based on previous error
    error.range = 0
    err.specificity =0
    while(error.range>120 || error.range <80 ||  err.specificity>80 ||err.specificity<20  ){
      err.specificity <- sample(0:(120-err), 1) #in %
      error.range= err.specificity + err
    }
    err = err.specificity 
    err.idx = sample(which(Y==0), round(err*length(Y[which(Y==0)])/100,0))
    spam.A[err.idx,t]=1
    
  }
  
  library(PresenceAbsence)
  
  all.A = cbind(expt.A,spam.A)
  
  entropy_input <- all.A[,1]
  write.table(entropy_input,paste(folder,"/entropy_input.txt", sep=""),row.names=FALSE,col.names=FALSE)
  
  evaluation = matrix(0,nrow=2,ncol=(expt_num+spam_num)) #auc and accuracy
  for (t in 1:(expt_num+spam_num)){
    evaluation[1,t] = auc(cbind(1:length(Y),Y, all.A[,t]),st.dev=FALSE)
    for (i in 1:length(Y)){
      if(all.A[i,t]==Y[i])    
        evaluation[2,t]= evaluation[2,t] + 1
    }
    evaluation[2,t]= evaluation[2,t]/length(Y)
  }
  
  Sens = matrix(data = NA, nrow = (expt_num+spam_num), ncol=1)
  Spec = matrix(data = NA, nrow = (expt_num+spam_num), ncol=1)
  t=1
  for (t in 1:(expt_num+spam_num)){
    xtab=table(all.A[,t],Y)
    Sens[t] = xtab[2,2]/(xtab[1,2]+xtab[2,2])
    Spec[t] = xtab[1,1]/(xtab[2,1]+xtab[1,1])
  }
  result2 = cbind((1-Spec),Sens)

##print the generated data    
  pdf(paste(folder,"SpammerScore.pdf",sep="/"),width=5,height=5,paper='special')
  header = paste(expt_num, " experts and ", spam_num," spammers", " (",filename," dataset)", sep = "")
  plot(1, type="n", xlab="1-Specificity", ylab="Sensitivity", xlim=c(0, 1), ylim=c(0, 1))
  lines(x = c(0,1), y = c(0,1),col="green")
  points(result2[1:expt_num,],col="blue",pch=1)
  points(result2[(expt_num+1):100,],col="red",pch=2)
  text(0.05, 0.05, "biased",cex = 0.8)
  text(0.95, 0.95, "biased",cex = 0.8)
  text(0.05, 0.95, "experts",cex = 0.8)
  text(0.5, 0.5, "spammers",cex = 0.8)
  text(0.95, 0.05, "malicious",cex = 0.8)
  dev.off()

###########################
##define the entropy distribution reference by applying EM-Bayesian, see step 1  
  Y=all.A
  Y=as.matrix(Y)
  colnames(Y)=seq(1,ncol(Y),1)
  
  annot.exp=matrix(0,nrow=ncol(Y),ncol=ncol(Y))
  annot=1:ncol(Y)
  Y.corr=matrix(0,nrow=2,ncol=ncol(Y))
  Y.corr[1,]=annot
  Y.corr[2,]=1:ncol(Y.corr)
  colnames(Y)=seq(1,ncol(Y),1)
  
  good.a=1:ncol(Y)
  
  library(polynom) 
  
  # Initialisation de mu
  R = ncol(Y)  #the number of annotators
  mu = NULL
  for (i in (1:nrow(X.vali)))
  {
    mu[i] = (1/R)*sum(Y[i,])    #Algorithm 9. Initialization of \mu value, Wolley's dissertation page 131, line 6
  }
  
  eta  = 1
  g = matrix(0,nrow=nrow(X.vali),ncol=ncol(X.vali))
  a = NULL
  b = NULL
  p = NULL
  diff1= 10
  H=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
  
  #calculate the Alpha1 and Beta1 :
  alpha = matrix(0,nrow=R,ncol=10000)
  beta = matrix(0,nrow=R,ncol=10000)
  c=1   #
  library(MASS)
  
  for (j in 1:R)
  {
    alpha[j,c] = (sum(mu*Y[,j]))/(sum(mu))
    beta[j,c] = (sum((1-mu)*(1-Y[,j])))/(sum(1-mu))
  }  
  alpha[,1]
  beta[,1]
  
  
  #### Calculate w, the parameter that has to be estimated by using Newton-Raphson method ######
  y2=NULL
  for ( i in 1:nrow(X.vali))
  {
    if (mu[i] < 0.5){
      y2[i] = 0
    }else {
      y2[i]=1
    }
  }
  w = matrix(0,nrow=ncol(X.vali),ncol=10000)
  w[,c] = ginv(X.vali)%*%y2
  
  
  # Calculate p, the probability fo giving true value, z=1
  p = 1/(1+exp(-(t(w[,c])%*%t(X.vali))))
  length(p)
  
  #Calculate g, gradient vector:
  for (i in 1:nrow(X.vali)){
    g[i,] = (mu[i]-p[i])*X.vali[i,]
  }
  
  g_new=NULL
  g_new=colSums(g)
  
  # Calcul H^-1, H is a Hessian Matrix:
  H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
  dim(H)
  for (i in 1:nrow(X.vali)){
    H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
  }
  
  H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
  H_new=H[1,,]
  for (i in 2:nrow(X.vali)){
    H_new=H_new+H[i,,]
  }
  H_new=-H_new
  dim(H_new)
  
  # Calculate the new w
  w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new
  
  # Calculate the differences of previous parameters with the current one
  diff1=sum(sqrt((w[,c+1]-w[,c])^2))
  diff1
  
########### E-STEP ############
  # Calculate a and b: a is the probability of true positive, and b is the probability of true negative
  ##set the convergence threshold
  convergence.Limit = 0.005
  while (diff1 > convergence.Limit)
  {
    for (i in 1:nrow(X.vali))
    {
      a[i] = prod((alpha[,c]^(Y[i,]))*((1-alpha[,c])^(1-Y[i,])))
      b[i] = prod((beta[,c]^(1-Y[i,]))*((1-beta[,c])^(Y[i,])))
      
    }
    
    for (i in 1:nrow(X.vali)){
      if (a[i]==0) a[i]=10^{-10}
      if (b[i]==0) b[i]=10^{-10}
    }
    
    # Calculate p:
    p = 1/(1+exp(-(t(w[,c+1])%*%t(X.vali))))
    
    # Calculate the new mu
    for (i in 1:nrow(X.vali))
    {
      mu[i] = (a[i]*p[i])/(a[i]*p[i]+b[i]*(1-p[i]))
    }
    
    #iteration to store new alpha and beta values
    c=c+1
    
    for (j in 1:R)
    {
      alpha[j,c] = (sum(mu*Y[,j]))/(sum(mu))
      beta[j,c] = (sum((1-mu)*(1-Y[,j])))/(sum(1-mu))
    }
    
    alpha[,c]
    beta[,c]
    
    #### Calculate w ######
    
    # Calculate g:
    for (i in 1:nrow(X.vali))
    {
      g[i,] = (mu[i]-p[i])*X.vali[i,]
    }
    
    g_new=NULL
    g_new=colSums(g)
    
    # Calcul de H^-1:
    H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
    dim(H)
    for (i in 1:nrow(X.vali)){
      H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
    }
    
    H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
    H_new=H[1,,]
    for (i in 2:nrow(X.vali)){
      H_new=H_new+H[i,,]
    }
    H_new=-H_new
    
    # Calculate w
    w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new
    
    #calculate the difference
    diff1=sum(sqrt((w[,c+1]-w[,c])^2))
    diff1
    
    #if the iteration reach the limit of allocated dimension, then break
    if(c==9999) diff1=0
    
  }
  c
  ##################################################
  if(is.nan(beta[,c]))
    c=c-1
  ##################################################
  alpha[,c]
  beta[,c]
  w[,c]
  
  
  ####Calculate the score, higher is better
  S=NULL
  for (j in 1:ncol(Y)) S[j]=(alpha[j,c]+beta[j,c]-1)^2
  S
  
  c3=1
  annot.exp[c3,1:length(good.a)]=good.a
  c3=c3+1
  
  #select only the highest value
  info= paste(info,"  EM: ",which(S==max(S)), sep="", collapse =NULL)
  #entropy reference is selected
  entropy_input <- Y[,which(S==max(S))]
  write.table(entropy_input,paste(folder,"/entropy_input_EM.txt", sep=""),row.names=FALSE,col.names=FALSE)
  
##EM is done##
  
######select entropy distribution by considering Majority voting
  R = ncol(Y)
  mu = NULL
  for (i in (1:length(Yb.vali)))
  {
    mu[i] = (1/R)*sum(Y[i,])
  }
  mu
  
  mu2 = as.numeric(mu)
  for (i in 1:nrow(Y)){
    if(mu2[i]>0.5)
      mu2[i]=1
    else
      mu2[i]=0
  }
  
  #entropy is selected
  entropy_input <- mu2
  write.table(entropy_input,paste(folder,"/entropy_input_MV.txt", sep=""),row.names=FALSE,col.names=FALSE)
  time.result= proc.time() - ptm2

####select entropy distribution Randomly
  Random.Entropy = sample(1:ncol(Y), 1)
  info= paste(info,"  Random: ",Random.Entropy, "  time: ",time.result["elapsed"], sep="", collapse =NULL)
  entropy_input <- Y[,Random.Entropy]
  write.table(entropy_input,paste(folder,"/entropy_input_Random.txt", sep=""),row.names=FALSE,col.names=FALSE)
  
####write the results to the files:
  write.table(info,paste(folder,"/info_entropy.txt", sep=""),row.names=FALSE,col.names=FALSE)
  write.table(all.A,paste(folder,"/Y.txt", sep=""),row.names=FALSE,col.names=FALSE)
  write.table(expt.A,paste(folder,"/expert.txt", sep=""),row.names=FALSE,col.names=FALSE)
  write.table(spam.A,paste(folder,"/spam.txt", sep=""),row.names=FALSE,col.names=FALSE)
  write.table(evaluation,paste(folder,"/eval_auc_accuracy.txt", sep=""),row.names=FALSE,col.names=FALSE)
}