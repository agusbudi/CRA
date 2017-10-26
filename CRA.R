##preparation, define the filename
rm(list=ls(all=TRUE))

#select the number of reliable annotator exist
expnum=2
filename= "diabetes"
folder = paste(filename,expnum,sep="/")
convergence.Limit = 0.001

ptm2=NULL

nb.S=1
ptm2 <- proc.time()
X.vali=read.table(paste(folder,"/X-vali.txt", sep=""),header=FALSE)
X.vali=as.matrix(X.vali)
Yb.vali=read.table(paste(folder,"/Yb-vali.txt", sep=""),header=FALSE)
Y.entropy = as.matrix(read.table(paste(folder,"/entropy_input_EM.txt", sep=""),header=FALSE))
Yb.vali=Yb.vali[,1]

Yname= paste(folder,"/Y.txt", sep="")
Y=read.table(Yname,header=FALSE)
Y=as.matrix(Y)
colnames(Y)=seq(1,ncol(Y),1)
Y[1,]

##define the spammer score threshold, see (Raykar, 2012)
threshold = 0.6
epsilon.S = ((1- threshold) + (1- threshold) -1)^2

##### calculate entropy: see (Zighed et al, 2010)
f=matrix(0,nrow=2,ncol=ncol(Y))
lamb=matrix(0,nrow=2,ncol=ncol(Y))
w=NULL
res.ent=NULL

for (t in 1:ncol(Y)){
  for (i in 1:2){
    w[i]=length(which(Y.entropy==(i-1)))/length(Yb.vali)

    f[i,t]=length(which(Y[,t]==(i-1)))/length(Yb.vali) #frequency of class 0 and 1 prediction for each annotator
    lamb[i,t]=(nrow(X.vali)*f[i,t]+1)/(nrow(X.vali)+2) #define lambda, 6.13
  }
  res.ent[t]=sum((lamb[,t]*(1-lamb[,t]))/((-2*w+1)*lamb[,t]+w^2)) #calculate entropy, 6.12
}


##################################################################################################
##choose with the default threshold (0.6 of sensitivity and specificity based on the reference distribution)
thres=NULL
for (i in 1:2){
  thres[i]=threshold*length(which(Yb.vali==(i-1)))
  f[i,t]=thres[i]/length(Yb.vali) #frequency of class 0 and 1 prediction for each annotator
  lamb[i,t]=(nrow(X.vali)*f[i,t]+1)/(nrow(X.vali)+2) #define lambda, 6.13
}
thres.res.ent=sum((lamb[,t]*(1-lamb[,t]))/((-2*w+1)*lamb[,t]+w^2)) #calculate entropy, 6.12
nb.annot = length(res.ent[which(res.ent>=thres.res.ent)])


###recent editing, in order to limit the group number, if top-K has >50% of annotators, select only top 30%
if(nb.annot > 0.5*ncol(Y))
  nb.annot=0.3*nb.annot

##choose the 30% from all annotators based on the entropy, higher is better
ordre.ent=order(res.ent, na.last = TRUE, decreasing = TRUE)
annot=ordre.ent[1:nb.annot]

#print the EC order based on S rank
ordre.ent
res.ent[ordre.ent]
annot

annot.exp=matrix(0,nrow=ncol(Y),ncol=ncol(Y))

#initiate Y with only the selected annotators
Y=Y[,annot]
Y.corr=matrix(0,nrow=2,ncol=ncol(Y))
Y.corr[1,]=annot
Y.corr[2,]=1:ncol(Y.corr)
colnames(Y)=seq(1,ncol(Y),1)

good.a=1:ncol(Y) #only the order

Y[1:2,]

library(polynom)

# Initialization of mu
R = ncol(Y)  #the number of annotators
mu = NULL
for (i in (1:nrow(X.vali)))
{
  mu[i] = (1/R)*sum(Y[i,])    #Algorithm 9. Initialization of \mu value, Wolley's dissertation page 131, line 6
}
mu

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

# Calculate the difference
diff1=sum(sqrt((w[,c+1]-w[,c])^2))
diff1


########### E-STEP ###################
# Calculate a and b: a is the probability of true positive, and b is the probability of true negative

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

  a
  b

  # Calcul de p:
  p = 1/(1+exp(-(t(w[,c+1])%*%t(X.vali))))

  p

  # Calcul du nouveau mu
  for (i in 1:nrow(X.vali))
  {
    mu[i] = (a[i]*p[i])/(a[i]*p[i]+b[i]*(1-p[i]))
  }

  mu

  c=c+1
  for (j in 1:R)
  {
    alpha[j,c] = (sum(mu*Y[,j]))/(sum(mu))
    beta[j,c] = (sum((1-mu)*(1-Y[,j])))/(sum(1-mu))
  }

  alpha[,c]
  beta[,c]

  #### Calculate w:

  # Calculate g:
  for (i in 1:nrow(X.vali))
  {
    g[i,] = (mu[i]-p[i])*X.vali[i,]
  }

  g_new=NULL
  g_new=colSums(g)

  # Calculate H^-1:
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

  # Calculate the new w
  w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new

  # Calculate the difference
  diff1=sum(sqrt((w[,c+1]-w[,c])^2))
  diff1
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
length(S)

c3=1
annot.exp[c3,1:length(good.a)]=good.a
c3=c3+1

#inspired from Condorcet Jury's theorem
if(length(which(S>epsilon.S)) > 0){
  annot.exp[c3,1:length(which(S>epsilon.S))]=which(S>epsilon.S)
  annot.exp[c3,1:length(which(S>epsilon.S))]=Y.corr[1,which(S>epsilon.S)]
}else{
  ordre.S=order(S, na.last = TRUE, decreasing = TRUE)
  temp.epsilon.S=S[ordre.S[nb.S]]
  annot.exp[c3,1:length(which(S>=temp.epsilon.S))]=which(S>=temp.epsilon.S)
  annot.exp[c3,1:length(which(S>=temp.epsilon.S))]=Y.corr[1,which(S>=temp.epsilon.S)]
  
#    annot.exp[c3,1:length(which(S==max(S)))]=which(S==max(S))
#    annot.exp[c3,1:length(which(S==max(S)))]=Y.corr[1,which(S==max(S))]
}

###iterate until convergent
while (sum((annot.exp[c3,]-annot.exp[c3-1,]))!=0){
  ##merge previous potential reliable annotator with the new group
  annot=ordre.ent[(((c3-1)*(nb.annot)+1):(c3*nb.annot))]
  annot=na.omit(annot)


  ###re-initiate the annotators only with the selected group

  Y=read.table(Yname,header=FALSE)
  Y=as.matrix(Y)
  colnames(Y)=seq(1,ncol(Y),1)
  Y[1,]

  Y=Y[,c(annot,annot.exp[c3,(1:length(which(S>epsilon.S)))])]

  #if there is still available group
  if(!is.null(ncol(Y))){
    Y.corr=matrix(0,nrow=2,ncol=ncol(Y))
    Y.corr[1,]=c(annot,annot.exp[c3,(1:length(which(S>epsilon.S)))])
    Y.corr[2,]=1:ncol(Y.corr)
    colnames(Y)=seq(1,ncol(Y),1)
    Y[1:2,]


    library(polynom)

    good.a=1:ncol(Y)


    # Initialization of mu:
    R = ncol(Y)
    mu = NULL
    for (i in (1:nrow(X.vali)))
    {
      mu[i] = (1/R)*sum(Y[i,])
    }
    mu

    eta  = 1
    g = matrix(0,nrow=nrow(X.vali),ncol=ncol(X.vali))
    a = NULL
    b = NULL
    p = NULL
    diff1= 10
    H=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))


    # Calculate alpha1 & beta1, alpha1 is P(y=1|z=1), beta1 is P(y=0|z=0)
    alpha = matrix(0,nrow=R,ncol=10000)
    beta = matrix(0,nrow=R,ncol=10000)
    c=1
    library(MASS)

    for (j in 1:R)
    {
      alpha[j,c] = (sum(mu*Y[,j]))/(sum(mu))
      beta[j,c] = (sum((1-mu)*(1-Y[,j])))/(sum(1-mu))
    }

    alpha[,1]
    beta[,1]


    #Calculate w
    y2=NULL
    for ( i in 1:nrow(X.vali))
    {
      if (mu[i] < 0.5)
      {
        y2[i] = 0
      }else
      {
        y2[i]=1
      }
    }
    w = matrix(0,nrow=ncol(X.vali),ncol=30000)
    w[,c] = ginv(X.vali)%*%y2


    # Calculate p:
    p = 1/(1+exp(-(t(w[,c])%*%t(X.vali))))
    length(p)
    #Calculate g:
    for (i in 1:nrow(X.vali))
    {
      g[i,] = (mu[i]-p[i])*X.vali[i,]
    }

    g_new=NULL
    g_new=colSums(g)

    # Calculate H^-1:
    H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
    dim(H)
    for (i in 1:nrow(X.vali))
    {
      H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
    }

    H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
    H_new=H[1,,]
    for (i in 2:nrow(X.vali))
    {
      H_new=H_new+H[i,,]
    }
    H_new=-H_new
    dim(H_new)

    # Calculate the new w
    w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new

    # Calculate the difference
    diff1=sum(sqrt((w[,c+1]-w[,c])^2))
    diff1

    ########### E-STEP ###################
    # Calculate a & b:
    while (diff1 > convergence.Limit){
      for (i in 1:nrow(X.vali))
      {
        a[i] = prod((alpha[,c]^(Y[i,]))*((1-alpha[,c])^(1-Y[i,])))
        b[i] = prod((beta[,c]^(1-Y[i,]))*((1-beta[,c])^(Y[i,])))
      }

      for (i in 1:nrow(X.vali)){
        if (a[i]==0) a[i]=10^{-1}
        if (b[i]==0) b[i]=10^{-1}
      }

      # Calculate p:
      p = 1/(1+exp(-(t(w[,c+1])%*%t(X.vali))))

      # Calculate the new mu
      for (i in 1:nrow(X.vali)){
        mu[i] = (a[i]*p[i])/(a[i]*p[i]+b[i]*(1-p[i]))
      }

      c=c+1
      for (j in 1:R)
      {
        alpha[j,c] = (sum(mu*Y[,j]))/(sum(mu))
        beta[j,c] = (sum((1-mu)*(1-Y[,j])))/(sum(1-mu))
      }

      alpha[,c]
      beta[,c]

      # Calculate g:
      for (i in 1:nrow(X.vali))
      {
        g[i,] = (mu[i]-p[i])*X.vali[i,]
      }

      g_new=NULL
      g_new=colSums(g)

      # Calculate H^-1:
      H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
      dim(H)
      for (i in 1:nrow(X.vali))
      {
        H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
      }

      H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
      H_new=H[1,,]
      for (i in 2:nrow(X.vali))
      {
        H_new=H_new+H[i,,]
      }
      H_new=-H_new

      # Calculate w
      w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new

      # Calculate the difference
      diff1=sum(sqrt((w[,c+1]-w[,c])^2))
      if(is.nan(diff1))
        diff1=0
      
      ##if the iteration reach the limit, the break
      if(c==9999) diff1=0
    }
    c
    alpha[,c]
    beta[,c]
    w[,c+1]


    ####calculate Spammer score, higher is better
    S=NULL
    for (j in 1:ncol(Y)) S[j]=(alpha[j,c]+beta[j,c]-1)^2
    
  ##if there is no more group
  } else{
    Y.corr=matrix(0,nrow=2,ncol=1)
    Y.corr[1,]=c(annot,annot.exp[c3,(1:length(which(S>epsilon.S)))])
    Y.corr[2,]=1:1

    library(polynom)

    good.a=1

    # Initialization of mu:
    R = 1
    mu = Y

    eta  = 1
    g = matrix(0,nrow=nrow(X.vali),ncol=ncol(X.vali))
    a = NULL
    b = NULL
    p = NULL
    diff1= 10
    H=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))


    # Calculate alpha1 & beta1:
    alpha = matrix(0,nrow=R,ncol=10000)
    beta = matrix(0,nrow=R,ncol=10000)
    c=1
    library(MASS)

    alpha[1,c] = (sum(mu*Y))/(sum(mu))
    beta[1,c] = (sum((1-mu)*(1-Y)))/(sum(1-mu))

    #Calculate w
    y2=NULL
    for ( i in 1:nrow(X.vali))
    {
      if (mu[i] < 0.5)
      {
        y2[i] = 0
      }else
      {
        y2[i]=1
      }
    }
    w = matrix(0,nrow=ncol(X.vali),ncol=30000)
    w[,c] = ginv(X.vali)%*%y2


    # Calculate p:
    p = 1/(1+exp(-(t(w[,c])%*%t(X.vali))))
    length(p)
    #Calculate g:
    for (i in 1:nrow(X.vali))
    {
      g[i,] = (mu[i]-p[i])*X.vali[i,]
    }

    g_new=NULL
    g_new=colSums(g)

    # Calculate H^-1:
    H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
    dim(H)
    for (i in 1:nrow(X.vali))
    {
      H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
    }

    H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
    H_new=H[1,,]
    for (i in 2:nrow(X.vali))
    {
      H_new=H_new+H[i,,]
    }
    H_new=-H_new
    dim(H_new)

    # Calculate the new w
    w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new

    # Calculate the difference
    diff1=sum(sqrt((w[,c+1]-w[,c])^2))
    diff1

    ########### E-STEP : ###################

    # Calculate a & b:
    while (diff1 > convergence.Limit){
      for (i in 1:nrow(X.vali))
      {
        a[i] = prod((alpha[,c]^(Y[i,]))*((1-alpha[,c])^(1-Y[i,])))
        b[i] = prod((beta[,c]^(1-Y[i,]))*((1-beta[,c])^(Y[i,])))

      }

      for (i in 1:nrow(X.vali)){
        if (a[i]==0) a[i]=10^{-1}
        if (b[i]==0) b[i]=10^{-1}
      }

      # Calculate p:
      p = 1/(1+exp(-(t(w[,c+1])%*%t(X.vali))))
      
      # Calculate the new mu
      for (i in 1:nrow(X.vali))
      {
        mu[i] = (a[i]*p[i])/(a[i]*p[i]+b[i]*(1-p[i]))
      }

      mu

      c=c+1
      alpha[1,c] = (sum(mu*Y))/(sum(mu))

      for (i in 1:nrow(X.vali))
      {
        g[i,] = (mu[i]-p[i])*X.vali[i,]
      }

      g_new=NULL
      g_new=colSums(g)

      H = array(0,dim=c(nrow(X.vali),ncol(X.vali),ncol(X.vali)))
      dim(H)
      for (i in 1:nrow(X.vali))
      {
        H[i,,] = p[i]*(1-p[i])*(X.vali[i,]%*%t(X.vali[i,]))
      }

      H_new=matrix(0,nrow=ncol(X.vali),ncol=ncol(X.vali))
      H_new=H[1,,]
      for (i in 2:nrow(X.vali))
      {
        H_new=H_new+H[i,,]
      }
      H_new=-H_new

      w[,c+1] = w[,c] - eta*ginv(H_new)%*%g_new

      diff1=sum(sqrt((w[,c+1]-w[,c])^2))
      if(is.nan(diff1))
        diff1=0
      if(c==9999) diff1=0
    }
    c
    alpha[,c]
    beta[,c]
    w[,c+1]

    #### Calculate spammer score

    S=NULL
    S=(alpha[1,c]+beta[1,c]-1)^2
  }
  S

  c3=c3+1
  #inspired from Condorcet Jury's theorem
  if(length(which(S>epsilon.S)) > 0){
    annot.exp[c3,1:length(which(S>epsilon.S))]=which(S>epsilon.S)
    annot.exp[c3,1:length(which(S>epsilon.S))]=Y.corr[1,which(S>epsilon.S)]
  }else{
  #  annot.exp[c3,1:length(which(S==max(S)))]=which(S==max(S))
  #  annot.exp[c3,1:length(which(S==max(S)))]=Y.corr[1,which(S==max(S))]
    
    ordre.S=order(S, na.last = TRUE, decreasing = TRUE)
    temp.epsilon.S=S[ordre.S[nb.S]]
    annot.exp[c3,1:length(which(S>=temp.epsilon.S))]=which(S>=temp.epsilon.S)
    annot.exp[c3,1:length(which(S>=temp.epsilon.S))]=Y.corr[1,which(S>=temp.epsilon.S)]
  }
}
###################the iteration is finished########


Y=read.table(Yname,header=FALSE)
Y=as.matrix(Y)
colnames(Y)=seq(1,ncol(Y),1)
c3

RC.f = NULL
##if there are several final reliable annotators
if(length(which(S>epsilon.S)) > 1){
  round(S,5)
  #rank:
  Maxconf= max(S)
  Maxaa=annot.exp[c3,which(S==max(S))]
  RC.f = annot.exp[c3,1:length(which(S>epsilon.S))]
  Y=Y[,annot.exp[c3,1:length(which(S>epsilon.S))]]
  
  R = ncol(Y)
  mu = NULL
  for (i in (1:nrow(X.vali)))
  {
    mu[i] = (1/R)*sum(Y[i,])  #define the final prediction by applying Majority Voting
  }
##if no reliable annotator, then select only one  
}else{
  RC.f = annot.exp[c3,1:length(which(S==max(S)))]
  Y=Y[,annot.exp[c3,1:length(which(S==max(S)))]]
  mu <- Y
}

#print the result
mu
write.table(mu,paste(folder,"/mu-CondexpertEM.txt",sep=""),row.names=FALSE,col.names=TRUE)

library(PresenceAbsence)
result2<-matrix(0,nrow=nrow(X.vali),ncol=3)
result2[,1]<-1:nrow(X.vali)
result2[,2]<-as.numeric(Yb.vali)
result2[,3]=mu
result2<-as.data.frame(result2)
dimnames(result2)[[2]]<-c("ID" ,"Observed","ExpertS")
dimnames(result2)[[2]]
is.data.frame(result2)

############################## List of reliable annotators
annot.exp[c3,1:length(which(S>epsilon.S))]

#########################computation time
proc.time() - ptm2
time.result= proc.time() - ptm2

header = paste("ROC curves Using Condorcet-ExpertS for ", filename ," data \n based on 100 Annotators", sep = "")
png(filename=paste(folder,"/ROC_Condorcet-ExpertSEM.png", sep=""))
auc.roc.plot(result2,col=TRUE,line.type=TRUE,main=header)
dev.off()

auc.result = auc(result2, st.dev=FALSE)
acc.result = 0
mu2 = as.numeric(mu)
for (i in 1:nrow(X.vali)){
  pred=0
  if(mu2[i]>0.5)
    pred=1
  
  if(pred==Yb.vali[i])    
    acc.result= acc.result + 1
}
acc.result= acc.result/nrow(X.vali)
expertS_result = paste(time.result["elapsed"], " ",auc.result, " ",acc.result )

write.table(RC.f,paste(folder,"/RC.f.txt", sep=""),row.names=FALSE,col.names=FALSE)

write.table(expertS_result,paste(folder,"/Condorcet-expertSEM_result.txt", sep=""),row.names=FALSE,col.names=FALSE)

###############################################################