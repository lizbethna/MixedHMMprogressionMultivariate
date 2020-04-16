##################################################
### 
### A mixed hidden Markov model for multivariate monotone disease processes in the presence of measurement errors
###  
### Lizbeth Naranjo (1), Emmanuel Lesaffre (2), Carlos J. Pérez (3).
###
### (1) Departamento de Matemáticas, Facultad de Ciencias, Universidad Nacional Autónoma de México, Ciudad de México, Mexico.
### (2) L-BioStat, School of Public Health, KU Leuven, Leuven, Belgium.
### (3) Departamento de Matemáticas, Facultad de Veterinaria, Universidad de Extremadura,  Cáceres, Spain.
### 
##################################################

##################################################
### R packages
###
### Instructions: 
### Load the R libraries
##################################################
library(rjags)
library(MCMCpack)
library(mnormt)
library(coda)
library(rjags)
library(dclone) # To run MCMC in
library(snow)   # several cores

##################################################
### Truncated Normal I[tra < Z]
rnormleft <- function(tra,mu,sig){
	rp <- pnorm(tra, mean=mu, sd=sig)
	u <- rp+(1-rp)*runif(1)
	qnorm(u, mean=mu, sd=sig)
}
##################################################

##################################################
### Instructions: 
### Change the address where the data and codes are located. 
### setwd("HERE")
##################################################
setwd("~/FileDirectory/")
setwd("~/Documents/Articulos/Emmanuel/4Codes/")
setwd("C:/Users/propietario 1/Documents/ArticuloEmmanuel/4Codes/")
getwd()
##################################################

##################################################
### Simulate mimic STdata
##################################################

##################################################
### Data and design matrices 

J = 3   # teeth 
N = 1000   # subject
K = 5   # time points

### p exogenous covariates associated with first examination of the jth tooth of the ith subject
p = 2    
X = array(runif(N*p,0,1),dim=c(N,p))

### q exogenous but possibly time-varying covariates associated with the jth tooth of the ith subject at time point t(i,k)
q = 2
Z = array(runif(N*K*q,0,1),dim=c(N,K,q))

### Q examiners
Q = 4
Xi = array(NA,dim=c(N,K))
for(k in 1:K){
	Xi[,k] = rep(1:Q,length.out=N)
}

##################################################

##################################################
### Parameters

BetaP = as.vector(rep(1.5,p))
GammaP = as.vector(rep(1.5,q)) 
BetaI = as.vector(rep(1.5,p)) 
GammaI = as.vector(rep(1.5,q)) 

Alpha = as.vector(rep(-4,K))

OmegaP = OmegaLP = matrix(1:(J*J),J,J)
OmegaI = OmegaLI = matrix(1:(J*J),J,J)
OmegaLP[lower.tri(OmegaLP)] = 0
OmegaLI[lower.tri(OmegaLI)] = 0
OmegaCovP = t(OmegaLP)%*%(OmegaLP)
OmegaCovI = t(OmegaLI)%*%(OmegaLI)
for(j1 in 1:J){
	for(j2 in j1:J){
		OmegaP[j1,j2] = OmegaP[j2,j1] = OmegaCovP[j1,j2]/(sqrt(OmegaCovP[j1,j1]*OmegaCovP[j2,j2]))
		OmegaI[j1,j2] = OmegaI[j2,j1] = OmegaCovI[j1,j2]/(sqrt(OmegaCovI[j1,j1]*OmegaCovI[j2,j2]))
	}
}
InvOmegaP = solve(OmegaP)
InvOmegaI = solve(OmegaI)
Sigma = as.vector((1:Q)/Q)
Sigma2 = Sigma^2

##################################################
### Simulate Covariates & Response

set.seed(12345)
	
UP = rmnorm(N, mean=rep(0,J), varcov=OmegaP)
UI = rmnorm(N, mean=rep(0,J), varcov=OmegaI)
eta = array(NA,dim=c(N,J,K))
for(j in 1:J){
	eta[,j,1] = Alpha[1] + X%*%BetaP + Z[,1,]%*%GammaP + UP[,j]
	for(k in 2:K){
		eta[,j,k] = Alpha[k] + X%*%BetaI + Z[,k,]%*%GammaI + UI[,j]
	}	
}
W = array(NA,dim=c(N,J,K))
Wmis = array(NA,dim=c(N,J,K))
for(i in 1:N){
	for(j in 1:J){
		W[i,j,1] =  rnorm(1,eta[i,j,1],1)
		Wmis[i,j,1] = rnorm(1,W[i,j,1], Sigma[Xi[i,1]])
		for(k in 2:K){
			W[i,j,k] = rnormleft(W[i,j,k-1], eta[i,j,k],1)
			Wmis[i,j,k] = rnorm(1,W[i,j,k], Sigma[Xi[i,k]])
		}	
	}	
}
Y = ifelse(W>0,1,0)
Ymis = ifelse(Wmis>0,1,0)

table(Y,Ymis)
for(k in 2:K){
	print(table(Y[,,k-1],Y[,,k]))
	print(table(Ymis[,,k-1],Ymis[,,k]))
}

##################################################

##################################################
### Parameters, Initials

param.mis <- c("Beta.P","Gamma.P", 
               "Beta.I","Gamma.I", 
               "Alpha", 
               "Sigma", 
               "Omega.P","Omega.I"
               )

inits.mis <- function(){	list(
  "Beta.P" = rep(0,p) , 
  "Beta.I" = rep(0,p) , 
  "Gamma.P" = rep(0,q) ,
  "Gamma.I" = rep(0,q) ,
  "Alpha" = rep(0,K) ,
  "Tau" = rep(1,Q) 
)	}

### Parameters of interest to check posterior distribution
monit_par <- c(
           "Beta.P[1]","Beta.P[2]",#"Beta.P[3]","Beta.P[4]",
		   "Gamma.P[1]","Gamma.P[2]",#"Gamma.P[3]",
		   "Beta.I[1]","Beta.I[2]",#"Beta.I[3]","Beta.I[4]",
		   "Gamma.I[1]","Gamma.I[2]",#"Gamma.I[3]",
		   "Alpha[1]","Alpha[2]","Alpha[3]","Alpha[4]",#"Alpha[5]","Alpha[6]",
            "Sigma[1]","Sigma[2]","Sigma[3]","Sigma[4]",#  "Sigma[5]","Sigma[6]","Sigma[7]","Sigma[8]",  
		   #"Sigma[9]","Sigma[10]","Sigma[11]","Sigma[12]",  "Sigma[13]","Sigma[14]","Sigma[15]","Sigma[16]",
		   "Omega.P[1,1]","Omega.P[1,2]","Omega.P[1,3]",#"Omega.P[1,4]",  
		   "Omega.P[2,1]","Omega.P[2,2]","Omega.P[2,3]",#"Omega.P[2,4]",  
		   "Omega.P[3,1]","Omega.P[3,2]","Omega.P[3,3]",#"Omega.P[3,4]",  
		   #"Omega.P[4,1]","Omega.P[4,2]","Omega.P[4,3]","Omega.P[4,4]"
		   "Omega.I[1,1]","Omega.I[1,2]","Omega.I[1,3]",#"Omega.I[1,4]",  
		   "Omega.I[2,1]","Omega.I[2,2]","Omega.I[2,3]",#"Omega.I[2,4]",  
		   "Omega.I[3,1]","Omega.I[3,2]","Omega.I[3,3]"#,"Omega.I[3,4]",  
		   #"Omega.I[4,1]","Omega.I[4,2]","Omega.I[4,3]","Omega.I[4,4]"
)

##################################################

##################################################
### Data
data.mis <- list(
  Yobs = Ymis ,
  N = N , 
  J = J ,
  K = K , 
  X = X , 
  Z = Z ,
  Q = Q ,
  Xi = Xi ,
  Uzeros = rep(0,J) ,
  p = p ,
  q = q
)
##################################################

##################################################
### Mixed HMM for misclassification model

fit.mis <- jags.model("MisMHMM.bug", data=data.mis, inits=inits.mis, n.chains=2)

update(fit.mis,5000)

sample.mis <- coda.samples(fit.mis, param.mis, n.iter=5000, thin=2)

post.mis <- summary(sample.mis)
post.mis
plot(sample.mis)

par(mfrow=c(4,4))
traceplot(sample.mis)

##################################################

##################################################
### Mixed HMM for misclassification model
### By using paralell in JAGS

### Creating clusters to run a chain in each core
detectCores()
cl <- makeCluster(3, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))

### Sampling from model
nburn <- 5000
niter <- 5000
nthin <- 2
thous <- 8

parJagsModel(cl, name="fit.mis", file="MisMHMM.jag", data=data.mis, inits=inits.mis, n.chains=2)
parUpdate(cl, "fit.mis", n.iter=nburn, thin=nthin)

for(part in 1:thous){
  fit.mis_s <- parCodaSamples(cl, model="fit.mis", variable.names=monit_par, n.iter=niter, thin=nthin)
  if(part==1){
    sampl_par <- dimnames(fit.mis_s[[1]])[[2]]
    for(chain in 1:2)
      for(par in sampl_par){
        fileout <- paste0("Chain",chain,par,".txt")
        cat(par, file=fileout, sep = "\n", append = FALSE)
      }
    }
  for(chain in 1:2){
    for(par in sampl_par){
      fileout <- paste0("Chain",chain,par,".txt")
      cat(fit.mis_s[[chain]][,par], file=fileout, sep = "\n", append = TRUE)
    }
  }
}
for(chain in 1:2){
  fileout <- paste0("Chain",chain,".txt")
  import <- NULL
  for(par in sampl_par){
    filein <- paste0("Chain",chain,par,".txt")
    import <- rbind(import, matrix(read.table(filein, header=T)[,1], ncol=1))
  }
  rows <- seq(nburn+1, nburn+niter*thous-nthin+1, by=nthin)
  write.table(cbind(rows,import), file=fileout,quote=FALSE, sep=" ", row.names=FALSE, col.names=FALSE)
}

index <- 0
indexM <- matrix(0, nrow=length(sampl_par), ncol=2)
for(ii in 1:length(sampl_par)){
  indexM[ii,] <- c(index+1,index+thous*niter/nthin)
  index <- index+thous*niter/nthin
}
rownames(indexM) <- sampl_par
write.table(indexM, file=paste0("CODAindex",".txt"), quote=FALSE, sep=" ", row.names=TRUE, col.names=FALSE)

fit.mis1 <- read.coda("Chain1.txt","CODAindex.txt")
fit.mis2 <- read.coda("Chain2.txt","CODAindex.txt")
###fit.mis3 <- read.coda("Chain3.txt","CODAindex.txt")

fit.mis_s <- mcmc.list(fit.mis1,fit.mis2)###,fit.mis3)
traceplot(fit.mis_s)
postSum1 <- summary(fit.mis_s[,monit_par[]])


table1 <- cbind(postSum1$statistics[,1], postSum1$quantiles[,"50%"], postSum1$statistics[,2], postSum1$quantiles[,c("2.5%","97.5%")])
colnames(table1)[1:3] <- c("Mean", "Median", "SD")

print(table1)
write.table(table1, file=paste0("table1",".txt"), quote=FALSE, sep=" ", row.names=TRUE)

##################################################
##################################################
