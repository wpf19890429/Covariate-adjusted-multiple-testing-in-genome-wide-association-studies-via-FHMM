bwfw.fhmm=function(x,pii,piii,A,J,f00,f01,f10,f11)
{ 
###################################################################

## USAGE
 # bwfw.fhmm(x,pii,piii,A,J,f00,f01,f10,f11)

## ARGUMENTS
 # x=(x[1], ..., x[m]): the observed data
 # pii=(pii[1], pii[2]): initial state distribution for main effects
 # piii=(piii[1], piii[2]): initial state distribution for covariate effects
 # A=(a00, a01 \\ a10 a11): transition matrix for main effects
 # J=(j00, j01 \\ j10 j11): transition matrix for covariate effects
 # f00: parameter set for the null distribution
 # f01,f10,f11: parameter set for the non-null distribution

## DETAILS
 # bwfw.fhmm calculates values for backward, forward variables, probabilities of hidden states, 
 # --the lfdr variables and etc. 
 # --using the forward-backward procedure (Baum et al.) 
 # --based on a sequence of observations for a given factorial hidden markov model 

## VALUES
 # alpha: rescaled backward variables
 # beta: rescaled forward variables
 # lfdr: lfdr variables
 # gamma: probabilities of hidden states
 # dgamma: rescaled transition variables


################################################################## 
## Initialize
      NUM=length(x)
## Densities
      f00x=dnorm(x, f00[1], f00[2])
      f01x=dnorm(x, f01[1], f01[2])
      f10x=dnorm(x, f10[1], f10[2])
      f11x=dnorm(x, f11[1], f11[2])
## the backward-forward procedure

# a. the backward variables
# --rescaled 
      alpha=matrix(rep(0, NUM*4), NUM, 4, byrow=T)
# scaling variable c0
      c0=rep(0, NUM)
      alpha[1, 1]=pii[1]*piii[1]*f00x[1]
      alpha[1, 2]=pii[2]*piii[1]*f10x[1]
      alpha[1, 3]=pii[1]*piii[2]*f01x[1]
      alpha[1, 4]=pii[2]*piii[2]*f11x[1]
# rescaling alpha
      c0[1]=1/sum(alpha[1, ])
      alpha[1, ]=c0[1]*alpha[1, ]
      for (k in 1:(NUM-1))
      { 
          alpha[k+1, 1]=(alpha[k, 1]*A[1, 1]*J[1,1]+alpha[k, 2]*A[2, 1]*J[1,1]
                        +alpha[k, 3]*A[1, 1]*J[2,1]+alpha[k, 4]*A[2, 1]*J[2,1])*f00x[k+1]
          alpha[k+1, 2]=(alpha[k, 1]*A[1, 2]*J[1,1]+alpha[k, 2]*A[2, 2]*J[1,1]
                        +alpha[k, 3]*A[1, 2]*J[2,1]+alpha[k, 4]*A[2, 2]*J[2,1])*f10x[k+1]
          alpha[k+1, 3]=(alpha[k, 1]*A[1, 1]*J[1,2]+alpha[k, 2]*A[2, 1]*J[1,2]
                        +alpha[k, 3]*A[1, 1]*J[2,2]+alpha[k, 4]*A[2, 1]*J[2,2])*f01x[k+1]
          alpha[k+1, 4]=(alpha[k, 1]*A[1, 2]*J[1,2]+alpha[k, 2]*A[2, 2]*J[1,2]
                        +alpha[k, 3]*A[1, 2]*J[2,2]+alpha[k, 4]*A[2, 2]*J[2,2])*f11x[k+1]
          # rescaling alpha
          c0[k+1]=1/sum(alpha[k+1, ])
          alpha[k+1, ]=c0[k+1]*alpha[k+1, ]
      }
# b. the forward variables
# --rescaled
      beta=matrix(rep(0, NUM*4), NUM, 4, byrow=T)
      beta[NUM, 1]=c0[NUM]
      beta[NUM, 2]=c0[NUM]
      beta[NUM, 3]=c0[NUM]
      beta[NUM, 4]=c0[NUM]
      for (k in (NUM-1):1)
      {
           beta[k, 1]=A[1, 1]*J[1, 1]*f00x[k+1]*beta[k+1, 1]+A[1, 2]*J[1, 1]*f10x[k+1]*beta[k+1, 2]+
                      A[1, 1]*J[1, 2]*f01x[k+1]*beta[k+1, 3]+A[1, 2]*J[1, 2]*f11x[k+1]*beta[k+1, 4]
           beta[k, 2]=A[2, 1]*J[1, 1]*f00x[k+1]*beta[k+1, 1]+A[2, 2]*J[1, 1]*f10x[k+1]*beta[k+1, 2]+
                      A[2, 1]*J[1, 2]*f01x[k+1]*beta[k+1, 3]+A[2, 2]*J[1, 2]*f11x[k+1]*beta[k+1, 4]
           beta[k, 3]=A[1, 1]*J[2, 1]*f00x[k+1]*beta[k+1, 1]+A[1, 2]*J[2, 1]*f10x[k+1]*beta[k+1, 2]+
                      A[1, 1]*J[2, 2]*f01x[k+1]*beta[k+1, 3]+A[1, 2]*J[2, 2]*f11x[k+1]*beta[k+1, 4]
           beta[k, 4]=A[2, 1]*J[2, 1]*f00x[k+1]*beta[k+1, 1]+A[2, 2]*J[2, 1]*f10x[k+1]*beta[k+1, 2]+
                      A[2, 1]*J[2, 2]*f01x[k+1]*beta[k+1, 3]+A[2, 2]*J[2, 2]*f11x[k+1]*beta[k+1, 4]
           # rescaling beta
           # using the same scaling factors as alpha 

           beta[k,  ]=c0[k]*beta[k, ]
       }
# c. lfdr variables
# --original
# --the same formulae hold for the rescaled alpha and beta
       lfdr=rep(0, NUM)
       for (k in 1:NUM)
       { 
           q1=alpha[k, 1]*beta[k, 1]+alpha[k, 3]*beta[k, 3]
           q2=alpha[k, 2]*beta[k, 2]+alpha[k, 4]*beta[k, 4]
           lfdr[k]=q1/(q1+q2)
       }
# d. probabilities of hidden states
# -- and transition variables
# -- both are rescaled
       gamma=matrix(rep(0, NUM*4),NUM,4,byrow=T)
       # initialize gamma
       gamma[NUM,1]=alpha[NUM, 1]*beta[NUM, 1]/(q1+q2)
       gamma[NUM,2]=alpha[NUM, 2]*beta[NUM, 2]/(q1+q2)
       gamma[NUM,3]=alpha[NUM, 3]*beta[NUM, 3]/(q1+q2)
       gamma[NUM,4]=alpha[NUM, 4]*beta[NUM, 4]/(q1+q2)
       dgamma=array(rep(0, (NUM-1)*16), c(2, 2, 2, 2, (NUM-1)))
       for (i in 1:(NUM-1))
       {
            denom0=0
            denom1=0
            denom2=0
            denom3=0
            for (t in 0:1)
            {   
                for(j in 0:1)
                { 
                   f0x=(1-j)*f00x[i+1]+j*f10x[i+1]
                   f1x=(1-j)*f01x[i+1]+j*f11x[i+1]
                   fx=c(f0x,f1x)
                   denom0=denom0+alpha[i, t+1]*A[t+1, j+1]*J[1,1]*fx[1]*beta[i+1, j+1]
                   denom1=denom1+alpha[i, t+3]*A[t+1, j+1]*J[2,1]*fx[1]*beta[i+1, j+1]
                   denom2=denom2+alpha[i, t+1]*A[t+1, j+1]*J[1,2]*fx[2]*beta[i+1, j+3]
                   denom3=denom3+alpha[i, t+3]*A[t+1, j+1]*J[2,2]*fx[2]*beta[i+1, j+3]
                 }    
             }
             for(t in 0:1)
             {
                 gamma[i,t+1]=0
                 gamma[i,t+3]=0
                 for (j in 0:1)
                 { 
                     f0x=(1-j)*f00x[i+1]+j*f10x[i+1]
                     f1x=(1-j)*f01x[i+1]+j*f11x[i+1]
                     fx=c(f0x,f1x)
                     dgamma[t+1,j+1,1,1,i]=alpha[i, t+1]*A[t+1, j+1]*J[1,1]*fx[1]*beta[i+1, j+1]/(denom0+denom1+denom2+denom3)
                     dgamma[t+1,j+1,2,1,i]=alpha[i, t+3]*A[t+1, j+1]*J[2,1]*fx[1]*beta[i+1, j+1]/(denom0+denom1+denom2+denom3)
                     dgamma[t+1,j+1,1,2,i]=alpha[i, t+1]*A[t+1, j+1]*J[1,2]*fx[2]*beta[i+1, j+3]/(denom0+denom1+denom2+denom3)
                     dgamma[t+1,j+1,2,2,i]=alpha[i, t+3]*A[t+1, j+1]*J[2,2]*fx[2]*beta[i+1, j+3]/(denom0+denom1+denom2+denom3)   
                  }
              }
              gamma[i,1]=dgamma[1,1,1,1,i]+dgamma[1,2,1,1,i]+dgamma[1,1,1,2,i]+dgamma[1,2,1,2,i]
              gamma[i,2]=dgamma[2,1,1,1,i]+dgamma[2,2,1,1,i]+dgamma[2,1,1,2,i]+dgamma[2,2,1,2,i]
              gamma[i,3]=dgamma[1,1,2,1,i]+dgamma[1,2,2,1,i]+dgamma[1,1,2,2,i]+dgamma[1,2,2,2,i]
              gamma[i,4]=dgamma[2,1,2,1,i]+dgamma[2,2,2,1,i]+dgamma[2,1,2,2,i]+dgamma[2,2,2,2,i]
       }
# f. return the results of the bwfw proc.
       bwfw.var<-list(bw=alpha, fw=beta, calis.or=lfdr, pr=gamma, ts=dgamma,c0=c0)
       return(bwfw.var)
 }