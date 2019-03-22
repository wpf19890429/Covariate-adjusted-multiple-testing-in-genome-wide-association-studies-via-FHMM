rdata.fhmm=function(NUM,pii,piii,A,J,f00,f01,f10,f11)
{    

## USAGE
 # rdata.fhmm(NUM,pii,piii,A,J,...)

## ARGUEMENTS
 # NUM: sample size
 # pii=(pii[1], pii[2]): initial state distribution for main effects
 # piii=(piii[1], piii[2]): initial state distribution for covariate effects
 # A=(a00, a01 \\ a10 a11): transition matrix for main effects
 # J=(j00, j01 \\ j10 j11): transition matrix for covariate effects
 # f00: parameter set for the null distribution
 # f01,f10,f11: parameter set for the non-null distribution

## VALUES
 # rdata.fhmm generates random variables 
 # from a four-component normal mixture via a factorial hidden markov model
 # x: continuous observed data
 # theta: binary unobserved states for main effects
 # gamma: binary unobserved states for covariate effects
     theta=rep(0,NUM) 
     x=rep(0,NUM)
     gamma=rep(0,NUM)
## generating the states
 # initial state
     theta[1]=rbinom(1,1,pii[2])
     gamma[1]=rbinom(1,1,piii[2])
 # other states
     for(i in 2:NUM)
     {
         if(theta[i-1]==0)
           theta[i]=rbinom(1,1,A[1, 2])
         else
           theta[i]=rbinom(1,1,A[2, 2])
             
         if(gamma[i-1]==0)
            gamma[i]=rbinom(1,1,J[1,2])
          else
            gamma[i]=rbinom(1,1,J[2,2])
     }
## generating the observations
     for(i in 1:NUM)
     { 
         if(theta[i]==0&&gamma[i]==0)
         {
            x[i]=rnorm(1,mean=f00[1],sd=f00[2])
         }
         else
         {
            if(theta[i]==0&&gamma[i]==1)
            {
               x[i]=rnorm(1,mean=f01[1],sd=f01[2])
            }
            else
            {
               if(theta[i]==1&&gamma[i]==0)
               {
                  x[i]=rnorm(1,mean=f10[1],sd=f10[2])
               }
               else
               {
                  x[i]=rnorm(1,mean=f11[1],sd=f11[2])
               }
    
            }
          }
       }
       data=list(s=theta,o=x)
       return(data)
}