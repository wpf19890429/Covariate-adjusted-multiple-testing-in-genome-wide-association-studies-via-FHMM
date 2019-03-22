em.fhmm=function(x, maxiter)
{
#####################################################################################

## USAGE
 # em.fhmm(x, maxiter)

## ARGUMENTS
 # x: the observed data
 # maxiter: the maximum number of iterations

## DETAILS
 # em.fhmm calculates the parameters of factorial hidden Markov models via the EM algrithm.

## VALUES
 # pii.new: the update of initial state distribution for main effects
 # piii.new: the update of initial state distribution for covariate effects
 # A.new=(a00, a01 \\ a10 a11): the update of transition matrix for main effects
 # J.new=(j00, j01 \\ j10 j11): the update of transition matrix for covariate effects
 # f00: parameter set for the null distribution
 # f01.new,f10.new,f11.new: the update of parameter set for the non-null distribution 
 # lfdr: the lfdr variables
 # niter: number of iterations

####################################################################################
   NUM=length(x)
# precision tolerance level
   ptol=1e-3
   niter=0 
### initializing model parameters
   pii.new=c(0.6, 0.4)
   piii.new=c(0.5,0.5)
   A.new=matrix(c(0.95, 0.05, 0.10, 0.90), 2, 2, byrow=T)
   J.new=matrix(c(0.90, 0.10, 0.10, 0.90), 2, 2, byrow=T)
   f00=c(0,1)
   f01.new=c(0, 1)
   f10.new=c(1, 1)
   f11.new=c(2, 1)
   diff=1
### The E-M Algorithm
   while(diff>ptol && niter<maxiter)
   {
         niter=niter+1
         pii.old=pii.new
         piii.old=piii.new
         A.old=A.new
         J.old=J.new
         f01.old=f01.new
         f10.old=f10.new
         f11.old=f11.new
## updating the weights and probabilities of hidden states
         bwfw.fhmm.res=bwfw.fhmm(x,pii.old,piii.old,A.old,J.old,f00,f01.old,f10.old,f11.old)
         gamma=bwfw.fhmm.res$pr
         xi0=bwfw.fhmm.res$ts
         xi=matrix(rep(0,(NUM-1)*16),(NUM-1),16,byrow=T)
         xi[ ,1]=xi0[1,1,1,1, ]
         xi[ ,2]=xi0[2,1,1,1, ]
         xi[ ,3]=xi0[1,2,1,1, ]
         xi[ ,4]=xi0[2,2,1,1, ]
         xi[ ,5]=xi0[1,1,2,1, ]
         xi[ ,6]=xi0[2,1,2,1, ]
         xi[ ,7]=xi0[1,2,2,1, ]
         xi[ ,8]=xi0[2,2,2,1, ]
         xi[ ,9]=xi0[1,1,1,2, ]
         xi[ ,10]=xi0[2,1,1,2, ]
         xi[ ,11]=xi0[1,2,1,2, ]
         xi[ ,12]=xi0[2,2,1,2, ]
         xi[ ,13]=xi0[1,1,2,2, ]
         xi[ ,14]=xi0[2,1,2,2, ]
         xi[ ,15]=xi0[1,2,2,2, ]
         xi[ ,16]=xi0[2,2,2,2, ]
         q00=sum((xi[1:(NUM-1),1]+xi[1:(NUM-1),5]+xi[1:(NUM-1),9]+xi[1:(NUM-1),13]))
         p00=sum((gamma[1:(NUM-1),1]+gamma[1:(NUM-1),3]))
         q01=sum((xi[1:(NUM-1),3]+xi[1:(NUM-1),7]+xi[1:(NUM-1),11]+xi[1:(NUM-1),15]))
         p01=sum((gamma[1:(NUM-1),1]+gamma[1:(NUM-1),3]))
         q10=sum((xi[1:(NUM-1),2]+xi[1:(NUM-1),6]+xi[1:(NUM-1),10]+xi[1:(NUM-1),14]))
         p10=sum((gamma[1:(NUM-1),2]+gamma[1:(NUM-1),4]))
         q11=sum((xi[1:(NUM-1),4]+xi[1:(NUM-1),8]+xi[1:(NUM-1),12]+xi[1:(NUM-1),16]))
         p11=sum((gamma[1:(NUM-1),2]+gamma[1:(NUM-1),4]))
         j00=sum((xi[1:(NUM-1),1]+xi[1:(NUM-1),2]+xi[1:(NUM-1),3]+xi[1:(NUM-1),4]))
         k00=sum((gamma[1:(NUM-1),1]+gamma[1:(NUM-1),2]))
         j01=sum((xi[1:(NUM-1),9]+xi[1:(NUM-1),10]+xi[1:(NUM-1),11]+xi[1:(NUM-1),12]))
         k01=sum((gamma[1:(NUM-1),1]+gamma[1:(NUM-1),2]))
         j10=sum((xi[1:(NUM-1),5]+xi[1:(NUM-1),6]+xi[1:(NUM-1),7]+xi[1:(NUM-1),8]))
         k10=sum((gamma[1:(NUM-1),3]+gamma[1:(NUM-1),4]))
         j11=sum((xi[1:(NUM-1),13]+xi[1:(NUM-1),14]+xi[1:(NUM-1),15]+xi[1:(NUM-1),16]))
         k11=sum((gamma[1:(NUM-1),3]+gamma[1:(NUM-1),4]))
## updating the parameter estimates
# a. initial state distribution
         pii.new[1]=gamma[1, 1]+gamma[1, 3]
         pii.new[2]=gamma[1, 2]+gamma[1, 4]
         piii.new[1]=gamma[1, 1]+gamma[1, 2]
         piii.new[2]=gamma[1, 3]+gamma[1, 4]
# b. transition matrix
         A.new[1,1]=q00/p00
         A.new[1,2]=q01/p01
         A.new[2,1]=q10/p10
         A.new[2,2]=q11/p11
         J.new[1,1]=j00/k00
         J.new[1,2]=j01/k01
         J.new[2,1]=j10/k10
         J.new[2,2]=j11/k11
# c. non-null distribution 
         q1_01=sum(gamma[, 3])
         q2_01=sum(gamma[, 3]*x)
         mu01=q2_01/q1_01
         q3_01=sum(gamma[, 3]*(x-mu01)*(x-mu01))
         sd01=sqrt(q3_01/q1_01)
         f01.new=c(mu01, sd01)
         q1_10=sum(gamma[, 2])
         q2_10=sum(gamma[, 2]*x)
         mu10=q2_10/q1_10
         q3_10=sum(gamma[, 2]*(x-mu10)*(x-mu10))
         sd10=sqrt(q3_10/q1_10)
         f10.new=c(mu10, sd10)
         q1_11=sum(gamma[, 4])
         q2_11=sum(gamma[, 4]*x)
         mu11=q2_11/q1_11
         q3_11=sum(gamma[, 4]*(x-mu11)*(x-mu11))
         sd11=sqrt(q3_11/q1_11)
         f11.new=c(mu11, sd11)

         df1=abs(A.old-A.new)
         df1=abs(J.old-J.new)
         df2=abs(f01.old-f01.new)
         df3=abs(f10.old-f10.new)
         df4=abs(f11.old-f11.new)
         diff=max(df1,df2,df3,df4)
}
# f. the final local fdr statistic
lfdr=gamma[, 1]+gamma[, 3]
# g. return the results of the E-M algorithm
em.var=list(pii=pii.new,piii=piii.new,A=A.new,J=J.new,f01=f01.new,f10=f10.new,f11=f11.new,calis=lfdr,ni=niter,diff=diff)
return (em.var)
}




