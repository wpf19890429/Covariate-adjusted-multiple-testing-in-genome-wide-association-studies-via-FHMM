mt.fhmm=function(lfdr,q)
{  

## USAGE
 # mt.fhmm(lfdr,q)
 
## ARGUMENTS
 # lfdr: local false discovery rate sequence
 # q: the FDR level

## DETAILS
 # mt.fhmm gives a multiple testing rule in a factorial hidden markov model
 # --that controls the FDR at level q, based on sequence of lfdr

## VALUES
 # nr: the number of rejected hypotheses
 # th: the threshold
 # re: the rejected hypotheses
 # ac: the accepted hypotheses
 # de: the decision rule

        m=length(lfdr)
        st.lfdr=sort(lfdr)
        hps=rep(0,m)
        if(min(lfdr)>q)
         {
             k=0
             threshold=1
             reject=NULL
             accept=1:m

          }
          else
          {
              k=1
              while(k<m && (1/k)*sum(st.lfdr[1:k])<q){
                 k=k+1   
               }
               k=k-1
               threshold=st.lfdr[k]
               reject=which(lfdr<=threshold)
               accept=which(lfdr>threshold)
               hps[reject]=1
       
           }
           y=list(nr=k,th=threshold,re=reject,ac=accept,de=hps)
           return(y)
}