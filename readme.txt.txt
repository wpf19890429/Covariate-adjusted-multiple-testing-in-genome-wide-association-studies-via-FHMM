
# This document contains four main files: 
  # function 'rdata.fhmm' (rdata.fhmm) is the factorial HMM data generator
  # function 'bwfw.fhmm' (bwfw.fhmm) realizes the forward-backward procedure
  # function 'em.fhmm' (em.fhmm) realizes the E-M algorithm
  # function 'mt.fhmm' realizes the new multiple testing procedures 
    # (given the estimates of factorial HMM parameters)

# Examples on how to use these functions are given below. 

# More detailed instructions are given in each separate file. 

###############################
######   EXAMPLES   ###########
###############################

source("rdata.fhmm.R")
source("bwfw.fhmm.R")
source("em.fhmm.R")
source("mt.fhmm.R")

 # the number of observations
NUM=3000
 # the prespecified FDR level
q=0.10

####################################
# Example 1: the HMM data generator
####################################

 # the initial state distribution
pii=c(0.95,0.05)
piii=c(0.8,0.2)	
 # the transition matrix
A=matrix(c(0.95,0.05,0.1,0.9),2,2,byrow=T)
J=matrix(c(0.9,0.1,0.05,0.95),2,2,byrow=T)
 # the null distribution
f00=c(0,1)
 # the alternative distribution
f01=c(-1, 1)
f10=c(1, 1)
f11=c(3, 1)
 # the factorial HMM data
set.seed(123456)
rdata1=rdata.fhmm(NUM, pii, piii, A, J, f00, f01, f10, f11)
 # the observed values
x1=rdata1$o
 # the unobserved states
theta1=rdata1$s


############################################
# Example 2: the forward-backward procedure
############################################

x1=rdata1$o
fb.res1=bwfw.fhmm(x1, pii, piii, A, J, f00, f01, f10, f11)
# the backward variable
backward.var=fb.res1$bw
# the backward variable
forward.var=fb.res1$fw
# the oracle calis (calis.or) variable
calis.or.var=fb.res1$calis.or


##############################################################################
# Example 3: the E-M Algorithm for calculating parameters of the factorial HMM 
##############################################################################

 # the EM algorithm
em.res1=em.fhmm(x1, maxiter=200)
 # the estimates for factorial HMM parameters
em.res1$A
em.res1$J
em.res1$f01
em.res1$f10
em.res1$f11
 # the data-driven calis (calis.dd) variables 
em.res1$calis
 # the number of interations 
em.res1$n1

#################################################
# Example 4: The CALIS.or and CALIS.dd procedures
#################################################

## (4.a) the CALIS.or procedure
 
 # the CALIS.or values
calis.or=fb.res1$calis.or
calis.or.pi=mt.fhmm(calis.or, q)
 # the decision rule
calis.or.de=calis.or.pi$de

## (4.b) the CALIS.dd procedure

 # the CALIS.dd variables 
calis.dd=em.res1$calis
calis.dd.pi=mt.hmm(calis.dd, q)
 # the decision rule
calis.dd.de=calis.dd.pi$de

####################################################################
###### Example 5: the analysis of the Bipolar Disorder data ########
####################################################################

###calculate the plug-in CALIS statistic for individual chromosomes
## (5.a) The data from chromosome 2
 # the Bipolar Disorder data of chromosome 2
BD.data.chr2=read.table("BD_data_chr2.txt")
BD.data.chr2=as.matrix(BD.data.chr2)

 # the CALIS values
BD.chr2.res=em.fhmm(BD.data.chr2, maxiter=200)
BD.chr2.calis=BD.chr2.res$calis
BD.chr2.calis=as.matrix(BD.chr2.calis)
write.table(BD.chr2.calis,"BD.chr2.calis.txt")

## (5.b) The data from chromosome 3
 # the Bipolar Disorder data of chromosome 3
BD.data.chr3=read.table("BD_data_chr3.txt")
BD.data.chr3=as.matrix(BD.data.chr3)

 # the CALIS values
BD.chr3.res=em.fhmm(BD.data.chr3, maxiter=200)
BD.chr3.calis=BD.chr3.res$calis
BD.chr3.calis=as.matrix(BD.chr3.calis)
write.table(BD.chr3.calis,"BD.chr3.calis.txt")

## (5.c) The data from chromosome 4
 # the Bipolar Disorder data of chromosome 4
BD.data.chr4=read.table("BD_data_chr4.txt")
BD.data.chr4=as.matrix(BD.data.chr4)

 # the CALIS values
BD.chr4.res=em.fhmm(BD.data.chr4, maxiter=200)
BD.chr4.calis=BD.chr4.res$calis
BD.chr4.calis=as.matrix(BD.chr4.calis)
write.table(BD.chr4.calis,"BD.chr4.calis.txt")

## (5.d) The data from chromosome 6
 # the Bipolar Disorder data of chromosome 6
BD.data.chr6=read.table("BD_data_chr6.txt")
BD.data.chr6=as.matrix(BD.data.chr6)

 # the CALIS values
BD.chr6.res=em.fhmm(BD.data.chr6, maxiter=200)
BD.chr6.calis=BD.chr6.res$calis
BD.chr6.calis=as.matrix(BD.chr6.calis)
write.table(BD.chr6.calis,"BD.chr6.calis.txt")


## (5.e) The data from chromosome 8
 # the Bipolar Disorder data of chromosome 8
BD.data.chr8=read.table("BD_data_chr8.txt")
BD.data.chr8=as.matrix(BD.data.chr8)

 # the CALIS values
BD.chr8.res=em.fhmm(BD.data.chr8, maxiter=200)
BD.chr8.calis=BD.chr8.res$calis
BD.chr8.calis=as.matrix(BD.chr8.calis)
write.table(BD.chr8.calis,"BD.chr8.calis.txt")

## (5.f) The data from chromosome 9
 # the Bipolar Disorder data of chromosome 9
BD.data.chr9=read.table("BD_data_chr9.txt")
BD.data.chr9=as.matrix(BD.data.chr9)

 # the CALIS values
BD.chr9.res=em.fhmm(BD.data.chr9, maxiter=200)
BD.chr9.calis=BD.chr9.res$calis
BD.chr9.calis=as.matrix(BD.chr9.calis)
write.table(BD.chr9.calis,"BD.chr9.calis.txt")

## (5.g) The data from chromosome 14
 # the Bipolar Disorder data of chromosome 14
BD.data.chr14=read.table("BD_data_chr14.txt")
BD.data.chr14=as.matrix(BD.data.chr14)

 # the CALIS values
BD.chr14.res=em.fhmm(BD.data.chr14, maxiter=200)
BD.chr14.calis=BD.chr14.res$calis
BD.chr14.calis=as.matrix(BD.chr14.calis)
write.table(BD.chr14.calis,"BD.chr14.calis.txt")

## (5.h) The data from chromosome 16
 # the Bipolar Disorder data of chromosome 16
BD.data.chr16=read.table("BD_data_chr16.txt")
BD.data.chr16=as.matrix(BD.data.chr16)

 # the CALIS values
BD.chr16.res=em.fhmm(BD.data.chr16, maxiter=200)
BD.chr16.calis=BD.chr16.res$calis
BD.chr16.calis=as.matrix(BD.chr16.calis)
write.table(BD.chr16.calis,"BD.chr16.calis.txt")

## (5.i) The data from chromosome 20
 # the Bipolar Disorder data of chromosome 20
BD.data.chr20=read.table("BD_data_chr20.txt")
BD.data.chr20=as.matrix(BD.data.chr20)

 # the CALIS values
BD.chr20.res=em.fhmm(BD.data.chr20, maxiter=200)
BD.chr20.calis=BD.chr20.res$calis
BD.chr20.calis=as.matrix(BD.chr20.calis)
write.table(BD.chr20.calis,"BD.chr20.calis.txt")

## (5.j) The data from chromosome 22
 # the Bipolar Disorder data of chromosome 22
BD.data.chr22=read.table("BD_data_chr22.txt")
BD.data.chr22=as.matrix(BD.data.chr22)

 # the CALIS values
BD.chr22.res=em.fhmm(BD.data.chr22, maxiter=200)
BD.chr22.calis=BD.chr22.res$calis
BD.chr22.calis=as.matrix(BD.chr22.calis)
write.table(BD.chr22.calis,"BD.chr22.calis.txt")

############################################################################
###combine and rank the plug-in CALIS statistic from all chromosomes
 
 # combine the plug-in CALIS statistic from the ten chromosomes

BD.calis.all=rbind(BD.chr2.calis,BD.chr3.calis,BD.chr4.calis,
                   BD.chr6.calis,BD.chr8.calis,BD.chr9.calis,
                   BD.chr14.calis,BD.chr16.calis,BD.chr20.calis,
                   BD.chr22.calis)
 # the FDR level
q=1e-07

 # the CALIS procedure
BD.calis.pi=mt.fhmm(BD.calis.all, q)
 # the threshold for calis
BD.calis.pi$th
 # the number of rejections
BD.calis.pi$nr
 # the indices of rejected hypotheses
BD.calis.pi$re
