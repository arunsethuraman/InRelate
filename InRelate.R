# This program estimates pairwise genetic relatedness by the MC2013 and MC2013WI methods.
# To run, source in R. Make sure that you install Rsolnp. Also, to run this program successfully, you will have to edit the STR file
# The standard STRUCTURE (Pritchard et al. 2000) format does not have a separate column name for Individuals, and Populations.
# So add the words "Indiv	Pop	" to the beginning of the STR file before running this program.
# This program assumes that you have obtained the INDIVQ and PKLA files from running MULTICLUST or STRUCTURE.
# Written by Arun Sethuraman (arun@temple.edu). 
# Modified on October 28, 2013.

# Initialize
#library(Rsolnp);
# Here are some examples: You have to read these files first before proceeding to run the script.
#hs<-read.table("example.str1",header=TRUE);
#etaik<-read.table("example.str.indivq");
#pkla<-read.table("example.str.pkla");

# Log Likelihood function that is to be maximized. Note the negative return
# For maximization

loglikibd<-function(ibds) {
# array has to be of dim=c(1, number of loci)
loglik<-array(0,dim=c(1,300))
for(l in 1:300){
sumibdloci<-0
for(i in 1:9){
sumibdloci<-sumibdloci+predictors[l,i]*ibds[i]
}

loglik[l]<-log(sumibdloci)
}
return(-sum(loglik))
}


# Equality constraint, sum of deltas = 1

eqn1<-function(ibds){
 sum=ibds[1]+ibds[2]+ibds[3]+ibds[4]+ibds[5]+ibds[6]+ibds[7]+ibds[8]+ibds[9]
 return(sum)
 }
 
 
# Inequality constraints: all should be >0, <1 
ineq1<-function(ibds){
 z1=ibds[1]
 z2=ibds[2]
 z3=ibds[3]
 z4=ibds[4]
 z5=ibds[5]
 z6=ibds[6]
 z7=ibds[7]
 z8=ibds[8]
 z9=ibds[9]
 return(c(z1,z2,z3,z4,z5,z6,z7,z8,z9))
 }

# Set value of K, number of subpopulations. This can be decided a priori based on sampling info, or 
# by using STRUCTURE (Pritchard 2000), MULTICLUST (Sethuraman et al.), and the methods of Evanno et al. (2005)

K<-3

# Array that holds the final relatedness results. First column has the pair, second has the MC2013WI relatedness, third has MC2013 relatedness
# Change dimensions according to the number of pairs. Here I have 50 pairs => total number of pairwise comparisons= 50*49/2 = 1225
relat<-array(0,dim=c(1225,3))

# Arrays for storing the etaiks. Should be of size K.

indiv1.etaik<-array(dim=c(K))
indiv2.etaik<-array(dim=c(K))

# Predictor array should be of size = Number of Loci x 9

predictors<-array(0,dim=c(300,9))

# Allele frequency arrays for each of four alleles (assume diploid genotypes $A_{i}A_{j}$ and $A_{k}A_{l}$, two individuals).
pklallelei<-array(dim=c(K))
pklallelej<-array(dim=c(K))
pklallelek<-array(dim=c(K))
pklallelel<-array(dim=c(K))

# Loop control variables - don't change!

g<-1
#decrementi<-0
#decrementj<-1

# Compute over pairs of individuals
# Here loop has to be changed according to the STRUCTURE format. I have here 100 individuals, with two lines for each genotype. So each pair has 4 lines.
# So loop over number of dyads x 4, increments of 4
# TODO: AS - loop over all pairs of individuals instead! DONE: 2/25/2014


for(i in 1:50) {
	for (j in i:50) {
		for(x in 1:K)	{
			indiv1.etaik[x]<-etaik[i,x+5]
		}
		for(x in 1:K)	{
			indiv2.etaik[x]<-etaik[j,x+5]
		}
	
	# Storing the individual genotypes
	indiv1a<-hs[i*2-1,]
	indiv1b<-hs[i*2,]
	indiv2a<-hs[j*2-1,]
	indiv2b<-hs[j*2,]

	
	#indiv1a<-hs[i,]
	#indiv1b<-hs[i+1,]
	#indiv2a<-hs[j,]
	#indiv2b<-hs[j+1,]

	# Compute across loci. Here I have 300 loci, so skipping first two lines, run loop from 3 to 300+2=302

	for(l in 3:302)	{
		# Indentify IBS mode.
	
		if(indiv1a[l]==indiv1b[l] && indiv1b[l]==indiv2a[l] && indiv2a[l] ==indiv2b[l])
			IBS<-1
		if(indiv1a[l]==indiv1b[l] && indiv2a[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l])
			IBS<-2
		if((indiv1a[l]==indiv1b[l] && indiv1a[l]==indiv2a[l] && indiv2a[l]!=indiv2b[l]) || (indiv1a[l]==indiv1b[l] && indiv1a[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l]))
			IBS<-3
		if((indiv1a[l]==indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]!=indiv2a[l] && indiv1a[l]!=indiv2b[l]))
			IBS<-4
		if((indiv1a[l]!=indiv1b[l] && indiv1a[l]==indiv2a[l] && indiv2a[l]==indiv2b[l]) || (indiv1a[l]!=indiv1b[l] && indiv1b[l]==indiv2a[l] && indiv2a[l]==indiv2b[l]))
			IBS<-5
		if((indiv2a[l]==indiv2b[l] && indiv1a[l]!=indiv1b[l] && indiv1a[l]!=indiv2a[l] && indiv1b[l]!=indiv2b[l]))
			IBS<-6
		if((indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]==indiv2a[l] && indiv1b[l]==indiv2b[l]) || (indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]==indiv2b[l] && indiv1b[l]==indiv2a[l]))
			IBS<-7
		if(indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && ((indiv1a[l]==indiv2a[l] && indiv1b[l]!=indiv2b[l]) || (indiv1a[l]==indiv2b[l] && indiv1b[l]!=indiv2a[l]) || (indiv1b[l]==indiv2a[l] && indiv1a[l]!=indiv2b[l]) || (indiv1b[l]==indiv2b[l] && indiv1a[l]!=indiv2a[l])))
			IBS<-8
		if((indiv1a[l]!=indiv1b[l] && indiv2a[l]!=indiv2b[l] && indiv1a[l]!=indiv2a[l] && indiv1a[l]!=indiv2b[l] && indiv1b[l]!=indiv2a[l] && indiv1b[l]!=indiv2b[l]))
			IBS<-9

		# Number of uniquealleles in locus 1

		length(which(pkla$V1==l-3))

		# Save locus1 as a subset

		locus1<-subset(pkla, pkla$V1==l-3)

		# Calculate coefficients under S1:

		if(IBS==1) {
			a<-indiv1a[l]
			x<-which(locus1$V2==a[,])
			allelei<-locus1[x,]
	
			for(x in 1:K)	{
				pklallelei[x]<-allelei[1,x+2]
			}

			zi1=0
			zi2=0

			for(x in 1:K){
				zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
				zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			}

			#Under S1:
			d1<-(zi1+zi2)/2
			d2<-zi1*zi2
			d3<-zi1*zi2
			d4<-(zi1*zi2*zi2+zi1*zi1*zi2)/2
			d5<-zi1*zi2
			d6<-(zi1*zi2*zi2+zi1*zi1*zi2)/2
			d7<-zi1*zi2
			d8<-(zi1*zi1*zi2+zi2*zi2*zi1)/2
			d9<-(zi1*zi1*zi2*zi2)

		}

		if(IBS==2)	{

			allelei<-array(dim=c(K))
			allelej<-array(dim=c(K))

			a<-indiv1a[l]
			b<-indiv2a[l]
			x<-which(locus1$V2==a[,])
			y<-which(locus1$V2==b[,])

			allelei<-locus1[x,]
			allelej<-locus1[y,]


			for(x in 1:K)	{
				pklallelei[x]<-allelei[1,x+2]
				pklallelej[x]<-allelej[1,x+2]
			}

			zi1=0
			zi2=0
			zj1=0
			zj2=0

			for(x in 1:K)	{
				zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
				zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
				zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
				zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			}

			# Under S2
			d1<-0
			d2<-(zi1*zj2+zj1*zi2)/2
			d3<-0
			d4<-(zi1*zj2*zj2+zj1*zi2*zi2)/2
			d5<-0
			d6<-(zi1*zj2*zj2+zj1*zi2*zi2)/2
			d7<-0
			d8<-0
			d9<-(zi1*zi1*zj2*zj2+zj1*zj1*zi2*zi2)/2

		}

		if(IBS==3)	{

			pklallelei<-array(dim=c(K))
			pklallelej<-array(dim=c(K))

			a<-indiv1a[l]
			b<-indiv2b[l]
			x<-which(locus1$V2==a[,])
			y<-which(locus1$V2==b[,])

			allelei<-locus1[x,]
			allelej<-locus1[y,]

			for(x in 1:K)	{
				pklallelei[x]<-allelei[1,x+2]
				pklallelej[x]<-allelej[1,x+2]
			}
			zi1=0
			zi2=0
			zj1=0
			zj2=0

			for(x in 1:K)	{
				zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
				zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
				zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
				zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			}

			# Under S3
			d1<-0
			d2<-0
			d3<-(zi1*zj2+zj1*zi2)/2
			d4<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
			d5<-0
			d6<-0
			d7<-0
			d8<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
			d9<-(zi1*zi1*zi2*zj2+zj1*zj1*zj2*zi2)/2

		}

	if(IBS==4)	{

		a<-indiv1a[l]
		b<-indiv2a[l]
		c<-indiv2b[l]
		x<-which(locus1$V2==a[,])
		y<-which(locus1$V2==b[,])
		z<-which(locus1$V2==c[,])

		allelei<-locus1[x,]
		allelej<-locus1[y,]
		allelek<-locus1[z,]

		pklallelei<-array(dim=c(K))
		pklallelej<-array(dim=c(K))
		pklallelek<-array(dim=c(K))

		for(x in 1:K)	{
			pklallelei[x]<-allelei[1,x+2]
			pklallelej[x]<-allelej[1,x+2]
			pklallelek[x]<-allelek[1,x+2]
		}

		zi1=0
		zi2=0
		zj1=0
		zj2=0
		zk1=0
		zk2=0

		for(x in 1:K)	{
			zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
			zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
			zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
			zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
		}

		# Under S4
		d1<-0
		d2<-0
		d3<-0
		d4<-(zi1*zj2*zk2+zi2*zj1*zk1)/2
		d5<-0
		d6<-0
		d7<-0
		d8<-0
		d9<-(zi1*zi1*zj2*zk2+zi2*zi2*zj1*zk1)/2
	}

	if(IBS==5)	{

		a<-indiv2b[l]
		b<-indiv1a[l]
		x<-which(locus1$V2==a[,])
		y<-which(locus1$V2==b[,])

		allelej<-locus1[x,]
		allelei<-locus1[y,]
		pklallelei<-array(dim=c(K))
		pklallelej<-array(dim=c(K))

		for(x in 1:K)	{
			pklallelei[x]<-allelei[1,x+2]
			pklallelej[x]<-allelej[1,x+2]
		}

#need to define zi1,zi2,zj2,zj1

		zi1=0
		zi2=0
		zj1=0
		zj2=0

		for(x in 1:K)	{
			zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
			zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
			zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
		}

		# Now to get conditional probabilities:
		#Under S5
		d1<-0
		d2<-0
		d3<-0
		d4<-0
		d5<-(zi1*zj2+zj1*zi2)/2
		d6<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
		d7<-0
		d8<-(zi1*zi2*zj2+zj1*zi2*zj2)/2
		d9<-(zi1*zi1*zi2*zj2+zj1*zj1*zj2*zi2)/2
	}

	if(IBS==6)	{

		a<-indiv1a[l]
		b<-indiv1b[l]
		c<-indiv2b[l]
		x<-which(locus1$V2==a[,])
		y<-which(locus1$V2==b[,])
		z<-which(locus1$V2==c[,])

		allelei<-locus1[z,]
		allelej<-locus1[x,]
		allelek<-locus1[y,]


		pklallelei<-array(dim=c(K))
		pklallelej<-array(dim=c(K))
		pklallelek<-array(dim=c(K))

		for(x in 1:K)	{
			pklallelei[x]<-allelei[1,x+2]
			pklallelej[x]<-allelej[1,x+2]
			pklallelek[x]<-allelek[1,x+2]
		}
		zi1=0
		zi2=0
		zj1=0
		zj2=0
		zk1=0
		zk2=0

		for(x in 1:K)	{
			zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
			zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
			zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
			zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
		}

		#Under S6
		d1<-0
		d2<-0
		d3<-0
		d4<-0
		d5<-0
		d6<-(zi1*zj2*zk2+zi2*zj1*zk1)/2
		d7<-0
		d8<-0
		d9<-(zi1*zi1*zj2*zk2+zi2*zi2*zj1*zk1)/2

		}

		if(IBS==7) {

			a<-indiv1a[l]
			b<-indiv1b[l]

			x<-which(locus1$V2==a[,])
			y<-which(locus1$V2==b[,])

			allelei<-locus1[x,]
			allelej<-locus1[y,]
			pklallelei<-array(dim=c(K))
			pklallelej<-array(dim=c(K))

			for(x in 1:K)	{
				pklallelei[x]<-allelei[1,x+2]
				pklallelej[x]<-allelej[1,x+2]
			}

		zi1=0
		zi2=0
		zj1=0
		zj2=0

	for(x in 1:K)	{
		zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
		zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
		zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
		zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
	}

		#Under S7
		d1<-0
		d2<-0
		d3<-0
		d4<-0
		d5<-0
		d6<-0
		d7<-(zi1*zj2+zj1*zi2)/2
		d8<-(zj1*zj2*(zi1+zi2*0.5)+(zi1*zi2*(zj1+zj2)*0.5))/2
		d9<-zi1*zi2*zj1*zj2
	}

	if(IBS==8) {
		a<-indiv1a[l]
		b<-indiv1b[l]
		c<-indiv2b[l]

		x<-which(locus1$V2==a[,])
		y<-which(locus1$V2==b[,])
		z<-which(locus1$V2==c[,])

		allelei<-locus1[x,]
		allelej<-locus1[y,]
		allelek<-locus1[z,]

		pklallelei<-array(dim=c(K))
		pklallelej<-array(dim=c(K))
		pklallelek<-array(dim=c(K))

		for(x in 1:K)	{
			pklallelei[x]<-allelei[1,x+2]
			pklallelej[x]<-allelej[1,x+2]
			pklallelek[x]<-allelek[1,x+2]
		}

		zi1=0
		zi2=0
		zj1=0
		zj2=0
		zk1=0
		zk2=0

		for(x in 1:K)	{
			zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
			zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
			zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
			zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
		}

		#Under S8
	d1<-0
	d2<-0
	d3<-0
	d4<-0
	d5<-0
	d6<-0
	d7<-0
	d8<-0.5*(((zi1+zi2)/2)*(zj1*zk2+zk1*zj2))
	d9<-0.5*((zi1*zi2)*(zj1*zk2+zj2*zk1))
	}

	if(IBS==9)	{
		a<-indiv1a[l]
		b<-indiv1b[l]
		c<-indiv2a[l]
		d<-indiv2b[l]

		x<-which(locus1$V2==a[,])
		y<-which(locus1$V2==b[,])
		z<-which(locus1$V2==c[,])
		z1<-which(locus1$V2==d[,])

		allelei<-locus1[x,]
		allelej<-locus1[y,]
		allelek<-locus1[z,]
		allelel<-locus1[z1,]

		pklallelei<-array(dim=c(K))
		pklallelej<-array(dim=c(K))
		pklallelek<-array(dim=c(K))
		pklallelel<-array(dim=c(K))

		for(x in 1:K)	{
			pklallelei[x]<-allelei[1,x+2]
			pklallelej[x]<-allelej[1,x+2]
			pklallelek[x]<-allelek[1,x+2]
			pklallelel[x]<-allelel[1,x+2]
		}

		zi1=0
		zi2=0
		zj1=0
		zj2=0
		zk1=0
		zk2=0
		zl1=0
		zl2=0

		for(x in 1:K)	{
			zi1<-zi1+pklallelei[x]*indiv1.etaik[x]
			zi2<-zi2+pklallelei[x]*indiv2.etaik[x]
			zj1<-zj1+pklallelej[x]*indiv1.etaik[x]
			zj2<-zj2+pklallelej[x]*indiv2.etaik[x]
			zk1<-zk1+pklallelek[x]*indiv1.etaik[x]
			zk2<-zk2+pklallelek[x]*indiv2.etaik[x]
			zl1<-zl1+pklallelel[x]*indiv1.etaik[x]
			zl2<-zl2+pklallelel[x]*indiv2.etaik[x]
		}

		#Under S9
		d1<-0
		d2<-0
		d3<-0
		d4<-0
		d5<-0
		d6<-0
		d7<-0
		d8<-0
		d9<-(zi1*zj1*zk2*zl2+zk1*zl1*zi2*zj2)/2
	}

#trying flexmix - might be old code...
predictors[l-2,1]<-d1
predictors[l-2,2]<-d2
predictors[l-2,3]<-d3
predictors[l-2,4]<-d4
predictors[l-2,5]<-d5
predictors[l-2,6]<-d6
predictors[l-2,7]<-d7
predictors[l-2,8]<-d8
predictors[l-2,9]<-d9

}
# Setting initial values for optimization.
# You can change this if you like, if you know some of the related pairs. Else leave it be. 
ibds<-c(0.0,0.0,0.0,0.0,0.0,0.0,0.25,0.5,0.25)

# Do optimization using Rsolnp package
mle=solnp(ibds,fun=loglikibd,eqfun=eqn1,eqB=c(1),ineqfun=ineq1,ineqLB=c(0,0,0,0,0,0,0,0,0),ineqUB=c(1,1,1,1,1,1,1,1,1))

# Calculate the relatedness for MC2013WI method
relatedness1<-2*(mle$pars[1]+0.5*(mle$pars[3]+mle$pars[5]+mle$pars[7])+0.25*(mle$pars[8]))

# Append delta values to an output file - you can change name here as you like.
# adding bootstrap code to this
deltafile<-paste(filename1,".deltas",sep="")
cat(c(mle$pars[1],mle$pars[2],mle$pars[3],mle$pars[4],mle$pars[5],mle$pars[6],mle$pars[7],mle$pars[8],mle$pars[9]),"\n",file=deltafile,append=TRUE)

# Calculate relatedness for MC2013
relatedness2<-2*(mle$pars[7]*0.5+0.25*mle$pars[8])

# Loop controls - don't change!
#r<-i-decrementi
#s<-i+2-decrementj

# Build relatedness table
relat[g,1]<-sprintf("%d_%d",i,j)
relat[g,2]<-relatedness1
relat[g,3]<-relatedness2

# Loop controls - don't change!
#decrementj<-decrementj+2
#decrementi<-decrementi+2
g<-g+1
}
}
# Write final relatedness output to a file. You can change name as you like.
relatfile<-paste(filename1,".relat",sep="")
write.table(relat,relatfile)

