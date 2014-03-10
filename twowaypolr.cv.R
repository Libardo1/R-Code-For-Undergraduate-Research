polr2cv=function(x, y, n=200, method=c("logistic", "probit", "cloglog", "multi"), flevels=levels(y))
{
	if (missing(method)==TRUE)
	{
	method="logistic"
	}
	# THIS FUNCTION DOES TWO-WAY CROSS-VALIDATION n TIMES TO ASSESS PREDICTIVE ACCURACY BETWEEN TWO SUBSETS OF EQUAL SIZE 
	# OF A PROPORTIONAL ODDS MODEL OR A MULTINOMIAL LOGIT MODEL. THIS DOES THIS FOR A FACTORED RESPONSE y 
	# AND A SINGLE PREDICTOR x. THE LEVELS OF THE FACTOR y CAN OPTIONALLY BE SPECIFIED, 
	# IF IT DOES NOT CORRESPOND WITH DEFAULT.
	require(MASS); require(nnet)
	# SIMPLE FUNCTION THAT COMPUTES THE LENGTH OF THE LEVELS OF THE ARGUMENT IT TAKES
	ll=function(x) return(length(levels(factor(x))))
	predictive_acc=NULL; counts=NULL
	#CHECKS TO SEE IF LENGTH OF X AND Y ARE EQUAL
	if (length(x)!=length(y))
	{
		stop("lengths of arguments are not equal")	
	}
	# FACTORS AN UNFACTORED ARGUMENT
	if (is.factor(y)==FALSE)
	{
		flevels=levels(factor(y))
	}
	lex=length(x)
	# DOES TWO-WAY CROSS-VALIDATION ON THE RESPECTIVE DATA SUBSETS
	for (i in 1:n)
	{
		k=(2*i-1); j=2*i
		checkdex=NULL
		checknot=NULL
		#THIS SETS TO ZERO THE LENGTH OF THE LEVELS OF YDEX (LYD) AND YNOT (LYN)
		ylevels=ll(y); lyd=0; lyn=0; count=0
		while (lyd<ylevels | lyn<ylevels)
		{
			dex=sample(lex, lex/2)
			xdex=x[dex]; ydex=y[dex]
			xnot=x[-dex]; ynot=y[-dex]
			lyd=ll(ydex); lyn=ll(ynot)
			count=count+1
			if (count > 10)
			{
			break	
			}
		}
		counts[i]=count
		if (method=="multi")
		{
			# FACTORS Y-VALUES
			ydex=factor(ydex, ordered=FALSE)
			ynot=factor(ynot, ordered=FALSE)
			# FITS MULTINOMIAL LOGIT MODELS
			mod_dex=multinom(ydex~xdex)
			mod_not=multinom(ynot~xnot)
			# PREDICTS CATEGORIES FOR OTHER HALF OF XVALUES NOT USED IN RESPECTIVE MODELS
			checkdex=predict(mod_dex, newdata=data.frame(xdex=xnot))
			checknot=predict(mod_not, newdata=data.frame(xnot=xdex))
			# FACTORS PREDICTED Y VALUES
			checkdex=factor(checkdex, levels=flevels) 
			checknot=factor(checknot, levels=flevels)
		}
		else 
		{
			# FACTORS Y-VALUES
			ydex=ordered(ydex)
			ynot=ordered(ynot)
			# FITS CUMULATIVE LOGIT (PROPORTIONAL ODDS) MODELS
			mod_dex=polr(ydex~xdex, method=method)
			mod_not=polr(ynot~xnot, method=method)
			# PREDICTS CATEGORIES FOR OTHER HALF OF XVALUES NOT USED IN RESPECTIVE MODELS
			checkdex=predict(mod_dex, newdata=data.frame(xdex=xnot))
			checknot=predict(mod_not, newdata=data.frame(xnot=xdex))
			# FACTORS PREDICTED Y VALUES
			checkdex=ordered(checkdex, levels=flevels)
			checknot=ordered(checknot, levels=flevels)
		}
		# RETURNS THE PREDICTIVE ACCURACY (PERCENT OF CATEGORIES CORRECT) AND RETURNS A VECTOR OF THESE 2n VALUES
		predictive_acc[k]=mean(checkdex==ynot)
		predictive_acc[j]=mean(checknot==ydex)
	}
	list("acc"=predictive_acc, "counts"=counts)
}


riditize=function(x, latent=FALSE, generate=FALSE)
{
	# THIS FUNCTION DOES RIDIT SCORING (AGRESTI, 2010) FOR A SET OF ORDERED FACTORS.
	# THE SCORES ASSIGNED ARE BASED ON CUMULATIVE PROBABILITIES. 
	# THAT IS, THE PROBABILITY EXHAUSTED BY THE FIRST K CATEGORIES PLUS HALF OF THE PROBABILITY OF K+1. 
	# IT ALSO CAN REPLACE THE VALUES OF X WITH THESE SCORES. 
	# ALTERNATIVELY IT CAN RETURN THE RESPECTIVE QUANTILES OF THE STANDARD NORMAL CDF THAT THE RIDIT VALUES CORRESPOND TO.
	# CHECKS TO SEE IF FUNCTION TAKES ORDERED FACTOR AS ARGUMENT
	if (is.ordered(x)==FALSE)
	{
		stop("must take ordered factor as argument")	
	}
	xlevels=levels(x); length_levels=length(xlevels)
	
	# CHECKS TO SEE IF FUNCTION IS AT LEAST THREE CATEGORIES.
	if (length_levels<3)
	{
		stop("ridit scoring requires 3 or more categories")	
	}
	
	# COMPUTES RIDIT SCORES
	lex=length(x)
	counts=NULL
	for (i in 1:length_levels)
	{
		sump=0
		for (j in 1:lex)
		{
			if (xlevels[i]==x[j])
			{
			sump=sump+1		
			}
		}
		counts[i]=sump
	}
	pk=counts/sum(counts)
	ridits=NULL
	for (i in 1:length_levels)
	{
		if (i==1)
		{
			ridits[i]=(pk[i]/2)
		}
		else
		{
			ridits[i]=sum(pk[1:(i-1)])+(pk[i]/2)
		}
	}
	# PRODUCES VALUES FROM NORMAL CDF FROM RIDITS.
	latentnormalridits=qnorm(ridits)
	# REPLACES ORIGINAL VECTOR VALUES WITH RIDIT SCORES
	if (generate==TRUE)
	{
		riditx=NULL
		for (i in 1:lex)
		{
			for (j in 1:length_levels)
			{
				if (xlevels[j]==x[i] && latent==FALSE)
				{
					riditx[i]=ridits[j]
				}
				else if (xlevels[j]==x[i] && latent==TRUE)
				{
					riditx[i]=latentnormalridits[j]
				}
			}
		}
		return(riditx)
	}
	# RETURNS RIDIT VALUES.
	else
	{
		if (latent==TRUE)
		{	
			return(latentnormalridits)	
		}
		else
		{
			return(ridits)	
		}
	}
}

# THIS FUNCTION GENERATES A 3-LEVEL CATEGORICAL VARIABLE FROM THE PROPORTIONAL ODDS MODEL ACCORDING TO SOME VECTOR OF X
# VALUES, EITHER USING SPECIFIED PARAMETER VALUES OR RANDOMLY GENERATING FROM AN (ARBITRARILY CHOSEN) STANDARD LOGISTIC
# DISTRIBUTION AND ORDERING THEM. NOTE ONLY THIS IS ONLY FOR A MODEL WITH A SINGLE PARAMETER.
generatefrompolr3=function(x, b, a1, a2, random=FALSE)
{
	require(MASS)
	if (random==TRUE)
	{
		pine=rlogis(3)
		cone=sample(pine, 2)
		a2=max(cone); a1=min(cone)
		b=setdiff(pine, cone)
	}
	proresp=NULL
	ypro=NULL
	for (i in 1:length(x))
	{
		ypro[i]=rlogis(1, location=x[i]*b)
		if (ypro[i] <= a1)
		{
			proresp[i]=1
		} 
		else if (ypro[i] > a1 && ypro[i] <= a2)
		{
			proresp[i]=2
		} 
		else 
		{
			proresp[i]=3
		}
	}
	cat("A1: ", a1, "\n", "A2: ", a2, "\n", "B: ", b, "\n")
	return(proresp)
}

# THIS FUNCTION GENERATES A 3-LEVEL CATEGORICAL RESPONSE VARIABLE FROM THE MULTINOMIAL LOGIT MODEL ACCORDING TO SOME VECTOR
# X VALUES, EITHER USING SPECIFIED PARAMETER VALUES OR RANDOMLY GENERATING FROM AN (ARBITRARILY CHOSEN) STANDARD NORMAL
# DISTRIBUTION AND ORDERING THEM. NOTE ONLY THIS IS ONLY FOR A MODEL WITH A SINGLE PARAMETER.
generatefrommulti3=function(x, a2, a3, b2, b3, random=FALSE)
{
	require(nnet)
	if (random==TRUE)
	{
		pine=rnorm(4)
		b2=pine[1]; b3=pine[2]; a2=pine[3]; a3=pine[4]
	}
	denom_i=NULL
	mulresp=NULL
	for (i in 1:length(x))
	{
		denom_i[i]=1+exp(a2+x[i]*b2)+exp(a3+x[i]*b3)
		prob=c(1/(denom_i[i]), exp(a2+x[i]*b2), exp(a3+x[i]*b3))
		jhuh=rmultinom(1, 1, prob=prob)
		for (j in 1:length(jhuh))
		{
			#SINCE THE RANDOM MULTINOMIAL ONLY GENERATES A SINGLE VALUE IN ANY OF THE THREE CATEGORIES
			#THE BELOW CODE DETERMINES WHICH ONE OF THE J IT IS IN.
			if (jhuh[j]>0)
			{
			mulresp[i]=j
			}
		}
	}
	cat("A2: ", a2, "\n", "A3: ", a3, "\n", "B2: ", b2, "\n", "B3: ", b3, "\n")
	return(mulresp)
}
