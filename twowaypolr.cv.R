polr2cv=function(x, y, n=200, method=c("logistic", "probit", "cloglog", "multi"), flevels=levels(y))
{
	#THIS FUNCTION DOES TWO-WAY CROSS-VALIDATION n TIMES TO ASSESS PREDICTIVE ACCURACY BETWEEN TWO SUBSETS OF EQUAL SIZE 	#OF A PROPORTIONAL ODDS MODEL OR A MULTINOMIAL LOGIT MODEL. THIS DOES THIS FOR A FACTORED RESPONSE y AND A SINGLE 		#PREDICTOR x. THE LEVELS OF THE FACTOR y CAN OPTIONALLY BE SPECIFIED, IF IT DOES NOT CORRESPOND WITH DEFAULT.
	require(MASS); require(nnet)
	predictive_acc=NULL
	#CHECKS TO SEE IF LENGTH OF X AND Y ARE EQUAL
	if (length(x)!=length(y))
	{
		stop("lengths of arguments are not equal")	
	}
	#FACTORS AN UNFACTORED ARGUMENT
	if (is.factor(y)==FALSE)
	{
		flevels=levels(factor(y))
	}
	lex=length(x)
	#DOES TWO-WAY CROSS-VALIDATION ON THE RESPECTIVE DATA SUBSETS
	for (i in 1:n)
	{
		k=(2*i-1); j=2*i
		checkdex=NULL
		checknot=NULL
		dex=sample(lex, lex/2)
		xdex=x[dex]; ydex=y[dex]
		xnot=x[-dex]; ynot=y[-dex]
		if (method=="multi")
		{
			#FACTORS Y-VALUES
			ydex=factor(ydex, ordered=FALSE)
			ynot=factor(ynot, ordered=FALSE)
			#FITS MULTINOMIAL LOGIT MODELS
			mod_dex=multinom(ydex~xdex)
			mod_not=multinom(ynot~xnot)
			#PREDICTS CATEGORIES FOR OTHER HALF OF XVALUES NOT USED IN RESPECTIVE MODELS
			checkdex=predict(mod_dex, newdata=data.frame(xdex=xnot))
			checknot=predict(mod_not, newdata=data.frame(xnot=xdex))
			#FACTORS PREDICTED Y VALUES
			checkdex=factor(checkdex, levels=flevels) 
			checknot=factor(checknot, levels=flevels)
		}
		else 
		{
			#FACTORS Y-VALUES
			ydex=ordered(ydex)
			ynot=ordered(ynot)
			#FITS CUMULATIVE LOGIT (PROPORTIONAL ODDS) MODELS
			mod_dex=polr(ydex~xdex, method=method)
			mod_not=polr(ynot~xnot, method=method)
			#PREDICTS CATEGORIES FOR OTHER HALF OF XVALUES NOT USED IN RESPECTIVE MODELS
			checkdex=predict(mod_dex, newdata=data.frame(xdex=xnot))
			checknot=predict(mod_not, newdata=data.frame(xnot=xdex))
			#FACTORS PREDICTED Y VALUES
			checkdex=ordered(checkdex, levels=flevels)
			checknot=ordered(checknot, levels=flevels)
		}
		#RETURNS THE PREDICTIVE ACCURACY (PERCENT OF CATEGORIES CORRECT) AND RETURNS A VECTOR OF THESE 2n VALUES
		predictive_acc[k]=mean(checkdex==ynot)
		predictive_acc[j]=mean(checknot==ydex)
	}
	return(prediction_errors)
}


riditize=function(x, latent=FALSE, generate=FALSE)
{
	#THIS FUNCTION DOES RIDIT SCORING (AGRESTI, 2010) FOR A SET OF ORDERED FACTORS
	#THE SCORES ASSIGNED ARE BASED ON CUMULATIVE PROBABILITIES. THAT IS, THE PROBABILITY EXHAUSTED
	#BY THE FIRST K CATEGORIES PLUS HALF OF THE PROBABILITY OF K+1. IT ALSO CAN REPLACE THE VALUES
	# OF X WITH THESE SCORES. IT CAN ALSO RETURN THE RESPECTIVE QUANTILES OF THE STANDARD NORMAL 
	#CDF THAT THE RIDIT VALUES CORRESPOND TO.
	
	#CHECKS TO SEE IF FUNCTION TAKES ORDERED FACTOR AS ARGUMENT
	if (is.ordered(x)==FALSE)
	{
		stop("must take ordered factor as argument")	
	}
	xlevels=levels(x); length_levels=length(xlevels)
	
	#CHECKS TO SEE IF FUNCTION IS AT LEAST THREE CATEGORIES.
	if (length_levels<3)
	{
		stop("ridit scoring requires 3 or more categories")	
	}
	
	#COMPUTES RIDIT SCORES
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
	#PRODUCES VALUES FROM NORMAL CDF FROM RIDITS.
	latentnormalridits=qnorm(ridits)
	#REPLACES ORIGINAL VECTOR VALUES WITH RIDIT SCORES
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
	#RETURNS RIDIT VALUES.
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