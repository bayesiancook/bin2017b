
data <- read.table("data10000.txt",header=TRUE)
y <- as.vector(data)
l <- length(y[1,])
countmatrix <- y[,2:l]

npos <- 1000
nind <- 14

# m_i(a)
basecount <- array(1,dim=c(npos,4))

# n_ij^tot
indcount <- array(0,dim=c(npos,nind))

for (i in 1:npos)	{
	index <- 1
	for (j in 1:nind)	{
		for (a in 1:4)	{
			basecount[i,a] <- basecount[i,a] + countmatrix[i,index]
			indcount[i,j] <- indcount[i,j] + countmatrix[i,index]
			index <- index+1
		}
	}
}

pi <- array(0,dim=c(npos,4))
for (i in 1:npos)	{
	pi[i,] <- basecount[i,] / sum(basecount[i,])	
}


# define array for total log likelihood for 
# error rate 
eps <- 0.01

# define arrays for:

# best genotype ...
genotype <- array(0,dim=c(npos,nind))

# ... and corrresponding posterior probability score
score <- array(0,dim=c(npos,nind))

for (i in 1:npos)	{
	for (j in 1:nind)	{

		# store prior and likelihood for the 10 genotypes
		logprior <- array(0,dim=10)
		loglikelihood <- array(0,dim=10)

		# g is a running index, from 1 to 10
		g <- 1

		# double loop over all possible diploid genotypes
		for (a in 1:4)	{
			for (b in a:4)	{

				# compute prior[g] and likelihood[g]
				if (a == b)	{
					logprior[g] <- 2*log(pi[i,a])
					loglikelihood[g] <- countmatrix[i,(j-1)*4+a]*log(1-eps) + (indtot[i,j]-countmatrix[i,(j-1)*4+a])*log(eps/3);
				}
				else	{
					logprior[g] <- log(2*pi[i,a]*pi[i,b])
					totab <- countmatrix[i,(j-1)*4+a]+countmatrix[i,(j-1)*4+b]
					loglikelihood[g] <- log(0.5*(1-epsilon) + 0.5*epsilon/3)*totab  + log(epsilon/3)*(indcount[i,j]-totab)
				}

				# increment index running over 10 genotypes
				g <- g+1
			}
		}

		# joint probability is product of prior and likelihood
		# thus, sum in log
		logjointprob <- logprior + loglikelihood

		# get the maximum joint prob over 10 genotypes
		maxlogprob <- max(logjointprob);

		# offset the exponential to avoid numerical errors
		jointprob <- exp(logjointprob-maxlogprob)

		# now make the sum
		totprob <- sum(jointprob)

		# and normalize, to get posterior probabilities
		posterior <- jointprob / totprob;

		# log likelihood
		logl <- logl + log(totprob) + maxlogprob

		# get best genotype
		genotype[i,j] <- which.max(posterior)
		score[i,j] <- max(posterior)
	}
}



