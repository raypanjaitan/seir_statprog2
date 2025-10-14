n <- 10000 ## number of population
hmax <- 5 ## maximum household size

h <- c() ## initiate the result variable, h
houseID <- 1 ##id

while (length(h) < n) { ## function run while the size of the result is less than number of population
  r <- n - length(h) ## counter of the remaining population that has not been processed
  
  ## condition so the size won't be greater than the remaining
  if (r >= hmax) { 
    sz <- sample(1:hmax, 1) ## if true take sample from 1 to hmax
  } else {
    sz <- sample(1:r, 1) ## if false take sample from 1 remaining
  }
  
  h <- c(h, rep(houseID, sz)) ## add household id to indices of the size of sz variable
  houseID <- houseID + 1 ## add 1 value to houseID so it can process the next iteration
}

h <- sample(h) ## make the h variable values to be randomized

get.net=function(beta, h, nc=15)
{
  #Creating a matrix that stores the network link between i-th and j-th person,
  #which is used to create the contact network model. Wherever people are
  #from the same household, their link is set to 0. This caters to the same 
  #person having the link with himself to be zero as well. Else, the network is
  #created using the probability formula. The Bernoulli distribution takes this
  #probability as an input and uses it to create a link between people i & j. 
  #1 denotes link and 0 denotes no link. Finally, creating a n-dimension list
  #that stores the indices wherever a 1 is potted in the matrix across rows.
  n=length(beta) #Population size initialization in the function
  links=matrix(data=NA, nrow=n, ncol=n) #Matrix to store the links between 
  #people i & j
  beta_bar=mean(beta) #Mean of beta vector for the probability formula
  for(i in 1:n)
  { #Looping through the rows (Identifying Person i) in the matrix
    for (j in i:n)
    { #Looping through the rows (Identifying Person j) in the matrix
      if(h[i]==h[j])
      { #Setting the link probability=0 where people i & j belong to the same
        #household
        links[i,j]=links[j,i]=0
      }
      else
      {
        #sum(runif(n)<p)/n
        p=(nc*beta[i]*beta[j])/((beta_bar^2)*(n-1)) #Creating a probability acc.
        #to the sociability factor
        links[i,j]=links[j,i]=rbinom(1, 1, p) #Assigning 1/0 link between i and
        #j according to the probability calculated
      }
    }
  }
  alink=apply(links, 1, function(row){ which(row==1) }) #For each person i, 
  #noting the index if the sociability link is established, i.e., if the value 
  #of a particular column corresponding to the i-th person's row is 1
  return(alink)
}

beta=runif(n, min=0, max=1) #Drawing the sociability parameter from a uniform 
#distribution since the probability of a person catching the disease is variable
alink=get.net(beta, h, nc=15) 



nseir <- function(beta, h, alink, alpha = c(0.1, 0.01, 0.01), 
                  delta = 0.2, gamma = 0.4, nc = 15, nt = 100, pinf = 0.005){
  
  n = length(beta)
  beta_bar = mean(beta)
  
  x = rep(0, n)
  ni = n * pinf # infecting a proportion of initial population
  
  init_inf = sample(1:n, ni) # sampling out of the whole population so we get their indexes
  x[init_inf] = 2
  
  S <- E <- I <- R <- rep(0, nt) # initializing the states
  
  time <- 1:nt
  for(i in 2:nt){
    
    u = runif(n)
    
    ## Step 1: I → R (Recovery)
    x[x == 2 & u < delta] <- 3
    
    ## Step 2: E → I (Becoming infectious)
    x[x == 1 & u < gamma] <- 2
    
    ## Step 3: S → E (New exposures due to infection)
    infectious <- which(x == 2)
    
    if(length(infectious) > 0){
      for (infector in infectious){
        inf_hous = h[infector]
        susceptible <- which(x == 0)
        
        
        for (person in susceptible) {
          prob_infect <- 0
          
          ##next comes the 3 independent events which im a bit confused about how to code
        }
      }
    }
    ## Update counts for each state
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  return(data.frame(time, S, E, I, R))
}
