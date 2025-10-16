# Aditya Sreekumar Achary - s2844915 - Part 3
# Sanjoi Sethi - s2891732 - Part 2 and 5
# Trisno Raynaldy Panjaitan - s2779061 - Part 1 and 4

# Repo Link: https://github.com/raypanjaitan/seir_statprog2

n <- 1000 ## number of population
hmax <- 5 ## maximum household size
h <- c() ## initiate the result variable, h
houseID <- 1 ## id of which house a person is in

while (length(h) < n) {
## function runs while the size of the result is less than number of population
## function assigns each person to a household, labeled by houseID
  
  r <- n - length(h) ## counter of the remaining population that has not been processed
  
  ## condition so the size of the household won't be greater than the remaining number of people
  if (r >= hmax) { 
    sz <- sample(1:hmax, 1) ## if true take sample from 1 to hmax
  } else {
    sz <- sample(1:r, 1) ## if false take sample from 1 to remaining
  }
  
  h <- c(h, rep(houseID, sz)) ## add household id to indices of the size of sz sample variable
  houseID <- houseID + 1 ## add 1 value to houseID so it can process the next iteration
}

h <- sample(h) ## randomize the h variable values

#Function to create the regular contacts network
get.net=function(beta, h, nc=15)
{
  #1. Matrix initialisation to store the network link between i-th and j-th 
  #person for the contact network model
  #2. If Household networks are encountered, they are set to 0 in the matrix.
  # (Catering to ij link=0 where i=j)
  #3 Else, the network is created using the sociability formula
  #4 Using the probability to sample from a Bernoulli distribution to create a 
  #link between people i & j. (1 denotes link and 0 denotes no link)
  #5 Creating a n-dimension list that stores the indices wherever a 1 is logged
  #in the matrix across rows.
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

nseir <- function(beta, h, alink, alpha = c(0.1, 0.01, 0.01), 
                  delta = 0.2, gamma = 0.4, nc = 15, nt = 100, pinf = 0.005){
  
  n = length(beta)
  beta_bar = mean(beta)
  
  ## Initial population state: 0=S, 1=E, 2=I, 3=R
  x = rep(0, n)
  ni = round(n * pinf) # infecting a proportion of initial population
  init_inf = sample(1:n, ni) # sampling out of the whole population so we get their indexes
  x[init_inf] = 2 # assign initial infectors
  
  S <- E <- I <- R <- rep(0, nt) # initializing the states
  time <- 1:nt
  
  ## Initial counts (day 1)
  S[1] <- n - ni
  E[1] <- 0
  I[1] <- ni
  R[1] <- 0
  
  
  for(i in 2:nt){
    
    u = runif(n)
    
    ## Step 1: I → R (recovery)
    x[x == 2 & u < delta] <- 3
    
    ## Step 2: E → I (becoming infectious)
    x[x == 1 & u < gamma] <- 2
    
    ## Step 3: S → E (new exposures due to infection)
    infectious <- which(x == 2)
    
    if(length(infectious) > 0){
      for (infector in infectious){
        inf_house = h[infector] # House of the infector
        susceptible <- which(x == 0) # Returns only true indices for susceptible
        
        inf_house = h[infector] # House of the infector
        susceptible <- which(x == 0) # Returns only true indices for susceptible
        
        #Transmitted within household
        prob_house <- rep(0, length(susceptible))
        prob_house[h[susceptible] == inf_house] <- alpha[1]
        
        #Transmitted from regular networks
        prob_reg <- rep(0, length(susceptible))
        prob_reg[susceptible %in% alink[[infector]]] <- alpha[2]
        
        #Random mixing
        prob_random <- alpha[3] * nc * beta[infector] * beta[susceptible] / 
          (beta_bar^2 * (n - 1))
        
        #Calculating the Probability for infection considering all the 3 independent events
        prob_infect = 1 - (1 - prob_house) * (1 - prob_reg) * (1 - prob_random)
        
        # Simulating the infection for all susceptible people 
        infection_sim = rbinom(length(susceptible), 1, prob_infect)
        infected = susceptible[infection_sim==1]
        
        #Assigning new Infected state
        x[infected]<-1
      }
    }
    ## Update counts for each state
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  return(list(
    S = S,   # Susceptible count per day
    E = E,   # Exposed count per day
    I = I,   # Infectious count per day
    R = R,   # Recovered count per day
    t = time # Time (days)
  ))
}

## plot the dynamics of the population by states
seirPlot <-function(epi, title){
  plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N",col=1,main=title) ## set
  #the maximum size of graph, label, put Susceptible data to the plot (black)
  points(epi$E,col=4) ## put Exposed data into the graph (blue)
  points(epi$I,col=2) ## put Infected data into the graph (red)
  points(epi$R,col=3) ## put Recovered data into the graph (green)
  legend(70,800,legend = c("Susceptible", "Exposed", "Infected", "Recovered"),
         col=c(1,4,2,3),pch=1) #Legend for more info on the graph
}

beta=runif(n, min=0, max=1) #Drawing the sociability parameter from a uniform 
#distribution since the probability of a person catching the disease is variable
alink=get.net(beta, h, nc=15) #Contact model
nseirResult1=nseir(beta, h, alink) #Calling the function and storing the result
#with standard parameters
nseirResult2=nseir(beta, h, alink, alpha = c(0, 0, 0.04)) #Storing the result
#only with the random mixing model
beta_new=rep(mean(beta), n) #Calculating the new beta vector with the mean of 
#the values of the previous beta
alink=get.net(beta_new, h, nc=15) #Contact model with new beta
nseirResult3=nseir(beta_new, h, alink) #Storing the result using the new beta 
#vector and standard parameters
nseirResult4=nseir(beta_new, h, alink, alpha = c(0, 0, 0.04)) #Storing the 
#result only with the random mixing model

#Plotting all the 4 results
par(mfrow = c(2,2), mar=c(4,4,1,1)) # Set plot window up for multiple plots
seirPlot(nseirResult1, "Result: 1")
seirPlot(nseirResult2, "Result: 2")
seirPlot(nseirResult3, "Result: 3")
seirPlot(nseirResult4, "Result: 4")

#Result Comments:
#In the first model where we are considering all 3 sociability parameters, it is
#observed that the population exposed and infected gradually over time. Whereas
#in the second model, we can see a higher peak for the infected and exposed
#states. From this, we can conclude that the epidemic is longer in the first 
#model as compared to the second one because the household and regular contact 
#network introduces the factor of local clusters, because of the which the
#number of people each individual can infect gets limited.
#When the beta value is put to a constant in the third and the fourth model, the 
#variability gets removed. Because of this, the severity of the epidemic gets
#overestimated, resulting in more number of people being exposed and infected in
#the population, which can be observed when we compare the 1st model vs 3rd and
#2nd model vs 4th.
#Out of the 4 models, the epidemic is largest when the beta becomes constant 
#and only random mixing is a sociability factor.