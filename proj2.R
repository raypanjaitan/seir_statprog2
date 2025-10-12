n <- 100 ## number of population
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