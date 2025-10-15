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


## Task 3 dummy (delete this later)
# Day 1-10
S_vals <- c(99, 99, 98, 95, 92, 89, 85, 82, 78, 75,
            # Day 11-20
            71, 68, 63, 58, 53, 49, 45, 42, 38, 35,
            # Day 21-30
            32, 28, 26, 24, 22, 20, 19, 18, 17, 16,
            # Day 31-40
            15, 15, 14, 14, 14, 13, 13, 13, 13, 13,
            # Day 41-50
            13, 13, 12, 12, 12, 12, 12, 12, 12, 12,
            # Continue pattern to day 100...
            rep(12, 50))  # Days 51-100 stay at 12

E_vals <- c(0, 0, 1, 3, 5, 6, 8, 9, 10, 11,
            12, 13, 14, 15, 15, 14, 13, 12, 11, 10,
            9, 8, 7, 6, 5, 4, 3, 3, 2, 2,
            1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
            rep(0, 60))

I_vals <- c(1, 1, 1, 2, 3, 5, 7, 9, 12, 14,
            17, 19, 23, 27, 32, 37, 42, 46, 51, 55,
            59, 64, 67, 70, 73, 76, 78, 79, 81, 82,
            84, 84, 85, 86, 86, 87, 87, 87, 87, 87,
            rep(0, 60))

R_vals <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            87, 87, 88, 88, 88, 88, 88, 88, 88, 88,
            rep(88, 50))


beta <- rnorm(100, 1, 0.2)

# Combine into matrix
epi <- list(S = S_vals, 
            E = E_vals, 
            I = I_vals, 
            R = R_vals,
            beta = beta)
## end of task 3

## Task 4
par(mfcol=c(2,3),mar=c(4,4,1,1)) ## set plot window up for multiple plots
epi <- seir(bmu=7e-5,bsc=1e-7) ## run simulation
hist(epi$beta,xlab="beta",main="") ## beta distribution
plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N") ## S black
points(epi$E,col=4) ## E (blue)
points(epi$I,col=2) ## I (red)
points(epi$R,col=3) ## R (green)