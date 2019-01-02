library(polynom)


# E (days before expiration)
days <- 30
# r (annual interest rate)
rate <- 0.01
# S (stock price at time 0)
stock <- 100.0
# h0, b0, b1, b2, c
h0 <- 0.010469
b0 <- 0.000006575
b1 <- 0.9
b2 <- 0.04
ccc <- 0.0
# X (strike price)
strike <- 100.0
# n1 (number of partitions per day)
n1 <- 2
# n2 (number of variances per node)
n2 <- 2


get_eta <- function(h2, rate, gamma, n1){
  eta <- ceiling(h2^0.5/gamma)
  while (TRUE){
    low_bound <- abs(rate - (h2/2)) / (2*eta*gamma*(n1^0.5))
    up_bound <- min(1 - low_bound, 0.5)
    mid <- h2 / (2*eta*eta*gamma*gamma)
    
    if (low_bound <= mid && mid <= up_bound) return(eta)
    else if (mid < low_bound){
      print("Valid eta does not exist!")
      return(0)
    }
    eta <- eta + 1
  }
}


get_next_h2 <- function(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n){
  eps <- (l*eta*gamma_n-rate+h2/2) / (h2^0.5)
  return(b0+b1*h2+b2*h2*((eps-ccc)^2))
}


get_pu_pm_pd <- function(h2, eta, gamma, n1, rate){
  tmp1 <- h2 / (eta*eta*gamma*gamma)
  tmp2 <- (rate-h2/2) / (2*eta*gamma*(n1^0.5))
  return(c(tmp1*0.5+tmp2, 1-tmp1, tmp1*0.5-tmp2))
}


# Daily interest rate
rate <- rate/365
# Well, you know :)
gamma <- h0
gamma_n <- h0 / (n1^0.5)
# Eta (jump parameter) tree
eta_tree <- vector("list", days)
# Variance tree
h2_tree <- vector("list", days+1)
h2_tree[[1]][["0"]] <- rep(h0^2, n2)
# names(h2_tree[1]) <- "1"
# Probability tree
p_tree <- vector("list", days)


for(i in 1:days){
  for(j in paste(sort(as.numeric(names(h2_tree[[i]]))))){
    eta_tree[[i]][[j]] <- rep(0, n2)
    p_tree[[i]][[j]] <- matrix(rep(0, (2*n1+1)*n2), ncol=2*n1+1, nrow=n2)
    for (k in 1:n2){
      eta <- get_eta(h2_tree[[i]][[j]][k], rate, gamma, n1)
      if (eta == 0) next
      
      eta_tree[[i]][[j]][k] <- eta
      
      pu_pm_pd <- get_pu_pm_pd(h2_tree[[i]][[j]][k], eta, gamma, n1, rate)
      pu <- pu_pm_pd[1]
      pm <- pu_pm_pd[2]
      pd <- pu_pm_pd[3]
      poly <- polynomial(coef=c(pd, pm, pu))^n1
      coefs <- coef(poly)
      for (l in -n1:n1) p_tree[[i]][[j]][k, l+n1+1] = coefs[n1+l+1]
    }
  }
  
  for (j in paste(sort(as.numeric(names(h2_tree[[i]]))))){
    for (k in 1:n2){
      eta <- eta_tree[[i]][[j]][k]
      if (eta == 0) next
      h2 <- h2_tree[[i]][[j]][k]
      
      for (l in -n1:n1){
        next_j <- paste(as.numeric(j) + eta*l)
        next_h2 <- get_next_h2(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n)
        if (!(next_j %in% paste(sort(as.numeric(names(h2_tree[[i+1]])))))) {
          h2_tree[[i+1]][[next_j]]=rep(next_h2, n2)}
        else{
          min_ <- h2_tree[[i+1]][[next_j]][1]
          h2_tree[[i+1]][[next_j]][1] <- min(next_h2, min_)
          max_ <- h2_tree[[i+1]][[next_j]][2]
          h2_tree[[i+1]][[next_j]][2] <- max(next_h2, max_)
        }
      }
    }
  }
  
  # Interpolation of variance (for n2 > 2)
  for (next_j in paste(sort(as.numeric(names(h2_tree[[i+1]]))))){
    min_ <- h2_tree[[i + 1]][[next_j]][1]
    max_ <- h2_tree[[i + 1]][[next_j]][2]
    for (k in 1:n2)
      h2_tree[[i + 1]][[next_j]][k] <- min_ + (k-1)*(max_ - min_) / (n2-1)
  }
}



# Pricing at the last day
put_tree <- vector("list", days+1)
for (j in paste(sort(as.numeric(names(tail(h2_tree, n=1)[[1]]))))){
  # put = max(stock * np.exp(h0 * j) - strike, 0.)
  put <- max(strike - stock * exp(gamma_n * as.numeric(j)), 0)
  put_tree[[length(put_tree)]][[j]] <- rep(put, n2)
}


# Backward induction
for (i in days:1){
  for (j in paste(sort(as.numeric(names(h2_tree[[i]]))))){
    put_tree[[i]][[j]] <- rep(put, n2)
    for (k in 1:n2){
      eta <- eta_tree[[i]][[j]][k]
      if (eta == 0) next
      h2 <- h2_tree[[i]][[j]][k]
      
      put = 0
      for (l in -n1:n1){
        next_j <- paste(as.numeric(j) + eta * l)
        next_h2<- get_next_h2(b0, b1, b2, ccc, h2, rate, l, eta, gamma_n)
        
        # Find the next (k, k+1) interval bounding next_h2.
        for (next_k in 1:(n2 - 1)){
          low <- h2_tree[[i + 1]][[next_j]][next_k]
          up <- h2_tree[[i + 1]][[next_j]][next_k + 1]
          if (low <= next_h2 && next_h2 <= up) break
        } 
        
        x <- ifelse((low-up) != 0, (next_h2-up)/(low-up), 0)
        put_ <- x * put_tree[[i + 1]][[next_j]][next_k] + (1-x) * put_tree[[i + 1]][[next_j]][next_k + 1]
        put <- put + p_tree[[i]][[j]][k, l+n1+1]*put_
      }
      # American put pricing
      exercise <- strike - stock * exp(gamma_n * as.numeric(j))
      put_tree[[i]][[j]][k] <- max(put/exp(rate), exercise)
    }
  }
}

cat("Price: ",  put_tree[[1]][["0"]][1])

