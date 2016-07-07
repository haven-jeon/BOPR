

prior_mean <- function(prior_prob, beta, num_features){
  bias_mean <- qnorm(prior_prob) * (beta ** 2 + num_features)
  return(bias_mean)
}



clip <- function(x, a, b) {
  a + (x-a > 0)*(x-a) - (x-b > 0)*(x-b)
}


