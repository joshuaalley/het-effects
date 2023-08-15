// Joshua Alley
// Bush and Prather Reanalysis 

data {
  int< lower = 0> N;
  real<lower = 0, upper = 4> y[N];
  int< lower = 0> T; // number of treatment groups
  int<lower = 0, upper = T> treat[N]; // treatment group indicator
  int<lower = 0> L; // number of het effects variables
  matrix[T, L] M; // het effects vars matrix
}


parameters {
  real alpha; // overall intercept
  real<lower = 0> sigma; // outcome variance
  vector[L] lambda; // het effects params
  vector[T] theta_std; // stdized het effect param
  real<lower = 0> sigma_theta; // het effect var
}

transformed parameters{
  
    vector[T] theta; // stdized treatment param
    vector[T] eta; // mean impact of treatment for each group

    // last mu_theta is zero
    eta = M * lambda;

    theta =  eta + sigma_theta * theta_std;

}


model {
  
    // define priors
  alpha ~ std_normal(); 
  lambda ~ normal(0, .5);
  theta_std ~ std_normal();
  sigma_theta ~ normal(0, .5);
  
  y ~ normal(alpha + theta[treat], sigma);

}
