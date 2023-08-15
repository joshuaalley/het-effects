// Joshua Alley 
// example heterogeneous effects model

data {
  int< lower = 0> N;
  int<lower = 0, upper = 1> y[N];
  int< lower = 0> C; // number of combat treatments
  int<lower = 0, upper = C> combat[N]; // combat group indicator
  int<lower = 0> L; // number of mil variables
  matrix[C, L] M; // military experience vars matrix
  int<lower = 0> P; // number of mil variables
  matrix[N, P] X; // ctontrols matrix
}


parameters {
  real alpha; // overall intercept
  vector[L] lambda; // het effects params
  vector[P] beta; // control effects params
  vector[C] theta_std; // stdized combat exp param
  real<lower = 0> sigma_theta; // combat exp impact var
}

transformed parameters{
  
    vector[C] theta; // stdized combat exp param
    vector[C] mu_theta; // impact of combat

    // last mu_theta is zero
    mu_theta = M * lambda;

    theta =  mu_theta + sigma_theta * theta_std;

}


model {
  
    // define priors
  alpha ~ std_normal(); 
  lambda ~ std_normal();
  beta ~ std_normal();
  theta_std ~ std_normal();
  sigma_theta ~ normal(0, .5);
  
  y ~ binomial_logit(1, alpha +
  theta[combat] + X * beta);

}

