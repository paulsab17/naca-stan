//Stan model for sea ice test

data{
  int < lower = 1 > N; // Sample size
  vector[N] year; // Predictor
  vector[N] extent; // Outcome
}

parameters{
  real alpha; // Intercept
  real beta; // Slope (regression coefficients)
  real < lower = 0 > sigma; // Error SD
}

model {
  beta ~ normal(0,10); //set prior distribution
  extent ~ normal(alpha + year * beta , sigma);
}
generated quantities {
  real ext_rep[N];

 for (n in 1:N) {
   ext_rep[n] = normal_rng(year[n] * beta + alpha, sigma);
 }
 //generates predicted values for each data point with rng

}
// The posterior predictive distribution"
