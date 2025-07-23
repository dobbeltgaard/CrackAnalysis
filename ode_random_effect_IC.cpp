
#include <TMB.hpp>



template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_IVECTOR(ID);
  DATA_VECTOR(dTime);
  PARAMETER(log_sigma_Y);
  //PARAMETER_VECTOR(u);
  PARAMETER(log_u_mu);
  //PARAMETER(log_u_sigma);
  PARAMETER_VECTOR(log_theta);
  
  Type sigma_Y = exp(log_sigma_Y);
  //vector<Type> u = exp(log_u); //transform to ensure positivity; u ~ logN(mu, sigma^2)
  //Type sigma_u = exp(log_u_sigma);
  vector<Type> theta = exp(log_theta);
  
  Type nll = 0.0;
  Type dt = 0.01; 
  int prev_id = -1;
  Type max_X = 140.0; 
  Type min_X = 0.0001;
  Type X = min_X;
  Type u_mu = exp(log_u_mu); 
  
  for(int i = 0; i < Y.size(); i++){
    if (ID(i) != prev_id) {
      X = u_mu;
      //X = u(ID(i));
      prev_id = ID(i);
    }
    int n_steps = CppAD::Integer(dTime(i) / dt);
    for(int k = 0; k < n_steps; k++){
      X += dt * theta(0) * pow(-X + theta(1), theta(2));
      //X += dt * theta(0) * pow(X, theta(1));
      //X += dt * (theta(0) * (theta(1) - X));
      //X += dt * (theta(0) * pow(theta(1) - X, theta(2)));
      if(X > max_X){ X = max_X; break; }
      if(X < 0){ X = min_X; break; }
    }
    nll += dnorm(Y(i), X, sigma_Y, true); 
  }
  //nll += sum(dnorm(u, u_mu, sigma_u, true));
  
  if (CppAD::isnan(nll)) {nll = -1e10;}
  
  return -nll;
}


