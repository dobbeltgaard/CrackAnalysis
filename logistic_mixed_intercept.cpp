#include <TMB.hpp>

template <class Type>
Type log_sigmoid(Type x) {
  if (asDouble(x) > 0) {
    return -log(1 + exp(-x));
  } else {
    return x - log(1 + exp(x));
  }
}

// Random intercept

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  DATA_IVECTOR(ID);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(u); 
  PARAMETER(log_sigma_u); 
  
  int n = Y.size();
  Type foo = 0.0;
  Type nll = 0.0;
  Type sigma_u = exp(log_sigma_u);
  
  vector<Type> eta = X * beta;
  for(int i = 0; i < n; i++){
    foo = eta(i);
    foo += u(ID(i));
    nll += Y(i) * log_sigmoid(foo) + (1 - Y(i)) * log_sigmoid(-foo);
    //nll += Y(i)*foo - log(1+exp(foo));
  }
  nll+=sum(dnorm(u, Type(0.0), sigma_u, true));
  return -nll;
}
