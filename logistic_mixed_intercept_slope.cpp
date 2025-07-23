#include <TMB.hpp>

template <class Type>
Type log_sigmoid(Type x) {
  if (asDouble(x) > 0) {
    return -log(1 + exp(-x));
  } else {
    return x - log(1 + exp(x));
  }
}

// Random intercept + slopes

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  DATA_IVECTOR(ID);
  PARAMETER_VECTOR(beta);
  PARAMETER_MATRIX(u);
  PARAMETER_VECTOR(log_sigma_u); 
  
  int n = Y.size();
  int m = X.cols();
  Type foo = 0.0;
  Type nll = 0.0;
  vector<Type> sigma_u = exp(log_sigma_u);
  
  vector<Type> eta = X * beta;
  for(int i = 0; i < n; i++){
    foo = eta(i);
    for(int j = 0; j < m; j++){foo += u(ID(i),j) * X(i, j);}
    nll += Y(i) * log_sigmoid(foo) + (1 - Y(i)) * log_sigmoid(-foo);
  }
  for(int g = 0; g < u.rows(); g++){
    for(int j = 0; j < m; j++){ nll += dnorm(u(g, j), Type(0.0), sigma_u(j), true); }
  }
  
  return -nll;
}

