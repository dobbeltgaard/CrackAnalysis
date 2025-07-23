
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  PARAMETER_VECTOR(beta);
  int n = Y.size();
  Type nll = 0.0;
  vector<Type> eta = X * beta;
  for(int i = 0; i < n; i++){nll += Y(i)*eta(i) - log(1+exp(eta(i))); }
  return -nll;
}
