
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


// template<class Type>
// Type objective_function<Type>::operator() ()
// {
//   DATA_VECTOR(Y);
//   DATA_MATRIX(X);
//   PARAMETER_VECTOR(beta);
//   int n = Y.size();
//   int m = X.cols();
//   Type foo = 0.0;
//   Type nll = 0.0;
//   
//   for(int i = 0; i < n; i++){
//     foo = 0.0;
//     for(int j = 0; j < m; j++){foo += beta(j)*X(i,j) ;}
//     nll += Y(i)*foo - log(1+exp(foo));
//   }
//   return -nll;
// }

//*** EXPLICIT MODEL STRUCTURE ***//

// template<class Type>
// Type sigmoid(Type x) {
//   if (x >= 0) return 1 / (1 + exp(-x));
//   else        return exp(x) / (1 + exp(x));
// }
// 
// 
// 
// template<class Type>
// Type objective_function<Type>::operator() ()
// {
//   DATA_VECTOR(Y);
//   DATA_MATRIX(X);
//   PARAMETER_VECTOR(beta);
//   int n = Y.size();
//   int m = X.cols();
//   Type foo = 0.0; 
//   Type nll = 0.0;
//   Type p;
//   
//   for(int i = 0; i < n; i++){
//     foo = 0.0;
//     for(int j = 0; j < m; j++){foo += beta(j)*X(i,j) ;}
//     p = sigmoid(foo) + 1e-12;
//     nll += Y(i)*log(p) + (1-Y(i))*log(1-p);
//   }
//   return -nll;
// }
