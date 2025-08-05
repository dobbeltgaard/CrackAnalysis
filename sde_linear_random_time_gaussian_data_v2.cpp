#include <TMB.hpp>


template <class Type>
Type softplus(Type x) {
  return log(1+exp(x));
}


template<class Type>
Type maintenance(Type X, Type alpha, Type M) {
  X += alpha * M;
  //if (X < min_X) X = min_X;
  //if (X > max_X) X = max_X;
  return X;
}

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Y);     
  DATA_IVECTOR(ID);
  //DATA_VECTOR(M);
  DATA_VECTOR(dTime); 
  
  PARAMETER(obs_noise_trans); 
  PARAMETER_VECTOR(theta_trans);
  PARAMETER_VECTOR(X);
  PARAMETER_VECTOR(t);
  PARAMETER_VECTOR(dist_par_trans);
  PARAMETER(sigma_x_trans);
  Type sigma_x = exp(sigma_x_trans);
  //PARAMETER(maint_slope); 
  
  Type dist_par1 = (dist_par_trans(0));
  Type dist_par2 = exp(dist_par_trans(1));
  vector<Type> w = exp(dist_par1 + dist_par2*t); 
  //vector<Type> v = pnorm(t,Type(0),Type(1));  // Uniformly distributed variables (on [0,1])
  //vector<Type> w = qgamma(v,dist_par1, dist_par2);
  
  Type obs_noise = exp(obs_noise_trans);
  vector<Type> theta = exp(theta_trans);
  
  Type nll = 0.0; 
  int prev_id = -1;
  Type time =0.0;
  Type X0 = 0.0001; 
  
  nll += sum(dnorm(t, Type(0.0), Type(1.0), true)); 
  
  for(int i = 0; i < X.size() - 1; i++){
    if (ID(i) != prev_id) { // if new crack appears 
      prev_id = ID(i);
      time = w(ID(i));
      Type mu = theta(1) + (X0 - theta(1))*exp(-theta(0)*time); 
      Type var = (sigma_x * sigma_x) / (2 * theta(0)) * (1 - exp(-2 * theta(0) * time));
      nll += dnorm(X(i), mu, sqrt(var), true );
    } else if( ID(i) == ID(i+1) ){ // crack continues
      time = dTime(i);
      Type mu = theta(1) + (X(i) - theta(1))*exp(-theta(0)*time); 
      Type var = (sigma_x * sigma_x) / (2 * theta(0)) * (1 - exp(-2 * theta(0) * time));
      nll += dnorm(X(i+1), mu, sqrt(var), true );
    }
  }
  
  for(int i = 0; i < Y.size(); i++){ 
    if(!isNA(Y(i)) && Y(i) > 0){ 
      nll += dnorm(Y(i), X(i), obs_noise, true); 
    } 
  }
  
  ADREPORT(X);
  ADREPORT(w);
  
  return -nll;
}



















