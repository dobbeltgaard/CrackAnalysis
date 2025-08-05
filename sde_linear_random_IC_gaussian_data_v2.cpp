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
  PARAMETER_VECTOR(dist_par_trans);
  PARAMETER(sigma_x_trans);
  Type sigma_x = softplus(sigma_x_trans);
  //PARAMETER(maint_slope); 
  
  Type dist_par1 = softplus(dist_par_trans(0));
  Type dist_par2 = softplus(dist_par_trans(1));
  //vector<Type> w = exp(dist_par1 + dist_par2*u); 
  //vector<Type> v = pnorm(u,Type(0),Type(1));  // Uniformly distributed variables (on [0,1])
  //vector<Type> Z = qgamma(v,dist_par1, dist_par2);
  
  Type obs_noise = softplus(obs_noise_trans);
  vector<Type> theta = softplus(theta_trans);
  
  Type nll = 0.0; 
  int prev_id = -1;
  Type time =0.0;

  for(int i = 0; i < X.size() - 1; i++){
    if (ID(i) != prev_id) { // if new crack appears 
      prev_id = ID(i);
      nll += dnorm(X(i), dist_par1, dist_par2, true);
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
  
  return -nll;
}

//Type X_maint = maintenance(X(i), maint_slope, M(i));



// #include <TMB.hpp>
// 
// 
// template <class Type>
// Type softplus(Type x) {
//   return log(1+exp(x));
// }
// 
// 
// template<class Type>
// Type maintenance(Type X, Type alpha, Type M) {
//   X += alpha * M;
//   //if (X < min_X) X = min_X;
//   //if (X > max_X) X = max_X;
//   return X;
// }
// 
// template<class Type>
// bool isNA(Type x){
//   return R_IsNA(asDouble(x));
// }
// 
// 
// template<class Type>
// Type objective_function<Type>::operator() ()
// {
//   DATA_VECTOR(Y);     
//   DATA_IVECTOR(ID);
//   DATA_VECTOR(M);
//   DATA_VECTOR(dTime); 
//   
//   PARAMETER(obs_noise_trans); 
//   PARAMETER_VECTOR(theta_trans);
//   //PARAMETER_VECTOR(u);
//   PARAMETER_VECTOR(X);
//   //PARAMETER_VECTOR(dist_par_trans);
//   PARAMETER(sigma_x_trans);
//   Type sigma_x = softplus(sigma_x_trans);
//   PARAMETER(maint_slope); 
//   
//   //Type dist_par1 = softplus(dist_par_trans(0));
//   //Type dist_par2 = softplus(dist_par_trans(1));
//   //vector<Type> w = exp(dist_par1 + dist_par2*u); 
//   //vector<Type> v = pnorm(u,Type(0),Type(1));  // Uniformly distributed variables (on [0,1])
//   //vector<Type> Z = qgamma(v,dist_par1, dist_par2);
//   
//   Type obs_noise = softplus(obs_noise_trans);
//   vector<Type> theta = softplus(theta_trans);
//   vector<Type> X_sol(Y.size());
//   
//   Type nll = 0.0;
//   int prev_id = -1;
//   int n = X.size(); 
//   Type dt = sum(dTime) / n; 
//   
//   Type crack_time = 0.0;
//   int j = int(0);
//   
//   
//   for(int i=0; i<n-1; i++){
//     if (ID(j) != prev_id) { // if new crack appears 
//       crack_time = 0.0; 
//       prev_id = ID(j); 
//     } else { //propagate forward 
//       nll += dnorm(X(i+1),X(i) + theta(0)*(theta(1) - X(i))*dt,sigma_x*sqrt(dt),true); 
//     }
//     crack_time += dt;
//     if(crack_time >= dTime(j)){
//       X(i) = maintenance(X(i), maint_slope, M(j));
//       if(!isNA(Y(j)) && Y(j) > 0){ nll += dnorm(Y(j), X(i), obs_noise, true); }
//       X_sol(j) = X(i);
//       j += 1; 
//       crack_time = 0.0;
//     }
//   }
//   
//   // for(int i=0; i<n-1; i++){
//   //   if (ID(j) != prev_id) { // if new crack appears 
//   //     crack_time = 0.0; 
//   //     prev_id = ID(j); 
//   //     } else { //propagate forward 
//   //       nll += dnorm(X(i+1),X(i) + theta(0)*(theta(1) - X(i))*dt,sigma_x*sqrt(dt),true); 
//   //     }
//   //   crack_time += dt;
//   //   if(crack_time >= dTime(j)){
//   //     X(i) = maintenance(X(i), maint_slope, M(j));
//   //     if(!isNA(Y(j)) && Y(j) > 0){ nll += dnorm(Y(j), X(i), obs_noise, true); }
//   //     X_sol(j) = X(i);
//   //     j += 1; 
//   //     crack_time = 0.0;
//   //   }
//   // }
//   
//   ADREPORT(X_sol);
// 
//   return -nll;
// }
// 
// 
// 
// //Type max_X = 140.0;  
// //Type min_X = 0.0001; 
// // 
// // Type time =0.0;
// // Type E_X = 0.0;
// // Type V_X = 0.0;
// // Type X = 0.0;
// // vector<Type> E_sol(Y.size());
// // vector<Type> V_sol(Y.size());
// // 
// // nll += sum(dnorm(u, Type(0.0), Type(1.0), true)); 
// 
// // for(int i = 0; i < Y.size(); i++){
// //   if (ID(i) != prev_id) {
// //     X = w(ID(i)); // initialize
// //     prev_id = ID(i);
// //     time = 0.0;
// //   } else {
// //     
// //     time += dTime(i);
// //     X = theta(1) + (X - theta(1)) * exp(-theta(0) * dTime(i));
// //   }
// //   X = maintenance(X, maint_slope, M(i), min_X, max_X);
// //   E_X = X;
// //   V_X = (sigma_x * sigma_x) / (2 * theta(0)) * (1 - exp(-2 * theta(0) * time));
// //   
// //   if(!isNA(Y(i)) && Y(i) > 0){
// //     Type total_var = V_X + obs_noise * obs_noise;
// //     Type total_sd  = sqrt(total_var);
// //     nll += dnorm(Y(i), E_X, total_sd, true);
// //   }
// //   
// //   E_sol(i) = E_X;
// //   V_sol(i) = V_X;
// // }
// 
