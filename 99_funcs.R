

sigmoid = function(x){
  return(1/(exp(-x) +1))
}

softplus = function(x){
  return(log(exp(x)+1))
}

softplus_inv <- function(y) {
  log(exp(y) - 1)
}



maintenance = function(X, alpha, M, min_X, max_X){
  res = X + alpha*M
  if(res < min_X) res = min_X
  if(res > max_X) res = max_X
  return(res)
}

simulate_process_w_maintenance <- function(
    n_ids = 500,
    obs_per_id = 10,
    theta = c(2, 15),
    obs_dist = c("normal", "lognormal", "gamma"),
    obs_noise_trans = -1, # Interpreted as: SD for normal, SDlog for lognormal, CV for gamma
    dist_par = c(1.3, 0.7),
    dist = c("normal", "lognormal","weibull", "gamma"),
    dt = 0.001,
    max_X = 140,
    min_X = 0.0001,
    drift_type = c("linear1","linear2", "logistic", "power_down", "power_up"),
    init_mode = c("X0", "t_tilde"),
    fixed = NULL,
    seed = 123, 
    n_mnt = 500, 
    maint_slope = 0.5, 
    type = "ode"
) {
  set.seed(seed)
  obs_dist <- match.arg(obs_dist)
  dist <- match.arg(dist)
  drift_type <- match.arg(drift_type)
  init_mode <- match.arg(init_mode)
  N <- n_ids * obs_per_id
  
  # -- Latent Gaussian transformed to Weibull/Gamma --
  u <- rnorm(n_ids)
  v <- pnorm(u)
  w <- switch(
    dist,
    weibull   = qweibull(v, shape = dist_par[1], scale = dist_par[2]),
    gamma     = qgamma(v, shape = dist_par[1], scale = dist_par[2]),
    normal    = qnorm(v, mean = dist_par[1], sd = dist_par[2]),
    lognormal = qlnorm(v, meanlog = dist_par[1], sdlog = dist_par[2])
  )
  if (!is.null(fixed)) { if (length(fixed) == 1L) {w <- rep(fixed, n_ids);} else if (length(fixed) == n_ids) { w <- fixed; } }
  
  if(type == "ode"){ }
  else if (type == "sde"){}
  
  # -- ID and observation times --
  ID <- rep(1:n_ids, each = obs_per_id)
  dTime <- rep(NA_real_, N)
  for (id in 1:n_ids) {
    dTime[ID == id] <- c(0, runif(obs_per_id - 1, 0.1, 1))
  }
  
  # -- Map drift_type to integer code --
  ode_type <- switch(
    drift_type,
    linear1     = 1L,
    linear2     = 2L,
    power_down = 3L,
    power_up   = 4L,
    logistic   = 5L
  )
  
  # -- Drift function selection --
  drift <- switch(
    drift_type,
    linear1     = function(x, th) th[1] * x,
    linear2     = function(x, th) th[1] * (-x + th[2]),
    logistic   = function(x, th) th[1] * x * (1 - x / th[2]),
    power_down = function(x, th) th[1] * pmax(th[2] - x, 0)^th[3],
    power_up   = function(x, th) th[1] * pmax(x, 0)^th[2]*(1-x/max_X)
  )
  
  
  # -- Simulate process --
  X <- numeric(N)
  Y <- numeric(N)
  M = numeric(N)
  obs_noise <- exp(obs_noise_trans)
  
  idx_mnt = sample(1:N, n_mnt, TRUE)
  M[idx_mnt] = 1
  Y[idx_mnt] = NA_real_
  
  for (id in 1:n_ids) {
    idx <- which(ID == id)
    if (init_mode == "X0") { Xi <- w[id]}
    if (init_mode =="t_tilde") {
      Xi <- min_X
      warmup_steps <- floor(w[id] / dt)
      for (k in 1:warmup_steps) {
        Xi <- Xi + dt * drift(Xi, theta)
        if (Xi > max_X) { Xi <- max_X; break }
        if (Xi < 0)     { Xi <- min_X; break }
      }
    }
    
    for (j in seq_along(idx)) {
      i <- idx[j]
      
      if (dTime[i] > 0) {
        n_steps <- floor(dTime[i] / dt)
        for (k in 1:n_steps) {
          Xi <- Xi + dt * drift(Xi, theta)
          if (Xi > max_X) { Xi <- max_X; break }
          if (Xi < 0)     { Xi <- min_X; break }
        }
      }
      
      Xi = maintenance(Xi, maint_slope, M[i], min_X, max_X)
      X[i] <- Xi
      
      if(is.na(Y[i])) next; 
      
      if (obs_dist == "normal") {
        Y[i] <- rnorm(1, Xi, obs_noise)
      } else if (obs_dist == "lognormal") {
        Y[i] <- rlnorm(1, meanlog = log(Xi), sdlog = obs_noise)
      } else if (obs_dist == "gamma") {
        shape_Y <- 1 / (obs_noise^2)             # Variance = mean^2 / shape
        scale_Y <- Xi * obs_noise^2             # So that mean = Xi
        Y[i] <- rgamma(1, shape = shape_Y, scale = scale_Y)
      }
    }
  }
  
  sim_data <- list(
    Y = Y,
    M = M,
    ID = as.integer(ID - 1),
    dTime = dTime,
    X0 = w,
    ode_type = ode_type
  )
  
  
  
  return(list(
    sim_data = sim_data,
    #init = init,
    X = X,
    Y = Y,
    M = M,
    latent = list(u = u, v = v, w = w),
    params = list(theta = theta, drift_type = drift_type, init_mode = init_mode,maint_slope = maint_slope)
  ))
}



# --------  WITHOUT MAINTENANCE ------------
simulate_process <- function(
    n_ids = 500,
    obs_per_id = 10,
    theta = c(2, 15),
    obs_dist = c("normal", "lognormal", "gamma"),
    obs_noise_trans = -1, # Interpreted as: SD for normal, SDlog for lognormal, CV for gamma
    dist_par = c(1.3, 0.7),
    dist = c("normal", "lognormal","weibull", "gamma"),
    dt = 0.001,
    max_X = 140,
    min_X = 0.0001,
    drift_type = c("linear", "logistic", "power_down", "power_up"),
    init_mode = c("X0", "t_tilde"),
    fixed = NULL,
    seed = 123 
) {
  set.seed(seed)
  obs_dist <- match.arg(obs_dist)
  dist <- match.arg(dist)
  drift_type <- match.arg(drift_type)
  init_mode <- match.arg(init_mode)
  N <- n_ids * obs_per_id
  
  # -- Latent Gaussian transformed to Weibull/Gamma --
  u <- rnorm(n_ids)
  v <- pnorm(u)
  w <- switch(
    dist,
    weibull   = qweibull(v, shape = dist_par[1], scale = dist_par[2]),
    gamma     = qgamma(v, shape = dist_par[1], scale = dist_par[2]),
    normal    = qnorm(v, mean = dist_par[1], sd = dist_par[2]),
    lognormal = qlnorm(v, meanlog = dist_par[1], sdlog = dist_par[2])
  )
  if (!is.null(fixed)) { if (length(fixed) == 1L) {w <- rep(fixed, n_ids);} else if (length(fixed) == n_ids) { w <- fixed; } }
  
  
  # -- ID and observation times --
  ID <- rep(1:n_ids, each = obs_per_id)
  dTime <- rep(NA_real_, N)
  for (id in 1:n_ids) {
    dTime[ID == id] <- c(0, runif(obs_per_id - 1, 0.1, 1))
  }
  
  # -- Map drift_type to integer code --
  ode_type <- switch(
    drift_type,
    linear1     = 1L,
    linear2     = 2L,
    power_down = 3L,
    power_up   = 4L,
    logistic   = 5L
  )
  
  # -- Drift function selection --
  drift <- switch(
    drift_type,
    linear1     = function(x, th) th[1] * x,
    linear2     = function(x, th) th[1] * (-x + th[2]),
    logistic   = function(x, th) th[1] * x * (1 - x / th[2]),
    power_down = function(x, th) th[1] * pmax(th[2] - x, 0)^th[3],
    power_up   = function(x, th) th[1] * pmax(x, 0)^th[2]
  )
  
  
  # -- Simulate process --
  X <- numeric(N)
  Y <- numeric(N)
  obs_noise <- exp(obs_noise_trans)

  for (id in 1:n_ids) {
    idx <- which(ID == id)
    if (init_mode == "X0") {
      Xi <- w[id]
    } else {
      Xi <- min_X
      warmup_steps <- floor(w[id] / dt)
      for (k in 1:warmup_steps) {
        Xi <- Xi + dt * drift(Xi, theta)
        if (Xi > max_X) { Xi <- max_X; break }
        if (Xi < 0)     { Xi <- min_X; break }
      }
    }
    
    for (j in seq_along(idx)) {
      i <- idx[j]
      
      if (dTime[i] > 0) {
        n_steps <- floor(dTime[i] / dt)
        for (k in 1:n_steps) {
          Xi <- Xi + dt * drift(Xi, theta)
          if (Xi > max_X) { Xi <- max_X; break }
          if (Xi < 0)     { Xi <- min_X; break }
        }
      }
      
      X[i] <- Xi
      if (obs_dist == "normal") {
        Y[i] <- rnorm(1, Xi, obs_noise)
      } else if (obs_dist == "lognormal") {
        Y[i] <- rlnorm(1, meanlog = log(Xi), sdlog = obs_noise)
      } else if (obs_dist == "gamma") {
        shape_Y <- 1 / (obs_noise^2)             # Variance = mean^2 / shape
        scale_Y <- Xi * obs_noise^2             # So that mean = Xi
        Y[i] <- rgamma(1, shape = shape_Y, scale = scale_Y)
      }
    }
  }
  
  sim_data <- list(
    Y = Y,
    ID = as.integer(ID - 1),
    dTime = dTime,
    ode_type = ode_type
  )
  
  # if(is.null(fixed)){
  #   init <- list(
  #     theta_trans = log(exp(theta)+1) + rnorm(length(theta)),
  #     dist_par = dist_par + rnorm(length(dist_par), mean = rep(0,length(dist_par)), sd = rep(0.05,length(dist_par))), 
  #     obs_noise_trans = obs_noise_trans,
  #   )  
  # } else{
  #   init <- list(
  #     theta_trans = log(exp(theta)+1) + rnorm(length(theta)),
  #     obs_noise_trans = obs_noise_trans + rnorm(length(fixed)), 
  #     X0_trans = fixed + rnorm(length(fixed))
  #   )
  # }
  
  return(list(
    sim_data = sim_data,
    init = init,
    X = X,
    Y = Y,
    latent = list(u = u, v = v, w = w),
    params = list(theta = theta, drift_type = drift_type, init_mode = init_mode)
  ))
}
