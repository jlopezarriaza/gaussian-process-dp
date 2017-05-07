################################################################################
#' @title 
#' @description 
#' @param X A number.
#' @param Y
#' @param mu

#' @return The sum of \code{x} and \code{y}.
LogPosterior<-function(params,X,Y,mu){
  r=abs(min(X)-max(X))
  phi=params[1:dim(X)[2]]
  tau=params[dim(X)[2]+1]
  var_GP=params[dim(X)[2]+2]
  
  if(tau >= var(Y) |var_GP <1e-5|tau <1e-5| var_GP>=var(Y) |max(phi<0)){
    LP = 10000
  }else{
    
    K_obs = CalculateCovariance(X,X,phi,tau);
    tau.prior = dbeta(tau / var(Y),1.1,1.1);
    var.gp.prior = dbeta(var_GP / var(Y),1.1,1.1);
    phi_sig = r * sqrt(pi)/sqrt(2);
    phi.prior = (sqrt(2)/(phi_sig*sqrt(pi)))*exp(-(phi)^2/(2*phi_sig^2));
    eye = diag(1,dim(K_obs))
    LP = -0.5 * log(2 * pi) - 
      0.5 * log(det(K_obs + var_GP * eye)) - 
      0.5 * t(Y - mu) %*% solve(K_obs + var_GP * eye) %*% (Y - mu) +
      log(tau.prior * var.gp.prior * prod(phi.prior))
    LP=-LP}
  
  return(LP)
}
################################################################################
#' @title 
#' @description 
#' @param X A number.
#' @param Y
#' @param l
#' @param tau
#' @return The sum of \code{x} and \code{y}.
CalculateCovariance <- function(X,Y,l,tau) {
  K = 1
  for (u in 1:ncol(X)) {
    K = K * exp(-.5 * sapply(1:nrow(Y), function(i)
      abs(X[,u] - Y[i,u]) ^ 2) * l[u])
  }
  K = K * tau
  return(K)
}

################################################################################
#' @title 
#' @description 
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.

jet.colors <-
  colorRampPalette(
    c(
      "#00007F", "blue", "#007FFF", "cyan",
      "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
    )
  )
jet.cols = jet.colors(10)
################################################################################
#' @title 
#' @description 
#' @param x A number.
#' @param r_ricker
#' @param k_ricker
#' @param competition.parameters
#' @return The sum of \code{x} and \code{y}.

generating.models <- function(x,
                              r_ricker,
                              k_ricker,
                              competition.parameters = 0) {
  N_species = length(x)
  if (N_species == 3) {
    competition.matrix = matrix(runif(n = N_species ^ 2,
                                      min = 0,
                                      max = 1), 
                                ncol = N_species,
                                nrow = N_species)
    competition.matrix[diag(N_species) == 1] = 0
  }
  if (N_species == 2) {
    competition.matrix = matrix(c(0, competition.parameters[1],
                                  competition.parameters[2], 0),
                                nrow = 2, byrow = TRUE)
  }
  if (N_species == 1) {
    competition.matrix = 0
  }
  X = matrix(NA,nrow = length(x))
  for (i in 1 : length(x)) {
    X[i] = x[i] * exp(r_ricker[i] * (1 - x[i] / k_ricker[i]))
  }
  X = X * exp(competition.matrix %*% X)
  return(X)
}
################################################################################
#' @title 
#' @description 
#' @param input.data
#' @return The sum of \code{x} and \code{y}.
################################################################################
fit.gaussian.process <- function(input.data) {
  gp.fit.outputs = list()
  number.species = dim(input.data)[2]
  grid.min = apply(input.data, MARGIN = 2, FUN = min) * 0 + 1e-5
  grid.max = apply(input.data, MARGIN = 2, FUN = max) * 1.25
  
  
  if (number.species == 1)
    prediction.grid.size = round(150 ^ (1 / number.species))
  if(number.species == 2)
    prediction.grid.size = round(15)
  
  x = matrix(ncol = dim(input.data)[2], nrow = prediction.grid.size)
  for (i in 1 : dim(input.data)[2]) {
    x[,i] = t(seq(from = grid.min[i],
                  to = grid.max[i],
                  length.out = prediction.grid.size))
  }
  state.grid = expand.grid(split(x,col(x)))
  state_grid_size = dim(state.grid)[1]
  state.grid = matrix(unlist(state.grid),
                      nrow = nrow(state.grid),
                      ncol = ncol(state.grid))
  embedding.dimension = 1
  
  tau = 1
  TS_length = dim(input.data)[1];
  LagMatrix = array(dim = c(TS_length - (embedding.dimension) * tau,
                            embedding.dimension + 1,
                            number.species))
  for (d in seq_len(number.species)) {
    for (i in 1 : (embedding.dimension + 1)) {
      lag.index <- (1 + (embedding.dimension - i + 1)) : (TS_length - (i - 1))
      LagMatrix[,i,d] = input.data[lag.index, d];
    }
  }
  obs = matrix(log(LagMatrix[,1,]),
               ncol = number.species)
  hist = c()
  for (d in 2 : (embedding.dimension + 1)) {
    hist = cbind(hist,LagMatrix[,d,])
  }
  prior.mean = log(hist)
  optimal.parameters = list()
  phi = list()
  tau_var = list()
  var_GP = list()
  K_obs = list()
  eye = list()
  for (d in seq_len(number.species)) {
    var_GP_initial = 0.005
    phi_initial = 0.25 * rep(1, embedding.dimension * number.species)
    tau_initial = 0.005
    initial.parameters = c(phi_initial, tau_initial, var_GP_initial)
    optimizer = optim(par = initial.parameters,
                      fn = LogPosterior,
                      X = hist,
                      Y = matrix(obs[,d],ncol = 1),
                      mu = prior.mean[,d],
                      method = "Nelder-Mead",
                      control = list(maxit = 1500))
    optimal.parameters[[d]] = optimizer$par
    phi[[d]] = optimal.parameters[[d]][1:(embedding.dimension * number.species)]
    tau_var[[d]] = optimal.parameters[[d]][(embedding.dimension * number.species) + 1]
    var_GP[[d]] = optimal.parameters[[d]][(embedding.dimension * number.species) + 2]
    K_obs[[d]] = CalculateCovariance(hist, hist, phi[[d]], tau_var[[d]])
    eye[[d]] = diag(1, dim(K_obs[[d]])) * var_GP[[d]];
  }
  
  harvest.dimensions = number.species
  if (number.species == 1)
    harvest.grid.size = 100
  if(number.species == 2)
    harvest.grid.size = 10
  harv_temp = matrix(ncol = harvest.dimensions,
                     nrow = harvest.grid.size)
  for (i in 1 : dim(input.data)[2]) {
    harv_temp[,i] = t(seq(from = 0, 
                          to = 1 - 1e-4,
                          length.out = harvest.grid.size))
  }
  harvest.grid = expand.grid(split(harv_temp, col(harv_temp)))
  harvest.grid.size = dim(harvest.grid)[1]
  forward.state.grid = list()
  for (h in seq_len(harvest.grid.size)) {
    forward.state.grid[[h]] = sweep(matrix(state.grid[,1 : harvest.dimensions],
                                           ncol = number.species),
                                    2,
                                    unlist(1 - harvest.grid[h,]),
                                    '*')
    if (dim(state.grid)[2] > harvest.dimensions) {
      forward.state.grid[[h]] = cbind(forward.state.grid[[h]], 
                                      state.grid[,(harvest.grid + 1) : dim(state.grid)[2]])
    }
  }
  MM = list()
  K_c = list()
  MU_temp = list()
  CC2 = list()
  T_i = list()
  CC = list()
  for (d in seq_len(number.species)) {
    CC[[d]] = matrix(NA, nrow = state_grid_size, ncol = harvest.grid.size)
    MM[[d]] = matrix(NA, nrow = state_grid_size, ncol = harvest.grid.size)
    K_c[[d]] = list()
    T_i[[d]] = solve(K_obs[[d]] + eye[[d]])
    MU_temp[[d]] = matrix(sapply(seq_len(harvest.grid.size),
                                 function(x) forward.state.grid[[x]][,d]),
                          nrow = nrow(state.grid))
    for (h in seq_len(harvest.grid.size)) {
      K_c[[d]][[h]] = CalculateCovariance(forward.state.grid[[h]],
                                          hist,
                                          phi[[d]],
                                          tau_var[[d]])
      MM[[d]][,h] = log(MU_temp[[d]][,h]) + K_c[[d]][[h]] %*% T_i[[d]] %*%
        (obs[,d] - prior.mean[,d])
      CC[[d]][,h] = diag(abs(tau_var[[d]]
                             - K_c[[d]][[h]] %*% T_i[[d]] %*% t(K_c[[d]][[h]])))
    }
  }
  
  transition.prob.matrix = array(NA, dim = c(state_grid_size, 
                                             state_grid_size, 
                                             harvest.grid.size))
  for (h in seq_len(harvest.grid.size)) {
    for (j in seq_len(state_grid_size)) {
      transition.prob.matrix[j,,h] = apply(X = 
                                             sapply(seq_len(number.species), 
                                                    function(x) 
                                                      dlnorm(state.grid[,x],
                                                             MM[[x]][j,h],
                                                             CC[[x]][j,h] ^ 0.5)),
                                           MARGIN = 1,
                                           FUN = prod)
      if (max(transition.prob.matrix[j,,h]) == 0) {
        dis = apply(sapply(seq_len(number.species),function(x)
          abs(log(
            state.grid[,x]
          ) - MM[[x]][j,h]) ^ 2),1,sum)
        
        ind = which(dis == min(dis))
        transition.prob.matrix[j,,h] = 0
        transition.prob.matrix[j,ind,h] = 1
      }
      
      transition.prob.matrix[j,,h] = transition.prob.matrix[j,,h] / sum(transition.prob.matrix[j,,h])
    }
  }
  
  gp.fit.outputs$obs = exp(obs)
  gp.fit.outputs$hist = hist
  gp.fit.outputs$mu = MM
  gp.fit.outputs$sig = CalculateCovariance(forward.state.grid[[1]], 
                                           forward.state.grid[[1]],
                                           phi[[1]],tau_var[[1]]) -
    K_c[[1]][[1]] %*% T_i[[1]] %*% t(K_c[[1]][[1]])
  gp.fit.outputs$state.grid = state.grid
  gp.fit.outputs$harvest.grid = harvest.grid
  
  gp.fit.outputs$phi = phi
  gp.fit.outputs$tau = tau_var
  gp.fit.outputs$var = var_GP
  gp.fit.outputs$transition = transition.prob.matrix
  gp.fit.outputs$harvest.dimensions = harvest.dimensions
  gp.fit.outputs$MM = MM
  gp.fit.outputs$forward.state.grid = forward.state.grid
  return(gp.fit.outputs)
}
################################################################################
#' @title 
#' @description 
#' @param gp.fit.output
#' @return The sum of \code{x} and \code{y}.
################################################################################
markov.decision.process <- function(gp.fit.output) {
  markov.decicion.process = list()
  forward.state.grid = gp.fit.output$forward.state.grid
  MM = gp.fit.output$MM
  harvest.grid = gp.fit.output$harvest.grid
  state.grid = gp.fit.output$state.grid
  number.species = dim(state.grid)[2]
  transition.prob.matrix = gp.fit.output$transition
  harvest.dimensions = gp.fit.output$harvest.dimensions
  R = list()
  harvest.grid.size = dim(harvest.grid)[1]
  state_grid_size = dim(state.grid)[1]
  for (d in 1 : harvest.dimensions) {
    R[[d]] = matrix(NA,nrow = state_grid_size,ncol = harvest.grid.size)
    for (h in seq_len(harvest.grid.size)) {
      R[[d]][,h] = state.grid[,d] * harvest.grid[h,d]
    }
  }
  
  
  epsilon = 0.1;
  delt = 0.925;
  tol = epsilon * (1 - delt) / delt;
  Max_Iter = 1000;
  N_weights = 1
  w_vals = matrix(runif(N_weights * number.species),ncol = number.species)
  w_vals = matrix(1,ncol = number.species)
  w_vals = w_vals / apply(w_vals,1,sum)
  
  for (w_ind in 1 : N_weights) {
    weight = w_vals[w_ind,]
    iter = 2
    Vp = matrix(NA,state_grid_size,Max_Iter)
    Vp[,1] = delt / (1 - delt) * matrix(1,state_grid_size,1)
    Vp_part = list()
    for (h in 1 : harvest.dimensions) {
      Vp_part[[h]] = matrix(NA,state_grid_size,Max_Iter)
      Vp_part[[h]][,1] = delt / (1 - delt) * matrix(1,state_grid_size,1)
    }
    DD = matrix(1,state_grid_size,Max_Iter)
    check_value_function = 1;
    check_harvest_function = 1;
    while (iter < Max_Iter && check_value_function > tol) {
      V_part = list()
      for (m in 1 : harvest.dimensions) {
        V_part[[m]] = matrix(NA,state_grid_size,harvest.grid.size)
        for (h in seq_len(harvest.grid.size)) {
          V_part[[m]][,h] = transition.prob.matrix[,,h] %*% Vp_part[[m]][,iter - 1]
        }
        V_part[[m]] = R[[m]] + V_part[[m]] * delt
      }
      V = 0
      for (m in 1 : harvest.dimensions) {
        V = V + V_part[[m]] * weight[m]
      }
      Vp[,iter] = apply(V,1,FUN = max)
      DD[,iter] = max.col(V)
      for (m in 1 : harvest.dimensions) {
        for (j in seq_len(state_grid_size)) {
          Vp_part[[m]][j,iter] = V_part[[m]][j,DD[j,iter]]
        }
      }
      check_value_function = max(abs(Vp[,iter] - Vp[,iter - 1]));
      iter = iter + 1;
    }
  }
  policy = list()
  dynamics = list()
  if (number.species == 2) {
    for (i in seq_len(number.species)) {
      policy[[i]] <-
        akima::interp(
          state.grid[,1], state.grid[,2], (harvest.grid[DD[,iter - 1],i]) * state.grid[,i],
          xo = seq(1e-2,max(state.grid[,1]),length.out = 100),
          yo = seq(1e-2,max(state.grid[,2]),length.out = 100),duplicate =
            number.time.steps
        )
      dynamics[[i]] <-
        akima::interp(
          state.grid[,1], state.grid[,2], transition.prob.matrix[,,1] %*% state.grid[,i],
          xo = seq(1e-2,max(state.grid[,1]),length.out =
                     100),
          yo = seq(1e-2,max(state.grid[,2]),length.out =
                     100),duplicate = number.time.steps
        )
    }
    markov.decicion.process$interpolated = list(policy = policy,dynamics = dynamics)
  }
  markov.decicion.process$harvest = harvest.grid[DD[,iter - 1],]
  markov.decicion.process$state.grid = state.grid
  return(markov.decicion.process)
}
################################################################################
#' @title 
#' @description 
#' @param mdp.results
#' @return The sum of \code{x} and \code{y}.
################################################################################

forward.iterations <- function(mdp.results){
  transition <- function(x) {
    return(which.min(apply((x - mdp.results$state.grid) ^ 2, 1 , sum)))
  }
  number.time.steps = 50
  r_ricker = mdp.results$R
  k_ricker = mdp.results$K
  V = mdp.results$V
  competition.parameters = mdp.results$competition.parameters
  num = mdp.results$num
  number.repetitions = 50
  X = array(dim = c(number.time.steps, num,number.repetitions))
  H = array(dim = c(number.time.steps, num,number.repetitions))
  harvest = matrix(data = NA, nrow = number.time.steps, ncol = num)
  harv = matrix(unlist(mdp.results$harvest),ncol=num)
  for(reps in 1 : number.repetitions){
    x = matrix(data = NA, nrow = number.time.steps + 1, ncol = num)
    x[1,] = k_ricker * runif(num, min = 1, max = 3)
    for (t in 1 : number.time.steps) {
      harvest[t,] = harv[transition(x[t,]),] * x[t,]
      x[t,] = x[t,] - harvest[t,]
      x[t + 1,] = generating.models(x[t,], 
                                    r_ricker, 
                                    k_ricker, 
                                    competition.parameters) * 
        exp(rnorm(num, 0, sd = V ^ 0.5))
    }
    X[,,reps] = x[1 : t,]
    H[,,reps] = harvest
  }
  
  return(list(X = X,
              H = H))
}
################################################################################
#' @title 
#' @description 
#' @param mdp.results
#' @return The sum of \code{x} and \code{y}.
################################################################################

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}