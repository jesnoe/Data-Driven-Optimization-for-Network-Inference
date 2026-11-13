# setwd("C:/Users/User/Documents/R")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(Matrix)
library(sf)
find_state_index <- function(d_var) {
  if (grepl("w", d_var)) {
    comma_index <- which(str_split(d_var, "")[[1]] == ",")
    dot_index <- which(str_split(d_var, "")[[1]] == ".")
    return(c(substr(d_var, 2, comma_index-1), substr(d_var, comma_index+1, dot_index-1)) %>% as.numeric)
  }
  
  if (grepl("S", d_var)) {
    return(rep(substr(d_var, 2, str_length(d_var)), 2) %>% as.numeric)
  }
}
find_state_index <- Vectorize(find_state_index)
states_data <- read.csv("Cocaine Network Optimization/states_data.csv") %>% as_tibble
price_t_test_summary <- read.csv("Cocaine Network Optimization/price t_test_summary greater (2009-2014).csv") %>% as_tibble

{# functions to create decision variables
  d_var_names <- function(letter, edges_) {
    if (edges_ == "N") result <- paste0(letter, "_", 1:N)
    if (edges_ == "U") result <- paste0(letter, "_", undir_edges$i, ",", undir_edges$j)
    if (edges_ == "E") result <- paste0(letter, "_", edges$from, ",", edges$to)
    return(result)
  }
  d_var_names <- Vectorize(d_var_names)
  
  create_d_vars <- function(d_var_idx_) {
    letters <- names(d_var_idx_)
    n_vars <- sapply(d_var_idx_, length)
    edges_<- ifelse(n_vars == N, "N",
                    ifelse(n_vars == K, "U", "E"))
    
    result <- d_var_names(letters, edges_) %>% unlist
    names(result) <- NULL
    return(result)
  }
}
{# parameters from the data
  states_data <- states_data %>% 
    filter(!is.na(med_death))
  N <- nrow(states_data) # Total number of nodes
  states_data$states_index <- 1:N
  
  seizure <- states_data$max_weight
  log_seizure <- log(seizure)
  O <- states_data$med_death
  
  # price p-value matrix an element in [i,j] is the p-value of t test between i and j
  price_pval_mat <- sparseMatrix(i = price_t_test_summary$state_index,
                                 j = price_t_test_summary$bordering_state_index,
                                 x = price_t_test_summary$p_value,
                                 dims = c(N, N))
  price_pval_vec <- price_t_test_summary$p_value
  
  # all possible edges (directed) of U.S. state-level network
  edges <- tibble(from = price_t_test_summary$state_index, to = price_t_test_summary$bordering_state_index) %>% mutate(varnames = paste0("w_", from, ",", to))
  E <- nrow(edges) # number of all directed edges
  lookup <- setNames(
    seq_len(nrow(edges)),
    edges$varnames
  )
  
  # all possible edges (undirected) of U.S. state-level network
  undir_edges <- edges %>%
    mutate(
      i = pmin(from, to),
      j = pmax(from, to)
    ) %>%
    distinct(i, j) %>%
    mutate(varnames = paste0("w_", i, ",", j),
           varnames_rev = paste0("w_", j, ",", i),
           edges_idx = lookup[varnames],
           edges_idx_rev = lookup[varnames_rev])
  K <- nrow(undir_edges)  # number of all undirected edges
}

### functions to create a constraint matrix (LHS starts with A_ and RHS starts with b_)
  # consumption (overdose deaths)
A_consumption_create <- function(d_var_idx_) {
  rows <- c()
  cols <- c()
  vals <- c()
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  for(e in 1:E){
    from <- edges$from[e]
    to   <- edges$to[e]
    
    # outgoing flow
    rows <- c(rows, from)
    cols <- c(cols, d_var_idx_$w[e])
    vals <- c(vals, -1)
    
    # incoming flow
    rows <- c(rows, to)
    cols <- c(cols, d_var_idx_$w[e])
    vals <- c(vals, 1)
  }
  
  # S variables
  rows <- c(rows, 1:N)
  cols <- c(cols, d_var_idx_$s)
  vals <- c(vals, rep(1, N))
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(N, n_d_vars_)
  )
  
  return(result)
}
b_consumption <- 100*O/sum(O)

  # lower bound by maximum seizure
A_seizure_create <- function(d_var_idx_) {
  rows <- c()
  cols <- c()
  vals <- c()
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  for(e in 1:E){
    from <- edges$from[e]
    to   <- edges$to[e]
    
    # incoming flow
    rows <- c(rows, to)
    cols <- c(cols, d_var_idx_$w[e])
    vals <- c(vals, 1)
  }
  
  # S variables
  rows <- c(rows, 1:N)
  cols <- c(cols, d_var_idx_$s)
  vals <- c(vals, rep(1, N))
  
  # epsilon variables
  rows <- c(rows, 1:N)
  cols <- c(cols, d_var_idx_$epsilon)
  vals <- c(vals, -100 * seizure/sum(seizure))
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(N, n_d_vars_)
  )
  
  return(result)
}
b_seizure <- rep(0, N)

  # constraints for flow indicator "t"
A_flow_indicator_create <- function(d_var_idx_) {
  rows <- c(rep(1:E, 2),
            rep((E+1):(2*E), 2))
  cols <- c(d_var_idx_$w, d_var_idx_$x,
            d_var_idx_$w, d_var_idx_$x)
  vals <- c(rep(-1, E), rep(100, E),
            rep(1, E), rep(-0.0001, E))
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(2*E, n_d_vars_)
  )
  
  return(result)
}
b_flow_indicator <- rep(0, 2*E)

  # source limit: total source income is 100
A_source_limit_create <- function(d_var_idx_) {
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  result <- sparseMatrix(
    i = rep(1, N),
    j = d_var_idx_$s,
    x = rep(1, N),
    dims = c(1, n_d_vars_)
  )
  
  return(result)
}
b_source_limit <- 100

  # source indicator: s_i = 1 if i imports cocaine outside of the U.s., and 0 otherwise
A_source_indicator_create <- function(d_var_idx_) {
  rows <- c(rep(1:N, 2),
            rep((N+1):(2*N), 2))
  cols <- c(d_var_idx_$s, d_var_idx_$t,
            d_var_idx_$s, d_var_idx_$t)
  vals <- c(rep(-1, N), rep(100, N),
            rep(1, N), rep(-0.0001, N))
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(2*N, n_d_vars_)
  )
  
  return(result)
}
b_source_indicator <- rep(0, 2*N)

  # x_ij = 1 if w_ij > 0, and 0 otherwise
A_w_complement_create <- function(d_var_idx_) {
  
  rows <- c(rep(1:K, 2),
            rep((K+1):(2*K), 2))
  cols <- c(d_var_idx_$w[undir_edges$edges_idx], d_var_idx_$z,
            d_var_idx_$w[undir_edges$edges_idx_rev], d_var_idx_$z)
  vals <- c(rep(-1, K), rep(100, K),
            rep(-1, K), rep(-100, K))
  n_d_vars_ <- d_var_idx_[[length(d_var_idx_)]] %>% max
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(2*K, n_d_vars_)
  )
  
  return(result)
}
b_w_complement <- c(rep(0, K), rep(-100, K))


# index list of decision variable letters (e.g., w, x, z, s, t, ...)
n_E_var <- 2
d_var_idx <- list(
  w = 1:E,
  x = (E+1):(2*E),
  z = (2*E+1):(2*E+K),
  epsilon = (n_E_var*E+K+1):(n_E_var*E+K+N),
  s = (n_E_var*E+K+N+1):(n_E_var*E+K+2*N),
  t = (n_E_var*E+K+2*N+1):(n_E_var*E+K+3*N)
)

d_vars <- create_d_vars(d_var_idx)
n_d_vars <- length(d_vars)

{ # constraints (LHS and RHS)
A_source_limit <- A_source_limit_create(d_var_idx)
A_consumption <- A_consumption_create(d_var_idx)
A_seizure <- A_seizure_create(d_var_idx)
A_flow_indicator <- A_flow_indicator_create(d_var_idx)
A_w_complement <- A_w_complement_create(d_var_idx)
A_source_indicator <- A_source_indicator_create(d_var_idx)


A_ <- rbind(A_source_limit, A_consumption, A_seizure, A_flow_indicator, A_w_complement, A_source_indicator)
rhs_ <- c(b_source_limit, b_consumption, b_seizure, b_flow_indicator, b_w_complement, b_source_indicator)
}

sense_ <- c(rep("=", 1+N), rep(">", nrow(A_)-N-1)) # set the inequality of constraints

price_dir_cost_ <- 1 - price_pval_vec # statistical penalty for price (1 - p-value)

{# set variable types and lower/upper bounds
vtype_ <- rep("C", n_d_vars)
vtype_[d_var_idx$x] <- "B"
vtype_[d_var_idx$z] <- "B"
vtype_[d_var_idx$t] <- "B"

lb_ <- integer(n_d_vars)
ub_ <- rep(Inf, n_d_vars)
ub_[d_var_idx$epsilon] <- 1
ub_[d_var_idx$x] <- 1
ub_[d_var_idx$z] <- 1
ub_[d_var_idx$t] <- 1
}

solve_gurobi <- function(obj, A, sense, vtype, rhs, lb=NULL, ub=NULL) {
  # solve the optimization problem with given constraints and parameters
  model_list <- list()
  model_list$obj        <- obj
  model_list$modelsense <- "min"
  model_list$A          <- A
  model_list$sense      <- sense
  model_list$vtype      <- vtype
  
  if (!is.null(lb)) model_list$lb <- lb
  if (!is.null(ub)) model_list$ub <- ub
  
  model_list$rhs        <- rhs
  
  model_result          <- gurobi(model_list)
  
  return(model_result)
}

price_only_optimal_network_map <- function(A_, rhs_, sense_, vtype_, lb_, ub_, flow_cost, price_dir_cost, import_cost) {
  # get optimal network with statistical price penalty
  
  obj_ <- integer(n_d_vars)
  obj_[d_var_idx$w] <- flow_cost
  obj_[d_var_idx$x] <- price_dir_cost
  obj_[d_var_idx$epsilon] <- -1
  obj_[d_var_idx$t] <- import_cost
  
  model_result <- solve_gurobi(obj = obj_, A = A_, sense = sense_, vtype = vtype_, rhs = rhs_, lb= lb_, ub=ub_)
  
  opt_sol <- tibble(deicision_vars=d_vars, optimal_sols=model_result$x)
  
  opt_network <- opt_sol[c(d_var_idx$w, d_var_idx$s),]
  max_flow <- opt_sol[d_var_idx$w,]$optimal_sols %>% max
  
  states_sf <- get_urbn_map(map = "states", sf = TRUE)
  opt_network_map <- left_join(states_sf %>%
                                 rename(state=state_name) %>%
                                 filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                                 mutate(centroid = st_centroid(geometry)),
                               states_data %>% select(state, med_death, states_index),
                               by="state") %>% 
    left_join(data.frame(states_index = 1:N, epsilon = opt_sol[d_var_idx$epsilon,]$optimal_sols), by = "states_index") %>% 
    arrange(states_index)
  
  
  opt_network_map_coords <- st_coordinates(opt_network_map$centroid)
  coords <- tibble(long_i = opt_network_map_coords[,1][c(edges$from, 1:N)],
                   lat_i = opt_network_map_coords[,2][c(edges$from, 1:N)],
                   long_j = opt_network_map_coords[,1][c(edges$to, 1:N)],
                   lat_j = opt_network_map_coords[,2][c(edges$to, 1:N)])
  
  opt_network <- bind_cols(opt_network, coords)
  result_flow_index <- 1:E
  result_source_index <- (E+1):nrow(opt_network)
  non_zero_flows <- opt_network[result_flow_index,]
  non_zero_flows$price_cost <- price_dir_cost
  non_zero_flows <- non_zero_flows %>% filter(optimal_sols > 0)
  non_zero_imports <- opt_network[result_source_index,] %>% filter(optimal_sols > 0)
  
  opt_network_map %>% ggplot() +
    geom_sf(aes(fill = -epsilon),
            color="black") +
    expand_limits(x=opt_network_map$long, y=opt_network_map$lat) +
    labs(x=NULL, y=NULL, title=NULL) +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_fill_continuous(palette = c("white", "grey40"),
                          limits = c(-1, 0)) +
    geom_sf(data = opt_network_map$centroid, color = "black", size = 2) +
    geom_segment(data=non_zero_flows,
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j,
                     color = price_cost),
                 linewidth = 0.2+2*non_zero_flows$optimal_sols/max_flow,
                 arrow=arrow(angle=10,
                             length=unit(0.5, "cm"),
                             type="closed")
    ) +
    geom_text(data=non_zero_imports,
              aes(x=long_i + 10000, y=lat_i -50000),
              label=non_zero_imports$optimal_sols %>% round(1),
              color="black") +
    scale_color_continuous(palette = c("#197C83", "#C00000"),
                           limits = c(0, 1)) +
    labs(x=NULL, y=NULL, title=NULL, fill="seizure\ncost", color="price\ncost") -> optimal_network_map
  
  result <- list(optimal_network_map = optimal_network_map,
                 optimal_solution = opt_sol)
  return(result)
}

# solve the model with flow_cost (alpha) = 1 and import_cost (beta) = 1
optimal_network_map1 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, import_cost = 1, price_dir_cost = price_dir_cost_)
optimal_network_map1$optimal_network_map

# create an image file of the optimal network
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 1 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map1$optimal_network_map, scale=1)

