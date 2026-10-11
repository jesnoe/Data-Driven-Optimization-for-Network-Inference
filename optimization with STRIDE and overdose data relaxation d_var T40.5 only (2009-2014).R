# setwd("C:/Users/User/Documents/R")
# Sys.setenv(LANG = "en_US.UTF-8")
# Sys.setenv(GRB_LICENSE_FILE = "C:/gurobi1203/win64/licenses/gurobi.lic")
library(readxl)
library(gurobi)
library(stringi)
library(urbnmapr)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(Matrix)
library(sf)

{
  stride <- read.csv("Cocaine Network Optimization/STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month) & !(State %in% c("AK", "HI")))
  stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
  stride$Potency <- as.numeric(stride$Potency)
  stride$Post.Price <- as.numeric(stride$Post.Price)
  # for (i in c(1,2,6,7)) {
  #   stride[,i] <- as.factor(stride[,i]) 
  # }
  stride$State <- ifelse(stride$State == "NB", "NE", stride$State)
  
  neighbor_xlsx <- read_xlsx("Cocaine Network Optimization/states neighborhood.xlsx")
  neighbor <- list()
  for (i in 1:length(neighbor_xlsx$State)) {
    state_i <- neighbor_xlsx$State[i]
    neighbor[[state_i]] <- strsplit(neighbor_xlsx$Bordering_States[i], ", ")[[1]]
    if (state_i == "Maryland") neighbor[[state_i]] <- c(neighbor[[state_i]], "District of Columbia")
    if (state_i == "Virginia") neighbor[[state_i]] <- c(neighbor[[state_i]], "District of Columbia")
  }
  neighbor$`District of Columbia` <- c("Maryland", "Virginia")
  neighbor$Michigan <- neighbor$Michigan[-3] # remove Minnesota from Michigan's neighbors
  neighbor$Minnesota <- neighbor$Minnesota[-2] # remove Michigan from Minnesota's neighbors
  
  # Cocaine Data
  cocaine <- stride %>%
    filter(Drug %in% c("COCAINE", "COCAINE HYDROCHLORIDE") & Seize.Year != "NULL" & Seize.Month != "NULL" & State != "NULL") %>% 
    mutate(state=as.factor(State),
           MethAcq=as.factor(MethAcq),
           Drug=as.factor(Drug),
           Potency=ifelse(Potency > 100, Potency/10, Potency),
           Potency=ifelse(Potency > 0 & Potency < 1, Potency*100, Potency),
           Seize.Year=as.numeric(Seize.Year),
           Seize.Month=as.numeric(Seize.Month),
           adjusted_price=Post.Price/(Nt.Wt*Potency/100)) %>% 
    mutate(adjusted_price=ifelse(Potency == 0, NA, adjusted_price)) %>% 
    filter(Seize.Year > 1999) %>% 
    left_join(states %>% 
                rename(state=state_abbv) %>% 
                select(state, state_name) %>% 
                unique,
              by="state") %>% 
    select(-State, -state) %>% 
    rename(state=state_name) %>% 
    relocate(state)
  
  overdose <- read_xlsx("Cocaine Network Optimization/Overdose deaths T40.5 only 1999-2020.xlsx")
  VSRR <- read.csv("Cocaine Network Optimization/VSRR_Provisional_Drug_Overdose_Death_Counts (2015-2023).csv") %>%
    as_tibble %>% 
    mutate(state=State,
           Data.Value=as.numeric(Data.Value)) %>% 
    left_join(states %>% 
                rename(state=state_abbv) %>% 
                select(state, state_name) %>% 
                unique,
              by="state") %>% 
    select(-State, -state) %>% 
    rename(state=state_name) %>% 
    relocate(state)
  VSRR_cocaine <- VSRR %>% 
    filter(Indicator == "Cocaine (T40.5)") %>% 
    group_by(state, Year) %>% 
    summarize(overdose_death=sum(Data.Value, na.rm=T))
  
  years1 <- 2000:2003
  years2 <- 2004:2008
  years3 <- 2009:2014
  
  price_t_test_summary <- read.csv("Cocaine Network Optimization/price t_test_summary greater (2009-2014).csv") %>% as_tibble
  purity_t_test_summary <- read.csv("Cocaine Network Optimization/purity t_test_summary less (2009-2014).csv") %>% as_tibble
}

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

# no D.C. data
period <- years3
price_period <- cocaine %>%
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
  group_by(state) %>% 
  summarise(avg_price=mean(adjusted_price, na.rm=T),
            med_price=median(adjusted_price, na.rm=T))

cocaine_prices <- cocaine %>% 
  filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000)

seizure_period <- cocaine %>%
  filter(!is.na(state) & Seize.Year %in% period) %>%
  group_by(state) %>% 
  summarise(max_weight=max(Nt.Wt, na.rm=T),
            max_cocaine_weight=max(Nt.Wt*Potency/100, na.rm=T))


names(overdose)[3] <- "deaths"
overdose_period <- overdose %>% 
  filter(!is.na(state) & year %in% period) %>%
  group_by(state) %>% 
  summarise(avg_death=mean(deaths, na.rm=T),
            med_death=median(deaths, na.rm=T))


states_data <- full_join(price_period, seizure_period, by="state") %>% 
  left_join(overdose_period, by="state") %>% 
  arrange(state) %>% 
  right_join(states %>% 
               rename(state=state_name) %>% 
               filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
               group_by(state) %>% 
               summarise(long=mean(long), lat=mean(lat)), by="state") %>% 
  mutate(states_index=1:length(state))

states_data$med_death <- ifelse(is.na(states_data$med_death), 0, states_data$med_death)  
states_data$max_weight <- ifelse(is.na(states_data$max_weight), 0, states_data$max_weight)
# write.csv(states_data, "Cocaine Network Optimization/states_data.csv", row.names = F)

# importing states based on NDTA 2011: CA, TX, FL, NY
states_data %>% filter(state %in% c("California", "Texas", "Florida", "New York")) # 4, 9, 31, 42

# Optimization Model - Baseline (no relaxation on consumption and seizure, price_pval = 0.05)
# model with price onlyh

{
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
  
  states_data <- states_data %>% 
    filter(!is.na(med_death))
  N <- nrow(states_data) # Total number of nodes
  states_data$states_index <- 1:N
  
  seizure <- states_data$max_weight
  log_seizure <- log(seizure)
  O <- states_data$med_death
  
  price_pval_mat <- sparseMatrix(i = price_t_test_summary$state_index,
                                 j = price_t_test_summary$bordering_state_index,
                                 x = price_t_test_summary$p_value,
                                 dims = c(N, N))
  purity_pval_mat <- sparseMatrix(i = purity_t_test_summary$state_index,
                                  j = purity_t_test_summary$bordering_state_index,
                                  x = purity_t_test_summary$p_value,
                                  dims = c(N, N))
  
  price_pval_vec <- price_t_test_summary$p_value
  purity_pval_vec <- purity_t_test_summary$p_value
  
  edges <- tibble(from = price_t_test_summary$state_index, to = price_t_test_summary$bordering_state_index) %>% mutate(varnames = paste0("w_", from, ",", to))
  E <- nrow(edges)
  lookup <- setNames(
    seq_len(nrow(edges)),
    edges$varnames
  )
  
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
  K <- nrow(undir_edges)
  
  PO <- 100 # overdose proportional parameter
  # PS <- 5 # seizure proportionality parameter
  
}

{
  constraint_mat <- function(i, j, x, nrow, ncol) {
    sparseMatrix(
      i = i,
      j = j,
      x = x,
      dims = c(nrow, ncol)
    )
  }
  
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
  # b_seizure <- 100 * seizure/sum(seizure)
  
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
  
  solve_gurobi <- function(obj, A, sense, vtype, rhs, lb=NULL, ub=NULL) {
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
}

###### use model with z far below!!!!
## price d_var without w_complement
n_E_var <- 2
d_var_idx <- list(
  w = 1:E,
  x = (E+1):(2*E),
  epsilon = (n_E_var*E+1):(n_E_var*E+N),
  s = (n_E_var*E+N+1):(n_E_var*E+2*N),
  t = (n_E_var*E+2*N+1):(n_E_var*E+3*N)
)

d_vars <- create_d_vars(d_var_idx)
n_d_vars <- length(d_vars)

A_source_limit <- A_source_limit_create(d_var_idx)
A_consumption <- A_consumption_create(d_var_idx)
A_seizure <- A_seizure_create(d_var_idx)
A_flow_indicator <- A_flow_indicator_create(d_var_idx)
A_source_indicator <- A_source_indicator_create(d_var_idx)

A_ <- rbind(A_source_limit, A_consumption, A_seizure, A_flow_indicator,  A_source_indicator)
rhs_ <- c(b_source_limit, b_consumption, b_seizure, b_flow_indicator, b_source_indicator)
sense_ <- c(rep("=", 1+N), rep(">", nrow(A_)-N-1))
# sense_ <- c("=", rep(">", nrow(A_)-1))

{
  # price_pval_threshold <- 0.8
  price_violated_idx <- which(price_t_test_summary$mean_state <= price_t_test_summary$mean_bordering_state)
  price_pval_vec_filtered <- price_pval_vec
  price_pval_vec_filtered[price_violated_idx] <- 1
  obj_ <- integer(n_d_vars)
  # obj_[d_var_idx$epsilon] <- (1-price_pval_vec)
  obj_[d_var_idx$x[price_violated_idx]] <-  (1-price_pval_vec)[price_violated_idx]
  
  vtype_ <- rep("C", n_d_vars)
  vtype_[d_var_idx$x] <- "B"
  # vtype_[d_var_idx$z] <- "B"
  vtype_[d_var_idx$t] <- "B"
  
  lb_ <- integer(n_d_vars)
  ub_ <- rep(Inf, n_d_vars)
  ub_[d_var_idx$epsilon] <- 1
  ub_[d_var_idx$x] <- 1
  ub_[d_var_idx$t] <- 1
}


obj_[d_var_idx$w] <- 1
obj_[d_var_idx$epsilon] <- 1
obj_[d_var_idx$t] <- 1
price_dvar_model_b1 <- solve_gurobi(obj = obj_, A = A_, sense = sense_, vtype = vtype_, rhs = rhs_)

obj_[d_var_idx$w] <- 1
obj_[d_var_idx$t] <- 5
price_dvar_model_b5 <- solve_gurobi(obj = obj_, A = A_, sense = sense_, vtype = vtype_, rhs = rhs_)

obj_[d_var_idx$w] <- 1
obj_[d_var_idx$t] <- 10
price_dvar_model_b10 <- solve_gurobi(obj = obj_, A = A_, sense = sense_, vtype = vtype_, rhs = rhs_)

price_dvar_model_b1_sol <- tibble(variable = d_vars, opt_val = price_dvar_model_b1$x)
price_dvar_model_b5_sol <- tibble(variable = d_vars, opt_val = price_dvar_model_b5$x)
price_dvar_model_b10_sol <- tibble(variable = d_vars, opt_val = price_dvar_model_b10$x)

price_dvar_model_b1_sol[d_var_idx$epsilon,]
price_dvar_model_b5_sol[d_var_idx$epsilon,]
price_dvar_model_b10_sol[d_var_idx$epsilon,]


# (price_t_test_summary %>% filter(state %in% c("Alabama", "Florida", "Georgia")) %>% select(state, bordering_state, p_value))[c(1,2,5,6,7,8),]
# states_data %>% filter(state %in% c("Alabama", "Florida", "Georgia"))
# opt_sol$optimal_sols[lookup[c("w_1,9", "w_10,9")]]
# opt_sol$optimal_sols[lookup[c("w_9,1", "w_9,10")]]

# opt_sol$optimal_sols[d_var_idx$w] %>% sum
# opt_sol[d_var_idx$epsilon,]
# opt_sol[d_var_idx$epsilon[edges$from %in% c(1,9)],]
# opt_sol[d_var_idx$epsilon,] %>% filter(optimal_sols > 0)

price_only_optimal_network_map <- function(A_, rhs_, sense_, vtype_, lb_, ub_, flow_cost, price_dir_cost, import_cost) {
  # flow_cost = 1; price_dir_cost = price_dir_cost_; import_cost = 5
  
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

# price_dir_cost_ <- (1-price_pval_vec)
# price_dir_cost_ <- (1-price_pval_vec_filtered)
optimal_network_map1 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 1)
optimal_network_map1$optimal_network_map
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 1 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map1$optimal_network_map, scale=1)

optimal_network_map5 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 5)
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 5 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map5$optimal_network_map, scale=1)

optimal_network_map10 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 10)
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 10 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map10$optimal_network_map, scale=1)

## with w_complement
## price d_var
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

A_source_limit <- A_source_limit_create(d_var_idx)
A_consumption <- A_consumption_create(d_var_idx)
A_seizure <- A_seizure_create(d_var_idx)
A_flow_indicator <- A_flow_indicator_create(d_var_idx)
A_w_complement <- A_w_complement_create(d_var_idx)
A_source_indicator <- A_source_indicator_create(d_var_idx)

A_ <- rbind(A_source_limit, A_consumption, A_seizure, A_flow_indicator, A_w_complement, A_source_indicator)
rhs_ <- c(b_source_limit, b_consumption, b_seizure, b_flow_indicator, b_w_complement, b_source_indicator)
sense_ <- c(rep("=", 1+N), rep(">", nrow(A_)-N-1))
# sense_ <- c("=", rep(">", nrow(A_)-1))
price_dir_cost_ <- 1 - price_pval_vec
{
  # price_pval_threshold <- 0.8
  # price_violated_idx <- which(price_t_test_summary$mean_state <= price_t_test_summary$mean_bordering_state)
  # price_pval_vec_filtered <- price_pval_vec
  # price_pval_vec_filtered[price_violated_idx] <- 1
  obj_ <- integer(n_d_vars)
  # obj_[d_var_idx$epsilon] <- (1-price_pval_vec)
  # obj_[d_var_idx$x[price_violated_idx]] <-  (1-price_pval_vec)[price_violated_idx]
  
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

optimal_network_map1 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 1)
optimal_network_map1$optimal_network_map
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 1 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map1$optimal_network_map, scale=1)

optimal_network_map5 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 5)
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 5 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map5$optimal_network_map, scale=1)

optimal_network_map5$optimal_solution$optimal_sols[d_var_idx$w] %>% sum
optimal_network_map5$optimal_solution$optimal_sols[d_var_idx$t] %>% sum
-(optimal_network_map5$optimal_solution$optimal_sols[d_var_idx$epsilon]) %>% sum
(optimal_network_map5$optimal_solution$optimal_sols[d_var_idx$z]*price_dir_cost_) %>% sum


optimal_network_map10 <- price_only_optimal_network_map(A_, rhs_, sense_, vtype_, lb = lb_, ub = ub_, flow_cost = 1, price_dir_cost = price_dir_cost_, import_cost = 10)
ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks seizure=inflows import_cost = 10 (", period[1], "-", period[length(period)], ").png"),
       optimal_network_map10$optimal_network_map, scale=1)

optimal_network_map10$optimal_solution$optimal_sols[d_var_idx$w] %>% sum
optimal_network_map10$optimal_solution$optimal_sols[d_var_idx$t] %>% sum
-(optimal_network_map10$optimal_solution$optimal_sols[d_var_idx$epsilon]) %>% sum
(optimal_network_map10$optimal_solution$optimal_sols[d_var_idx$z]*price_dir_cost_) %>% sum

A_seizure_in_and_out_create <- function(d_var_idx_) {
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
    vals <- c(vals, 1)
    
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
A_seizure_in_and_out <- A_seizure_in_and_out_create(d_var_idx)

A_with_seizure_in_and_out <- rbind(A_source_limit, A_consumption, A_seizure_in_and_out, A_price, A_w_complement, A_source_indicator)
optimal_network_map2 <- price_only_optimal_network_map(A_with_seizure_in_and_out, rhs_, sense_, vtype_, flow_cost = 1, price_dir_cost = exp(price_pval_vec), import_cost = 1)
optimal_network_map2$optimal_network_map

# Cocaine seizure effectiveness rate (2012~2014): 42%, 3.7%, and 3.9% (3.93% on average)
A_seizure_ratio_est_create <- function(d_var_idx_) {
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
  
  result <- sparseMatrix(
    i = rows,
    j = cols,
    x = vals,
    dims = c(N, n_d_vars_)
  )
  
  return(result)
}

total_cocaine_flow <- sum(seizure)/0.0393
seizure_ratio_est <- seizure / total_cocaine_flow
b_source_limit_with_seizure_ratio_est <- 100 - sum(seizure_ratio_est)

A_seizure_ratio_est <- A_seizure_ratio_est_create(d_var_idx)

A_with_seizure_in_and_out <- rbind(A_source_limit, A_consumption, A_seizure_ratio_est, A_price, A_w_complement, A_source_indicator)
rhs_with_seizure_in_and_out <- c(b_source_limit_with_seizure_ratio_est, b_consumption, b_seizure, b_price, b_w_complement, b_source_indicator)
optimal_network_map3 <- price_only_optimal_network_map(A_with_seizure_in_and_out, rhs_with_seizure_in_and_out, sense_relaxed_consumption, vtype_, flow_cost = 1, price_dir_cost = exp(price_pval_vec), import_cost = 1)
optimal_network_map3$optimal_network_map
optimal_network_map3$optimal_solution

# Fixed importing states: CA, TX, FL, NY
# {
#   fixed_importing_states_index <- which(d_vars %in% c("x4", "x9", "x31", "x42"))
#   n_fixed_importing_states <- length(fixed_importing_states_index)
#   A_fixed_importing_states <- matrix(0, n_fixed_importing_states, n_d_vars)
#   for (i in 1:n_fixed_importing_states) {
#     A_fixed_importing_states[i, fixed_importing_states_index[i]] <- 1
#   }
#   b_fixed_importing_states <- rep(1, n_fixed_importing_states)
# }

## results with b_seizure_ratio = 100*seizure/sum(seizure)




# model with price and purity
{
  states_data <- states_data %>% 
    filter(!is.na(med_death))
  N <- nrow(states_data) # Total number of nodes
  states_data$states_index <- 1:N
  
  price_pval_mat <- matrix(0, N, N)
  purity_pval_mat <- matrix(0, N, N)
  
  for (i in 1:nrow(price_t_test_summary)) {
    state_index_i <- price_t_test_summary$state_index[i]
    botering_state_index_i <- price_t_test_summary$bordering_state_index[i]
    price_pval_i <- price_t_test_summary$p_value[i]
    price_pval_mat[state_index_i, botering_state_index_i] <- price_pval_i
    
    state_index_i <- purity_t_test_summary$state_index[i]
    botering_state_index_i <- purity_t_test_summary$bordering_state_index[i]
    purity_pval_i <- purity_t_test_summary$p_value[i]
    purity_pval_mat[state_index_i, botering_state_index_i] <- purity_pval_i
  }
  
  price_pval_vec <- c()
  purity_pval_vec <- c()
  d_vars <- c()
  for (i in 1:N) {
    bordering_states_index <- states_data %>%
      filter(state %in% neighbor[[states_data$state[i]]]) %>% pull(states_index)
    d_vars <- c(d_vars, paste0("w", i, ",", bordering_states_index, "."))
    
    for (j in bordering_states_index) {
      price_pval_vec <- c(price_pval_vec, price_pval_mat[i, j])
      purity_pval_vec <- c(purity_pval_vec, purity_pval_mat[i, j])
    }
    
  }
  
  x_ij <- gsub("w", "x", d_vars)
  y_ij <- gsub("w", "y", d_vars)
  d_vars <- c(d_vars, x_ij, y_ij)
  d_vars <- c(d_vars, paste0("S", 1:N))
  d_vars <- c(d_vars, paste0("t", 1:N))
  n_d_vars <- length(d_vars)
  n_w_ij <- (n_d_vars - 2*N) / 3
  
  PO <- 100 # overdose proportional parameter
  # PS <- 5 # seizure proportionality parameter
  epsilon_price <- 0.05 # price p-value to satisfy price directionality. The higher epsilon_price, the more relaxation
  epsilon_purity <- 0.05 # price p-value to satisfy price directionality. The higher epsilon_price, the more relaxation
}

{
  # seizure <- states_data$total_weight
  seizure <- states_data$max_weight
  log_seizure <- log(seizure)
  O <- states_data$med_death
  price <- states_data$med_price
  
  # Constraints
  # Conservation of flow/consumption
  A_consumption <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_consumption[i, which(str_match(d_vars[1:n_w_ij], paste0("w", i, ","))[,1] == paste0("w", i, ","))] <- -1
    A_consumption[i, which(str_match(d_vars[1:n_w_ij], paste0(",", i, "."))[,1] == paste0(",", i, "."))] <- 1
    A_consumption[i, which(d_vars == paste0("S", i))] <- 1
  }
  b_consumption <- 100*O/sum(O)
  
  # State seizure
  A_seizure <- matrix(0, N, n_d_vars)
  for (i in 1:N) {
    A_seizure[i, which(str_match(d_vars, paste0(",", i, "."))[,1] == paste0(",", i, "."))] <- 1
    A_seizure[i, which(d_vars == paste0("S", i))] <- 1
  }
  b_seizure <- 100*seizure/sum(seizure)
  b_seizure_log_ratio <- 100*log_seizure/sum(log_seizure)
  
  # Price indicator
  A_price_indicator <- matrix(0, n_w_ij, n_d_vars)
  for (i in 1:n_w_ij) {
    A_price_indicator[i, i] <- -1
    A_price_indicator[i, i + n_w_ij] <- 100
  }
  b_price_indicator <- rep(0, n_w_ij)
  
  # Purity indicator
  A_purity_indicator <- matrix(0, n_w_ij, n_d_vars)
  for (i in 1:n_w_ij) {
    A_purity_indicator[i, i] <- -1
    A_purity_indicator[i, i + n_w_ij] <- 100
  }
  b_purity_indicator <- rep(0, n_w_ij)
  
  # Source indicator
  A_source_indicator <- matrix(0, 2*N, n_d_vars)
  for (i in 1:N) {
    A_source_indicator[i, which(d_vars %in% c(paste0("S", i), paste0("t", i)))] <- c(-1, 100)
    A_source_indicator[i+N, which(d_vars %in% c(paste0("S", i), paste0("t", i)))] <- c(1, -0.0001)
  }
  b_source_indicator <- rep(0, 2*N)
  
  A_source_limit <- matrix(0, 1, n_d_vars)
  A_source_limit[1, grep("S", d_vars)] <- 1
  b_source_limit <- 100
}

# Fixed importing states: CA, TX, FL, NY
# {
#   fixed_importing_states_index <- which(d_vars %in% c("x4", "x9", "x31", "x42"))
#   n_fixed_importing_states <- length(fixed_importing_states_index)
#   A_fixed_importing_states <- matrix(0, n_fixed_importing_states, n_d_vars)
#   for (i in 1:n_fixed_importing_states) {
#     A_fixed_importing_states[i, fixed_importing_states_index[i]] <- 1
#   }
#   b_fixed_importing_states <- rep(1, n_fixed_importing_states)
# }

opt_fixed_pval_price_purity <- function(epsilon_price, epsilon_purity) {
  
  A_price_condition <- matrix(0, 2*n_w_ij, n_d_vars)
  for (i in 1:n_w_ij) {
    A_price_condition[i, i+n_w_ij] <- 1.1
    A_price_condition[i+n_w_ij, i+n_w_ij] <- -1.1
  }
  b_price_condition <- c((epsilon_price - price_pval_vec), -(epsilon_price - price_pval_vec + 1.1))
  
  A_purity_condition <- matrix(0, 2*n_w_ij, n_d_vars)
  for (i in 1:n_w_ij) {
    A_purity_condition[i, i+n_w_ij] <- 1.1
    A_purity_condition[i+n_w_ij, i+n_w_ij] <- -1.1
  }
  b_purity_condition <- c((epsilon_purity - purity_pval_vec), -(epsilon_purity - purity_pval_vec + 1.1))
  
  model_list <- list()
  model_list$obj        <- c(rep(0.01, n_w_ij), rep(0, 2*n_w_ij), rep(0, N), rep(1/49, N))
  # model_list$obj        <- c(rep(0.01, n_w_ij), rep(0, 2*n_w_ij), rep(0, N), rep(0, N))
  model_list$modelsense <- "min"
  model_list$A          <- rbind(A_consumption, A_seizure, 
                                 A_price_indicator, A_price_condition,
                                 A_purity_indicator, A_purity_condition,
                                 A_source_indicator, A_source_limit)
  
  model_list$sense      <- rep(">", nrow(model_list$A)) # initialize inequality >= for all
  model_list$sense[nrow(model_list$A)] <- "="
  
  model_list$vtype      <- c(rep("C", n_w_ij), rep("B", 2*n_w_ij), rep("C", N), rep("B", N)) # Continuous and binary
  model_list$lb         <- c(rep(0, n_d_vars)) # lower bounds: 0
  model_list$ub         <- c(rep(Inf, n_w_ij), rep(1, 2*n_w_ij), rep(100, N), rep(1, N)) # upper bounds
  model_list$rhs        <- c(b_consumption, b_seizure,
                             b_price_indicator, b_price_condition,
                             b_purity_indicator, b_purity_condition,
                             b_source_indicator, b_source_limit)
  model_result          <- gurobi(model_list, list(Method=-1))
  
  return(model_result)
}


epsilon_price_purity_results <- data.frame()
for (epsilon_price_ in seq(0.05, 1, by=0.05)) {
  for (epsilon_purity_ in seq(0.05, 1, by=0.05)) {
    result_epsilon_sum <- opt_fixed_pval_price_purity(epsilon_price_, epsilon_purity_)
    objvals <- c()
    for (j in 1:length(result_epsilon_sum$pool)) {
      objvals <- c(objvals, result_epsilon_sum$pool[[j]]$objval)
    }
    epsilon_price_purity_results <- rbind(epsilon_price_purity_results, c(epsilon_price_, epsilon_purity_, result_epsilon_sum$status, sum(objvals == min(objvals))))
  }
}
names(epsilon_price_purity_results) <- c("epsilon_price", "epsilon_purity", "feasibility", "n_opt_sols")
epsilon_price_purity_results %>% filter(feasibility == "OPTIMAL")
epsilon_price_purity_results %>% filter(epsilon_purity == 1)


alpha_results <- data.frame()
for (alpha in seq(1, 0, by=-0.1)) {
  model_alpha <- list()
  model_alpha$obj        <- c(rep(1-alpha, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha, N))
  model_alpha$modelsense <- "min"
  model_alpha$A          <- rbind(A_consumption[,-zero_d_vars_index], A_source_indicator[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_alpha$sense      <- rep(">", nrow(model_alpha$A))
  model_alpha$sense[1:N] <- "="
  model_alpha$sense[nrow(model_alpha$A)] <- "="
  model_alpha$vtype      <- c(rep("C", length(model_alpha$obj) - N), rep("B", N)) # Continuous and binary
  model_alpha$rhs        <- c(b_consumption, b_source_indicator, b_source_limit)
  result_alpha_sum       <- gurobi(model_alpha, list(Method=-1))
  objvals <- c()
  for (j in 1:length(result_alpha_sum$pool)) {
    objvals <- c(objvals, result_alpha_sum$pool[[j]]$objval)
  }
  alpha_results <- rbind(alpha_results, c(alpha, result_alpha_sum$status, sum(objvals == min(objvals))))
}
names(alpha_results) <- c("alpha", "feasibility", "n_opt_sols")
alpha_results


for (alpha in seq(0.99, 0.91, by=-0.01)) {
  model_alpha <- list()
  model_alpha$obj        <- c(rep(1-alpha, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha, N))
  model_alpha$modelsense <- "min"
  model_alpha$A          <- rbind(A_consumption[,-zero_d_vars_index], A_source_indicator[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_alpha$sense      <- rep(">", nrow(model_alpha$A))
  model_alpha$sense[1:N] <- "="
  model_alpha$sense[nrow(model_alpha$A)] <- "="
  model_alpha$vtype      <- c(rep("C", length(model_alpha$obj) - N), rep("B", N)) # Continuous and binary
  model_alpha$rhs        <- c(b_consumption, b_source_indicator, b_source_limit)
  result_alpha_sum       <- gurobi(model_alpha, list(Method=-1))
  
  opt_sol <- data.frame(deicision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_alpha_sum$pool[[1]]$xn) %>%
    filter(optimal_sols > 0 & !grepl("x", deicision_vars))
  result_flow_index <- grep("w", opt_sol$deicision_vars)
  result_source_index <- grep("S", opt_sol$deicision_vars)
  max_flow <- max(opt_sol$optimal_sols[result_flow_index])
  
  opt_sol <- cbind(opt_sol, find_state_index(opt_sol$deicision_vars) %>% t) %>% 
    rename(i="1", j="2")
  
  opt_sol$long_i <- states_data$long[opt_sol$i]
  opt_sol$lat_i <- states_data$lat[opt_sol$i]
  opt_sol$long_j <- states_data$long[opt_sol$j]
  opt_sol$lat_j <- states_data$lat[opt_sol$j]
  
  overdose_map_year <- left_join(states %>%
                                   rename(state=state_name) %>% 
                                   filter(!(state %in% c("Alaska", "Hawaii"))),
                                 states_data %>% select(state, med_death),
                                 by="state")
  overdose_map_year %>% ggplot() +
    geom_polygon(aes(x=long,
                     y=lat,
                     group=group,
                     fill=med_death),
                 color="black") +
    scale_fill_viridis_c(na.value="white") +
    expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
    coord_quickmap() +
    labs(x="", y="", title=bquote(alpha == .(alpha))) +
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          axis.line =  element_blank(),
          axis.text = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(data=states_data,
               aes(x=long, y=lat)) +
    geom_segment(data=opt_sol[result_flow_index,],
                 aes(x=long_i, 
                     y=lat_i, 
                     xend=long_j,
                     yend=lat_j),
                 linewidth = 0.2+2*opt_sol$optimal_sols[result_flow_index]/max_flow,
                 arrow=arrow(angle=10,
                             length=unit(0.3, "cm"),
                             type="closed")
    ) +
    geom_text(data=opt_sol[result_source_index,],
              aes(x=long_i, y=lat_i),
              label=opt_sol[result_source_index,]$optimal_sols %>% round(1),
              color="red",
              nudge_x = c(0.1, 0), nudge_y = c(-0.5, 0)) -> optimal_network_map
  ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks ",
                period[1], "-", period[length(period)],
                " (alpha=",
                alpha,
                ").png"),
         optimal_network_map, scale=1.5)
}


# epsilon_seizure test
alpha_i <- 0.9
epsilon_si <- 0.44
opt_val_i <- 18.96
epsilon_s_results <- data.frame()
for (epsilon_s in seq(0, 1, by=0.02)) {
  model_epsilon_s <- list()
  model_epsilon_s$obj        <- c(rep(1-alpha_i, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha_i, N))
  model_epsilon_s$modelsense <- "min"
  model_epsilon_s$A          <- rbind(A_consumption[,-zero_d_vars_index], A_seizure[,-zero_d_vars_index], A_source_indicator[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
  model_epsilon_s$sense      <- rep(">", nrow(model_epsilon_s$A))
  model_epsilon_s$sense[1:N] <- "="
  model_epsilon_s$sense[nrow(model_epsilon_s$A)] <- "="
  model_epsilon_s$vtype      <- c(rep("C", length(model_epsilon_s$obj) - N), rep("B", N)) # Continuous and binary
  model_epsilon_s$rhs        <- c(b_consumption, epsilon_s*b_seizure_ratio, b_source_indicator, b_source_limit)
  result_epsilon_s_sum       <- gurobi(model_epsilon_s, list(Method=-1))
  objvals <- c()
  for (j in 1:length(result_epsilon_s_sum$pool)) {
    objvals <- c(objvals, result_epsilon_s_sum$pool[[j]]$objval)
  }
  epsilon_s_results <- rbind(epsilon_s_results, c(epsilon_s, result_epsilon_s_sum$status, sum(objvals == min(objvals))))
}
names(epsilon_s_results) <- c("epsilon_s", "feasibility", "n_opt_sols")
epsilon_s_results

for (alpha_i in seq(1, 0, by=-.1)) {
  for (epsilon_si in seq(.4, .1, by=-.1)) {
    model_epsilon_s <- list()
    model_epsilon_s$obj        <- c(rep(1-alpha_i, n_d_vars - n_zero_d_vars - 2*N), rep(0, N), rep(alpha_i, N))
    model_epsilon_s$modelsense <- "min"
    model_epsilon_s$A          <- rbind(A_consumption[,-zero_d_vars_index], A_seizure[,-zero_d_vars_index], A_source_indicator[,-zero_d_vars_index], A_source_limit[,-zero_d_vars_index])
    model_epsilon_s$sense      <- rep(">", nrow(model_epsilon_s$A))
    model_epsilon_s$sense[1:N] <- "="
    model_epsilon_s$sense[nrow(model_epsilon_s$A)] <- "="
    model_epsilon_s$vtype      <- c(rep("C", length(model_epsilon_s$obj) - N), rep("B", N)) # Continuous and binary
    model_epsilon_s$rhs        <- c(b_consumption, epsilon_si*b_seizure_ratio, b_source_indicator, b_source_limit)
    result_epsilon_s_sum       <- gurobi(model_epsilon_s, list(Method=-1))
    
    opt_val_i <- round(result_epsilon_s_sum$objval, 2)
    
    opt_sol <- data.frame(decision_vars=d_vars[-zero_d_vars_index], optimal_sols=result_epsilon_s_sum$pool[[1]]$xn) %>%
      filter(optimal_sols > 0 & !grepl("x", decision_vars))
    result_flow_index <- grep("w", opt_sol$decision_vars)
    result_source_index <- grep("S", opt_sol$decision_vars)
    max_flow <- max(opt_sol$optimal_sols[result_flow_index])
    
    opt_sol %>% write.csv(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal solution ",
                                 period[1], "-", period[length(period)],
                                 " (alpha=",
                                 alpha_i,
                                 ", epsilton_s=",
                                 epsilon_si,
                                 ", z=",
                                 opt_val_i,
                                 ").csv"),
                          row.names=F)
    
    # opt_sol <- read.csv("Cocaine Network Optimization/Results with relaxation parameters with Overdose deaths T40.5 1999-2016/optimal solution 2009-2014 (alpha=0.9, epsilton_s=0.44, z=18.96).csv") %>%
    #   rename(decision_vars = deicision_vars)
    # result_flow_index <- grep("w", opt_sol$decision_vars)
    # result_source_index <- grep("S", opt_sol$decision_vars)
    # max_flow <- max(opt_sol$optimal_sols[result_flow_index])
    opt_sol <- cbind(opt_sol, find_state_index(opt_sol$decision_vars) %>% t) %>% 
      rename(i="1", j="2")
    
    states_sf <- get_urbn_map(map = "states", sf = TRUE)
    overdose_map_year <- left_join(states_sf %>%
                                     rename(state=state_name) %>%
                                     filter(!(state %in% c("Alaska", "Hawaii"))) %>% 
                                     mutate(centroid = st_centroid(geometry)),
                                   states_data %>% select(state, med_death, states_index),
                                   by="state") %>% 
      arrange(states_index)
    
    overdose_map_year_coords <- st_coordinates(overdose_map_year$centroid)
    opt_sol$long_i <- overdose_map_year_coords[,1][opt_sol$i]
    opt_sol$lat_i <- overdose_map_year_coords[,2][opt_sol$i]
    opt_sol$long_j <- overdose_map_year_coords[,1][opt_sol$j]
    opt_sol$lat_j <- overdose_map_year_coords[,2][opt_sol$j]
    
    opt_sol_Si <- opt_sol %>% filter(grepl("S", decision_vars)) %>% rename(states_index = i)
    opt_sol_Si_sf <- left_join(overdose_map_year %>% filter(states_index %in% opt_sol_Si$states_index), opt_sol_Si %>% select(states_index, optimal_sols), by="states_index")
    
    overdose_map_year %>% ggplot() +
      geom_sf(aes(fill = med_death),
              color="black") +
      # geom_polygon(aes(x=long,
      #                  y=lat,
      #                  group=group,
      #                  fill=med_death),
      #              color="black") +
      scale_fill_viridis_c(na.value="white") +
      expand_limits(x=overdose_map_year$long, y=overdose_map_year$lat) +
      # coord_quickmap() +
      labs(x="", y="", fill="") +#title=bquote(paste(epsilon[s], "=", .(epsilon_si), ", ", alpha, "=", .(alpha_i), ", opt_val=", .(opt_val_i)) )) +
      theme_bw() + 
      theme(axis.ticks = element_blank(),
            axis.line =  element_blank(),
            axis.text = element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_sf(data = overdose_map_year$centroid, color = "black", size = 2) +
      geom_segment(data=opt_sol[result_flow_index,],
                   aes(x=long_i, 
                       y=lat_i, 
                       xend=long_j,
                       yend=lat_j),
                   color="grey60",
                   linewidth = 0.2+2*opt_sol$optimal_sols[result_flow_index]/max_flow,
                   arrow=arrow(angle=10,
                               length=unit(0.3, "cm"),
                               type="closed")
      ) +
      geom_sf_text(data = opt_sol_Si_sf, aes(label=round(optimal_sols, 1)), color = "red", size =3) -> optimal_network_map
    ggsave(paste0("Cocaine Network Optimization/Results with relaxation parameters/optimal networks epsilton_s ",
                  period[1], "-", period[length(period)],
                  " (alpha=",
                  alpha_i,
                  ", epsilton_s=",
                  epsilon_si,
                  ", z=",
                  opt_val_i,
                  ").png"),
           optimal_network_map, scale=1)#width=20, height=10, unit="cm")
  }
}
