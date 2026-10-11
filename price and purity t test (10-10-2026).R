library(readxl)
library(stringi)
library(urbnmapr) # install by running -> devtools::install_github("UrbanInstitute/urbnmapr")
library(tidyverse)
library(gridExtra)
library(lubridate)

{
  stride <- read.csv("Cocaine Network Optimization/STRIDE_Raw.csv") %>% as_tibble %>% filter(!is.na(Seize.Year) & !is.na(Seize.Month) & !(State %in% c("AK", "HI")))
  stride$Nt.Wt <- as.numeric(stride$Nt.Wt)
  stride$Potency <- as.numeric(stride$Potency)
  stride$Post.Price <- as.numeric(stride$Post.Price)
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
  
  years1 <- 2000:2003
  years2 <- 2004:2008
  years3 <- 2009:2014
  
  period <- years3
  price_period <- cocaine %>%
    filter(!is.na(state) & Seize.Year %in% period & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000) %>%
    group_by(state) %>% 
    summarise(avg_price=mean(adjusted_price, na.rm=T),
              med_price=median(adjusted_price, na.rm=T))
  
  cocaine_period <- cocaine %>% 
    filter(!is.na(state) & Seize.Year %in% c(years2, years3) & adjusted_price > 0 & Nt.Wt >= 5 & Nt.Wt <= 1000)
  cocaine_period %>% filter(state %in% c("Vermont", "South Dakota", "Oregon"))
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

t_test_summary <- tibble()
for (i in 1:49) {
  state_i <- states_data$state[i]
  bordering_states <- neighbor[[state_i]]
  
  state_prices <- cocaine_period %>%
    filter(state == state_i) %>%
    pull(adjusted_price)
  
  if (length(state_prices) < 3) next
  
  for(bordering_state in bordering_states) {
    bordering_state_prices <- cocaine_period %>%
      filter(state == bordering_state) %>%
      pull(adjusted_price)
    
    if (length(bordering_state_prices) < 3) next
    price_t_test <- t.test(state_prices, bordering_state_prices, alternative = "greater")
    t_test_summary <- rbind(t_test_summary,
                            tibble(state_index=i,
                                   bordering_state_index=states_data$states_index[which(states_data$state == bordering_state)],
                                   state=state_i,
                                   bordering_state=bordering_state,
                                   mean_state=price_t_test$estimate[1],
                                   mean_bordering_state=price_t_test$estimate[2],
                                   p_value=price_t_test$p.value))
  }
}
t_test_summary
write.csv(t_test_summary, paste0("Cocaine Network Optimization/price t_test_summary greater (", period[1], "-", period[length(period)], ").csv"), row.names=F)