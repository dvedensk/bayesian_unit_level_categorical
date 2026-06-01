#takes in an n x m matrix
#returns an array of n/block_size matrices each of dim. block_size x m
mat_split <- function(input.matrix, block.size){
  n.rows <- nrow(input.matrix)
  n.cols <- ncol(input.matrix)
  if(n.rows %% block.size != 0){stop("Not a valid block size")}

  n.blocks <- n.rows/block.size
  output.array <- array(NA, dim=c(n.blocks, block.size, n.cols))
  for(i in 1:n.blocks){
    start_ind <- (i-1)*block.size +1
    end_ind <- i*block.size
    output.array[i,,]  <- input.matrix[start_ind:end_ind,]
  }
  return(output.array)
}

#takes an N x K matrix of probabilities and returns an N x (K+1) matrix
stick <- function(Theta){
    K <- ncol(Theta)+1
    pred.probs <- matrix(NA, nrow=nrow(Theta), ncol=K)
    pred.probs[,1] <- Theta[,1]
    for(kk in 1:(K-2)){
      pred.probs[ , kk + 1] <- Theta[ , kk + 1]/Theta[ , kk] * (1 - Theta[ , kk]) * pred.probs[,kk]
    }
    pred.probs[ , K] <- (1 - Theta[ , (K-1)])/Theta[,(K-1)] * pred.probs[, (K-1)]
    return(pred.probs)
}

#function for reconstructing stick-breaking probabilities output by our model
#input is an array of size n_pop x (K-1) x n_iter
#output flattens this array into a matrix of dimension  (n_pop*n_iter) x K
array_stick <- function(pred.probs, n_pop, n_iter, K){
    Theta <- aperm(pred.probs[,,], c(1,3,2))
    dim(Theta) <- c(n_pop*n_iter, K - 1)
    pred.probs <- matrix(NA, nrow=n_pop*n_iter, ncol=K)
    pred.probs[,1] <- Theta[,1]
    for(kk in 1:(K-2)){
      pred.probs[ , kk + 1] <- Theta[ , kk + 1]/Theta[ , kk] * (1 - Theta[ , kk]) * pred.probs[,kk]
    }
    pred.probs[ , K] <- (1 - Theta[ , (K-1)])/Theta[,(K-1)] * pred.probs[, (K-1)]
    return(pred.probs)
}

sample_MVN <- function(mu, prec){
  U <- base::chol(prec)
  tmp_norm <- rnorm(nrow(prec))
  ret_mat <- backsolve(U, backsolve(U, mu, transpose=TRUE) + tmp_norm)
  return(ret_mat)    
}

get_sample <- function(HPS_pop_long=HPS_pop_long, 
                       HPS_pop_wide=HPS_pop_wide,
                       sample_size=2500,
                       n_weeks,
                       n_areas){
#add a size variable for weights
   individual_means <- rowMeans(
                           sapply(HPS_pop_wide[ , c("RESPONSE_1", "RESPONSE_2", "RESPONSE_3")],
                                  as.numeric),
                           na.rm=T)

  HPS_pop_wide$INDIV_MEAN <- individual_means
  HPS_pop_wide <- HPS_pop_wide 

  HPS_pop_wide$SIZE_VAR <- as.numeric(exp(.1*scale(log(HPS_pop_wide$ORIG_WEIGHT)) +
                                          .2*scale(HPS_pop_wide$INDIV_MEAN)))

  ###These are *true* probabilities of selection
  inclusion_probs <- inclusionprobabilities(HPS_pop_wide$SIZE_VAR, sample_size)
  inclusion_probs <- inclusion_probs/sum(inclusion_probs) * sample_size
  
  HPS_pop_wide$PWEIGHT <- 1/inclusion_probs #These are true weights
  HPS_pop_wide$PROB_SELECTION <- inclusion_probs #These are true prob. selection
  indices <- UPpoisson(inclusion_probs)
  SCRAM_to_sample <- HPS_pop_wide$SCRAM[indices==1]
  
  HPS_sample_wide <- filter(HPS_pop_wide, SCRAM %in% SCRAM_to_sample)
  #Take the sample in long format
  #Need to join first because HPS_pop_wide now has the sampling weights we'll 
  #need later
  HPS_sample_long <- HPS_pop_long %>% 
                        filter(SCRAM %in% SCRAM_to_sample)

  reduced_wide <- select(HPS_sample_wide, SCRAM, PWEIGHT, PROB_SELECTION)
  HPS_sample_long <- merge(HPS_sample_long, reduced_wide, 
                           all=TRUE, by="SCRAM") %>% tibble
    
  return(list(HPS_sample_wide=HPS_sample_wide,
              HPS_sample_long=HPS_sample_long))
}

get_y1_y2 <- function(HPS_sample_long){
  #need to scale weights so that at each week
  #they sum to that week's sample size
  HPS_sample_long <- HPS_sample_long %>% group_by(WEEK) %>% 
                                    mutate(SCALE_WEIGHT=
                                              PWEIGHT*n()/sum(PWEIGHT)) %>% 
                                   ungroup() %>%
                                   arrange(WEEK)

  ##### First, get first-time respondents
  y_1.df <- filter(HPS_sample_long, !IS_FOLLOWUP)
  
  #get the format we need. i.e. list of N_{t,a} to track indices
  N_1 <- y_1.df %>% arrange(WEEK) %>% 
                 group_by(WEEK) %>% 
                 summarize(N_1=n()) %>%
                 select(N_1) %>% ungroup()
  N_1 <- N_1$N_1
  
  covar_names <-  y_1.df %>% select(starts_with("COVAR")) %>% 
                             names() %>% 
                             paste(collapse="+") 
  covar_formula <- as.formula(paste("~", covar_names))
            
  X_1 <- model.matrix(covar_formula, data=y_1.df)

  ###Get repeat respondents in a N_2 x 2 matrix
  ##for each week from 2 to T, check if there was a response in the previous week, 
  ##which we store as entry [i,1]
  #first handle respondents who answered twice
  y_2.df <- filter(HPS_sample_long, IS_FOLLOWUP)
  
  N_2 <- y_2.df %>% arrange(WEEK) %>% 
                 group_by(WEEK) %>% 
                 summarize(N_2=n()) %>%
                 select(N_2)

  N_2 <- c(0, N_2$N_2) # pad with 0 so we can treat it as a list of dim n_weeks
                       # obviously there are no repeat respondents in week 1.
  
  X_2 <- model.matrix(covar_formula, data=y_2.df)

  return(list(y_1.df=y_1.df, N_1=N_1, X_1=X_1,
              y_2.df=y_2.df, N_2=N_2, X_2=X_2))
}

get_direct_estimates <- function(HPS_sample_long, population_counts_by_time,
                                 n_areas, n_weeks){

  HPS_sample_long <-  HPS_sample_long %>% mutate(WEEK=factor(WEEK, levels=1:n_weeks))
  y_1 <- filter(HPS_sample_long, !IS_FOLLOWUP)
  y_2 <- filter(HPS_sample_long, IS_FOLLOWUP)
  #Direct estimate does not use scaled weights PWEIGHT instead
  # (sum w_i*y_i)/N_{population,t} where N_population is # unique respondents at time t
  response_1 <- y_1 %>% select(AREA, WEEK, RESPONSE, PWEIGHT) %>% 
                               mutate(INCLUSION_PROB=1/PWEIGHT)
  response_2 <- y_2  %>% select(AREA, WEEK, RESPONSE, PWEIGHT) %>% 
                               mutate(INCLUSION_PROB=1/PWEIGHT)
  response_by_week_area <- bind_rows(response_1, response_2)

  response_by_week_area <- data.frame(WEEK=factor(1:n_weeks), 
                                      N=population_counts_by_time) %>% 
                               right_join(response_by_week_area, by="WEEK") 

  response_by_week_area <- response_by_week_area %>%           
                            group_by(WEEK) %>%
                            mutate(SCALE_WEIGHT=
                                        PWEIGHT*N/sum(PWEIGHT)) %>%
                            mutate(INCLUSION_PROB=1/SCALE_WEIGHT) %>%
                            ungroup() 

  response_by_week_area <- response_by_week_area %>%
                             group_by(AREA, WEEK) 
  
  ##Add in NA for WEEK/AREA not represented in sample?
  tmp <- expand.grid(AREA=levels(HPS_sample_long$AREA),
                     WEEK=levels(HPS_sample_long$WEEK),
                     RESPONSE=levels(HPS_sample_long$RESPONSE)) #cover all possibilities regardless of sample
  dir_est_table <- HPS_sample_long %>% group_by(AREA, WEEK, RESPONSE) %>% summarize(N=n())
  dir_est_table <- left_join(tmp, dir_est_table, by=c("AREA", "WEEK","RESPONSE"))
  
  J <- nrow(dir_est_table)
  for(j in 1:J){
    if(j %% 50 == 0){print(paste0("Direct estimate #", j, "/", J))}
    curr_area <- dir_est_table[[j,"AREA"]]
    curr_week <-  dir_est_table[[j,"WEEK"]]
    curr_cat <- dir_est_table[[j,"RESPONSE"]]
    curr_resps <- filter(response_by_week_area, AREA==curr_area &
                                                WEEK==curr_week)
    vals <- as.integer(curr_resps$RESPONSE==curr_cat)                                      
    #try LinHB, if that doesn't work try LinHH
#    dir_est_table[which(complete.cases(dir_est_table)),]
    HT <- horvitzThompson(y=vals,
                          pi=curr_resps$INCLUSION_PROB, 
                          var_est=T, var_method = "LinHTSRS") 
    dir_est_table[j,"MEAN"] = HT$pop_mean
    dir_est_table[j,"VAR"] = HT$pop_mean_var
    dir_est_table[j,"CI_LOWER"] = HT$pop_mean - 1.96*sqrt(HT$pop_mean_var)
    dir_est_table[j,"CI_UPPER"] = HT$pop_mean + 1.96*sqrt(HT$pop_mean_var)
  }
  
  dir_est_table$WEEK <- as.double(dir_est_table$WEEK)
  return(dir_est_table)  
}


get_svy_direst <- function(HPS_pop_long, HPS_sample_long) {
  tst <- HPS_sample_long %>%
      group_by(WEEK) %>%
      mutate(SCALE_WEIGHT_POP = PWEIGHT/sum(PWEIGHT)*population_counts_by_time[WEEK])

  samp.design <- svydesign(ids = ~1, weights=~SCALE_WEIGHT_POP, data=tst)
  svy_de <- svyby(
      ~RESPONSE,
      ~AREA + WEEK,
      samp.design,
      svymean,
      na.rm = TRUE,
      vartype = "ci",
      keep.names = FALSE
  )

  svy_de_long <- svy_de %>%
      pivot_longer(
          cols = -c(AREA,WEEK),
          names_to = c("stat", "CATEGORY"),
          names_pattern = "(RESPONSE|ci_l\\.RESPONSE|ci_u\\.RESPONSE)([0-9]+)",
          values_to = "value"
      ) %>%
      mutate(
          CATEGORY = as.integer(CATEGORY),
          stat = recode(
              stat,
              "RESPONSE" = "point_est",
              "ci_l.RESPONSE" = "ci_lower",
              "ci_u.RESPONSE" = "ci_upper"
          )
      ) %>%
      pivot_wider(
          names_from = stat,
          values_from = value
      )
}

make_areal_plot_model <- function(plot_df, filepath=NULL){
 states <- tolower(levels(plot_df$AREA))
 state_poly <- map_data("state") %>% filter(region %in%  states)
 plot_df$region<- tolower(plot_df$AREA)
 plot_df <- merge(plot_df, state_poly, by="region")

 plot_df %>%
   ggplot(mapping = aes(x=long, y=lat,
                        fill=prop, group=group)) +
     geom_polygon(color="black", linewidth=0.15) +
     viridis::scale_fill_viridis(name="Estimated proportion", discrete=F)+
     guides(pattern = "none", #guide_legend(override.aes = list(pattern = "stripe")),        
            fill=guide_colorbar()) +
     coord_map() +
     facet_wrap(~WEEK) +
     ggthemes::theme_map(base_size=10) +
     theme(
          legend.position="right",
          legend.title = element_text(size = 12),
          legend.text  = element_text(size = 12),
          legend.key.height = grid::unit(0.8, "cm"),
          legend.key.width  = grid::unit(0.8, "cm"),
          legend.box.margin = margin(0, 8, 0, 0),
          legend.margin     = margin(8, 8, 8, 8),
          strip.text = element_text(size = 14, face = "bold"),
          plot.margin=grid::unit(c(0,0,0,0), "mm")
      )

   if(!is.null(filepath)){ 
       ggsave(filepath,             
              dpi=600,
              width = 14,
              height = 8,
              units = "in",
              limitsize=FALSE
              )
       system2(command = "pdfcrop",
               args    = c(filepath,
                           filepath))
   }
}

make_areal_plot_DE <- function(plot_df, filepath=NULL){
  plot_df <- plot_df %>%
      mutate(pattern_flag = ifelse(is.na(prop), "stripe", "none"))

 states <- tolower(levels(plot_df$AREA))
 state_poly <- map_data("state") %>% filter(region %in%  states)
 plot_df$region<- tolower(plot_df$AREA)
 plot_df <- merge(plot_df, state_poly, by="region")

 plot_df %>% ggplot() +
     geom_polygon_pattern(
         aes(
             x = long, y = lat,
             pattern = pattern_flag,
             fill = prop,
             group = group
         ),
         pattern_density = 0.1,
         pattern_angle = 45,
         color = "black",
         size = 0.15
     ) +
     coord_map() +
     facet_wrap(~WEEK) +
     ggthemes::theme_map(base_size=10) +
     theme(
         legend.position="right",
         legend.title = element_text(size = 12),
         legend.text  = element_text(size = 12),
         legend.key.height = grid::unit(0.8, "cm"),
         legend.key.width  = grid::unit(0.8, "cm"),
         legend.box.margin = margin(0, 8, 0, 0),
         legend.margin     = margin(8, 8, 8, 8),
         strip.text = element_text(size = 14, face = "bold"),
         plot.margin=grid::unit(c(0,0,0,0), "mm")
     ) +
     scale_pattern_manual(values = c(stripe = "stripe",
                                     none = "none"),
                          guide="none") +
     viridis::scale_fill_viridis(name="Estimated proportion", discrete=F) +
     guides(pattern = "none", #guide_legend(override.aes = list(pattern = "stripe")),        
            fill=guide_colorbar())

  if(!is.null(filepath)){ 
      ggsave(filepath,
             dpi = 600,
             width = 14,
             height = 8,
             units = "in",
             limitsize=FALSE)
      system2(command = "pdfcrop",
              args    = c(filepath,
                          filepath))
  }
}
