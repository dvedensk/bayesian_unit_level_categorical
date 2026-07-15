library(readr)
library(purrr)
library(tidyr)
library(dplyr)

response_type <- "ordinal"
type <- "simulation"

data_dir <- file.path("data", response_type)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

#Data can be downloaded with a command like e.g.
# for i in {1..9}; do wget https://www2.census.gov/programs-surveys/demo/datasets/hhp/2020/wk$i/HPS_Week0$i'_PUF_CSV.zip'; done
# unfortunately inconvenient that they are inconsistent about leading zeros, so handle 10-12 separately.
# unzip HPS*zip
# then run the code below

#following https://stackoverflow.com/a/65001063
filenames <- list.files(pattern="pulse2020_puf_[0-9]+.csv$",recursive=T)
filenames <- filenames[1:12]
combined_data <- purrr::map_df(filenames, ~read_csv(.x) %>% mutate(filename = .x))
                                
n_weeks <- length(filenames)

#-99 = Question seen but category not selected -88 = Missing/did not report
#1 = Not at all, 2 = several days, 3 = more than half the days, 4 = nearly every day
HPS_df <- select(combined_data, SCRAM, WEEK, AREA=EST_ST, TBIRTH_YEAR,
                 ORIG_WEIGHT=PWEIGHT, EGENDER, RRACE, ANXIOUS, WORRY) %>%
    filter(ANXIOUS %in% c(1,2,3,4),
           WORRY %in% c(1,2,3,4)) %>%
    mutate(ANXIOUS = ANXIOUS - 1,
           WORRY=WORRY-1,
           RESPONSE = ANXIOUS+WORRY) #GAD-2 score
HPS_df$RESPONSE <- as.factor(HPS_df$RESPONSE)

HPS_df <- HPS_df %>% 
              mutate(COVAR.1=factor(EGENDER, levels=c(1,2),
                                             labels=c("MALE","FEMALE"))) %>%
              mutate(COVAR.2=as.integer(2020-TBIRTH_YEAR)) %>% 
              mutate(COVAR.3=factor(RRACE, levels=c(1,2,3,4),
                                             labels=c("White","Black", "Asian","Other"))) %>%
              select(-TBIRTH_YEAR, -EGENDER, -RRACE) 

age_breaks <- c(17, seq(25,65,5),100)
HPS_df <- mutate(HPS_df,
                 COVAR.2=cut(COVAR.2, breaks=age_breaks))

HPS_df$AREA <- recode_factor(HPS_df$AREA,
                    '01'='Alabama', '02'='Alaska', '04'='Arizona',
                    '05'='Arkansas','06'='California','08'='Colorado',
                    '09'='Connecticut', '10'='Delaware', '11'='District of Columbia',
                    '12'='Florida','13'='Georgia','15'='Hawaii',
                    '16'='Idaho','17'='Illinois','18'='Indiana',
                    '19'='Iowa','20'='Kansas','21'='Kentucky',
                    '22'='Louisiana','23'='Maine','24'='Maryland',
                    '25'='Massachusetts','26'='Michigan','27'='Minnesota',
                    '28'='Mississippi','29'='Missouri','30'='Montana',
                    '31'='Nebraska','32'='Nevada','33'='New Hampshire',
                    '34'='New Jersey','35'='New Mexico','36'='New York',
                    '37'='North Carolina','38'='North Dakota','39'='Ohio', 
                    '40'='Oklahoma','41'='Oregon','42'='Pennsylvania',
                    '44'='Rhode Island','45'='South Carolina','46'='South Dakota',
                    '47'='Tennessee','48'='Texas','49'='Utah',
                    '50'='Vermont','51'='Virginia','53'='Washington',
                    '54'='West Virginia','55'='Wisconsin','56'='Wyoming')

HPS_df_long <- HPS_df %>% arrange(WEEK) %>% 
                          group_by(SCRAM) %>% 
                          mutate(RESPONSE_NUMBER=1:n()) %>% 
                          #Take only the first WEIGHT, which we'll 
                          #use to make a SIZE_VAR later:
                          mutate(ORIG_WEIGHT=ORIG_WEIGHT[1]) %>%
                          arrange(SCRAM) %>% 
                          ungroup()

#Want to add a logical flag for whether a response is a follow-up
#and a PREV_RESPONSE column for response at time t-1 if follow-up (o.w. NA)
HPS_df_long <- HPS_df_long %>% 
                 mutate(IS_FOLLOWUP=ifelse(RESPONSE_NUMBER>1, TRUE, FALSE)) %>%
                 mutate(PREV_RESPONSE=ifelse(IS_FOLLOWUP, lag(RESPONSE), NA)) 

new_covar <- addNA(factor(HPS_df_long$PREV_RESPONSE))
num_lev <- length(unique(new_covar)) - 1
levels(new_covar) <- c(paste0("PREV_RESP_",1:num_lev), "PREV_RESP_NA")
new_covar <- relevel(new_covar, ref="PREV_RESP_NA")
HPS_df_long$COVAR.NEW <- new_covar
    
#Also want a wide data frame to sample from. 
HPS_df_wide <- HPS_df_long %>%  pivot_wider(values_from=c(RESPONSE, WEEK), 
                                            names_from=c(RESPONSE_NUMBER),
                                            id_cols=SCRAM)

#Join back to get covariate values
HPS_df_wide <- HPS_df_wide %>% left_join(HPS_df_long, 
                                         by=c("SCRAM"="SCRAM", 
                                              "RESPONSE_1"="RESPONSE", 
                                              "WEEK_1"="WEEK")) %>%
                               select(-RESPONSE_NUMBER)

save(HPS_df_long, HPS_df_wide,
         file=file.path(data_dir, "HPS_empirical_pop_df_GAD2.RData"))

