

# Description: Calculator for determining power with WHO tube experiments
# Author: Frank Mechan (frank.mechan@lstmed.ac.uk)


## Outstanding features: 

# if ever have compute power/time
# smaller (1%) increments, 10000 simulations
# greater range of variability


# Load libraries
library(lme4)
library(lmtest)
library(tidyverse)



############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################


WHO_Tube_Power_Calculator <- function(Nsims, Ndays, Ntubes, WithinSD, BetweenSD) {
  
####################################################################################################
#### Simulation Parameters ####
    
  if(missing(Nsims))      { nsims       <- 1000  }    else { nsims       <- Nsims       }   # Number of iterations of each variance/difference combination
  
  if(missing(Ndays))      { n_days      <- 1     }    else { n_days      <- Ndays       }   # Number of testing days

  if(missing(Ntubes))     { n_tubes     <- 4     }    else { n_tubes     <- Ntubes      }   # Replicates OF EACH GROUP per testing day
  
  if(missing(WithinSD))   { within.sd   <- 0.2   }    else { within.sd   <- WithinSD    }   # standard deviation for within-day variation
  
  if(missing(BetweenSD))  { between.var <- 0.8   }    else { between.var <- BetweenSD  }   # standard deviation for between-day variation
  

Mosq_Num      <- 25      # Number of mosquitoes per assay
control.var   <- 0       #### power inflated when thus value high (why?) ####
treatment.var <- 0       #### power inflated when thus value high (why?) ####


##################################################################################################
#### Set Up Output Dataframe ####


All.TubeDay.Combos.df   <- expand.grid(Tubes = n_tubes,  # dataframe of all requested combinations of tubes and days
                                       Days  = n_days)

All.Variance.Combos.df  <- expand.grid(WithinSD  = within.sd,
                                       BetweenSD = between.var)    # dataframe of all variance value Tube*Day combinations

# Create empty dataframe(these will be filled with outputs as the model runs)
All.Combos.Output.df <- as.numeric()  # empty vector that will accumulate each tube/day/var combination of each pass

#####################################################################
          ###### (LEVEL 1 loop) #####
          #### Begin loop over Tube*Day combinations (highest level loop of function) ####
          
          
          # loop begins
          # Highest level of loop. loops over all tube*day*var combinations
          # on each iteration, each unique combination of tube*day*var evaluated
          # Purpose of this level is to set up dataframe of starting parameters for subsequent stages
          
          for ( i in 1:nrow(All.TubeDay.Combos.df)) { 
            
            
            itube  <- All.TubeDay.Combos.df$Tubes[i] # holds value of tube no. on this iteration
            iday   <- All.TubeDay.Combos.df$Days[i]  # holds value of total testing days on this iteration
          
            ThisTubeDay.Combo.df <- as.numeric() # empty vector that will accumulate outputs across each variance value for this unique tube*day combination
          
            ThisTubeDay.MeanDiff  <- as.numeric()  # mean estimated difference in mortality for each true difference
            ThisTubeDay.Power     <- as.numeric()  # proportion of comparisons significant for each true difference
          
          
                ##################################################################################
                #### Logic to distribute tubes over testing days ####
                
                # Aims to distribute tubes as evenly as possible over days
                # Example: 5 tubes over 3 days = 2,2,1
                # Developer Note: Ideally would be much more efficient and compact but works (for loop)
                
                
                # Calculate how tubes are split up over days 
                #  Spread out as evenly as possible
                #  produces vector (length = number of days) where each value represents number of tubes on that day 
                if (iday > 1) {
                  cut.vec <- tabulate(cut(1:itube , iday)) 
                } else {
                  cut.vec <- NA }
          
          
                    ## Create a vector (n = total number of tubes) of day numbers (representing each tube)
                    if (itube%%iday == 0) {
                      day.seq <- rep(1:iday, each=(itube*2)/iday, length.out=itube*2)
                      
                        }  else if  ((itube%%iday != 0) & (iday==2)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2))
                          
                        }  else if  ((itube%%iday != 0) & (iday==3)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2), rep(3, each=cut.vec[3]*2))
                        
                        }  else if  ((itube%%iday != 0) & (iday==4)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2), rep(3, each=cut.vec[3]*2), rep(4, each=cut.vec[4]*2)) 
                          
                        }  else if  ((itube%%iday != 0) & (iday==5)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2), rep(3, each=cut.vec[3]*2), rep(4, each=cut.vec[4]*2), rep(5, each=cut.vec[5]*2)) 
                          
                        }  else if  ((itube%%iday != 0) & (iday==6)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2), rep(3, each=cut.vec[3]*2), rep(4, each=cut.vec[4]*2), rep(5, each=cut.vec[5]*2), rep(6, each=cut.vec[6]*2)) 
                          
                        }  else if  ((itube%%iday != 0) & (iday==7)) {
                          day.seq <- c(rep(1, each=2*cut.vec[1]), rep(2, each=cut.vec[2]*2), rep(3, each=cut.vec[3]*2), rep(4, each=cut.vec[4]*2), rep(5, each=cut.vec[5]*2), rep(6, each=cut.vec[6]*2), rep(7, each=cut.vec[7]*2)) 
                          
                                } else if (iday>itube) {
                                next 
                    }  

                            ### Produces a dataset that distributes tubes of each treatment over testing days
                            # Tubes and days are numbered appropriately to allow simulation of variance as appropriate 
                            input.df     <- data.frame(treatment=rep(c("treatment", "control"), times=itube),
                                                       day=day.seq,
                                                       tube=rep(1:itube, each=2))


##############################################################################################################################################################################################

                                  ###### (LEVEL 2 loop) #####
                                  #### Loop over variance values for each Tube*Day combination ####
                    
                                    for (i in 1:nrow(All.Variance.Combos.df)) {
                                      tvar <- All.Variance.Combos.df$WithinSD[i]
                                      dvar <- All.Variance.Combos.df$BetweenSD[i]
                                      
                                      # dataframe of effect sizes to investigate (increments of 2.5% mortality)
                                      EffSiz.ThisCombo.TEMP.df <- expand.grid(Control   = 0.90, 
                                                                              Treatment = 0.92)
                                      
                                      tube_var = tvar      # sd between tubes on same day 
                                      day_var  = dvar      # sd between tubes on different day 
                                      
                              
#########################################################################################################################################################################################################                              
                                                    ###### (LEVEL 3 loop) #####
                                                    #### Loop over effect size values for each Tube*Day*Var combination  
                                                  
                                                    for(i in 1:nrow(EffSiz.ThisCombo.TEMP.df)) {
                                                    
                                                    codds  = EffSiz.ThisCombo.TEMP.df$Control[i]/(1 - EffSiz.ThisCombo.TEMP.df$Control[i])      # odds of dying in control arm  
                                                    todds  = EffSiz.ThisCombo.TEMP.df$Treatment[i]/(1 - EffSiz.ThisCombo.TEMP.df$Treatment[i])  # odds of dying in treatment arm 
                                                    
                                                    b0 <- log(codds)        # log odds of dying in control arm
                                                    b1 <- log(todds/codds)  # log odds of dying in treatment arm
                                                    
                                                    Runs_Mean       <- numeric() # empty vector that will store mean mortality for each Tube*Day*Var*EffSize combination
                                                    #Runs_Lwr        <- numeric() # empty vector that will store lower limit of CI 
                                                    #Runs_Upr        <- numeric() # empty vector that will store upper limit of CI
                                                    
                                                    Effect_Size_p <- numeric() # empty vector that will store power for each Tube*Day*Var*EffSize combination
      
  
####################################################################################################################################################################################################
  
                                                                          ###### (LEVEL 4 loop) #####
                                                                          #### Simulate data for each Tube*Day*Var*EffSiz combination  
                                                                          # Runs nsims number of times
                                                                          # loop that runs model of mortality for given control and treatment odds (runs n number of times, random effects is stchastic)
                                                                                  
                                                                                  for(i in 1:nsims) {    # runs nsims number of times
                                                                                    
                                                                                    sim.df = input.df # create local version on input dataframe (probably unnecessary)
                                                                                    
                                                                                    day.df  <- data.frame(day=unique(sim.df$day)) # a dataframe of unique day values
                                                                                    day.df$day_eff <- rnorm(n = nrow(day.df),  mean = 0, sd = day_var) # simulates a distribution of values (with mean 0) that represent variation due to between-day effect
                                                                                    sim.df.DayEff <- left_join(sim.df, day.df, by= "day") # joins day variance values with main input dataset
                                                                                    
                                                                                    tube.df <- data.frame(tube=unique(sim.df$tube)) # as above for variation between days
                                                                                    tube.df$tube_eff <- rnorm(n = nrow(tube.df),  mean = 0, sd = tube_var)
                                                                                    sim.Eff.TubeEff <- left_join(sim.df.DayEff, tube.df, by= "tube")
                                                                                    
                                                                                    Control.df           <- data.frame(treatment="control")
                                                                                    Treatment.df         <- data.frame(treatment="treatment")
                                                                                    Control.df$treat_eff   <- rnorm(n = nrow(Control.df),    mean = 0, sd = sqrt(control.var))
                                                                                    Treatment.df$treat_eff <- rnorm(n = nrow(Treatment.df),  mean = 0, sd = sqrt(treatment.var))
                                                                                    ConTreatVar.df <- rbind(Control.df, Treatment.df)
                                                                                    sim.Eff.df <- left_join(sim.Eff.TubeEff, ConTreatVar.df, by= "treatment")
                                                                                    
                                                                                  # calculate log odds of death for this iteration (based on randomised sd values obtained above )
                                                                                    sim.Eff.df$log_odds = with(sim.Eff.df, b0 + b1*(treatment == "treatment") + sim.Eff.df$day_eff + sim.Eff.df$tube_eff + sim.Eff.df$treat_eff)
                                                                                    sim.Eff.df$prop = plogis(sim.Eff.df$log_odds) # convert log odds to probability
                                                                                    sim.Eff.df$TotalMosq = Mosq_Num               # Number of mosquitoes in each tube
                                                                                    
                                                                                  # calculate number of dead in each tube. No stochasticisty in this step, simply uses proportions simulated above 
                                                                                    sim.Eff.df$Dead = rbinom(n = nrow(sim.Eff.df), size = sim.Eff.df$TotalMosq, prob = sim.Eff.df$prop)
                                                                            
                                                                              
                                                                                    
                                                                                ###########################################################################################################
                                                                                    
                                                                                  #### Perform GLM ####
                                                                                    
                                                                                  ### Assess if there was a difference in mortality between treatments
                                                                                        sim.Eff.df$treatment <- factor(sim.Eff.df$treatment, levels = c("control","treatment")) # convert 'treatment' variable to factor to allow assessment in GLM
                                                             
                                                                                        mod = glm(cbind(Dead, TotalMosq-Dead) ~ 1 + treatment, data = sim.Eff.df,
                                                                                                    family = quasibinomial) # uses quasi-binomial link function to permit variance greater than standard binomial GLM (i.e. to counter overdispersion)
                                                                                    
                                                                                        
                                                                                        mod.null <- glm(cbind(Dead, TotalMosq-Dead) ~ 1, data = sim.Eff.df,
                                                                                                        family = quasibinomial) # uses quasi-binomial link function to permit variance greater than standard binomial GLM (i.e. to counter overdispersion)
                                                                                        
                                                                         
                                                                                        mod.an <- anova(mod, test="Chisq") # extract p-value for 'treatment' for this iteration
                                                                                        my_out <- mod.an$`Pr(>Chi)`[2]
                                                                                        Effect_Size_p <- c(Effect_Size_p, my_out) # store this p-value in a vector
                                                                                  
                                                                                  
                                                                                  ### Obtain predicted mean mortality for each treatment group
                                                                                  
                                                                                        new.df <- data.frame(treatment=c("treatment", "control"))
                                                                                            new.df$pred.mean <- predict(mod, newdata=new.df, type="response")
                                                                                            abs.diff         <- new.df$pred.mean[1]- new.df$pred.mean[2]
                                                                                        Runs_Mean   <- c(Runs_Mean, abs.diff) # store absolute difference between means in a vector
                                                                                  
                                                                                  }
                                          
###################################################################################################################################################################################################                                            
                                                      
                                                      ### calculate outputs for the current Tube*Day*Var*EffSize combination ###
                                                      
                                                      # Effect Size
                                                      Mean_EachContValue    <- as.numeric(mean(Runs_Mean))
                                                      ThisTubeDay.MeanDiff  <- c(ThisTubeDay.MeanDiff, Mean_EachContValue) # add estimated mortality difference value on to string of true EffSizes
                                          
                                                      # Power
                                                      power              <- sum(Effect_Size_p < 0.05)/length(Effect_Size_p) # calculate proportion of comparisons that were statistically significant (i.e. 'power')
                                                      ThisTubeDay.Power  <- c(ThisTubeDay.Power, power)
                                                      
                                          }
  
###################################################################################################################################################################################################                                            
                              
                                  ### Store outputs (mean difference and power) from each effect size comparison performed above (for this specific Tube/Day/Variance combination) in a temporary dataframe 
                                      
                                  EffSiz.ThisCombo.TEMP.df$MeanDiff <- ThisTubeDay.MeanDiff
                                  EffSiz.ThisCombo.TEMP.df$Power    <- ThisTubeDay.Power
                                  
                                  Effect.sizes.ThisVar         <- EffSiz.ThisCombo.TEMP.df
                                  Effect.sizes.ThisVar$TubeVar <- tvar
                                  Effect.sizes.ThisVar$DayVar  <- dvar
                                  
                                  # Add the contents of this temporary dataframe (for this specific Tube/Day/Variance combination) to the cumulative dataframe of all Tube/Day/Variance combinations
                                  ThisTubeDay.Combo.df <- rbind(ThisTubeDay.Combo.df, Effect.sizes.ThisVar)
                                  
                                  # Clear out contents of temporary dataframe to prepare for next loop
                                  Effect.sizes.ThisVar        <- as.numeric()
                                  ThisTubeDay.MeanDiff        <- as.numeric()
                                  EffSiz.ThisCombo.TEMP.df    <- as.numeric()
                                  ThisTubeDay.Power           <- as.numeric()
                          
                    }
        
                      ThisTubeDay.Combo.df$Ntube <- itube
                      ThisTubeDay.Combo.df$Nday  <- iday
                      
                      
                      ThisTubeDay.Combo.df$TrueDiff      <- ThisTubeDay.Combo.df$Treatment-ThisTubeDay.Combo.df$Control
                      ThisTubeDay.Combo.df$TrueDiff.Perc <- ThisTubeDay.Combo.df$TrueDiff*100
                      All.Combos.Output.df <- rbind(All.Combos.Output.df, ThisTubeDay.Combo.df)
                      ThisTubeDay.Combo.df <- as.numeric()
                      
                      # move on to next tube/day combination
      
      } 

      # End Tube/Day/Variance Cycle

        # Tidy up outputs
        All.Combos.Output.df$Design  <- as.factor(as.character(All.Combos.Output.df$Ntube))
        All.Combos.Output.df$Design  <- recode_factor(All.Combos.Output.df$Design, '1'= "1x1",'2'="2x2",'3'="3x3",'4'="4x4",'5'="5x5",'6'="6x6", '7'="7x7")
        All.Combos.Output.df$Power.Perc <- All.Combos.Output.df$Power*100
              return(All.Combos.Output.df)


} # End Function





############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################

# 5000 not sufficient for fully smooth curves



#Power.10over5.4to7.day1  <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(4,5,6,7), Ndays=1, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
#Power.10over5.8to10.day1 <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(8,9,10),  Ndays=1, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))

#Power.10over5.4to7.day2  <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(4,5,6,7), Ndays=2, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
#Power.10over5.8to10.day2 <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(8,9,10),  Ndays=2, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))

test  <- WHO_Tube_Power_Calculator(Nsims=100, Ntubes=4, Ndays=1, WithinSD=0.25, BetweenSD=0.60)
test

print(Sys.time())
start_time <- Sys.time()
Power.10over5.4to7.day3  <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(4,5,6,7), Ndays=3, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
Power.10over5.8to10.day3 <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(8,9,10),  Ndays=3, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
print(Sys.time())

Power.10over5.4to7.day4  <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(4,5,6,7), Ndays=4, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
Power.10over5.8to10.day4 <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(8,9,10),  Ndays=4, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))

Power.10over5.4to7.day5  <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(4,5,6,7), Ndays=5, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))
Power.10over5.8to10.day5 <- WHO_Tube_Power_Calculator(Nsims=8000, Ntubes=c(8,9,10),  Ndays=5, WithinSD=c(0.15,0.20,0.25,0.30), BetweenSD=c(0.40,0.50,0.60,0.70,0.80))


end_time <- Sys.time()
print(Sys.time())


end_time - start_time

####################################################################

#Power.7over4.Full <- rbind(Power.7over4.Day1, Power.7over4.day2, Power.7over4.day3,  Power.7over4.day4)



### Save output as excel sheet (commented out to prevent accidental overwrite)
#setwd("C:/Users/Frank.Mechan/Documents/LocalShinyTest/DesignApp/WHO_Power_DesignApp/data")


#write.csv(Power.10over5.4to7.day1,  "P:/i2i/Bioassay Thresholds/Data archive/Day1_4to7.csv")
#write.csv(Power.10over5.8to10.day1, "P:/i2i/Bioassay Thresholds/Data archive/Day1_5to10.csv")

#write.csv(Power.10over5.4to7.day2,  "P:/i2i/Bioassay Thresholds/Data archive/Day2_4to7.csv")
#write.csv(Power.10over5.8to10.day2, "P:/i2i/Bioassay Thresholds/Data archive/Day2_5to10.csv")

#write.csv(Power.10over5.4to7.day3,  "P:/i2i/Bioassay Thresholds/Data archive/Day3_4to7.csv")
#write.csv(Power.10over5.8to10.day3, "P:/i2i/Bioassay Thresholds/Data archive/Day3_5to10.csv")

#write.csv(Power.10over5.4to7.day4,  "P:/i2i/Bioassay Thresholds/Data archive/Day4_4to7.csv")
#write.csv(Power.10over5.8to10.day4, "P:/i2i/Bioassay Thresholds/Data archive/Day4_5to10.csv")

#write.csv(Power.10over5.4to7.day5,  "P:/i2i/Bioassay Thresholds/Data archive/Day5_4to7.csv")
#write.csv(Power.10over5.8to10.day5, "P:/i2i/Bioassay Thresholds/Data archive/Day5_5to10.csv")

#