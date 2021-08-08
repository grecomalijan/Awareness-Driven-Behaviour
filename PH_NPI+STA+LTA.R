#+eval=FALSE

#########################################################################################
# This code is divided into three sections.                                             #
# Section 1 shows the model and the optimisation algorithm I used to generate parameter #
# confidence intervals. It also shows the script used to generate model plots and to    #
# formally evaluate performance. Only the full NPI+STA+LTA model fitted to Philippine   #
# data is shown because the script is analogous. Code for other models fitted to other  #
# territories can be provided by request.                                               # 
# Section 2 generates the model performance comparison plots for all territories.       #
# Section 3 generates the heat maps for the sensitivity analysis performed for the full #
# NPI+STA+LTA model. The performance results generated for the sensitivity analysis     #
# applied the same code as in Section 1.                                                #
#########################################################################################

#### SETUP WORKSPACE AND LIBRARIES. ####
# Clear the work space.
rm(list=ls())

# Install packages.
install.packages("tidyverse")
install.packages("deSolve")
install.packages("lubridate")
install.packages("zoo")
install.packages("doParallel")

# Load libraries.
library(tidyverse)
library(deSolve)
library(lubridate)
library(zoo)
library(doParallel)

# Set working directory.
setwd("~/Documents/MSc IHTM/Awareness-Driven Behaviour/Philippines")

#########################################################################################
# SECTION 1.                                                                            #
# NON-PHARMACEUTICAL INTERVENTIONS + SHORT-TERM AWARENESS + LONG-TERM AWARENESS MODEL   #
# Fitted to Philippine COVID-19 Data                                                    #
#########################################################################################

#### PERFORM DATA WRANGLING. ####

## Extract case data. 
dataCases <- read_csv("PHcases.csv")%>%
  mutate(date=as.Date(DateRepConf,format="%m/%e/%y")) %>%
  select(date, dailyCases=Count) %>%
  drop_na() %>%
  arrange(date)

## Extract death data. 
dataDeaths <- read_csv("PHdeaths.csv") %>%
  mutate(date=as.Date(DateDied,format="%m/%e/%y")) %>%
  select(date, dailyDeaths=Count) %>%
  drop_na() %>%
  arrange(date) 

## Extract stringency index.
dataSI <- read_csv("PHSI.csv") %>%
  mutate(date=as.Date(Day, format="%m/%e%y")) %>%
  select(date, stringency_index) %>%
  drop_na() %>%
  arrange(date)

## Set data timeline.
dateEpoch <- min(min(dataCases$date), min(dataDeaths$date)) %>%
  as.Date()

dataDeaths <- dataDeaths %>%
  mutate(time=as.numeric(date-dateEpoch)) %>%
  left_join(dataSI, by='date')

## Apply stringency index.
str_ind <- approxfun(dataDeaths$stringency_index, rule=2)
str_ind(seq(from=0, to=300, by=1))

#### DEFINE FUNCTIONS TO BE USED IN THE CODE. ####

## Model ##
ratesFun <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),
       {
         dday <- gammaH*H
         N <- S+E+I+R+H
         SI <- str_ind(t)
         ld <- (1-(NPI*(SI/100)))
         lambda <- (beta*ld*I/N)/((1+(dday/dcrit)^k)+(D/dtotcrit)^k)
         dS <- -lambda*S
         dE <- lambda*S - mu*E
         dI <- mu*E - gamma*I
         dR <- gamma*I*(1-fd)
         dH <- fd*gamma*I-dday
         dD <- dday
         return(list(c(dS, dE, dI, dR, dH, dD)))
       }
  )
}

runModel <- function(state, ratesFun, parms, times) {
  out <- ode(y = state, times = times, func = ratesFun, parms = parms)
  
}

## Log, Logit ##
inv.log <- exp 
logit <- function(x) {
  log(x/(1-x))
}
inv.logit <- function(x) {
  exp(x)/(1+exp(x))
}

## Root Mean Square Error ## 
RMSE <- function(m,o) {
  sqrt(mean((m-o)^2))
}

#### EMPLOY GOODNESS OF FIT AND OPTIMISATION ALGORITHM. ####

gofFunc <- function(par){
  # Define model parameters.
  R0 = 2.5                  #basic reproduction number 
  gamma = 1/7               #recovery rate
  parms <- c(beta=R0*gamma, #calculate the value transmission rate based on R0 and gamma
             mu=1/7,        #transition to infectious compartment rate
             gamma=gamma,   #recovery rate
             gammaH=1/13,   #inverse of the mean delay from infection to fatality
             fd=0.017,      #infection fatality probability
             dcrit=100,     #half saturation constant for short-term awareness
             dtotcrit=5000, #half saturation constant for long-term awareness
             k=2,           #sharpness of change in the force of infection
             NPI=0.5)       #non-pharmaceutical intervention term value
  
  # Convert parameters for optimisation into values that can be passed on to the optimisation algorithm.
  parms[['k']] = par[['logk']] %>%
    inv.log
  parms[['dcrit']] = par[['logdcrit']] %>%
    inv.log
  parms[['dtotcrit']] = par[['logdtotcrit']] %>%
    inv.log
  parms[['NPI']] = par[['logitNPI']] %>%
    inv.logit
  
  # Define initial conditions.
  initN = 1.109e8
  initE = 16177
  initI = 0
  initR = 0
  initH = 0
  initS = initN - (initE+initI+initH+initR)
  state <- c(S = initS, E = initE, I = initI, R = initR, H = initH, D = 0)
  
  # Define model duration.
  times <- seq(0, 300, by=1)
  
  # Run model.
  model.out <- runModel(state, ratesFun, parms, times) 
  
  # Perform maximum likelihood estimation using COVID-19 weekly death data.
  model.out[,'D'] %>% 
    diff() %>%
    tibble(modelDeaths=.) %>%
    rowid_to_column("day") %>%
    mutate(date=dateEpoch+lubridate::days(day)) %>%
    inner_join(dataDeaths, by='date') %>%
    mutate(Week=week(date), month=month(date), year=year(date)) %>%
    filter(date >= "2020-03-15") %>%
    group_by(year,Week) %>%
    summarise(modelDeaths=mean(modelDeaths), dataDeaths=round(mean(dailyDeaths)), date=mean(date), .groups = "drop") %>%
    mutate(negloglikelihood=-dpois(x=dataDeaths, lambda=modelDeaths, log=T)) %>%
    pull(negloglikelihood) %>%
    sum() 
}

# Perform parameter optimisation.
result <- optim(#initialise parameters
                par=c(
                    logk=log(2), 
                       logdcrit=log(100), 
                       logdtotcrit=log(5000), 
                       logitNPI=logit(0.5)),
                #lower bounds
                 lower=c(logk=log(1), 
                         logdcrit=log(50), 
                         logdtotcrit=log(500), 
                         logitNPI=logit(0.01)),
                #upper bounds
                 upper=c(logk=log(5), 
                         logdcrit=log(500), 
                         logdtotcrit=log(20000), 
                         logitNPI=logit(1)),
                #select limited memory Broyden-Fletcher-Goldfarb-Shanno algorithm
                 method="L-BFGS-B", 
                #pass model function
                 fn=gofFunc,
                #display Hessian matrix
                 hessian=T) 
result

# Estimate parameter confidence interval based on the Hessian matrix.
hessian <- result$hessian
hessian.inv <- solve(hessian)
hessian.inv
parm.se <- sqrt(diag(hessian.inv))
parm.se
CI.matrix <- as.data.frame(matrix(NA, nrow=4, ncol=3))
CI.matrix[,1] <- result$par
CI.matrix[,2] <- result$par - 1.96*parm.se
CI.matrix[,3] <- result$par + 1.96*parm.se
names(CI.matrix) <- c("ML", "95% Lower bound", "95% Upper bound")
rownames(CI.matrix) <- c("logk", "logdcrit", "logdtotcrit", "logitNPI")
CI.matrix

# Convert Hessian matrix-derived parameter confidence interval to their natural units.
Cal.matrix <- as.data.frame(matrix(NA, nrow=4, ncol=3))
Cal.matrix[1,] <- inv.log(CI.matrix[1,])
Cal.matrix[2,] <- inv.log(CI.matrix[2,])
Cal.matrix[3,] <- inv.log(CI.matrix[3,])
Cal.matrix[4,] <- inv.logit(CI.matrix[4,])
names(Cal.matrix) <- c("ML", "95% Lower bound", "95% Upper bound")
rownames(Cal.matrix) <- c("k", "dcrit", "dtotcrit", "NPI")
Cal.matrix

# Print calibrated confidence intervals.
CI <- bind_rows(CI.matrix, Cal.matrix)
write.csv(CI, "CI.csv")

#### GENERATE MODEL CI AND RELEVANT PLOTS. ####

# Define model parameters.
parms <- c(R0=2.5,                      #basic reproduction number
           mu = 1/7,                    #transformation to infectious compartment rate
           gamma = 1/7,                 #recovery rate
           gammaH = 1/13,               #inverse of the mean delay from infection to fatality
           beta=2.5/7,                  #calculated transmission rate
           fd=0.017,                    #infection fatality probability
           k=Cal.matrix[1,1],           #optimised value for the sharpness of change in the force of infection
           dcrit=Cal.matrix[2,1],       #optimised value for the half-saturation constant for short-term awareness
           dtotcrit=Cal.matrix[3,1],    #optimised value for the half-saturation constant for long-term awareness
           NPI=Cal.matrix[4,1])         #optimised value for the non-pharmaceutical intervention term 

# Define initial conditions.
initN = 1.109e8
initE = 16177
initI = 0
initR = 0
initH = 0
initS = initN - (initE+initI+initH+initR)
state <- c(S = initS, E = initE, I = initI, R = initR, H = initH, D = 0)

# Define model duration.
times <- seq(0, 300, by=1)

# Call model function.
model.out <- runModel(state, ratesFun, parms, times) 

# Sample 10,0000 parameter sets from the confidence ranges generated by the Hessian matrix.
k <- runif(10000, Cal.matrix[1,2], Cal.matrix[1,3])
dcrit <- runif(10000, Cal.matrix[2,2], Cal.matrix[2,3])
dtotcrit <- runif(10000, Cal.matrix[3,2], Cal.matrix[3,3])
NPI <-runif(10000, Cal.matrix[4,2], Cal.matrix[4,3])

# Run model simulations 10,000 times using the 10,000 parameter sets.
cl <- makeCluster(8)
registerDoParallel(cl)
modelResults <- foreach(ind = 1:10000, .packages=c("deSolve")) %dopar% {
  parms[['k']] = k[ind]
  parms[['dcrit']] = dcrit[ind]
  parms[['dtotcrit']] = dtotcrit[ind]
  parms[['NPI']] = NPI[ind]
  model.out <- runModel(state, ratesFun, parms, times) 
  return(model.out[,'D'])
} 
modelResults.final <- list()
for(i in 1:length(modelResults[[1]])){
  modelResults.final[[i]]=sapply(modelResults,`[[`, i)
}

# Extract 5%, 50%, and 95% output from each time point.
modelResults.CI <- sapply(modelResults.final, function(deathsAtTimePoint) {
  quantile(deathsAtTimePoint, c(.05,.50,.95)) %>% as_tibble
}) 

# Convert results to dataframe for visualisation.
modelResults.CI = do.call(rbind,modelResults.CI)
rownames(modelResults.CI)=1:length(modelResults[[1]])
colnames(modelResults.CI) = c("CI0.05","CI0.50","CI0.95")
res <- modelResults.CI %>% 
  diff() %>% 
  as.data.frame %>%
  rowid_to_column("day") %>%
  mutate(date=dateEpoch+lubridate::days(day)) %>%
  inner_join(dataDeaths, by='date') %>%
  mutate(Week=week(date), month=month(date), year=year(date)) %>%
  filter(date >= "2020-03-15") %>%
  group_by(year,Week) %>%
  summarise(modelDeaths_0.05=mean(CI0.05), modelDeaths_0.50=mean(CI0.50), modelDeaths_0.95=mean(CI0.95), dataDeaths=round(mean(dailyDeaths)), date=mean(date), .groups = "drop")  
 
# Plot model output.
plot <- ggplot(data=res, aes(x=date)) +
  geom_point(aes(y=dataDeaths), color='red', size=5) +
  geom_line(aes(y=modelDeaths_0.50), color='dark blue', size=1) +
  geom_ribbon(aes(x=date, ymax=modelDeaths_0.95, ymin=modelDeaths_0.05), fill="blue", alpha=0.1) +
  ggtitle("PHILIPPINES: NPI, Short Term, and Long Term Awareness Model") + xlab("Time") + ylab("Reported Deaths") + ylim(0,100) +
  theme(title = element_text(family = 'Avenir', hjust = 0.5, face="bold", size = 16),
        panel.background = element_rect(fill = "#FFF9F2"),
        axis.text = element_text(family = 'Avenir', size = 14))
plot
ggsave(filename="ph_npi.sta.lta.png", plot=plot, dpi=120)


#### PERFORM FORMAL MODEL EVALUATION. ####

# Root mean square error
rmse0.50 <- RMSE(res$modelDeaths_0.50, res$dataDeaths)
rmse0.05 <- RMSE(res$modelDeaths_0.05, res$dataDeaths)
rmse0.95 <- RMSE(res$modelDeaths_0.95, res$dataDeaths)
rmse.matrix <- as.data.frame(matrix(NA, nrow=1, ncol=3))
rmse.matrix[,1] <- rmse0.50
rmse.matrix[,2] <- rmse0.05
rmse.matrix[,3] <- rmse0.95
names(rmse.matrix) <- c("ML", "Lower bound", "Upper bound")
rownames(rmse.matrix) <- c("RMSE")

# Akaike Information Criterion
AIC0.50 = (-2*log(min(dpois(lambda=res$modelDeaths_0.50, x=res$dataDeaths)))) + (2*4)
AIC0.05 = (-2*log(min(dpois(lambda=res$modelDeaths_0.05, x=res$dataDeaths)))) + (2*4)
AIC0.95 = (-2*log(min(dpois(lambda=res$modelDeaths_0.95, x=res$dataDeaths)))) + (2*4)
aic.matrix <- as.data.frame(matrix(NA, nrow=1, ncol=3))
aic.matrix[,1] <- AIC0.50
aic.matrix[,2] <- AIC0.05
aic.matrix[,3] <- AIC0.95
names(aic.matrix) <- c("ML", "Lower bound", "Upper bound")
rownames(aic.matrix) <- c("AIC")

# Bayesian Information Criterion
BIC0.50 = (-2*log(min(dpois(lambda=res$modelDeaths_0.50, x=res$dataDeaths)))) + 4*log(42)
BIC0.05 = (-2*log(min(dpois(lambda=res$modelDeaths_0.05, x=res$dataDeaths)))) + 4*log(42)
BIC0.95 = (-2*log(min(dpois(lambda=res$modelDeaths_0.95, x=res$dataDeaths)))) + 4*log(42)
bic.matrix <- as.data.frame(matrix(NA, nrow=1, ncol=3))
bic.matrix[,1] <- BIC0.50
bic.matrix[,2] <- BIC0.05
bic.matrix[,3] <- BIC0.95
names(bic.matrix) <- c("ML", "Lower bound", "Upper bound")
rownames(bic.matrix) <- c("BIC")

# Print results.
evaluation <- bind_rows(rmse.matrix, aic.matrix, bic.matrix)
write.csv(evaluation, "evaluation.csv")
evaluation

#########################################################################################
# SECTION 2.                                                                            #
# GENERATE MODEL COMPARISON PLOTS PER TERRITORY                                         #
#########################################################################################

#### PERFORM DATA WRANGLING. ####
mc <- read_csv("Model_Comparison.csv")
mc$Models <- factor(mc$Models, levels = c("NPI", "STA", "LTA", "STA + LTA", "NPI + STA", "NPI + LTA", "NPI + STA + LTA"))
mc$test <- factor(mc$test, levels = c("rmse", "aic", "bic", "points"))
ph <- filter(mc, territory == "ph")
ncr <- filter(mc, territory == "ncr")
r3 <- filter(mc, territory == "r3")
r4a <- filter(mc, territory == "r4a")
r6 <- filter(mc, territory == "r6")
r7 <- filter(mc, territory == "r7")

#### GENERATE PLOTS. ####

#Philippines
phplot <- ggplot(data = ph, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "PHILIPPINES - Comparison of Model Performance")
phplot
ggsave(filename="phplot.png", plot=phplot, dpi = 120)

#National Capital Region
ncrplot <- ggplot(data = ncr, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "NATIONAL CAPITAL REGION - Comparison of Model Performance")
ncrplot
ggsave(filename="ncrplot.png", plot=ncrplot, dpi = 120)

#Central Luzon (Region 3)
r3plot <- ggplot(data = r3, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "CENTRAL LUZON (REGION 3) - Comparison of Model Performance")
r3plot
ggsave(filename="r3plot.png", plot=r3plot, dpi = 120)

#CALABARZON (Region 4A)
r4aplot <- ggplot(data = r4a, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "CALABARZON (REGION 4A) - Comparison of Model Performance")
r4aplot
ggsave(filename="r4aplot.png", plot=r4aplot, dpi = 120)

#Western Visayas (Region 6)
r6plot <- ggplot(data = r6, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "WESTERN VISAYAS (REGION 6) - Comparison of Model Performance")
r6plot
ggsave(filename="r6plot.png", plot=r6plot, dpi = 120)

#Central Visayas (Region 7)
r7plot <- ggplot(data = r7, aes(x = test, y = val2)) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#FCF6EE"),
        panel.border = element_rect(fill = NA, colour = "#053276", size = 2, linetype = "solid"), 
        axis.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14),
        legend.background = element_rect(fill = "#FCF6EE"),
        legend.key = element_rect(fill = "white"),
        legend.box.background = element_rect(fill = NA, colour = "#053276", size = 1.5, linetype = "solid")) +
  geom_vline(xintercept = ph$test, col = "#D3E0EE") +
  geom_point(aes(colour=Models, shape=Models), size=8) +
  scale_color_manual(values = c("#34169F", "#FFDF00", "#185A37", "#3C1361", "#800000", "#FFA500", "#03A9F4")) +
  scale_shape_manual(values = c(16, 18, 15, 14, 12, 9, 8)) +
  xlab(label = "Method of Evaluation") +
  scale_x_discrete(label = c("RMSE", "AIC", "BIC", "Points")) +
  ylab(label = "Relative Performance") +
  ggtitle(label = "CENTRAL VISAYAS (REGION 7) - Comparison of Model Performance")
r7plot
ggsave(filename="r7plot.png", plot=r7plot, dpi = 120)

#########################################################################################
# SECTION 3.                                                                            #
# GENERATE SENSITIVITY ANALYSIS HEAT MAPS BY TERRITORY                                  #
#########################################################################################

#### PERFORM DATA WRANGLING. ####
Sensi <- read_csv("Sensi.csv")
Philippines <- filter(Sensi, territory == "Philippines")
NCR <- filter(Sensi, territory == "NCR")
R3 <- filter(Sensi, territory == "Region3")
R4A <- filter(Sensi, territory == "Region4A")
R6 <- filter(Sensi, territory == "Region6")
R7 <- filter(Sensi, territory == "Region7")

#### GENERATE PLOTS. ####

## PHILIPPINES ##
#Outcome: RMSE
ph.rmse <- ggplot(data = Philippines, mapping = aes(x = th,
                                                    y = inite,
                                                    fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(4044, 8089, 12133, 16177)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 6.396513) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Philippines, Root Mean Square Error")
ph.rmse
ggsave(filename="ph.rmse.png", plot=ph.rmse, dpi = 120)

#Outcome: AIC
ph.aic <- ggplot(data = Philippines, mapping = aes(x = th,
                                                   y = inite,
                                                   fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(4044, 8089, 12133, 16177)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 20.1728795) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Philippines, Akaike Information Criterion")
ph.aic
ggsave(filename="ph.aic.png", plot=ph.aic, dpi = 120)

#Outcome: BIC
ph.bic <- ggplot(data = Philippines, mapping = aes(x = th,
                                                   y = inite,
                                                   fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(4044, 8089, 12133, 16177)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 26.723227) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Philippines, Bayesian Information Criterion")
ph.bic
ggsave(filename="ph.bic.png", plot=ph.bic, dpi = 120)

#Outcome: Points
ph.points <- ggplot(data = Philippines, mapping = aes(x = th,
                                                      y = inite,
                                                      fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(4044, 8089, 12133, 16177)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.486842106) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Philippines, Proportion of Observations within Model CI")
ph.points
ggsave(filename="ph.points.png", plot=ph.points, dpi = 120)


## NATIONAL CAPITAL REGION ##
#Outcome: RMSE
ncr.rmse <- ggplot(data = NCR, mapping = aes(x = th,
                                             y = inite,
                                             fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(2941, 5883, 8824, 11765)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 6.6467945) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "National Capital Region, Root Mean Square Error")
ncr.rmse
ggsave(filename="ncr.rmse.png", plot=ncr.rmse, dpi = 120)

#Outcome: AIC
ncr.aic <- ggplot(data = NCR, mapping = aes(x = th,
                                            y = inite,
                                            fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(2941, 5883, 8824, 11765)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 26.446926) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "National Capital Region, Akaike Information Criterion")
ncr.aic
ggsave(filename="ncr.aic.png", plot=ncr.aic, dpi = 120)

#Outcome: BIC
ncr.bic <- ggplot(data = NCR, mapping = aes(x = th,
                                            y = inite,
                                            fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(2941, 5883, 8824, 11765)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 32.997271) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "National Capital Region, Bayesian Information Criterion")
ncr.bic
ggsave(filename="ncr.bic.png", plot=ncr.bic, dpi = 120)

#Outcome: Points
ncr.points <- ggplot(data = NCR, mapping = aes(x = th,
                                               y = inite,
                                               fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(2941, 5883, 8824, 11765)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.236842106) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "National Capital Region, Proportion of Observations within Model CI")
ncr.points
ggsave(filename="ncr.points.png", plot=ncr.points, dpi = 120)

## CENTRAL LUZON (REGION 3) ##
#Outcome: RMSE
r3.rmse <- ggplot(data = R3, mapping = aes(x = th,
                                           y = inite,
                                           fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 2.0646757) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Luzon, Root Mean Square Error")
r3.rmse
ggsave(filename="r3.rmse.png", plot=r3.rmse, dpi = 120)

#Outcome: AIC
r3.aic <- ggplot(data = R3, mapping = aes(x = th,
                                          y = inite,
                                          fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 15.35021555) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Luzon, Akaike Information Criterion")
r3.aic
ggsave(filename="r3.aic.png", plot=r3.aic, dpi = 120)

#Outcome: BIC
r3.bic <- ggplot(data = R3, mapping = aes(x = th,
                                          y = inite,
                                          fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 21.90055985) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Luzon, Bayesian Information Criterion")
r3.bic
ggsave(filename="r3.bic.png", plot=r3.bic, dpi = 120)

#Outcome: Points
r3.points <- ggplot(data = R3, mapping = aes(x = th,
                                             y = inite,
                                             fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.578947369) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Luzon, Proportion of Observations within Model CI")
r3.points
ggsave(filename="r3.points.png", plot=r3.points, dpi = 120)

## CALABARZON (REGION 4A) ##
#Outcome: RMSE
r4a.rmse <- ggplot(data = R4A, mapping = aes(x = th,
                                             y = inite,
                                             fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 1.843911) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "CALABARZON, Root Mean Square Error")
r4a.rmse
ggsave(filename="r4a.rmse.png", plot=r4a.rmse, dpi = 120)

#Outcome: AIC
r4a.aic <- ggplot(data = R4A, mapping = aes(x = th,
                                            y = inite,
                                            fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 14.485824) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "CALABARZON, Akaike Information Criterion")
r4a.aic
ggsave(filename="r4a.aic.png", plot=r4a.aic, dpi = 120)

#Outcome: BIC
r4a.bic <- ggplot(data = R4A, mapping = aes(x = th,
                                            y = inite,
                                            fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 21.241342) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "CALABARZON, Bayesian Information Criterion")
r4a.bic
ggsave(filename="r4a.bic.png", plot=r4a.bic, dpi = 120)

#Outcome: Points
r4a.points <- ggplot(data = R4A, mapping = aes(x = th,
                                               y = inite,
                                               fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(1103, 2206, 3309, 4412)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.355263158) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "CALABARZON, Proportion of Observations within Model CI")
r4a.points
ggsave(filename="r4a.points.png", plot=r4a.points, dpi = 120)

## WESTERN VISAYAS (REGION 6) ##
#Outcome: RMSE
r6.rmse <- ggplot(data = R6, mapping = aes(x = th,
                                           y = inite,
                                           fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 1.4958185) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Western Visayas, Root Mean Square Error")
r6.rmse
ggsave(filename="r6.rmse.png", plot=r6.rmse, dpi = 120)

#Outcome: AIC
r6.aic <- ggplot(data = R6, mapping = aes(x = th,
                                          y = inite,
                                          fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 13.604533) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Western Visayas, Akaike Information Criterion")
r6.aic
ggsave(filename="r6.aic.png", plot=r6.aic, dpi = 120)

#Outcome: BIC
r6.bic <- ggplot(data = R6, mapping = aes(x = th,
                                          y = inite,
                                          fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 18.933351) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Western Visayas, Bayesian Information Criterion")
r6.bic
ggsave(filename="r6.bic.png", plot=r6.bic, dpi = 120)

#Outcome: Points
r6.points <- ggplot(data = R6, mapping = aes(x = th,
                                             y = inite,
                                             fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.464285714) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Western Visayas, Proportion of Observations within Model CI")
r6.points
ggsave(filename="r6.points.png", plot=r6.points, dpi = 120)

## CENTRAL VISAYAS (REGION 7) ##
#Outcome: RMSE
r7.rmse <- ggplot(data = R7, mapping = aes(x = th,
                                           y = inite,
                                           fill = rmse)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "RMSE",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 3.353555) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Visayas, Root Mean Square Error")
r7.rmse
ggsave(filename="r7.rmse.png", plot=r7.rmse, dpi = 120)

#Outcome: AIC
r7.aic <- ggplot(data = R7, mapping = aes(x = th,
                                          y = inite,
                                          fill = aic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "AIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 19.341073) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Visayas, Akaike Information Criterion")
r7.aic
ggsave(filename="r7.aic.png", plot=r7.aic, dpi = 120)

#Outcome: BIC
r7.bic <- ggplot(data = R7, mapping = aes(x = th,
                                          y = inite,
                                          fill = bic)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "BIC",
                       low = "#008FF5",
                       mid = "white",
                       high = "#C2320A",
                       midpoint = 25.8914175) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Visayas, Bayesian Information Criterion")
r7.bic
ggsave(filename="r7.bic.png", plot=r7.bic, dpi = 120)

#Outcome: Points
r7.points <- ggplot(data = R7, mapping = aes(x = th,
                                             y = inite,
                                             fill = points)) +
  geom_tile() +
  xlab(label = "Mean Time from Infection to Fatality (days)") +
  scale_x_discrete(limits = c(7,14,21)) +
  ylab(label = "Initial Seed Value") +
  scale_y_discrete(limits = c(354, 709, 1063, 1417)) +
  facet_wrap(~ Basic.Reproduction.Number, strip.position = "top",
             labeller = labeller(Basic.Reproduction.Number =
                                   c("2.5" = "Basic Reproduction Number: 2.5",
                                     "3" = "Basic Reproduction Number: 3"))) +
  scale_fill_gradient2(name = "Points",
                       low = "#C2320A",
                       mid = "white",
                       high = "#008FF5",
                       midpoint = 0.342105263) +
  theme(title = element_text(family = 'Avenir', size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
        panel.background = element_rect(fill = "#E4F1FF"),
        strip.background = element_rect(fill = "#F2D6CE", color = "#FFFFFF"),
        axis.text = element_text(family = 'Avenir', size = 16),
        strip.text = element_text(family = 'Avenir', size = 16, face = "bold"),
        strip.placement = "outside",
        legend.key.size = unit(2, "line"),
        legend.text = element_text(family = 'Avenir', size = 14)) +
  ggtitle(label = "Central Visayas, Proportion of Observations within Model CI")
r7.points
ggsave(filename="r7.points.png", plot=r7.points, dpi = 120)


#### NOTHING FOLLOWS ####