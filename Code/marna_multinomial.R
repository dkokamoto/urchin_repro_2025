### some terrible hastily written code by DKO.  Annotations to follow
library(ggplot2)
library(mgcv)
library(rstan)
library(bridgesampling)
library(lubridate)
library(dplyr)
library(brms)
library(tidyr)

### load biometric data
data <- read.csv("Data/Biometrics_Final.csv") %>%
  mutate(nads = dry_gonad_plus_wb_mass_g - gonad_wb_mass_g,
         nads_wet = wet_gonad_mass_g*5) %>% 
  mutate(mean_temp = ifelse(temp_deg_c == "21-18",20,ifelse(temp_deg_c == "18-14",16,temp_deg_c)),
         treatment = ifelse(temp_deg_c == "21-18","El Nino",ifelse(temp_deg_c == "18-14","La Nina", "static")))
data$mean_temp= as.numeric(data$mean_temp)

histo_biodata <- data %>%
  #filter(treatment != "static")%>%
  select(pco2_h_l,treatment,tank_id,temp_deg_c,HISTO.CASETTE.NO, nads,nads_wet)

### read in histology data
histo_data <- read.csv("Data/Marna_histology_2021.csv") %>%
  mutate(temp= ifelse(temp_deg_c == "21-18",20, ifelse(temp_deg_c == "18-14",16,temp_deg_c)),
         stage = cat)

### summarize data for plotting
histo_summary <- histo_data %>%
  filter(pco2_h_l =="l")%>%
  data.frame()%>%
  group_by(tank_id,cat,treatment,sex,temp,pco2_h_l) %>%
  dplyr::summarize(N = length(cat))

histo_summary <- histo_summary %>%
  group_by(treatment,tank_id,sex,temp,pco2_h_l) %>%
  dplyr::summarize(TOT = sum(N)) %>%
  left_join(histo_summary)

histo_summary <- histo_data %>%
  select(-stage)%>%
  group_by(cat,treatment,sex,temp_deg_c,pco2_h_l,temp,tank_id)%>%
  dplyr::summarize(N= length(cat)) %>%
  data.frame()%>%
  pivot_wider(names_from = cat,values_from = N, values_fill= 0) 

histo_summary$N <- rowSums(histo_summary[,7:10])
names(histo_summary)[7:10] <- c("x1","x2","x3","x4")
histo_summary$y <- with(histo_summary,cbind(x1,x2,x3,x4))

### subset to only low pCO2 treatments
data =  subset(histo_data, pco2_h_l =="l")

### specify gradient for the Gaussian process
gpx = as.numeric(seq(10,20,by = 0.5))

### create treatment subset
data <- data %>%
  mutate(Treatment2 = factor(treatment):factor(sex))

### create projection datasets
newdata <- data %>% 
  select(temp,treatment,sex) %>%
  unique() %>%
  filter(treatment != "static") 

newdata <- expand.grid(temp= gpx, 
                       treatment= "static",
                       sex= c("M","F")) %>%
  rbind(newdata) %>%
  mutate(Treatment2 = factor(factor(treatment):factor(sex),levels= levels(data$Treatment2)))

### model 1: model with only treatment by sex with random mesocosm effects
X = model.matrix(~  Treatment2 - 1, data = data)
Z = model.matrix(~ factor(tank_id), data = data)

X_pred = model.matrix( ~ Treatment2-1, data = newdata)
Z_pred = matrix( 1 / 20, ncol = 20, nrow= nrow(newdata))

M_ID1 <- (1:nrow(data))[data$sex == "M"]
F_ID1 <- (1:nrow(data))[data$sex == "F"]
M_ID2 <- match(data$temp[M_ID1],gpx)
F_ID2 <-match(data$temp[F_ID1],gpx)

M_ID1_pred <- (1:nrow(newdata))[newdata$sex == "M"]
F_ID1_pred <- (1:nrow(newdata))[newdata$sex == "F"]
M_ID2_pred <- match(newdata$temp[M_ID1_pred],gpx)
F_ID2_pred <-match(newdata$temp[F_ID1_pred],gpx)

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  N_pred = nrow(newdata),
  X = X,
  Z = Z,
  gpx = gpx,
  Ngp = length(gpx),
  N_gp_ID = c(length(M_ID1),length(F_ID1)),
  N_gp_ID_pred = c(length(M_ID1_pred),length(F_ID1_pred)),
  ID_gp1 = rbind(M_ID1,M_ID2),
  ID_gp2 = rbind(F_ID1,F_ID2),
  ID_pred_gp1 = rbind(M_ID1_pred,M_ID2_pred),
  ID_pred_gp2 = rbind(F_ID1_pred,F_ID2_pred),
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

### model 1: model treatment by sex with random mesocosm effects and GP of temperature
mod_gp <- stan_model("Code/reg_multinom_GP.stan") 

fit_mult <- sampling(
  mod_gp,
  data = stan_data,
  iter = 12000,
  warmup = 1000,
  chains = 4,
  cores = 4
)

X = model.matrix(~  sex  - 1, data = data)
Z = model.matrix(~ factor(tank_id), data = data)

X_pred = model.matrix( ~ sex -1, data = newdata)
Z_pred = matrix( 1 / 20, ncol = 20, nrow= nrow(newdata))

### model 2: model with only sex with random mesocosm effects
stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  N_pred = nrow(newdata),
  X = X,
  Z = Z,
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

mod <- stan_model("Code/reg_multinom.stan") 

fit_mult2 <- sampling(
  mod,
  data = stan_data,
  iter = 12000,
  warmup = 1000,
  chains = 4,
  cores = 4
)


b1 <- bridge_sampler(fit_mult, silent = TRUE)
b2 <- bridge_sampler(fit_mult2, silent = TRUE)

### bayes factor comparison of models 1-4
bayes_factor(b1,b2)

### extract predictons from full model
p1 <- rstan::extract(fit_mult,"mu_pred")$mu_pred

### format prediction summaries
test <- apply(p1,c(2,3),quantile,probs= c(0.025,0.1,0.5,0.9,0.975))
test <- plyr::adply(test,c(2,3))
names(test)[3:7] <- c("CIL","CIL2","med","CIU2","CIU")
newdata2 <- cbind(test,newdata)
newdata2$temp <- as.numeric(newdata2$temp)
histo_summary$temp <- as.numeric(histo_summary$temp)
histo_summary$X2 <- factor(histo_summary$cat)

### plot model predictions and raw data
p1 <- ggplot(aes(temp,med), 
       data= subset(newdata2, treatment =="static"))+
  geom_ribbon(aes(ymin= CIL2,ymax =CIU2, fill= X2),alpha=0.5)+
  geom_ribbon(aes(ymin= CIL,ymax =CIU, fill= X2),alpha= 0.1)+
  geom_line(aes(y= med, colour= X2,linetype= X2),size= 1,alpha= 0.5)+
  facet_wrap(ifelse(X2==4,1,2)~sex)+
  geom_linerange(aes(x = temp - 0.4 + as.numeric(as.factor(X2))*0.1,
                     ymin= CIL,ymax =CIU, colour= X2),
                 alpha = 0.5,
                 subset(newdata2, treatment !="static"))+
  geom_linerange(aes(x = temp - 0.4 + as.numeric(as.factor(X2))*0.1,
                     ymin= CIL2,ymax =CIU2, colour= X2),
                 alpha = 0.3,
                 subset(newdata2, treatment !="static"),linewidth= 2)+
  geom_point(aes(y= N/TOT,x= jitter(temp,0.5),fill= X2,shape= treatment, size= treatment), 
             data= subset(histo_summary,pco2_h_l=="l"&treatment=="static"))+
  geom_point(aes(y= N/TOT,x= jitter(temp,0.5),colour= X2,shape= treatment, size= treatment), 
             data= subset(histo_summary,pco2_h_l=="l"&treatment!="static"))+
  facet_grid(ifelse(X2=="4","stage IV","stages I-III")~sex,scales= "free_y")+
  geom_point(aes( x = temp - 0.4 + as.numeric(as.factor(X2))*0.1,
                  shape= treatment,fill= X2),
             subset(newdata2, treatment !="static"),alpha=0.7,size= 3)+
  scale_y_continuous(breaks= c(0,0.2,0.4,0.6,0.8,1))+
  geom_hline(yintercept= 0)+
  scale_shape_manual(values=c(22,23,21),name= "TRT")+
  xlab(expression(paste("temperature (",degree,"C)")))+
  gg_options()+
  ylab("proportion")+
  theme(legend.position = c(0.2,0.4),
        legend.direction = "vertical",
        legend.key.width= unit(2,"cm"))+
  guides(guide_legend = list(fill = FALSE))+
  scale_size_manual(values= c(2,2,1.5),name= "TRT")+
  scale_linetype_manual(values= c(3,2,1,1),name= "TRT")+
  scale_fill_manual(values=c("red","grey","purple","blue"),name= "TRT")+
  scale_colour_manual(values=c("red","grey","purple","blue"),name= "TRT")+
  scale_x_continuous(breaks=c(10,13,16,17,18,20)); p1