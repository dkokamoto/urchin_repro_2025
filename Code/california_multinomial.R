library(dplyr)
library(tidyr)
library(gdata)
library(rstan)
library(ggplot2)

### read in data
data <- read.csv("Data/pop_comparison_histology.csv")

### format treatments for regressions
data <- data %>%
  mutate(treatment3 = factor(location):factor(sex):factor(treatment),
         treatment4 = factor(treatment):factor(sex),
         treatment5 = factor(sex):factor(location))

### create aggregated data for visualization
data_sum <- data %>%
  group_by(treatment,location,sex,cat) %>%
  dplyr::summarize(N= length(cat))

data_sum <- with(data_sum,
              expand.grid(treatment = unique(treatment),
                          location= unique(location),
                          sex = unique(sex),
                          cat = unique(cat)
              ))%>%
  full_join(data_sum)%>%
  mutate(N= ifelse(is.na(N),0,N))

data_sum <- data_sum %>%
  group_by(treatment,location,sex) %>%
  dplyr::summarize(N_tot= sum(N))%>%
  left_join(data_sum)

### create projection dataframes
newdata <- data %>% 
  select(treatment,location,sex,treatment3,treatment4,treatment5) %>%
  unique() %>%
  mutate(row= 1:nrow(.))

newdata2<- data %>% filter(treatment%in%c(10,20,"HW")) %>%
  select(treatment,sex,treatment4) %>%
  unique() %>%
  mutate(row= 1:nrow(.))

### compile Stan model
mod <- stan_model("Code/reg_multinom.stan") 

### model 1: full model with all interactions
### create model matrices for fixed (X) and random (Z) effects
X = model.matrix(~ treatment3-1, data = data)
Z = model.matrix(~ location, data = data)

### create projection model matrices for fixed (X) and random (Z) effects
X_pred = model.matrix(~ treatment3-1, data = newdata)
Z_pred = model.matrix(~ location, data = newdata)

### create the data for Stan
stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  X = X,
  Z = Z,
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

### sample the posterior
fit_mult <- sampling(
  mod,
  data = stan_data,
  iter = 12000,
  warmup = 2000,
  chains = 4,
  cores = 4
)

### run the bridgesampler and the LOO
b1 <- bridge_sampler(fit_mult, silent = TRUE)
loo1 <- loo(fit_mult)

### extract and create the dataframe for plotting
draws <- rstan::extract(fit_mult,"mu_pred")
pred <- draws$mu_pred
dquants <- plyr::adply(apply(draws$mu_pred,c(2,3),function(x) c(quantile(x,probs= c(0.025,0.1,0.5,0.9,0.975)),mean(x))),c(2,3))
names(dquants) <-  c("row1","cat","CIL2","CIL","M","CIU","CIU2","mu")
dquants <-cbind(dquants,newdata)

### model 2: model with only treatment and Sex with interactions and random locaton effects
X = model.matrix(~ treatment4-1, data = data)
Z = model.matrix(~ location, data = data)
X_pred = model.matrix(~ treatment4-1, data = newdata2)
Z_pred = matrix( 1 / 3, ncol = 3, nrow= nrow(newdata2))

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  X = X,
  Z = Z,
  re= 0,
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

fit_mult2 <- sampling(
  mod,
  data = stan_data,
  iter = 12000,
  warmup = 2000,
  chains = 4,
  cores = 4
)

b2 <- bridge_sampler(fit_mult2, silent = TRUE)

draws <- rstan::extract(fit_mult2,"mu_pred")
pred_trt <- draws$mu_pred
dquants2 <- plyr::adply(apply(draws$mu_pred,c(2,3),function(x) c(quantile(x,probs= c(0.025,0.1,0.5,0.9,0.975)),mean(x))),c(2,3))
names(dquants2) <-  c("row1","cat","CIL2","CIL","M","CIU","CIU2","mu")
dquants2 <-cbind(dquants2,newdata2)

### model 3: model with only location and sex with interactions and random locatoon effects
X = model.matrix(~ treatment5-1, data = data)
Z = model.matrix(~ location, data = data)
X_pred = model.matrix(~ treatment5-1, data = newdata)
Z_pred = matrix( 1 / 3, ncol = 3, nrow= nrow(newdata))

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  X = X,
  Z = Z,
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

fit_mult3 <- sampling(
  mod,
  data = stan_data,
  iter = 12000,
  warmup = 2000,
  chains = 4,
  cores = 4
)

b3 <- bridge_sampler(fit_mult3, silent = TRUE)

### model 4: model with only sex and random location effects
X = model.matrix(~ sex-1, data = data)
Z = model.matrix(~ location, data = data)
X_pred = model.matrix(~ sex-1, data = newdata2)
Z_pred = matrix( 1 / 3, ncol = 3, nrow= nrow(newdata2))

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  G = ncol(Z),
  X = X,
  Z = Z,
  Y = data$cat,
  ncat = 4,
  N_pred = nrow(X_pred),
  X_pred = X_pred,
  Z_pred = Z_pred
)

fit_mult4 <- sampling(
  mod,
  data = stan_data,
  iter = 12000,
  warmup = 2000,
  chains = 4,
  cores = 4
)

b4 <- bridge_sampler(fit_mult4, silent = TRUE)

### bayes factor comparison of models 1-4
bayes_factor(b1,b2)
bayes_factor(b1,b3)
bayes_factor(b2,b3)
bayes_factor(b2,b4)
post_prob(b1,b2,b3,b4)

### loo comparison of models 1-4
model_list <- list(fit_mult, fit_mult2, fit_mult3, fit_mult4)
log_lik_list <- lapply(model_list, loo::extract_log_lik)

r_eff_list <- lapply(model_list, function(x) {
  ll_array <- loo::extract_log_lik(x, merge_chains = FALSE)
  loo::relative_eff(exp(ll_array))
})

mod_weights <- loo_model_weights(
  log_lik_list,
  method = "stacking",
  r_eff_list = r_eff_list,
  optim_control = list(reltol=1e-10)
)

mod_weights

dquants$location <- factor(dquants$location, levels= c("SB","SD","SON"))
data_sum$location  <- factor(data_sum$location, levels= c("SB","SD","SON"))
data_sum2$location  <- factor(data_sum2$location, levels= c("SB","SD","SON"))

### plot posteriors for model 1 and 2
p1 <- ggplot(aes(y= M, 
                 x= as.numeric(factor(treatment))+(as.numeric(factor(location))-1.5)*0.1), 
             data= dquants)+
  geom_point(aes(x= as.numeric(factor(treatment))),
             size= 2.5,alpha= 0.5,data = dquants2)+
  geom_linerange(aes(x= as.numeric(factor(treatment)), ymin= CIL2, ymax = CIU2),
                 linewidth=.2,alpha=0.2,data = dquants2)+
  geom_linerange(aes(x= as.numeric(factor(treatment)),ymin= CIL, ymax = CIU),
                 linewidth=1,alpha=0.2,data = dquants2)+
  geom_line(aes(x= as.numeric(factor(treatment))),
            data = dquants2)+
  geom_point(aes(colour= factor(location),
                 group= factor(location),
                 shape= factor(location)),
             size= 1.5,alpha= 0.5)+
  geom_linerange(aes(ymin= CIL2, ymax = CIU2,
                     colour= factor(location),
                     group= factor(location)),
                 linewidth=.2,alpha=0.2)+
  geom_linerange(aes(ymin= CIL, ymax = CIU,
                     colour= factor(location),
                     group= factor(location)),
                 linewidth= .7,alpha=0.2)+
  geom_line(aes(colour= factor(location),
                group= factor(location),
                linetype= factor(location)),
            linewidth=0.2)+
  scale_colour_manual(values= c("red","grey40","blue"))+
  scale_fill_manual(values= c("red","grey40","blue"))+
  scale_shape_manual(values=c(21,22,23))+
  xlab("treatment")+
  ylab("Proportion")+
  geom_point(aes(y= N/N_tot,colour = factor(location),
                 x = as.numeric(factor(treatment)) + (as.numeric(factor(location))-1.5)*0.1),
             data= data_sum,
             shape= 8,
             size = 1)+
  #geom_text(aes(y= 0.82,colour = factor(location),
  #              label= N_tot,
  #              x = as.numeric(factor(treatment)) + (as.numeric(factor(location))-2.5)*0.2),
  #          data= data_sum2,
  #          size = 2)+
  facet_grid(cat~sex)+
  scale_x_continuous(breaks= c(1,2,3),labels= c("10", "20","HW"))+
  theme_bw()+
  scale_linetype_manual(values = c(2,3,4))+
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.background= element_blank(),
        legend.box = element_blank(),
        legend.position= c(0.4,0.67),
        legend.direction= "vertical",
        legend.title= element_blank(),
        legend.key.width = unit(1,"cm"))+
  scale_y_continuous(limits= c(0,0.88), 
                     breaks= c(0,0.2,0.4,0.6,0.8,1.0),
                     labels= c(0,0.2,0.4,0.6,0.8,""), 
                     expand= c(0,0));p1

#pdf(width = 5, height = 7, file= "Figures/pop_comparison.pdf")
#p1
#dev.off()
