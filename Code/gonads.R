### some terrible hastily written code by DKO.  Annotations to follow
library(ggplot2)
library(mgcv)
library(lubridate)
library(dplyr)
library(brms)
library(tidyr)

### read in biometric data
data <- read.csv("Data/Biometrics_Final.csv") %>%
  mutate(nads_wet = wet_gonad_mass_g*5) %>% 
  mutate(mean_temp = ifelse(temp_deg_c == "21-18",20,ifelse(temp_deg_c == "18-14",16,temp_deg_c)),
         treatment = ifelse(temp_deg_c == "21-18","El Nino",ifelse(temp_deg_c == "18-14","La Nina", "static")),
         mean_temp= as.numeric(mean_temp))

### set brms priors
priors2 <- c(set_prior("student_t(3,0,2.5)", class="Intercept"),
             set_prior("normal(0,1)", class="b"),
             set_prior("normal(0,1)", class="sd"),
             set_prior("cauchy(0,1)",class= "shape"))

### define and sample the model
fit_full <- brm((nads_wet)~ poly(mean_temp, 3) * sex + log(diam_mm) +
             treatment*sex+(1|tank_id), 
           prior = priors2,
           iter = 12000,
           warmup = 2000,
           data= subset(data, pco2_h_l =="l"&nads_wet>0),
           family=Gamma(link="log"),
           chains = 4, cores= 4,
           save_pars = save_pars(all = TRUE)) 

fit <-  brm((nads_wet)~poly(mean_temp, 3) + log(diam_mm) +
              treatment+(1|tank_id), 
            prior = priors2,
            iter = 12000,
            warmup = 2000,
            data= subset(data, pco2_h_l =="l"&nads_wet>0),
            family=Gamma(link="log"),
            chains = 4, cores= 4,
            save_pars = save_pars(all = TRUE)) 

### subset priors for null model
priors3 <-c(set_prior("student_t(3,0,2.5)", class="Intercept"),
            set_prior("normal(0,1)", class="sd"),
            set_prior("cauchy(0,1)",class= "shape"))

fit_null <- brm((nads_wet)~1+(1|tank_id) + log(diam_mm), 
                prior = priors3,
                iter = 12000,
                warmup = 2000,
                data= subset(data, pco2_h_l =="l"&nads_wet>0),
                family=Gamma(link="log"),
                chains = 4, cores= 4,
                save_pars = save_pars(all = TRUE))  

### conduct bridgesampling
b1 <- bridge_sampler(fit)
b2 <- bridge_sampler(fit_null)
b3 <- bridge_sampler(fit_full)

bayes_factor(b1, b2)
bayes_factor(b1, b3)

### plot predictions and se.
newdata <- data.frame(list(mean_temp= seq(10,20,by= 0.1), diam_mm = 50, treatment= "static"))
pred1 <- posterior_epred(fit,newdata= newdata, re.form = NA)
pred <- t(apply(pred1,2,quantile,probs= c(0.025,0.05,0.1,0.5,0.9,0.95,0.975)))
newdata <- cbind(newdata,pred)
temp <- unique(subset(subset(data, pco2_h_l =="l"&nads_wet>0),treatment!= "static")[,c("treatment","mean_temp")])
temp$diam_mm <- 50
pred2 <-posterior_epred(fit,newdata= temp, re.form = NA)
pred <- t(apply(pred2,2,quantile,probs=  c(0.025,0.05,0.1,0.5,0.9,0.95,0.975)))
newdata2 <- cbind(temp,pred)
newdata <- rbind(newdata,newdata2)

names(newdata)[4:10] <- c("CIL","CIL2","CIL3","Med","CIU3","CIU2","CIU")

# ggplot options for visual appeal:
gg_options <- function() theme_bw()+theme(
  panel.grid=element_blank(), # removes ugly grid lines
  plot.background= element_blank(), # removes ugly plot background
  strip.background=element_blank(), # removes ugly strip background
  panel.background= element_blank(), # removes ugly panel background
  legend.title=element_blank(), # not a fan of legend titles if not needed
  legend.background= element_blank(), # removes ugly panel background
  legend.box.background=element_blank(), # removes legend box background
  legend.key= element_blank()) #

p1 <-  ggplot(aes( x =  as.numeric(as.character(mean_temp))+
                     as.numeric(ifelse(treatment!="static",-0.2,0))), 
              data= subset(newdata, treatment=="static"))+ 
  geom_ribbon(aes(ymin = CIL/5,ymax= CIU/5),fill= "black", 
              alpha =0.25)+
  geom_ribbon(aes(ymin = CIL2/5,ymax= CIU2/5),fill= "black",
              alpha =0.25)+
  geom_ribbon(aes(ymin = CIL3/5,ymax= CIU3/5),fill= "black",
              alpha =0.25, guide= "none")+
  geom_line(aes(y= Med/5,colour= treatment),data= subset(newdata, treatment=="static"))+
  geom_linerange(aes(ymin = CIL/5,ymax= CIU/5,
                     col= treatment),
                 data=subset(newdata,treatment!="static"))+
  geom_linerange(aes(ymin = CIL3/5,ymax= CIU3/5,
                     col= treatment),size =2,alpha=0.75,
                 data=subset(newdata,treatment!="static"))+
  geom_point(aes(y= wet_gonad_mass_g*5/(diam_mm^2.36)*(50^2.36),x = jitter(mean_temp),
                 colour = treatment,fill= treatment, 
                 group = tank_id, shape= treatment),alpha= 0.75,
            data= subset(data, pco2_h_l =="l"&nads_wet>0), size= 0.5,alpha=0.1)+
  stat_summary(aes(y=wet_gonad_mass_g*5/(diam_mm^2.36)*(50^2.36),fill= treatment, shape= treatment),
               subset(data, pco2_h_l =="l"&nads_wet>0), geom= "point", size= 3)+   # geom_point(aes(y= mean/N,colour= treatment, group = tank_id),subset(data2, pco2_h_l =="l" &treatment!="static"))+
  scale_colour_manual(values= c("red","blue","black"),name= "")+
  scale_fill_manual(values= c("red","blue","black"),name= "")+
  scale_size_manual(values= c(1.5,1.5,0.5), name= "")+
  scale_shape_manual(values= c(23,22,21), name= "")+
  scale_y_continuous(expand= c(0,0),
                     breaks = c(0,5,10,15,20,25,30,35),
                     limits= c(0,22))+
  theme_bw()+
  theme(panel.grid.major.y = element_line(size= 0.1))+
  scale_x_continuous(breaks= c(10,13,16,17,18,20))+
  theme_bw()+
  gg_options()+
  ylab("gonad mass (g)")+
  xlab(expression(paste("Temperature (", degree, "C)")))+
  theme(legend.position= c(0.2,0.2));p1

#pdf("Figures/nads.pdf",width = 4,height= 5)
#p1
#dev.off()

