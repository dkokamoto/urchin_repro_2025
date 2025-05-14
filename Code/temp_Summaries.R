###########################################################################
###  (Sloppy) Script to generate plots and running means for temps
###  From Quadra, SIO, and Van Damme
###  Author D.K. Okamoto
###  Updated Sept 22, 2021
###########################################################################

library(ggplot2);library(dplyr);library(lubridate)

### load quadra time series
quadra <- read.csv("Data/temperature_time_series/Quadra_temp_data.csv", skip = 7)%>%
  mutate(date = as.Date(measurementTime, "%Y-%m-%d %H:%M:%S"))%>%
  dplyr::select(WaterTemp,date,year)%>%
  mutate(MONTH = month(date),
         DAY= day(date))%>%
  rename(temp = WaterTemp, YEAR = year)%>%
  mutate(site= "Quadra")%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))

### load Van Damme time series
VD <- read.csv("Data/temperature_time_series/VD_temp_data.csv")%>%
  mutate(Date = parse_date_time(Date, "%m/%d/%y"))%>%
  rename(temp= Temp,
         date= Date)%>%
  mutate(YEAR = year(date),
         MONTH= month(date),
         DAY= day(date))%>%
  mutate(site= "Van Damme")%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))


### load BML time series
BML <- read.csv("Data/temperature_time_series/BML_temp_data.csv")%>%
  mutate(Date = parse_date_time(Date, "%m/%d/%y"))%>%
  rename(temp= wtemp,
         date= Date)%>%
  mutate(YEAR = year(date),
         MONTH= month(date),
         DAY= day(date))%>%
  mutate(site= "Sonoma")%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))%>%
  group_by(date, YEAR,MONTH,DAY,YEAR2,date2,site)%>%
  summarize(temp= mean(temp))%>%
  data.frame()

### load SBC time series
SB<- read.csv("Data/temperature_time_series/SBC_temp_data.csv")%>%
  mutate(Date = parse_date_time(DATE_LOCAL, "%m/%d/%y"))%>%
  rename(temp= TEMP_C,
         date= Date)%>%
  mutate(YEAR = year(date),
         MONTH= month(date),
         DAY= day(date))%>%
  mutate(site= "Santa Barbara")%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))%>%
  group_by(date, YEAR,MONTH,DAY,YEAR2,date2,site) %>%
  summarize(temp= mean(temp))%>%
  data.frame()

### load Scripps pier time series
SIO <- read.csv("Data/temperature_time_series/SIO_temp_data.csv")%>%
  mutate(date = parse_date_time(paste(YEAR,MONTH,DAY),"%Y-%m-%d"))%>%
  dplyr::select(date,YEAR,DAY,MONTH,BOT_TEMP_C)%>%
  rename(temp = BOT_TEMP_C)%>%
  mutate(site= "San Diego")%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))

### load Amphitrite Point lighthouse data
AP <- read.csv("Data/temperature_time_series/Amphitrite_Point_temp_data.csv") %>%
  pivot_longer(cols= 2:13, names_to = "MONTH", values_to = "temp") %>%
  mutate(temp = ifelse(temp >= 999,NA,temp),
         DAY= 15,
         date = parse_date_time(paste(DAY,MONTH,YEAR, sep= "-"), orders= "%d-%B-%Y"),
         MONTH = month(date),
         site = "Amphitrite Point") %>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))

  

### join time series

temp_data <- full_join(SIO,SB)%>%
  full_join(quadra) %>%
  full_join(BML) %>%
  full_join(AP) %>%
  group_by(site,MONTH,YEAR,DAY) %>%
  dplyr::summarize(temp= mean(temp,na.rm=T))%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,DAY, sep= "-"), "%Y-%m-%d"))

  

### set plotting options
options <- function() theme_bw()+theme(strip.text.y = element_text(size =14),
  legend.position =  "top",
  axis.title.x= element_blank(),
  legend.direction="horizontal",
  legend.background=element_blank(),
  legend.key.height= unit(1,"lines"),
  legend.title = element_blank(),
  legend.text = element_text(size = 12),
  legend.title.align = 0,
  legend.key.width = unit(2.5,"lines"),
  legend.key = element_rect(colour=NA),
  panel.background=element_blank(),
  plot.background=element_blank())

### generate a figure
p1 <- ggplot(aes(as.Date(date2),temp), data = subset(SIO,YEAR>1990))+
  scale_x_date(date_labels = "%b")+
  geom_rect(aes(xmin= as.Date(parse_date_time("9-15-2000","%m-%d-%Y")),
                xmax=as.Date(parse_date_time("12-21-2000","%m-%d-%Y")),
                ymin= -Inf,ymax =Inf), data= SIO)+
  geom_point(alpha = 0.1,size = 0.2)+
  stat_smooth(method= "gam", se= F, formula = y~s(x, bs = "cc"),
              alpha = 0.5, aes(colour = "mean"),size = 1)+
  theme_bw() +
  options() +
  theme(legend.text = element_text(size= 10),
        legend.background=element_blank(),
        legend.key = element_blank())+
  scale_colour_manual(values = c("red","blue","black"))+
  stat_smooth(method = "gam", formula = y ~ s(x, k = 10), se = F,
              data= subset(SIO, MONTH%in%c(1:12)&YEAR2%in%c(2014,1997,1998,1999,2011)),aes(colour = ifelse(YEAR2%in%c(1999,1998,2011),"La Nina","El Nino")))+
  ylab(expression(paste("Bottom Temp", ~degree, " C")))+
  geom_hline(yintercept=c(10,13,16,17,18,20), linetype= 3)+
  scale_y_continuous(breaks= c(10,12,14,16,18,20,22),limits= c(9,24),expand= c(0,0))

#png(file= "Figures/temp_summaries_SCRIPPS.png",width = 500,height = 500)
#p1
#dev.off()

temp_sub <- subset(temp_data, MONTH%in%c(1:12)&YEAR>=2009)
temp_sub$site = factor(temp_sub$site, levels = c("San Diego","Santa Barbara","Sonoma","Quadra","Amphitrite Point"), 
                       labels= c("A) Scripps Pier (SD)","B) Arroyo Burro (SB)","C) Bodega Marine Lab (SON)","D) Quadra Island (BC)","E) Ucluelet (BC)"))

monthly_mean <- temp_sub %>%
  group_by(MONTH,YEAR,site) %>%
  dplyr::summarize(temp= mean(temp,na.rm=T))%>%
  mutate(YEAR2 = ifelse(MONTH%in%c(1:3),YEAR-1,YEAR))%>%
  mutate(date2 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,15, sep= "-"), "%Y-%m-%d"),
         date3 = parse_date_time(paste(ifelse(MONTH%in%c(1:3),"2001","2000"),MONTH,1, sep= "-"), "%Y-%m-%d"))

p1 <- ggplot(aes(as.Date(date2),temp), data = subset(temp_sub, YEAR>=2015 & YEAR<2021))+
  scale_x_date(date_labels = "%b",  date_breaks = '3 months', date_minor_breaks= '1 months', expand= c(0,0))+
  geom_point(alpha = 0.1,size=0.1)+
  geom_point(data = monthly_mean,alpha= 0.5,size= 1)+
  theme_bw() +
  options() +
  facet_grid(~site)+
  geom_segment(aes(y = 6, yend = 6.2, x=as.Date(date3)), data= monthly_mean)+
  theme(legend.text = element_text(size= 10),
        legend.background=element_blank(),
        axis.minor.ticks.length = unit(0.1,"cm"),
        panel.background = element_blank(),
        strip.background= element_blank(),
        panel.grid= element_blank())+
  ylab(expression(paste("Bottom Temp", ~degree, "C")))+
  geom_hline(yintercept=c(10,20), linetype= 3)+
  scale_y_continuous(breaks= c(6,8,10,12,14,16,18,20,22),
                     limits= c(6,24),expand= c(0,0));p1

#pdf(file= "Figures/temp_summaries_full.pdf",width = 10,height = 3)
#p1
#dev.off()