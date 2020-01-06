#Utkarsh Srivastava
#HW 2
#EBGN 590
install.packages(c("readxl","stargazer","ggplot2","corrgram"))
install.packages("tidyverse")
install.packages("plotly")
install.packages("tseries")
install.packages("readr")
library(readxl)
library(tidyverse)
library(stargazer)
library(ggplot2)
library(corrgram)
library(plotly)
library(gghighlight)
library(readxl)
library(stargazer)
library(tseries)
library(readr)

data <- read.csv(file="D:\\Academics\\CSM\\Spring 2019\\Econometrics\\Homeworks\\HW 3\\oil_gaso_data.csv",
         header=TRUE)
data$DATE <- as.Date(data$DATE, "%m/%d/%Y")
data_real <- mutate(data,
                    real_brent = Brent*(CPI2015/100),
                    real_gasoline = Gasoline*(CPI2015/100),
                    real_gasoline_barrel = real_gasoline*42
)

p <- ggplot(data=data_real) + 
  geom_line(aes(x=DATE,y=real_brent, colour="Brent Oil")) +
  geom_line(aes(x=DATE,y=real_gasoline*25, colour="Gasoline"))+ 
  theme_bw() + 
  scale_colour_manual(values = c("green", "red")) + 
  labs(y = "2015 $",
       x = "Date",
       colour = "Type of oil") + 
  theme(legend.position = c(0.2, 0.8))+
  scale_x_date(date_breaks = "2 year", date_labels = "%b-%Y")
p
summary(data_real$real_brent)
summary(data_real$real_gasoline_barrel)

adf.test(data_real$real_brent, k=1)
adf.test(data_real$real_gasoline, k=1)

stats::PP.test(data_real$real_brent)
stats::PP.test(data_real$real_gasoline)

data_real$diff_brent <- (lead(data_real$real_brent)-data_real$real_brent)/data_real$real_brent
data_real$diff_gasoline <- (lead(data_real$real_gasoline)-data_real$real_gasoline)/data_real$real_gasoline

p.difference <- ggplot(data=data_real) + 
  geom_line(aes(x = DATE, y =diff_brent,colour="Brent Oil"),size=1) +
  geom_line(aes(x = DATE, y =diff_gasoline,colour="Gasoline"),size=1) +
  scale_colour_manual(values = c("green", "red")) + 
  labs(y = "Difference",
       x = "Date",
       colour = "Oil type") +
  theme_bw()
p.difference

p.density <- ggplot(data=data_real) + 
  geom_density(aes(x = diff_brent, colour="Brent Oil"),size=1) +
  geom_density(aes(x = diff_gasoline, colour="Gasoline"),size=1) +
  scale_colour_manual(values = c("green", "red")) + 
  labs(y = "Density",
       x = "",
       colour = "Oil type") +
  theme_bw()
p.density


data_real <- data_real %>%
  mutate(moy = as.Date(DATE, format = "%Y-%m-%d"))

model1 <- lm(real_brent ~ real_gasoline + as.factor(moy), data = data_real)
model2 <- lm(real_brent ~ real_gasoline, data = data_real)
model1
model2
stargazer(model1, model2, type="text")

data_real <- mutate(data_real,
                    crack_spread = (real_brent/42) - real_gasoline
                    )

p.crack <- ggplot(data=data_real) + 
  geom_line(aes(x=DATE,y=crack_spread)) +
  theme_bw() + 
  labs(y = "2015 crack spread $",
       x = "Date")+
  scale_x_date(date_breaks = "2 year", date_labels = "%b-%Y")
p.crack

data_real$crack_diff <- (lead(data_real$crack_spread)-data_real$crack_spread)/data_real$crack_spread

p.difference2 <- ggplot(data=data_real) + 
  geom_line(aes(x = DATE, y =diff_brent, colour="Brent Oil"),size=1) +
  geom_line(aes(x = DATE, y =crack_diff, color ="Crack spread"),size=1) +
  scale_colour_manual(values = c("green", "red")) +
  labs(y = "Difference",
       x = "Date",
       colour = "Type")+
  theme_bw()
p.difference2

adf.test(data_real$crack_spread, k=1)
stats::PP.test(data_real$crack_spread)