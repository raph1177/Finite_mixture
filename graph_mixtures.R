#### Graph for the simulation exercises ####
library(data.table)
library(tidyr)
library(ggplot2)
library(nnet)
library(haven)

#### Directory to save the graphs ####
setwd("C:\\Users\\raph1\\OneDrive\\Bureau\\arxiv\\mixture")

### Set the inferior and superior limits for all graphs ###
quant_sup=0.95
quant_inf=0.05

#### Graph first simulation exercise ####
##### Normal, mu=0.25, sd=1,1, pi1=0.5 #####
norm_mu1_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_025_sd1_1_1.txt")))
norm_mu2_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_025_sd1_1_1.txt")))
norm_pi1_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_025_sd1_1_1.txt")))
norm_pi2_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_025_sd1_1_1.txt")))
norm_sd1_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_025_sd1_1_1.txt")))
norm_sd2_high_025_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_025_sd1_1_1.txt")))
nrow(norm_mu1_high_025_sd1[V1> 2])/1000
nrow(norm_mu1_high_025_sd1[V2> 2])/1000
nrow(norm_mu1_high_025_sd1[V3> 2])/1000
colnames(norm_mu1_high_025_sd1) = c("100","1000","10000")
norm_mu1_high_025_sd1 = pivot_longer(norm_mu1_high_025_sd1, c(1,2,3))
norm_mu1_high_025_sd1$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_025_sd1) = c("N","value","Coefficient")
colnames(norm_mu2_high_025_sd1) = c("100","1000","10000")
norm_mu2_high_025_sd1 = pivot_longer(norm_mu2_high_025_sd1, c(1,2,3))
norm_mu2_high_025_sd1$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_025_sd1) = c("N","value","Coefficient")
colnames(norm_pi1_high_025_sd1) = c("100","1000","10000")
norm_pi1_high_025_sd1 = pivot_longer(norm_pi1_high_025_sd1, c(1,2,3))
norm_pi1_high_025_sd1$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_025_sd1) = c("N","value","Coefficient")
colnames(norm_pi2_high_025_sd1) = c("100","1000","10000")
norm_pi2_high_025_sd1 = pivot_longer(norm_pi2_high_025_sd1, c(1,2,3))
norm_pi2_high_025_sd1$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_025_sd1) = c("N","value","Coefficient")
colnames(norm_sd1_high_025_sd1) = c("100","1000","10000")
norm_sd1_high_025_sd1 = pivot_longer(norm_sd1_high_025_sd1, c(1,2,3))
norm_sd1_high_025_sd1$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_025_sd1) = c("N","value","Coefficient")
colnames(norm_sd2_high_025_sd1) = c("100","1000","10000")
norm_sd2_high_025_sd1 = pivot_longer(norm_sd2_high_025_sd1, c(1,2,3))
norm_sd2_high_025_sd1$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_025_sd1) = c("N","value","Coefficient")

norm_high_025_sd1 = as.data.table(rbind(norm_mu1_high_025_sd1, norm_mu2_high_025_sd1, norm_pi1_high_025_sd1, norm_pi2_high_025_sd1,norm_sd1_high_025_sd1, norm_sd2_high_025_sd1))
norm_high_025_sd1_mean = norm_high_025_sd1[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_025_sd1,norm_mu1_high_025_sd1,norm_pi1_high_025_sd1,norm_pi2_high_025_sd1,norm_sd1_high_025_sd1,norm_sd2_high_025_sd1)

intercept = as.data.frame(unique(norm_high_025_sd1$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0.25,0.25,-0.25,-0.25,0.5,0.5,0.5,0.5,1,1,1,1)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_025_sd1_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig2_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(norm_high_025_sd1[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_025_sd1[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_025_sd1[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_025_sd1[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_025_sd1[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_025_sd1[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(-0.25,0.5,1)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0.25,0.5,1)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig2_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Normal, mu=0.75, sd=1,1, pi1=0.5 #####
norm_mu1_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_075_sd1_1_1.txt")))
norm_mu2_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_075_sd1_1_1.txt")))
norm_pi1_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_075_sd1_1_1.txt")))
norm_pi2_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_075_sd1_1_1.txt")))
norm_sd1_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_075_sd1_1_1.txt")))
norm_sd2_high_075_sd1 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_075_sd1_1_1.txt")))
nrow(norm_mu1_high_075_sd1[V1> 2])/1000
nrow(norm_mu1_high_075_sd1[V2> 2])/1000
nrow(norm_mu1_high_075_sd1[V3> 2])/1000
colnames(norm_mu1_high_075_sd1) = c("100","1000","10000")
norm_mu1_high_075_sd1 = pivot_longer(norm_mu1_high_075_sd1, c(1,2,3))
norm_mu1_high_075_sd1$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_075_sd1) = c("N","value","Coefficient")
colnames(norm_mu2_high_075_sd1) = c("100","1000","10000")
norm_mu2_high_075_sd1 = pivot_longer(norm_mu2_high_075_sd1, c(1,2,3))
norm_mu2_high_075_sd1$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_075_sd1) = c("N","value","Coefficient")
colnames(norm_pi1_high_075_sd1) = c("100","1000","10000")
norm_pi1_high_075_sd1 = pivot_longer(norm_pi1_high_075_sd1, c(1,2,3))
norm_pi1_high_075_sd1$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_075_sd1) = c("N","value","Coefficient")
colnames(norm_pi2_high_075_sd1) = c("100","1000","10000")
norm_pi2_high_075_sd1 = pivot_longer(norm_pi2_high_075_sd1, c(1,2,3))
norm_pi2_high_075_sd1$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_075_sd1) = c("N","value","Coefficient")
colnames(norm_sd1_high_075_sd1) = c("100","1000","10000")
norm_sd1_high_075_sd1 = pivot_longer(norm_sd1_high_075_sd1, c(1,2,3))
norm_sd1_high_075_sd1$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_075_sd1) = c("N","value","Coefficient")
colnames(norm_sd2_high_075_sd1) = c("100","1000","10000")
norm_sd2_high_075_sd1 = pivot_longer(norm_sd2_high_075_sd1, c(1,2,3))
norm_sd2_high_075_sd1$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_075_sd1) = c("N","value","Coefficient")
norm_high_075_sd1 = as.data.table(rbind(norm_mu1_high_075_sd1, norm_mu2_high_075_sd1, norm_pi1_high_075_sd1, norm_pi2_high_075_sd1,norm_sd1_high_075_sd1, norm_sd2_high_075_sd1))
norm_high_075_sd1_mean = norm_high_075_sd1[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_075_sd1,norm_mu1_high_075_sd1,norm_pi1_high_075_sd1,norm_pi2_high_075_sd1,norm_sd1_high_075_sd1,norm_sd2_high_075_sd1)

intercept = as.data.frame(unique(norm_high_075_sd1$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0.75,0.75,-0.75,-0.75,0.5,0.5,0.5,0.5,1,1,1,1)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_075_sd1_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig3_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(norm_high_075_sd1[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_075_sd1[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_075_sd1[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_075_sd1[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_075_sd1[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_075_sd1[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(-0.75,0.5,1)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0.75,0.5,1)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig3_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Normal, mu=0.25, sd=0.95,1.05, pi1=0.7 #####
norm_mu1_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_025_sd095_105_1.txt")))
norm_mu2_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_025_sd095_105_1.txt")))
norm_pi1_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_025_sd095_105_1.txt")))
norm_pi2_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_025_sd095_105_1.txt")))
norm_sd1_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_025_sd095_105_1.txt")))
norm_sd2_high_025_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_025_sd095_105_1.txt")))
nrow(norm_mu1_high_025_sd095[V1> 2])/1000
nrow(norm_mu1_high_025_sd095[V2> 2])/1000
nrow(norm_mu1_high_025_sd095[V3> 2])/1000
colnames(norm_mu1_high_025_sd095) = c("100","1000","10000")
norm_mu1_high_025_sd095 = pivot_longer(norm_mu1_high_025_sd095, c(1,2,3))
norm_mu1_high_025_sd095$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_025_sd095) = c("N","value","Coefficient")
colnames(norm_mu2_high_025_sd095) = c("100","1000","10000")
norm_mu2_high_025_sd095 = pivot_longer(norm_mu2_high_025_sd095, c(1,2,3))
norm_mu2_high_025_sd095$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_025_sd095) = c("N","value","Coefficient")
colnames(norm_pi1_high_025_sd095) = c("100","1000","10000")
norm_pi1_high_025_sd095 = pivot_longer(norm_pi1_high_025_sd095, c(1,2,3))
norm_pi1_high_025_sd095$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_025_sd095) = c("N","value","Coefficient")
colnames(norm_pi2_high_025_sd095) = c("100","1000","10000")
norm_pi2_high_025_sd095 = pivot_longer(norm_pi2_high_025_sd095, c(1,2,3))
norm_pi2_high_025_sd095$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_025_sd095) = c("N","value","Coefficient")
colnames(norm_sd1_high_025_sd095) = c("100","1000","10000")
norm_sd1_high_025_sd095 = pivot_longer(norm_sd1_high_025_sd095, c(1,2,3))
norm_sd1_high_025_sd095$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_025_sd095) = c("N","value","Coefficient")
colnames(norm_sd2_high_025_sd095) = c("100","1000","10000")
norm_sd2_high_025_sd095 = pivot_longer(norm_sd2_high_025_sd095, c(1,2,3))
norm_sd2_high_025_sd095$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_025_sd095) = c("N","value","Coefficient")
norm_high_025_sd095 = as.data.table(rbind(norm_mu1_high_025_sd095, norm_mu2_high_025_sd095, norm_pi1_high_025_sd095, norm_pi2_high_025_sd095,norm_sd1_high_025_sd095, norm_sd2_high_025_sd095))
norm_high_025_sd095_mean = norm_high_025_sd095[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_025_sd095,norm_mu1_high_025_sd095,norm_pi1_high_025_sd095,norm_pi2_high_025_sd095,norm_sd1_high_025_sd095,norm_sd2_high_025_sd095)

intercept = as.data.frame(unique(norm_high_025_sd095$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0.25,0.25,-0.25,-0.25,0.7,0.7,0.3,0.3,0.95,0.95,1.05,1.05)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_025_sd095_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c1_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(norm_high_025_sd095[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_025_sd095[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_025_sd095[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_025_sd095[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_025_sd095[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_025_sd095[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(-0.25,0.3,1.05)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0.25,0.7,0.95)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c1_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Normal, mu=0.75, sd=0.95,1.05, pi1=0.7 #####
norm_mu1_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_075_sd095_105_1.txt")))
norm_mu2_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_075_sd095_105_1.txt")))
norm_pi1_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_075_sd095_105_1.txt")))
norm_pi2_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_075_sd095_105_1.txt")))
norm_sd1_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_075_sd095_105_1.txt")))
norm_sd2_high_075_sd095 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_075_sd095_105_1.txt")))
nrow(norm_mu1_high_075_sd095[V1> 2])/1000
nrow(norm_mu1_high_075_sd095[V2> 2])/1000
nrow(norm_mu1_high_075_sd095[V3> 2])/1000
colnames(norm_mu1_high_075_sd095) = c("100","1000","10000")
norm_mu1_high_075_sd095 = pivot_longer(norm_mu1_high_075_sd095, c(1,2,3))
norm_mu1_high_075_sd095$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_075_sd095) = c("N","value","Coefficient")
colnames(norm_mu2_high_075_sd095) = c("100","1000","10000")
norm_mu2_high_075_sd095 = pivot_longer(norm_mu2_high_075_sd095, c(1,2,3))
norm_mu2_high_075_sd095$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_075_sd095) = c("N","value","Coefficient")
colnames(norm_pi1_high_075_sd095) = c("100","1000","10000")
norm_pi1_high_075_sd095 = pivot_longer(norm_pi1_high_075_sd095, c(1,2,3))
norm_pi1_high_075_sd095$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_075_sd095) = c("N","value","Coefficient")
colnames(norm_pi2_high_075_sd095) = c("100","1000","10000")
norm_pi2_high_075_sd095 = pivot_longer(norm_pi2_high_075_sd095, c(1,2,3))
norm_pi2_high_075_sd095$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_075_sd095) = c("N","value","Coefficient")
colnames(norm_sd1_high_075_sd095) = c("100","1000","10000")
norm_sd1_high_075_sd095 = pivot_longer(norm_sd1_high_075_sd095, c(1,2,3))
norm_sd1_high_075_sd095$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_075_sd095) = c("N","value","Coefficient")
colnames(norm_sd2_high_075_sd095) = c("100","1000","10000")
norm_sd2_high_075_sd095 = pivot_longer(norm_sd2_high_075_sd095, c(1,2,3))
norm_sd2_high_075_sd095$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_075_sd095) = c("N","value","Coefficient")
norm_high_075_sd095 = as.data.table(rbind(norm_mu1_high_075_sd095, norm_mu2_high_075_sd095, norm_pi1_high_075_sd095, norm_pi2_high_075_sd095,norm_sd1_high_075_sd095, norm_sd2_high_075_sd095))
norm_high_075_sd095_mean = norm_high_075_sd095[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_075_sd095,norm_mu1_high_075_sd095,norm_pi1_high_075_sd095,norm_pi2_high_075_sd095,norm_sd1_high_075_sd095,norm_sd2_high_075_sd095)

intercept = as.data.frame(unique(norm_high_075_sd095$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0.75,0.75,-0.75,-0.75,0.7,0.7,0.3,0.3,0.95,0.95,1.05,1.05)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_075_sd095_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c2_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(norm_high_075_sd095[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_075_sd095[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_075_sd095[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_075_sd095[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_075_sd095[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_075_sd095[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(-0.75,0.3,1.05)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0.75,0.7,0.95)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c2_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Normal, mu=0, sd=0.75,1.25, pi1=0.5 #####
norm_mu1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_0_sd_1.txt")))
norm_mu2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_0_sd_1.txt")))
norm_pi1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_0_sd_1.txt")))
norm_pi2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_0_sd_1.txt")))
norm_sd1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_0_sd_1.txt")))
norm_sd2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_0_sd_1.txt")))
nrow(norm_mu1_high_0_sd[V1> 2])/1000
nrow(norm_mu1_high_0_sd[V2> 2])/1000
nrow(norm_mu1_high_0_sd[V3> 2])/1000
colnames(norm_mu1_high_0_sd) = c("100","1000","10000")
norm_mu1_high_0_sd = pivot_longer(norm_mu1_high_0_sd, c(1,2,3))
norm_mu1_high_0_sd$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_mu2_high_0_sd) = c("100","1000","10000")
norm_mu2_high_0_sd = pivot_longer(norm_mu2_high_0_sd, c(1,2,3))
norm_mu2_high_0_sd$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_0_sd) = c("N","value","Coefficient")
colnames(norm_pi1_high_0_sd) = c("100","1000","10000")
norm_pi1_high_0_sd = pivot_longer(norm_pi1_high_0_sd, c(1,2,3))
norm_pi1_high_0_sd$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_pi2_high_0_sd) = c("100","1000","10000")
norm_pi2_high_0_sd = pivot_longer(norm_pi2_high_0_sd, c(1,2,3))
norm_pi2_high_0_sd$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_0_sd) = c("N","value","Coefficient")
colnames(norm_sd1_high_0_sd) = c("100","1000","10000")
norm_sd1_high_0_sd = pivot_longer(norm_sd1_high_0_sd, c(1,2,3))
norm_sd1_high_0_sd$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_sd2_high_0_sd) = c("100","1000","10000")
norm_sd2_high_0_sd = pivot_longer(norm_sd2_high_0_sd, c(1,2,3))
norm_sd2_high_0_sd$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_0_sd) = c("N","value","Coefficient")
norm_high_0_sd = as.data.table(rbind(norm_mu1_high_0_sd, norm_mu2_high_0_sd, norm_pi1_high_0_sd, norm_pi2_high_0_sd,norm_sd1_high_0_sd, norm_sd2_high_0_sd))
norm_high_0_sd_mean = norm_high_0_sd[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_0_sd,norm_mu1_high_0_sd,norm_pi1_high_0_sd,norm_pi2_high_0_sd,norm_sd1_high_0_sd,norm_sd2_high_0_sd)

intercept = as.data.frame(unique(norm_high_075_sd1$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0,0,0,0,0.5,0.5,0.5,0.5,0.75,0.75,1.25,1.25)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_0_sd_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c3_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_0_sd[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_0_sd[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_0_sd[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(0,0.5,1.25)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0,0.5,0.75)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c3_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Normal, mu=0, sd=0.75,1.25, pi1=0.7 (not in the paper) #####
norm_mu1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu1_high_0_sd_1_07.txt")))
norm_mu2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_mu2_high_0_sd_1_07.txt")))
norm_pi1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi1_high_0_sd_1_07.txt")))
norm_pi2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_pi2_high_0_sd_1_07.txt")))
norm_sd1_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd1_high_0_sd_1_07.txt")))
norm_sd2_high_0_sd = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\norm_sd2_high_0_sd_1_07.txt")))
nrow(norm_mu1_high_0_sd[V1> 2])/1000
nrow(norm_mu1_high_0_sd[V2> 2])/1000
nrow(norm_mu1_high_0_sd[V3> 2])/1000
colnames(norm_mu1_high_0_sd) = c("100","1000","10000")
norm_mu1_high_0_sd = pivot_longer(norm_mu1_high_0_sd, c(1,2,3))
norm_mu1_high_0_sd$Coefficient = "Mean, Group 1"
colnames(norm_mu1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_mu2_high_0_sd) = c("100","1000","10000")
norm_mu2_high_0_sd = pivot_longer(norm_mu2_high_0_sd, c(1,2,3))
norm_mu2_high_0_sd$Coefficient = "Mean, Group 2"
colnames(norm_mu2_high_0_sd) = c("N","value","Coefficient")
colnames(norm_pi1_high_0_sd) = c("100","1000","10000")
norm_pi1_high_0_sd = pivot_longer(norm_pi1_high_0_sd, c(1,2,3))
norm_pi1_high_0_sd$Coefficient = "Mixing Weight, Group 1"
colnames(norm_pi1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_pi2_high_0_sd) = c("100","1000","10000")
norm_pi2_high_0_sd = pivot_longer(norm_pi2_high_0_sd, c(1,2,3))
norm_pi2_high_0_sd$Coefficient = "Mixing Weight, Group 2"
colnames(norm_pi2_high_0_sd) = c("N","value","Coefficient")
colnames(norm_sd1_high_0_sd) = c("100","1000","10000")
norm_sd1_high_0_sd = pivot_longer(norm_sd1_high_0_sd, c(1,2,3))
norm_sd1_high_0_sd$Coefficient = "Standard Deviation, Group 1"
colnames(norm_sd1_high_0_sd) = c("N","value","Coefficient")
colnames(norm_sd2_high_0_sd) = c("100","1000","10000")
norm_sd2_high_0_sd = pivot_longer(norm_sd2_high_0_sd, c(1,2,3))
norm_sd2_high_0_sd$Coefficient = "Standard Deviation, Group 2"
colnames(norm_sd2_high_0_sd) = c("N","value","Coefficient")
norm_high_0_sd = as.data.table(rbind(norm_mu1_high_0_sd, norm_mu2_high_0_sd, norm_pi1_high_0_sd, norm_pi2_high_0_sd,norm_sd1_high_0_sd, norm_sd2_high_0_sd))
norm_high_0_sd_mean = norm_high_0_sd[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(norm_mu2_high_0_sd,norm_mu1_high_0_sd,norm_pi1_high_0_sd,norm_pi2_high_0_sd,norm_sd1_high_0_sd,norm_sd2_high_0_sd)

intercept = as.data.frame(unique(norm_high_075_sd1$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(0,0,0,0,0.7,0.7,0.3,0.3,0.75,0.75,1.25,1.25)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(norm_high_0_sd_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=3, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p

data_scatter = as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],norm_high_0_sd[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],norm_high_0_sd[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
data_scatter = rbind(data_scatter,as.data.frame(cbind(norm_high_0_sd[N=="10000"&Coefficient=="Standard Deviation, Group 1",.(V1=value, Coefficient="Standard Deviation")],norm_high_0_sd[N=="10000"&Coefficient=="Standard Deviation, Group 2",.(V2=value)])))

intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(0,0.3,1.25)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(0,0.7,0.75)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p

##### Poisson mean=5,4.5, pi1=0.5 #####
pois_mu1_high_5_45 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_mu1_high_5_45_1.txt")))
pois_mu2_high_5_45 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_mu2_high_5_45_1.txt")))
pois_pi1_high_5_45 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_pi1_high_5_45_1.txt")))
pois_pi2_high_5_45 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_pi2_high_5_45_1.txt")))
colnames(pois_mu1_high_5_45) = c("100","1000","10000")
pois_mu1_high_5_45 = pivot_longer(pois_mu1_high_5_45, c(1,2,3))
pois_mu1_high_5_45$Coefficient = "Mean, Group 1"
colnames(pois_mu1_high_5_45) = c("N","value","Coefficient")
colnames(pois_mu2_high_5_45) = c("100","1000","10000")
pois_mu2_high_5_45 = pivot_longer(pois_mu2_high_5_45, c(1,2,3))
pois_mu2_high_5_45$Coefficient = "Mean, Group 2"
colnames(pois_mu2_high_5_45) = c("N","value","Coefficient")
colnames(pois_pi1_high_5_45) = c("100","1000","10000")
pois_pi1_high_5_45 = pivot_longer(pois_pi1_high_5_45, c(1,2,3))
pois_pi1_high_5_45$Coefficient = "Mixing Weight, Group 1"
colnames(pois_pi1_high_5_45) = c("N","value","Coefficient")
colnames(pois_pi2_high_5_45) = c("100","1000","10000")
pois_pi2_high_5_45 = pivot_longer(pois_pi2_high_5_45, c(1,2,3))
pois_pi2_high_5_45$Coefficient = "Mixing Weight, Group 2"
colnames(pois_pi2_high_5_45) = c("N","value","Coefficient")
pois_high_5_45 = as.data.table(rbind(pois_mu1_high_5_45, pois_mu2_high_5_45, pois_pi1_high_5_45, pois_pi2_high_5_45))
pois_high_5_45_mean = pois_high_5_45[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(pois_mu2_high_5_45,pois_mu1_high_5_45,pois_pi1_high_5_45,pois_pi2_high_5_45)

intercept = as.data.frame(unique(pois_high_5_45$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(5,5,4.5,4.5,0.5,0.5,0.5,0.5)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(pois_high_5_45_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=2, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c4_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(pois_high_5_45[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],pois_high_5_45[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(pois_high_5_45[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],pois_high_5_45[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(4.5,0.5)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(5,0.5)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c4_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Poisson mean=5-3, pi1=0.5 #####
pois_mu1_high_5_3 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_mu1_high_5_3_1.txt")))
pois_mu2_high_5_3 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_mu2_high_5_3_1.txt")))
pois_pi1_high_5_3 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_pi1_high_5_3_1.txt")))
pois_pi2_high_5_3 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\pois_pi2_high_5_3_1.txt")))
colnames(pois_mu1_high_5_3) = c("100","1000","10000")
pois_mu1_high_5_3 = pivot_longer(pois_mu1_high_5_3, c(1,2,3))
pois_mu1_high_5_3$Coefficient = "Mean, Group 1"
colnames(pois_mu1_high_5_3) = c("N","value","Coefficient")
colnames(pois_mu2_high_5_3) = c("100","1000","10000")
pois_mu2_high_5_3 = pivot_longer(pois_mu2_high_5_3, c(1,2,3))
pois_mu2_high_5_3$Coefficient = "Mean, Group 2"
colnames(pois_mu2_high_5_3) = c("N","value","Coefficient")
colnames(pois_pi1_high_5_3) = c("100","1000","10000")
pois_pi1_high_5_3 = pivot_longer(pois_pi1_high_5_3, c(1,2,3))
pois_pi1_high_5_3$Coefficient = "Mixing Weight, Group 1"
colnames(pois_pi1_high_5_3) = c("N","value","Coefficient")
colnames(pois_pi2_high_5_3) = c("100","1000","10000")
pois_pi2_high_5_3 = pivot_longer(pois_pi2_high_5_3, c(1,2,3))
pois_pi2_high_5_3$Coefficient = "Mixing Weight, Group 2"
colnames(pois_pi2_high_5_3) = c("N","value","Coefficient")
pois_high_5_3 = as.data.table(rbind(pois_mu1_high_5_3, pois_mu2_high_5_3, pois_pi1_high_5_3, pois_pi2_high_5_3))
pois_high_5_3_mean = pois_high_5_3[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(pois_mu2_high_5_3,pois_mu1_high_5_3,pois_pi1_high_5_3,pois_pi2_high_5_3)

intercept = as.data.frame(unique(pois_high_5_3$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(5,5,3,3,0.5,0.5,0.5,0.5)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(pois_high_5_3_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=2, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c5_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(pois_high_5_3[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],pois_high_5_3[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(pois_high_5_3[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],pois_high_5_3[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(3,0.5)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(5,0.5)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c5_2.pdf",plot=p, width=8, height=2.5, dpi=300)

##### Exponential, mean=1,0.9, pi1=0.5 #####
expo_mu1_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu1_high_1_09_1.txt")))
expo_mu2_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu2_high_1_09_1.txt")))
expo_pi1_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi1_high_1_09_1.txt")))
expo_pi2_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi2_high_1_09_1.txt")))
colnames(expo_mu1_high_1_09) = c("100","1000","10000")
expo_mu1_high_1_09 = pivot_longer(expo_mu1_high_1_09, c(1,2,3))
expo_mu1_high_1_09$Coefficient = "Mean, Group 1"
colnames(expo_mu1_high_1_09) = c("N","value","Coefficient")
colnames(expo_mu2_high_1_09) = c("100","1000","10000")
expo_mu2_high_1_09 = pivot_longer(expo_mu2_high_1_09, c(1,2,3))
expo_mu2_high_1_09$Coefficient = "Mean, Group 2"
colnames(expo_mu2_high_1_09) = c("N","value","Coefficient")
colnames(expo_pi1_high_1_09) = c("100","1000","10000")
expo_pi1_high_1_09 = pivot_longer(expo_pi1_high_1_09, c(1,2,3))
expo_pi1_high_1_09$Coefficient = "Mixing Weight, Group 1"
colnames(expo_pi1_high_1_09) = c("N","value","Coefficient")
colnames(expo_pi2_high_1_09) = c("100","1000","10000")
expo_pi2_high_1_09 = pivot_longer(expo_pi2_high_1_09, c(1,2,3))
expo_pi2_high_1_09$Coefficient = "Mixing Weight, Group 2"
colnames(expo_pi2_high_1_09) = c("N","value","Coefficient")
expo_high_1_09 = as.data.table(rbind(expo_mu1_high_1_09, expo_mu2_high_1_09, expo_pi1_high_1_09, expo_pi2_high_1_09))
expo_high_1_09_mean = expo_high_1_09[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(expo_mu2_high_1_09,expo_mu1_high_1_09,expo_pi1_high_1_09,expo_pi2_high_1_09)

intercept = as.data.frame(unique(expo_high_1_09$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(1,1,0.9,0.9,0.5,0.5,0.5,0.5)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(expo_high_1_09_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=2, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c6_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(expo_high_1_09[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],expo_high_1_09[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(expo_high_1_09[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],expo_high_1_09[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(0.9,0.5)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(1,0.5)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c6_2.pdf",plot=p, width=8, height=2.5, dpi=300)


##### Exponential, mean=1,0.9, pi1=0.7 (not in the paper) #####
expo_mu1_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu1_high_1_09_73.txt")))
expo_mu2_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu2_high_1_09_73.txt")))
expo_pi1_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi1_high_1_09_73.txt")))
expo_pi2_high_1_09 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi2_high_1_09_73.txt")))
colnames(expo_mu1_high_1_09) = c("100","1000","10000")
expo_mu1_high_1_09 = pivot_longer(expo_mu1_high_1_09, c(1,2,3))
expo_mu1_high_1_09$Coefficient = "Mean, Group 1"
colnames(expo_mu1_high_1_09) = c("N","value","Coefficient")
colnames(expo_mu2_high_1_09) = c("100","1000","10000")
expo_mu2_high_1_09 = pivot_longer(expo_mu2_high_1_09, c(1,2,3))
expo_mu2_high_1_09$Coefficient = "Mean, Group 2"
colnames(expo_mu2_high_1_09) = c("N","value","Coefficient")
colnames(expo_pi1_high_1_09) = c("100","1000","10000")
expo_pi1_high_1_09 = pivot_longer(expo_pi1_high_1_09, c(1,2,3))
expo_pi1_high_1_09$Coefficient = "Mixing Weight, Group 1"
colnames(expo_pi1_high_1_09) = c("N","value","Coefficient")
colnames(expo_pi2_high_1_09) = c("100","1000","10000")
expo_pi2_high_1_09 = pivot_longer(expo_pi2_high_1_09, c(1,2,3))
expo_pi2_high_1_09$Coefficient = "Mixing Weight, Group 2"
colnames(expo_pi2_high_1_09) = c("N","value","Coefficient")
expo_high_1_09 = as.data.table(rbind(expo_mu1_high_1_09, expo_mu2_high_1_09, expo_pi1_high_1_09, expo_pi2_high_1_09))
expo_high_1_09_mean = expo_high_1_09[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(expo_mu2_high_1_09,expo_mu1_high_1_09,expo_pi1_high_1_09,expo_pi2_high_1_09)

intercept = as.data.frame(unique(expo_high_1_09$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(1,1,0.9,0.9,0.7,0.7,0.3,0.3)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(expo_high_1_09_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=2, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p

data_scatter = as.data.frame(cbind(expo_high_1_09[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],expo_high_1_09[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(expo_high_1_09[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],expo_high_1_09[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(0.9,0.3)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(1,0.7)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p


##### Exponential, mean=1,0.6, pi1=0.5 #####
expo_mu1_high_1_05 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu1_high_1_06_1.txt")))
expo_mu2_high_1_05 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_mu2_high_1_06_1.txt")))
expo_pi1_high_1_05 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi1_high_1_06_1.txt")))
expo_pi2_high_1_05 = as.data.table(t(fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\expo_pi2_high_1_06_1.txt")))
colnames(expo_mu1_high_1_05) = c("100","1000","10000")
expo_mu1_high_1_05 = pivot_longer(expo_mu1_high_1_05, c(1,2,3))
expo_mu1_high_1_05$Coefficient = "Mean, Group 1"
colnames(expo_mu1_high_1_05) = c("N","value","Coefficient")
colnames(expo_mu2_high_1_05) = c("100","1000","10000")
expo_mu2_high_1_05 = pivot_longer(expo_mu2_high_1_05, c(1,2,3))
expo_mu2_high_1_05$Coefficient = "Mean, Group 2"
colnames(expo_mu2_high_1_05) = c("N","value","Coefficient")
colnames(expo_pi1_high_1_05) = c("100","1000","10000")
expo_pi1_high_1_05 = pivot_longer(expo_pi1_high_1_05, c(1,2,3))
expo_pi1_high_1_05$Coefficient = "Mixing Weight, Group 1"
colnames(expo_pi1_high_1_05) = c("N","value","Coefficient")
colnames(expo_pi2_high_1_05) = c("100","1000","10000")
expo_pi2_high_1_05 = pivot_longer(expo_pi2_high_1_05, c(1,2,3))
expo_pi2_high_1_05$Coefficient = "Mixing Weight, Group 2"
colnames(expo_pi2_high_1_05) = c("N","value","Coefficient")
expo_high_1_05 = as.data.table(rbind(expo_mu1_high_1_05, expo_mu2_high_1_05, expo_pi1_high_1_05, expo_pi2_high_1_05))
expo_high_1_05_mean = expo_high_1_05[,.(Mean = mean(value), lim_inf = quantile(value, probs=quant_inf), lim_sup = quantile(value, probs=quant_sup)), by=c("Coefficient","N")]
rm(expo_mu2_high_1_05,expo_mu1_high_1_05,expo_pi1_high_1_05,expo_pi2_high_1_05)

intercept = as.data.frame(unique(expo_high_1_05$Coefficient))
intercept = rbind(intercept,intercept)
setorder(intercept)
intercept$int = c(1,1,0.6,0.6,0.5,0.5,0.5,0.5)
colnames(intercept) = c("Coefficient","Int")

p = ggplot(expo_high_1_05_mean, aes(x=N, y=Mean, group=1))+
  geom_point(colour="steelblue", alpha=1.1, size = 2)+
  facet_wrap(vars(Coefficient), nrow=2, scales = "free")+
  geom_hline(data=intercept, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup),colour="steelblue",width=.2,position=position_dodge(0.05))+
  xlab("Number of Observations per Replication")+ylab("Estimated Value")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c7_1.pdf",plot=p, width=8, height=5.5, dpi=300)

data_scatter = as.data.frame(cbind(expo_high_1_05[N=="10000"&Coefficient=="Mean, Group 1",.(V1=value, Coefficient="Mean Value")],expo_high_1_05[N=="10000"&Coefficient=="Mean, Group 2",.(V2=value)]))
data_scatter = rbind(data_scatter,as.data.frame(cbind(expo_high_1_05[N=="10000"&Coefficient=="Mixing Weight, Group 1",.(V1=value, Coefficient="Mixing Weight")],expo_high_1_05[N=="10000"&Coefficient=="Mixing Weight, Group 2",.(V2=value)])))
intercept_h = as.data.frame(unique(data_scatter$Coefficient))
intercept_h$int = c(0.6,0.5)
colnames(intercept_h) = c("Coefficient","Int")
intercept_v = as.data.frame(unique(data_scatter$Coefficient))
intercept_v$int = c(1,0.5)
colnames(intercept_v) = c("Coefficient","Int")

p = ggplot(data_scatter, aes(x=V1, y=V2))+
  geom_point(colour="black", alpha=0.7, size = 0.8)+
  xlab("Estimated Value, Group 1")+ylab("Estimated Value, Group 2")+
  geom_hline(data=intercept_h, aes(yintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  geom_vline(data=intercept_v, aes(xintercept = Int), linetype="dashed", colour="red", alpha=0.8)+
  facet_wrap(vars(Coefficient), nrow=1, scales = "free")+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"))
p
ggsave("fig_c7_2.pdf",plot=p, width=8, height=2.5, dpi=300)



#### Graph second simulation exercise ####
##### K=2, 1 covariate, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov1_K2_500.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov1_K2_500.txt")
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm,),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig4_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov1_K2_500.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov1_K2_500.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov1_K2_500.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov1_K2_500.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov1_K2_500.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov1_K2_500.txt")

data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig4_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig4_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=2, 5 covariates, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov5_K2_500.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov5_K2_500.txt") 
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig5_1.pdf",plot=p, width=8, height=4.5, dpi=300)

data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov5_K2_500.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov5_K2_500.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov5_K2_500.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov5_K2_500.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov5_K2_500.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov5_K2_500.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig5_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig5_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=2, 10 covariates, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov10_K2_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov10_K2_500_1.txt") 
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig6_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov10_K2_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov10_K2_500_1.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov10_K2_500_1.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov10_K2_500_1.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov10_K2_500_1.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov10_K2_500_1.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig6_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig6_3.pdf",plot=p, width=8, height=3, dpi=300)


##### K=2, 1 covariate, N=750, T=8 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov1_K2_750.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov1_K2_750.txt") 
data1 = data_EM[,1:10]
data2 = data_EM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:10]
data2 = data_CEM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c8_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov1_K2_750.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov1_K2_750.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov1_K2_750.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov1_K2_750.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov1_K2_750.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov1_K2_750.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c8_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c8_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=2, 5 covariates, N=750, T=8 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov5_K2_750.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov5_K2_750.txt") 
data1 = data_EM[,1:10]
data2 = data_EM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:10]
data2 = data_CEM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c9_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov5_K2_750.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov5_K2_750.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov5_K2_750.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov5_K2_750.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov5_K2_750.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov5_K2_750.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c9_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c9_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=2, 10 covariates, N=750, T=8 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov10_K2_750_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov10_K2_750_1.txt") 
data1 = data_EM[,1:10]
data2 = data_EM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_EM = rbind(data1, data2)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:10]
data2 = data_CEM[,11:20]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5","Time6","Time7","Time8") 
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Group = "Group 2"
data_CEM = rbind(data1, data2)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c10_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov10_K2_750_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov10_K2_750_1.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov10_K2_750_1.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov10_K2_750_1.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov10_K2_750_1.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov10_K2_750_1.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Mixing Weight, Group 1","Mixing Weight, Group 2","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c10_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c10_3.pdf",plot=p, width=8, height=3, dpi=300)



##### K=3, 1 covariates, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov1_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov1_K3_500_1.txt") 
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
data3 = data_EM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_EM = rbind(data1, data2,data3)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
data3 = data_CEM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_CEM = rbind(data1, data2, data3)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]


p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c11_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov1_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov1_K3_500_1.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov1_K3_500_1.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov1_K3_500_1.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov1_K3_500_1.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov1_K3_500_1.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Mixing Weight, Group3"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c11_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c11_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=3, 5 covariates, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov5_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov5_K3_500_1.txt") 
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
data3 = data_EM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_EM = rbind(data1, data2,data3)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
data3 = data_CEM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_CEM = rbind(data1, data2, data3)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]


p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c12_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov5_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov5_K3_500_1.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov5_K3_500_1.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov5_K3_500_1.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov5_K3_500_1.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov5_K3_500_1.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")

data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Mixing Weight, Group3"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c12_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c12_3.pdf",plot=p, width=8, height=3, dpi=300)

##### K=3, 10 covariates, N=500, T=5 #####
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_EM_cov10_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_beta_CEM_cov10_K3_500_1.txt") 
data1 = data_EM[,1:7]
data2 = data_EM[,8:14]
data3 = data_EM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_EM = rbind(data1, data2,data3)
data_EM$Algorith = "EM"

data1 = data_CEM[,1:7]
data2 = data_CEM[,8:14]
data3 = data_CEM[,15:21]
colnames(data1) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data2) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
colnames(data3) = c("Beta","Gamma","Time1","Time2","Time3","Time4","Time5") 
data1 = pivot_longer(data1, cols= seq(1,7,1))
data1$Group = "Group 1"
data2 = pivot_longer(data2, cols= seq(1,7,1))
data2$Group = "Group 2"
data3 = pivot_longer(data3, cols= seq(1,7,1))
data3$Group = "Group 3"
data_CEM = rbind(data1, data2, data3)
data_CEM$Algorith = "C-EM"
data = rbind(data_EM,data_CEM)
colnames(data) = c("Var","Value","Group","Algorithm")
data = as.data.table(data)
data_mean = data[,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Group","Algorithm","Var")]

p = ggplot(data_mean, aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  facet_grid(vars(Group))+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c13_1.pdf",plot=p, width=8, height=4.5, dpi=300)
data_EM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_EM_cov10_K3_500_1.txt")
data_CEM = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_var_CEM_cov10_K3_500_1.txt")
data_EM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_EM_cov10_K3_500_1.txt")
data_CEM_pi = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\bias_pi_CEM_cov10_K3_500_1.txt")
data_EM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_EM_cov10_K3_500_1.txt")
data_CEM_class = fread(file="C:\\Users\\raph1\\.spyder-py3\\Work\\class_error_CEM_cov10_K3_500_1.txt")
data1 = cbind(data_EM, data_EM_pi, data_EM_class)
data2 = cbind(data_CEM, data_CEM_pi, data_CEM_class)
colnames(data1) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")
colnames(data2) = c("Sigma2_e, Group 1","Sigma2_a, Group 1","Sigma2_e, Group 2","Sigma2_a, Group 2","Sigma2_e, Group3","Sigma2_a, Group3","Mixing Weight, Group 1","Mixing Weight, Group 2","Mixing Weight, Group3","Misclassification Rate")
test = as.data.table(rowSums(data_CEM_pi))
data1 = pivot_longer(data1, cols= seq(1,10,1))
data1$Algorithm = "EM"
data2 = pivot_longer(data2, cols= seq(1,10,1))
data2$Algorithm = "C-EM"
data = rbind(data1,data2)
colnames(data) = c("Var","Value","Algorithm")
data = as.data.table(data)
data_mean = data[Var=="Mixing Weight, Group 1"|Var=="Mixing Weight, Group 2"|Var=="Mixing Weight, Group3"|Var=="Misclassification Rate",Group:=2][is.na(Group),Group:=1][,.(Mean = mean(Value), lim_inf = quantile(Value, probs=c(quant_inf)),lim_sup = quantile(Value, probs=c(quant_sup))), by=c("Algorithm","Var","Group")]

p = ggplot(data_mean[Group==1], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), axis.title.x=element_blank(), legend.position="none")+
  ylab("Estimation Bias")
p
ggsave("fig_c13_2.pdf",plot=p, width=8, height=2.5, dpi=300)
p = ggplot(data_mean[Group==2], aes(y=Mean, x=Var, fill=Algorithm))+
  geom_point(aes(colour=Algorithm),size=2,position=position_dodge(0.4))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", alpha=0.6)+
  geom_errorbar(aes(ymin=lim_inf, ymax=lim_sup, colour=Algorithm),width=.2,position=position_dodge(0.4))+
  theme_bw()+
  theme(strip.background =element_rect(fill="lightgrey"), legend.position="bottom")+
  ylab("Estimation Bias")+
  xlab("Coefficient")
p
ggsave("fig_c13_3.pdf",plot=p, width=8, height=3, dpi=300)
