#############################################
### CWT models for the TransPlant network ###
#############################################

# J. Lynn 10/2025

# packages 
library(tidyverse); library(rjags); library(R2jags); library(broom.mixed); library(ggmcmc)
library(loo); library(reshape2); library(gridExtra)
options(max.print = 3000)

# load in data 
dat <- read_csv("MacroContextTransplantDat_cwtdat.csv")

# standardization function
standard <- function(x) (x - mean(x, na.rm=T)) / ((sd(x, na.rm=T)))

# visualize the fig
ggplot(dat, aes(y = (cwt_mat-mat_site), x = exp_time, colour = mat_site, dat$destSiteID, sep = "_") ) + geom_point() +
  facet_grid(~Treatment) + geom_smooth(method = "lm")+ theme(legend.position = "none")

### begin the model of how CWT changes through time in warm versus control plots

# dep var preparation
reg <- factor(dat$Region)
dat$exp_time_ln <- log(dat$exp_time+1)
exp_dur <- dat$exp_time_ln
regsit <- factor(paste(dat$Region, dat$destSiteID, sep = "_"))
N_regsit <- as.numeric(length(levels(regsit)))
warm <- factor(dat$Treatment)

# make the dependent variables
dep_var <- dat$cwt_mat - dat$mat_site
N <- as.numeric(length(dep_var))

# specify the model - can use for each analysis
mc_mod <- function(){
  # likelihood
  for(i in 1:N){
    dep_var[i]~dnorm(mu[i], prec)
    mu[i] <- beta0[regsit[i]] + beta1[regsit[i]]*warm[i] + beta 
    
    res[i] <- pow(dep_var[i]-mu[i], 2)
    dp.new[i] ~ dnorm(mu[i], prec)
    res.new[i] <- pow(dp.new[i]-mu[i],2)
  }
  #Priors
  for(j in 1:N_parm){ab[j]~dnorm(0,1.0E-6)}
  prec~dgamma(0.001,0.001)
  sig <- 1/sqrt(prec)
  
 # prec_reg~dgamma(0.001,0.001)
  #sig_reg <- 1/sqrt(prec_reg)
  
  # derived params
  rss <- sum(res[])
  rss.new <- sum(res.new[])
  
  # gelman R^2
  var_fit <- pow(sd(mu[]), 2)
  var_res <- pow(sd(res[]),2)
  r2 <- var_fit/(var_fit + var_res)
}

# assemble the data amd params for the model 
jags.data <- list( "dep_var", "treatmat","reg", "N", "N_parm", "N_reg")
jags.param <- c("ab", "sig", "sig_reg", "rss", "rss.new", "zlog.lik", "r2")

# run the model for mat differences 
cwmatmod <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000,model.file=mc_mod, n.thin=5, n.chains=3)

cwmatmod

# traceplots
s <- ggs(as.mcmc(cwmatmod))
ggs_traceplot(s, family="15") 

# posterior check
matdif.paramlist <- cwmatmod$BUGSoutput$sims.list

plot(matdif.paramlist$rss,matdif.paramlist$rss.new,main="Mat dif Cover",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(matdif.paramlist$rss>matdif.paramlist$rss.new) 

# assemble posteriors
matdif.post <- as_tibble(matdif.paramlist$ab)
names(matdif.post) <- colnames(treatmat)
boxplot(matdif.post)

write_csv(matdif.post, "CWMAT_mainPost_bysite_log.csv")

# need to obtain the individual intercepts and slopes
#int <- matdif.post$"(Intercept)"
##intNoWarm <- matdif.post %>% select(regsitCH_Calanda_Nes:regsitUS_Montana_Low)
#intNoWarm <- intNoWarm + int

# summarize the posterior
matdif.post.sum <- matdif.post %>%
  pivot_longer(cols = "regsitCH_Calanda_Cal":"regsitUS_Montana_Low:TreatmentWarm:exp_time_ln", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(matdif.post.sum, "CWMAT_mainPost_summarized_bysite_log.csv")

# r^2
mcmc <- cwmatmod$BUGSoutput$sims.matrix
coefs = mcmc[,1:190]
fit = coefs %*% t(treatmat)
resid = sweep(fit, 2, dep_var, "-")
var_f = apply(fit, 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
cwmat.newdata <- select(dat, Treatment, exp_time_ln, Region, mat) %>% 
  cbind(broom.mixed::tidy(as.mcmc(fit), conf.int = TRUE), dep_var, regsit) %>%
  mutate(Treatment = recode(Treatment, LocalControl = "control", Warm = "warm"))
#write_csv(cwmat.newdata, "cwmatDatForPlot_log.csv")
cwmat.newdata <- read_csv("cwmatDatForPlot_log.csv")

dat2 <- dat %>% cbind(dep_var) %>% rename(estimate = dep_var) %>%  
  mutate(Treatment = recode(Treatment, LocalControl = "control", Warm = "warm"))
dat2$regsit <- factor(paste(dat$Region, dat$destSiteID, sep = "_"))
cwmat.newdata <- left_join(cwmat.newdata, unique(select(dat2, regsit, mat_site)), by = "regsit")

matdif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=exp(exp_time_ln)-1, y=estimate, colour=factor(regsit)),
              width = 0.4, alpha=.35, size=0.5)+
  geom_hline(yintercept = 0)+
  geom_ribbon(data=cwmat.newdata, aes(ymin=conf.low, ymax=conf.high, x=exp(exp_time_ln)-1, 
                                   fill=factor(regsit)), alpha=0.15)+
  geom_line(data=cwmat.newdata, aes(y = estimate, x = exp(exp_time_ln)-1, colour=factor(regsit)))+
  facet_wrap(Region~Treatment, scales = "free_y")+
  #scale_colour_manual("Treatment", values=c("gray", "red")) + 
  #scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_continuous("Years") + 
  scale_y_continuous("CWMAT - Site MAT")+ labs(title = "a)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "none",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
matdif_plot

ggsave("CWMATMain_plot_bysite_wkey_log.jpeg", plot=matdif_plot, width = 20, height=12,
       units="cm", dpi=300)


matdif_plot<- ggplot()+ 
  geom_hline(yintercept = 0, color = "gray")+
  geom_jitter(data=dat2, aes(x=(exp(exp_time_ln)-1), y=estimate, shape=factor(regsit), colour = mat_site),
              width = 0.3, alpha=.35, size=0.5, shape = 19)+
  geom_ribbon(data=cwmat.newdata, aes(ymin=conf.low, ymax=conf.high, x=(exp(exp_time_ln)-1), shape = factor(regsit), fill = mat_site), alpha=0.15)+
  geom_line(data=cwmat.newdata, aes(y = estimate, x = exp(exp_time_ln)-1, shape=factor(regsit), colour = mat_site))+
  facet_wrap(~Treatment)+ scale_colour_gradient("Site MAT (ºC)", low = "cyan", high = "red") + 
  scale_fill_gradient("Site MAT (ºC)",  low = "cyan", high = "red") +  scale_x_continuous("Years") + 
  scale_y_continuous("CWMAT - Site MAT")+ labs(title = "a)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = c(.25,.16),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  #guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
matdif_plot

ggsave("CWMATMain_plot_bysite_log_siteMatCol.jpeg", plot=matdif_plot, width = 12, height=6,
       units="cm", dpi=300)

head(treatmat)

## Run a model where parms are estimated individually 
# set up the model matrix for warm data
datwarm <- dat %>% filter(Treatment == "Warm")
datwarm$regsit <- factor(paste(datwarm$Region, datwarm$destSiteID, sep = "_"))

# set up the mati
regsit <- datwarm$regsit
N_regsit <- as.numeric(length(levels(regsit)))
exp_dur <- datwarm$exp_time_ln

# get temperature for the submodel
datsubmodelwarm <-datwarm %>% select(regsit, mat_site) %>% unique()
temp <- datsubmodelwarm$mat_site

# make the dependent variables
dep_var <- datwarm$cwt_mat - datwarm$mat_site
N <- as.numeric(length(dep_var))

# specify the model - can use for each analysis
mc_mod <- function(){
  # likelihood
  for(i in 1:N){
    dep_var[i]~dnorm(mu[i], prec)
    mu[i] <- beta0[regsit[i]] + beta1[regsit[i]] * exp_dur[i]
    
    res[i] <- pow(dep_var[i]-mu[i], 2)
    dp.new[i] ~ dnorm(mu[i], prec)
    res.new[i] <- pow(dp.new[i]-mu[i],2)
  }
  #Priors
  for(j in 1:N_regsit){beta0[j]~dnorm(int_mu[j],prec_int)
    int_mu[j] <- beta0_int + beta1_int*temp[j]}
  beta0_int ~ dnorm(0,1.0E-6)
  beta1_int ~ dnorm(0,1.0E-6)
  prec_int~dgamma(0.001,0.001)
  sig_int<- 1/sqrt(prec_int)
  
  for(j in 1:N_regsit){beta1[j]~dnorm(slope_mu[j],prec_slope)
    slope_mu[j] <- beta0_slope + beta1_slope*temp[j]}
  beta0_slope ~ dnorm(0,1.0E-6)
  beta1_slope ~ dnorm(0,1.0E-6)
  prec_slope~dgamma(0.001,0.001)
  sig_slope<- 1/sqrt(prec_slope)
  
  prec~dgamma(0.001,0.001)
  sig <- 1/sqrt(prec)
  
  # derived params
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# assemble the data amd params for the model 
jags.data <- list( "dep_var", "regsit", "N", "N_regsit", "temp", "exp_dur")
jags.param <- c("beta0", "beta1","beta0_int", "beta1_int","beta0_slope", "beta1_slope", "sig", "sig_slope", "sig_int", "rss", "rss.new", "zlog.lik")

# run the model for mat differences 
cwmatmod_indwarm <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                                  n.iter=50000,model.file=mc_mod, n.thin=5, n.chains=3)

cwmatmod_indwarm

# posterior check
matdif.paramlistwarm <- cwmatmod_indwarm$BUGSoutput$sims.list

mean(matdif.paramlistwarm$beta1_slope>0)
mean(matdif.paramlistwarm$beta1_int>0)

# assemble posteriors
matdif.postindwarmint <- as_tibble(matdif.paramlistwarm$beta0)
base::names(matdif.postindwarmint) <- c(levels(regsit))
boxplot(matdif.postindwarmint)

matdif.postindwarmslope <- as_tibble(matdif.paramlistwarm$beta1)
base::names(matdif.postindwarmslope) <- c(levels(regsit))
boxplot(matdif.postindwarmslope)

#write_csv(matdif.post, "CWMAT_mainPost_bysite.csv")

# summarize the posterior
matdif.post.sumind_intwarm <- matdif.postindwarmint %>% 
  select("CH_Calanda_Cal":"US_Montana_Low")%>%
  pivot_longer(cols = "CH_Calanda_Cal":"US_Montana_Low", names_to = c("Region"), 
               values_to = "post") %>%
  mutate(Treatment = rep("warm", length(post))) %>%
  group_by(Region) %>% dplyr::summarize(mean_int = mean(post),
                                        sd_int = sd(post)) 
matdif.post.sumind_intwarm$Treatment <- rep("warm", length(matdif.post.sumind_intwarm$mean_int))

matdif.post.sumind_slopewarm <-  matdif.postindwarmslope %>% 
  select("CH_Calanda_Cal":"US_Montana_Low") %>%
  pivot_longer(cols = "CH_Calanda_Cal":"US_Montana_Low", 
               names_to = c("Region"),  values_to = "post")  %>%
  group_by(Region) %>% dplyr::summarize(mean_slope = mean(post), 
                                               sd_slope = sd(post)) 
matdif.post.sumind_slopewarm$Treatment <- rep("warm", length(matdif.post.sumind_slopewarm$mean_slope))

# put together the two datasets
slope_int_datawarm <- left_join(matdif.post.sumind_slopewarm, 
                                matdif.post.sumind_intwarm, by = c("Region", "Treatment")) %>% 
  left_join(., rename(datsubmodelwarm, Region = "regsit"), by = "Region")

# make data for 
mcmc <- cwmatmod_indwarm$BUGSoutput$sims.matrix
coefs = as.matrix(select(as_tibble(mcmc), c("beta0_int", "beta1_int")))
mat_graph <- seq(min(slope_int_datawarm$mat_site), max(slope_int_datawarm$mat_site), by =0.1)
fit_int = coefs %*% t(model.matrix(~ mat_graph))

coefs = as.matrix(select(as_tibble(mcmc), c("beta0_slope", "beta1_slope")))
fit_slope = coefs %*%t(model.matrix(~ mat_graph))

# make a graph
mat_line_warmint <- cbind(broom.mixed::tidy(as.mcmc(fit_int), conf.int = TRUE), mat_graph)
mat_line_warmslope <- cbind(broom.mixed::tidy(as.mcmc(fit_slope), conf.int = TRUE), mat_graph)
mat_line_warm <- bind_rows(list(int = mat_line_warmint, slope = mat_line_warmslope), .id = "parm")

#intmodel
ggplot() + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = mat_line_warmint, aes(x = mat_graph, y = estimate))+
  geom_ribbon(data=mat_line_warmint, aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  geom_point(data = slope_int_datawarm, aes(x= mat_site, y = mean_int, 
                                            colour = Region)) +
  geom_linerange(data = slope_int_datawarm, aes(x = mat_site, ymax=mean_int + sd_int,
                                    ymin = mean_int - sd_int, colour = Region)) +
  #facet_wrap(~Treatment, scales = "free")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Intercepts")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "right",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

#slope model
ggplot() + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = mat_line_warmslope, aes(x = mat_graph, y = estimate))+
  geom_ribbon(data=mat_line_warmslope, aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  geom_point(data = slope_int_datawarm, aes(x= mat_site, y = mean_slope, 
                                            colour = Region)) +
  geom_linerange(data = slope_int_datawarm, aes(x = mat_site, ymax=mean_slope + sd_slope,
                                                ymin = mean_slope - sd_slope, colour = Region)) +
  #facet_wrap(~Treatment, scales = "free")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Slopes")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "right",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

# set up the model matrix for control data
datcon <- dat %>% filter(Treatment == "LocalControl")
datcon$regsit <- factor(paste(datcon$Region, datcon$destSiteID, sep = "_"))

# set up the mati
regsit <- datcon$regsit
N_regsit <- as.numeric(length(levels(regsit)))
exp_dur <- datcon$exp_time_ln

# get temperature for the submodel
datsubmodelcon <-datcon %>% select(regsit, mat_site) %>% unique()
temp <- datsubmodelcon$mat_site

# make the dependent variables
dep_var <- datcon$cwt_mat - datcon$mat_site
N <- as.numeric(length(dep_var))

# specify the model - can use for each analysis
mc_mod <- function(){
  # likelihood
  for(i in 1:N){
    dep_var[i]~dnorm(mu[i], prec)
    mu[i] <- beta0[regsit[i]] + beta1[regsit[i]] * exp_dur[i]
    
    res[i] <- pow(dep_var[i]-mu[i], 2)
    dp.new[i] ~ dnorm(mu[i], prec)
    res.new[i] <- pow(dp.new[i]-mu[i],2)
  }
  #Priors
  for(j in 1:N_regsit){beta0[j]~dnorm(int_mu[j],prec_int)
    int_mu[j] <- beta0_int + beta1_int*temp[j]}
  beta0_int ~ dnorm(0,1.0E-6)
  beta1_int ~ dnorm(0,1.0E-6)
  prec_int~dgamma(0.001,0.001)
  sig_int<- 1/sqrt(prec_int)
  
  for(j in 1:N_regsit){beta1[j]~dnorm(slope_mu[j],prec_slope)
    slope_mu[j] <- beta0_slope + beta1_slope*temp[j]}
  beta0_slope ~ dnorm(0,1.0E-6)
  beta1_slope ~ dnorm(0,1.0E-6)
  prec_slope~dgamma(0.001,0.001)
  sig_slope<- 1/sqrt(prec_slope)
  
  prec~dgamma(0.001,0.001)
  sig <- 1/sqrt(prec)
  
  # derived params
  rss <- sum(res[])
  rss.new <- sum(res.new[])
}

# assemble the data amd params for the model 
jags.data <- list( "dep_var", "regsit", "N", "N_regsit", "temp", "exp_dur")
jags.param <- c("beta0", "beta1","beta0_int", "beta1_int","beta0_slope",
                "beta1_slope", "sig", "sig_slope", "sig_int", "rss", "rss.new", "zlog.lik")

# run the model for mat differences 
cwmatmod_indcon <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                                 n.iter=50000,model.file=mc_mod, n.thin=5, n.chains=3)

cwmatmod_indcon

# posterior check
matdif.paramlistcon <- cwmatmod_indcon$BUGSoutput$sims.list

mean(matdif.paramlistcon$beta1_slope>0)
mean(matdif.paramlistcon$beta1_int>0)

# assemble posteriors
matdif.postindconint <- as_tibble(matdif.paramlistcon$beta0)
base::names(matdif.postindconint) <- c(levels(regsit))
boxplot(matdif.postindconint)

matdif.postindconslope <- as_tibble(matdif.paramlistcon$beta1)
base::names(matdif.postindconslope) <- c(levels(regsit))
boxplot(matdif.postindconslope)

# summarize the posterior
matdif.post.sumind_intcon <- matdif.postindconint %>% 
  select("CH_Calanda_Cal":"US_Montana_Low")%>%
  pivot_longer(cols = "CH_Calanda_Cal":"US_Montana_Low", names_to = c("Region"), 
               values_to = "post") %>%
  mutate(Treatment = rep("control", length(post))) %>%
  group_by(Region) %>% dplyr::summarize(mean_int = mean(post),
                                        sd_int = sd(post)) 
matdif.post.sumind_intcon$Treatment <- rep("control", length(matdif.post.sumind_intcon$mean_int))

matdif.post.sumind_slopecon <-  matdif.postindconslope %>% 
  select("CH_Calanda_Cal":"US_Montana_Low") %>%
  pivot_longer(cols = "CH_Calanda_Cal":"US_Montana_Low", 
               names_to = c("Region"),  values_to = "post")  %>%
  group_by(Region) %>% dplyr::summarize(mean_slope = mean(post), 
                                        sd_slope = sd(post)) 
matdif.post.sumind_slopecon$Treatment <- rep("control", length(matdif.post.sumind_slopecon$mean_slope))

# put together the two datasets
slope_int_datacon <- left_join(matdif.post.sumind_slopecon, 
                               matdif.post.sumind_intcon, by = c("Region", "Treatment")) %>%
  left_join(., rename(datsubmodelcon, Region = "regsit"), by = "Region") %>%
  relocate(Treatment, .after = Region)

slope_int_data <- bind_rows(slope_int_datacon, slope_int_datawarm)

#write_csv(slope_int_data, "SlopIntDataCWMAT_forGraph_log_matmod.csv")
slope_int_data <- read_csv("SlopIntDataCWMAT_forGraph_log_matmod.csv")

# make data for 
mcmc <- cwmatmod_indcon$BUGSoutput$sims.matrix
coefs = as.matrix(select(as_tibble(mcmc), c("beta0_int", "beta1_int")))
mat_graph <- seq(min(slope_int_datacon$mat_site), max(slope_int_datacon$mat_site), by =0.1)
fit_int = coefs %*% t(model.matrix(~ mat_graph))

coefs = as.matrix(select(as_tibble(mcmc), c("beta0_slope", "beta1_slope")))
fit_slope = coefs %*%t(model.matrix(~ mat_graph))

# make a graph
mat_line_conint <- cbind(broom.mixed::tidy(as.mcmc(fit_int), conf.int = TRUE), mat_graph)
mat_line_conslope <- cbind(broom.mixed::tidy(as.mcmc(fit_slope), conf.int = TRUE), mat_graph)
mat_line_con <- bind_rows(list(int = mat_line_conint, slope = mat_line_conslope), .id = "parm")
mat_line <- bind_rows(list(warm = mat_line_warm, control = mat_line_con), .id = "Treatment")

#write_csv(mat_line, "SlopIntSubModDataCWMAT_forGraph_log_matmod.csv")
mat_line <- read_csv("SlopIntSubModDataCWMAT_forGraph_log_matmod.csv")

#intmodel
ggplot() + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = mat_line_conint, aes(x = mat_graph, y = estimate))+
  geom_ribbon(data=mat_line_conint, aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  geom_point(data = slope_int_datacon, aes(x= mat_site, y = mean_int, 
                                            colour = Region)) +
  geom_linerange(data = slope_int_datacon, aes(x = mat_site, ymax=mean_int + sd_int,
                                                ymin = mean_int - sd_int, colour = Region)) +
  #facet_wrap(~Treatment, scales = "free")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Intercepts")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "right",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

#slope model
ggplot() + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = mat_line_conslope, aes(x = mat_graph, y = estimate))+
  geom_ribbon(data=mat_line_conslope, aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  geom_point(data = slope_int_datacon, aes(x= mat_site, y = mean_slope, 
                                            colour = Region)) +
  geom_linerange(data = slope_int_datacon, aes(x = mat_site, ymax=mean_slope + sd_slope,
                                                ymin = mean_slope - sd_slope, colour = Region)) +
  #facet_wrap(~Treatment, scales = "free")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Slopes")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "right",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

# figures for pub
slopintgraph <- slope_int_data %>% ggplot(aes(x= mean_int, y = mean_slope, 
                                              colour = mat_site)) +
  geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_point() + geom_linerange(aes(x = mean_int, ymax=mean_slope + sd_slope,
                                    ymin = mean_slope - sd_slope)) +
  geom_linerange(aes(xmax=mean_int + sd_int, xmin = mean_int - sd_int)) + 
  facet_wrap(~Treatment, scales = "free")+
  #scale_shape_manual(values = 1:nlevels(factor(dat2$Region)))+
  scale_colour_gradient("Site MAT (ºC)", low = "cyan", high = "red") + 
  scale_x_continuous("Initial ΔCWMAT (ªC)") + 
  scale_y_continuous("Thermophilisation rate (ºC/year)")+ # labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = c(.05,.16),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  #guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

slopintgraph

summary(lm(mean_slope~mean_int*Treatment, data = slope_int_data))

ggsave("CWMAT_SlopInt_plot_bysite_wkey_log_TempColor.jpeg", plot=slopintgraph, width = 20,
       height=12, units="cm", dpi=300)

# part b
slopintgraph <- slope_int_data %>% #filter(mean_slope>-3) %>%
  ggplot(aes(x= mean_int, y = mean_slope, 
                                              colour = Region)) +
  geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_point() + geom_linerange(aes(x = mean_int, ymax=mean_slope + sd_slope,
                                    ymin = mean_slope - sd_slope)) +
  geom_linerange(aes(xmax=mean_int + sd_int, xmin = mean_int - sd_int)) + 
  facet_wrap(~Treatment, scales = "free_x")+
  scale_shape_manual(values = 1:nlevels(factor(dat2$Region)))+
  scale_x_continuous("Initial ΔCWMAT (ªC)") + scale_colour_discrete("Region by site") + 
  scale_fill_discrete("Region by site")+
  scale_y_continuous("Slopes")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = "right",
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

slopintgraph
ggsave("CWMAT_SlopInt_plot_bysite_log.jpeg", plot=slopintgraph, width = 12,
       height=6, units="cm", dpi=300)

# intercepts and slopes figure panesl

int_mat_graph <- ggplot(data = filter(mat_line, parm == "int")) + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = filter(mat_line, parm == "int"), aes(x = mat_graph, y = estimate))+
  geom_ribbon(data= filter(mat_line, parm == "int"), aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  facet_wrap(~Treatment, scales = "free")+
  geom_linerange(data = slope_int_data, aes(x = mat_site, ymax=mean_int + sd_int,
                                            ymin = mean_int - sd_int, colour = mat_site)) +
  geom_point(data = slope_int_data, aes(x= mat_site, y = mean_int, 
                                            colour = mat_site)) +
  scale_colour_gradient("Site MAT (ºC)", low = "cyan", high = "red") + 
  #scale_colour_discrete("Region by site") + scale_fill_discrete("Region by site")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Initial ΔCWMAT (ºC)")+  labs(title = "b)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = c(.6,.16),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  #guides(colour = guide_legend(ncol = 2))+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
int_mat_graph

# slope graph

slope_mat_graph <- ggplot(data = filter(mat_line, parm == "slope")) + geom_hline(aes(yintercept = 0), colour = "grey")+
  geom_line(data = filter(mat_line, parm == "slope"), aes(x = mat_graph, y = estimate, linetype = Treatment))+
  geom_ribbon(data= filter(mat_line, parm == "slope"), aes(ymin=conf.low, ymax=conf.high, x=mat_graph), alpha=0.15)+
  facet_wrap(~Treatment, scales = "free")+
  geom_linerange(data = slope_int_data, aes(x = mat_site, ymax=mean_slope + sd_slope,
                                            ymin = mean_slope - sd_slope, colour = mat_site)) +
  geom_point(data = slope_int_data, aes(x= mat_site, y = mean_slope, 
                                        colour = mat_site)) +
  scale_colour_gradient("Site MAT (ºC)", low = "cyan", high = "red") + 
 # scale_colour_discrete("Region by site") + scale_fill_discrete("Region by site")+
  #scale_shape_manual(values = 1:nlevels(factor(Region)))+
  scale_x_continuous("MAT (ºC)") + 
  scale_y_continuous("Thermophilisation rate (ºC/year)")+  labs(title = "c)")+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=8),
        legend.title=element_text(size=6), legend.text=element_text(size=5),
        legend.margin=margin(0,0,0,0),legend.position = c(.10,.16),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+
  guides(linetype = "none")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
slope_mat_graph


jpeg("fig4_transplant_log_wLegends_Matsite.jpeg", width =12, height=18, units="cm", res=300)
grid.arrange(matdif_plot,  int_mat_graph, slope_mat_graph, ncol = 1)
dev.off()

