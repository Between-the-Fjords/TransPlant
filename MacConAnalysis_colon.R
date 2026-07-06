#########################################################
### Macroecological context analyses for colonization ###
#########################################################

# J. Lynn 12.2024

# packages 
library(tidyverse); library(rjags); library(R2jags); library(broom.mixed); library(ggmcmc)
library(loo);library(car);library(reshape2)

# inverse logit function
invlogit<-function(x){a <- exp(x)/(1+exp(x))
a[is.nan(a)]=1
return(a)
}

# load in data 
dat <- read_csv("MacroContextTransplantDat_2.csv")

# clean data by getting rid of what rows without climate data and only keeping
## species range observations of 30 or mor 
dat <- dat %>% filter(!is.na(mat_meandiff)) %>% filter(n_obs_sub>=30) %>% 
  filter(mat_meandiff<20)  # exlude outlier

# look at the correlation of dependent data
cor(dplyr::select(dat, mat_meandiff, sum_meandiff, ap_meandiff, pet_meandiff))

## plot the data for a look 
ggplot(dat, aes(x=mat_meandiff, y = colonization, colour = Treatment))+ geom_jitter(height = 0.2)+
  facet_wrap(~Region)

## begin putting together the models 

## these variables sshould be stable
# random effects
region <- factor(dat$Region)
N_region <- as.numeric(length(levels(region)))
site <- factor(paste(dat$Region,dat$destSiteID, sep = "_"))
site2 <- dat %>% dplyr::select(Region,destSiteID)
site2 <- site2[!duplicated(site2),]
reg_nest <- factor(site2$Region)
N_site <- as.numeric(length(levels(site)))
# split spp and genus
spp_gen <- colsplit(dat$species, " ", c("genus", "species"))
spp <- factor(spp_gen$species)
N_spp <- as.numeric(length(levels(spp)))
genus <- factor(spp_gen$genus)
N_gen <- as.numeric(length(levels(genus)))
sppgen2 <- spp_gen[!duplicated(spp_gen),]
gen_nest <- factor(sppgen2$genus)
exp_dur <- factor(dat$exp_dur)
N_ed <-  as.numeric(length(levels(exp_dur)))

# dependent variables
col <- dat$colonization
N <- as.numeric(length(col))

# specify the model - can use for each analysis
mc_mod <- function(){
  # nested random effects
  for(j in 1:N_region){area[j]~dnorm(0, prec_reg)}
  for(j in 1:N_site){sted[j]~dnorm(area[reg_nest[j]], prec_site)}
  for(j in 1:N_gen){genera[j]~dnorm(0, prec_gen)}
  for(j in 1:N_spp){arter[j]~dnorm(genera[gen_nest[j]], prec_spp)}
  for(j in 1:N_ed){lenge[j]~dnorm(0, prec_ed)}
  # likelihood
  for(i in 1:N){
    col[i]~dbern(mu[i])
    logit(mu[i]) <- inprod(ab, treatmat[i,]) +sted[site[i]]+ area[region[i]]+arter[spp[i]]+
      lenge[exp_dur[i]]+genera[genus[i]]
    
    zlog.lik[i] <- logdensity.bern(col[i], mu[i])
    
    res[i] <- pow(col[i]-mu[i], 2)
    col.new[i] ~ dbern(mu[i])
    res.new[i] <- pow(col.new[i]-mu[i],2)
    
    # Fit assessments: Tjur R-Squared (fit statistic for logistic regression)
    col_true[i] <- ifelse(col[i]>0, mu[i], 0)
    col_false[i] <- ifelse(col[i] == 0, mu[i], 0)
  }
  # priors
  for(j in 1:N_parm){ab[j]~dnorm(0,1.0E-6)}
  
  prec_reg~dgamma(0.001,0.001)
  prec_site~dgamma(0.001,0.001)
  prec_gen~dgamma(0.001,0.001)
  prec_spp~dgamma(0.001,0.001)
  prec_ed~dgamma(0.001,0.001)
  sig_reg <- 1/sqrt(prec_reg)
  sig_site <- 1/sqrt(prec_site)
  sig_gen <- 1/sqrt(prec_gen)
  sig_spp <- 1/sqrt(prec_spp)
  sig_ed <- 1/sqrt(prec_ed)
  
  # derived params
  rss <- sum(res[])
  rss.new <- sum(res.new[])
  
  # Discrepancy measures for entire data set
  col_true.n <- sum(col[]>0)
  col_false.n <- sum(col[]==0)
  r2_tjur <- abs(sum(col_true[])/col_true.n - sum(col_false[])/col_false.n)
}

# parameters that change with each model
# independnt variables
treatmat <- model.matrix(~Treatment*mat_meandiff, data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
matdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

matdif

# traceplots
s <- ggs(as.mcmc(matdif))
ggs_traceplot(s, family="ab") 

# posterior check
matdif.paramlist <- matdif$BUGSoutput$sims.list

plot(matdif.paramlist$rss,matdif.paramlist$rss.new,main="MAT dif colonization",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(matdif.paramlist$rss>matdif.paramlist$rss.new) 

#Information criterion
mat.loglik <- matdif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
mat.waic <- waic(mat.loglik)
mat_eff <- relative_eff(exp(mat.loglik), chain_id=rep(1:1, each=15000))
mat.loo <- loo(mat.loglik, r_eff=mat_eff)

mat.waic
mat.loo

# assemble posteriors
matdif.post <- as_tibble(matdif.paramlist$ab)
names(matdif.post) <- c("intercept", "warm", "mat_dif", "warm_mat_dif")
boxplot(matdif.post)

write_csv(matdif.post, "matDif_mainPost_col.csv")

# summarize the posterior
matdif.post.sum <- matdif.post %>%
  pivot_longer(cols = "intercept":"warm_mat_dif", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                sd_post = sd(post),
                                prob_great_zero = mean(post>0))
write_csv(matdif.post.sum, "matDif_mainPost_summarized_col.csv")

# r^2
mcmc <- matdif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]", "ab[3]", "ab[4]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
mat.newdata <- expand.grid(mat_meandiff = seq(min(dat$mat_meandiff), max(dat$mat_meandiff), length = 100),
                           Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment*mat_meandiff, data= mat.newdata)
fit = coefs %*% t(xmat)
mat.newdata <- mat.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(mat.newdata, "matDiffDatForPlot_col.csv")

dat2 <- dat %>% rename(estimate = colonization)

matdif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=mat_meandiff, y=estimate, colour=factor(Treatment)),alpha=.15, 
              size=0.1,height = 0.1)+
  geom_ribbon(data=mat.newdata, aes(ymin=invlogit(conf.low), ymax=invlogit(conf.high), x=mat_meandiff, 
                                fill=factor(Treatment)), alpha=0.35)+
  geom_line(data=mat.newdata, aes(y = invlogit(estimate), x = mat_meandiff, colour=factor(Treatment)))+
  scale_colour_manual("Treatment", values=c("gray", "red")) + 
  scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_continuous("ΔΜΑΤ (ºC)") + 
  scale_y_continuous("Colonisation")+ 
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=5), legend.text=element_text(size=4),
        legend.margin=margin(0,0,0,0),legend.position = c(0.25,0.6),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
matdif_plot

ggsave("MatDifMain_plot_col.jpeg", plot=matdif_plot, width = 8, height=8, units="cm", dpi=300)

### ap differences
# parameters that change with each model
# independnt variables
treatmat <- model.matrix(~Treatment*ap_meandiff, data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
apdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

apdif

# traceplots
s <- ggs(as.mcmc(apdif))
ggs_traceplot(s, family="ab") 

# posterior check
apdif.paramlist <- apdif$BUGSoutput$sims.list

plot(apdif.paramlist$rss,apdif.paramlist$rss.new,main="AP dif Cover",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(apdif.paramlist$rss>apdif.paramlist$rss.new) 

#Information criterion
ap.loglik <- apdif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
ap.waic <- waic(ap.loglik)
ap_eff <- relative_eff(exp(ap.loglik), chain_id=rep(1:1, each=15000))
ap.loo <- loo(ap.loglik, r_eff=ap_eff)

ap.waic
ap.loo

# assemble posteriors
apdif.post <- as_tibble(apdif.paramlist$ab)
names(apdif.post) <- c("intercept", "warm", "ap_dif", "warm_ap_dif")
boxplot(apdif.post)

write_csv(apdif.post, "apDif_mainPost_col.csv")

# summarize the posterior
apdif.post.sum <- apdif.post %>%
  pivot_longer(cols = "intercept":"warm_ap_dif", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(apdif.post.sum, "apDif_mainPost_summarized_col.csv")

# r^2
mcmc <- apdif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]", "ab[3]", "ab[4]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
ap.newdata <- expand.grid(ap_meandiff = seq(min(dat$ap_meandiff), max(dat$ap_meandiff), length = 100),
                           Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment*ap_meandiff, data= ap.newdata)
fit = coefs %*% t(xmat)
ap.newdata <- ap.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(ap.newdata, "apDiffDatForPlot_col.csv")

apdif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=ap_meandiff, y=estimate, colour=factor(Treatment)),
             alpha=.15, size=0.1, height = 0.1)+
  geom_ribbon(data=ap.newdata, aes(ymin=invlogit(conf.low), ymax=invlogit(conf.high), x=ap_meandiff, 
                                    fill=factor(Treatment)), alpha=0.35)+
  geom_line(data=ap.newdata, aes(y = invlogit(estimate), x = ap_meandiff, colour=factor(Treatment)))+
  scale_colour_manual("Treatment", values=c("gray", "red")) + 
  scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_continuous("ΔAP (m)") + 
  scale_y_continuous("ΔRelative cover")+ 
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=5), legend.text=element_text(size=4),
        legend.margin=margin(0,0,0,0),legend.position = c(0.85,0.7),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

apdif_plot
ggsave("APDifMain_plot_col.jpeg", plot=apdif_plot, width = 8, height=8, units="cm", dpi=300)

### pet differecnes differences
# parameters that change with each model
# independnt variables
treatmat <- model.matrix(~Treatment*pet_meandiff, data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
petdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                       n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

petdif

# traceplots
s <- ggs(as.mcmc(petdif))
ggs_traceplot(s, family="ab") 

# posterior check
petdif.paramlist <- petdif$BUGSoutput$sims.list

plot(petdif.paramlist$rss,petdif.paramlist$rss.new,main="AP dif Extinciton",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(petdif.paramlist$rss>petdif.paramlist$rss.new) 

#Information criterion
pet.loglik <- petdif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
pet.waic <- waic(pet.loglik)
pet_eff <- relative_eff(exp(pet.loglik), chain_id=rep(1:1, each=15000))
pet.loo <- loo(pet.loglik, r_eff=pet_eff)

pet.waic
pet.loo

# assemble posteriors
petdif.post <- as_tibble(petdif.paramlist$ab)
names(petdif.post) <- c("intercept", "warm", "pet_dif", "warm_pet_dif")
boxplot(petdif.post)

write_csv(petdif.post, "petDif_mainPost_col.csv")

# summarize the posterior
petdif.post.sum <- petdif.post %>%
  pivot_longer(cols = "intercept":"warm_pet_dif", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(petdif.post.sum, "petDif_mainPost_summarized_col.csv")

# r^2
mcmc <- petdif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]", "ab[3]", "ab[4]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
pet.newdata <- expand.grid(pet_meandiff = seq(min(dat$pet_meandiff), max(dat$pet_meandiff), length = 100),
                          Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment*pet_meandiff, data= pet.newdata)
fit = coefs %*% t(xmat)
pet.newdata <- pet.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(pet.newdata, "petDiffDatForPlot_col.csv")

petdif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=pet_meandiff, y=estimate, colour=factor(Treatment)),
              alpha=.15, size=0.1, height = 0.1)+
  geom_ribbon(data=pet.newdata, aes(ymin= invlogit(conf.low), ymax= invlogit(conf.high),
                                    x=pet_meandiff, fill=factor(Treatment)), alpha=0.35)+
  geom_line(data=pet.newdata, aes(y = invlogit(estimate), x = pet_meandiff, colour=factor(Treatment)))+
  scale_colour_manual("Treatment", values=c("gray", "red")) + 
  scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_continuous("ΔPET (m)") + 
  scale_y_continuous("ΔRelative cover")+ 
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=5), legend.text=element_text(size=4),
        legend.margin=margin(0,0,0,0),legend.position = c(0.85,0.25),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())

petdif_plot
ggsave("PETDifMain_plot_col.jpeg", plot=petdif_plot, width = 8, height=8, units="cm", dpi=300)

# independnt variables mat interacting with AP
treatmat <- model.matrix(~Treatment*mat_meandiff*ap_meandiff, data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
matapdif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

matapdif

# traceplots
s <- ggs(as.mcmc(matapdif))
ggs_traceplot(s, family="ab") 

# posterior check
matapdif.paramlist <- matapdif$BUGSoutput$sims.list

plot(matapdif.paramlist$rss,matapdif.paramlist$rss.new,main="Mat dif Cover",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(matapdif.paramlist$rss>matapdif.paramlist$rss.new) 

#Information criterion
matap.loglik <- matapdif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
matap.waic <- waic(matap.loglik)
matap_eff <- relative_eff(exp(matap.loglik), chain_id=rep(1:1, each=15000))
matap.loo <- loo(matap.loglik, r_eff=matap_eff)

matap.waic
matap.loo

# assemble posteriors
matapdif.post <- as_tibble(matapdif.paramlist$ab)
names(matapdif.post) <- c("intercept", "warm", "mat_dif", "ap_dif", "warm_mat_dif", "warm_ap_dif", 
                        "mat_dif_ap_dif", "warm_mat_dif_ap_dif")
boxplot(matapdif.post)

write_csv(matapdif.post, "matapDif_mainPost_col.csv")

# summarize the posterior
matapdif.post.sum <- matapdif.post %>%
  pivot_longer(cols = "intercept":"warm_mat_dif_ap_dif", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(matapdif.post.sum, "matApDif_mainPost_summarized_col.csv")

# r^2
mcmc <- matapdif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]", "ab[3]", "ab[4]","ab[5]", "ab[6]", "ab[7]", "ab[8]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
matap.newdata <- expand.grid(mat_meandiff = seq(min(dat$mat_meandiff), max(dat$mat_meandiff), length = 50),
                             ap_meandiff = seq(min(dat$ap_meandiff), max(dat$ap_meandiff), length = 50),
                           Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment*mat_meandiff*ap_meandiff, data= matap.newdata)
fit = coefs %*% t(xmat)
matap.newdata <- matap.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(matap.newdata, "matApDiffDatForPlot_col.csv")

matapdif_plot<- ggplot()+ 
  geom_tile(data=matap.newdata, aes(y = ap_meandiff, x = mat_meandiff, fill=invlogit(estimate)))+
  geom_point(data=dat2, aes(x=mat_meandiff, y=ap_meandiff), size=0.1)+
  geom_vline(xintercept = 0, alpha = 0.5)+geom_hline(yintercept = 0,alpha = 0.5)+
  facet_grid(~Treatment)+scale_x_continuous("ΔΜΑΤ (ºC)") +scale_y_continuous("ΔAP (m)")+ 
  scale_fill_gradient2("Colonization", low="#833D91",mid="gray", high="#407F4A", midpoint = 0.35)+
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),legend.margin=margin(0,0,0,0),
        legend.title=element_text(size=5), legend.text=element_text(size=4),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
matapdif_plot

ggsave("MatAPDifMain_plot_col.jpeg", plot=matapdif_plot, width = 16, height=8, units="cm", dpi=300)

# parameters that change with each model mat squared
# independnt variables
treatmat <- model.matrix(~Treatment*poly(mat_meandiff, degree = 2, raw = TRUE), data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
mat2dif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                        n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

mat2dif

# traceplots
s <- ggs(as.mcmc(mat2dif))
ggs_traceplot(s, family="ab") 

# posterior check
mat2dif.paramlist <- mat2dif$BUGSoutput$sims.list

plot(mat2dif.paramlist$rss,mat2dif.paramlist$rss.new,main="Mat dif 2 ext",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(mat2dif.paramlist$rss>mat2dif.paramlist$rss.new) 

#Information criterion
mat2.loglik <- mat2dif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
mat2.waic <- waic(mat2.loglik)
mat2_eff <- relative_eff(exp(mat2.loglik), chain_id=rep(1:1, each=15000))
mat2.loo <- loo(mat2.loglik, r_eff=mat2_eff)

mat2.waic
mat2.loo

# assemble posteriors
mat2dif.post <- as_tibble(mat2dif.paramlist$ab)
names(mat2dif.post) <- c("intercept", "warm", "mat_dif", "mat_dif_2",  "warm_mat_dif",
                        "warm_mat_dif_2")
boxplot(mat2dif.post)

write_csv(mat2dif.post, "mat2Dif_mainPost_col.csv")

# summarize the posterior
mat2dif.post.sum <- mat2dif.post %>%
  pivot_longer(cols = "intercept":"warm_mat_dif_2", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(mat2dif.post.sum, "mat2Dif_mainPost_summarized_col.csv")

# r^2
mcmc <- mat2dif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]", "ab[3]", "ab[4]","ab[5]", "ab[6]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
mat2.newdata <- expand.grid(mat_meandiff = seq(min(dat$mat_meandiff), max(dat$mat_meandiff), length = 100),
                           Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment*poly(mat_meandiff, degree = 2, raw = TRUE), data= mat2.newdata)
fit = coefs %*% t(xmat)
mat2.newdata <- mat2.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(mat2.newdata, "mat2DiffDatForPlot_col.csv")

mat2dif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=mat_meandiff, y=estimate, colour=factor(Treatment)),
              height = 0.1, alpha=.15, size=0.1)+
  geom_ribbon(data=mat2.newdata, aes(ymin=invlogit(conf.low), ymax=invlogit(conf.high), x=mat_meandiff, 
                                    fill=factor(Treatment)), alpha=0.35)+
  geom_line(data=mat2.newdata, aes(y = invlogit(estimate), x = mat_meandiff, colour=factor(Treatment)))+
  scale_colour_manual("Treatment", values=c("gray", "red")) + 
  scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_continuous("ΔΜΑΤ (ºC)") + 
  scale_y_continuous("ΔRelative cover")+ 
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=5), legend.text=element_text(size=4),
        legend.margin=margin(0,0,0,0),legend.position = c(0.85,0.25),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
mat2dif_plot

ggsave("Mat2DifMain_plot_col.jpeg", plot=mat2dif_plot, width = 8, height=8, units="cm", dpi=300)

# parameters that change with each model null
# independnt variables
treatmat <- model.matrix(~Treatment, data = dat)
N_parm <- as.numeric(ncol(treatmat))

# assemble the data amd params for the model 
jags.data <- list("region", "site", "spp", "treatmat", "col", "reg_nest", "exp_dur", "genus",
                  "N", "N_region", "N_site", "N_spp", "N_parm", "N_gen", "N_ed", "gen_nest")
jags.param <- c("ab", "sig", "sig_reg", "sig_site", "sig_spp", "sig_ed", "sig_gen",
                "rss", "rss.new", "zlog.lik", "r2_tjur")

# run the model for mat differences 
nulldif <- jags.parallel(data=jags.data,inits=NULL, parameters.to.save=jags.param,
                         n.iter=50000 ,model.file=mc_mod, n.thin=5, n.chains=3)

nulldif

# traceplots
s <- ggs(as.mcmc(nulldif))
ggs_traceplot(s, family="ab") 

# posterior check
nulldif.paramlist <- nulldif$BUGSoutput$sims.list

plot(nulldif.paramlist$rss,nulldif.paramlist$rss.new,main="Mat dif 2 ext",
     xlab="SSQ observed", ylab="SSQ simulated")
abline(0,1)

mean(nulldif.paramlist$rss>nulldif.paramlist$rss.new) 

#Information criterion
null.loglik <- nulldif.paramlist$zlog.lik # loglikelihood value for each of n (where n = # saved) ith log-likelihood values  (n x i matrix)

# use waic() function from loo package
null.waic <- waic(null.loglik)
null_eff <- relative_eff(exp(null.loglik), chain_id=rep(1:1, each=15000))
null.loo <- loo(null.loglik, r_eff=null_eff)

null.waic
null.loo

# assemble posteriors
nulldif.post <- as_tibble(nulldif.paramlist$ab)
names(nulldif.post) <- c("intercept", "warm")
boxplot(nulldif.post)

write_csv(nulldif.post, "nullDif_mainPost_col.csv")

# summarize the posterior
nulldif.post.sum <- nulldif.post %>%
  pivot_longer(cols = "intercept":"warm", 
               names_to = "param", values_to = "post") %>%
  group_by(param) %>% dplyr::summarize(mean_post = mean(post),
                                       sd_post = sd(post),
                                       prob_great_zero = mean(post>0))
write_csv(nulldif.post.sum, "nullDif_mainPost_summarized_col.csv")

# r^2
mcmc <- nulldif$BUGSoutput$sims.matrix
coefs = mcmc[, c("ab[1]", "ab[2]")]
fit = coefs %*% t(treatmat)
resid = sweep(invlogit(fit), 2, col, "-")
var_f = apply(invlogit(fit), 1, var)
var_e = apply(resid, 1, var)
R2 = var_f/(var_f + var_e)
tidyMCMC(as.mcmc(R2), conf.int = TRUE, conf.method = "HPDinterval")

# make a graph
null.newdata <- expand.grid(Treatment = c("LocalControl", "Warm"))

xmat <- model.matrix(~Treatment, data= null.newdata)
fit = coefs %*% t(xmat)
null.newdata <- null.newdata %>% cbind(tidyMCMC(fit, conf.int = TRUE))
write_csv(null.newdata, "nullDiffDatForPlot_col.csv")

nulldif_plot<- ggplot()+ 
  geom_jitter(data=dat2, aes(x=factor(Treatment), y=estimate, colour=factor(Treatment)), height = 0.1, alpha=.15, size=0.1)+
  geom_point(data=null.newdata, aes(y = invlogit(estimate), x = factor(Treatment), colour=factor(Treatment)))+
  geom_errorbar(data=null.newdata, aes(y = invlogit(estimate), x = factor(Treatment), colour=factor(Treatment),
                                       ymin = invlogit(conf.low), ymax = invlogit(conf.high)))+
  scale_colour_manual("Treatment", values=c("gray", "red")) + 
  scale_fill_manual("Treatment", values=c("dark gray", "red")) + 
  scale_x_discrete("Treatment")+
  scale_y_continuous("Colonization")+ 
  theme(axis.text.y = element_text(size=7,colour= "black"),
        axis.text.x= element_text(size=7, colour="black"), 
        axis.title=element_text(size=7),strip.text=element_text(size=5),
        plot.title=element_text(size=7),
        legend.title=element_text(size=5), legend.text=element_text(size=4), 
        legend.margin=margin(0,0,0,0),legend.position = c(0.8,0.3),
        legend.box.margin=margin(-10,-2,-10,-5),legend.justification="left",
        legend.key.size = unit(0.15, "cm"))+ #labs(title="Colonization probability")+
  theme(panel.background = element_rect(fill='white', colour='black'))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
nulldif_plot

ggsave("nullDifMain_plot_col.jpeg", plot=mat2dif_plot, width = 8, height=8, units="cm", dpi=300)

