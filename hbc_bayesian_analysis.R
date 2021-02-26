### 
# Working on the individual call analysis 
#
#
#
### 
library(arm)
library(stringi)
library(stringr)
d <- read.csv('one_call_per_row_2020-12-17.csv')
# create a single/multi-bat categorical variable
d['multi_bat'] <- factor(d['num_bats']>1)
d['annot_id'] <- stri_sub(d$video_annot_id,8,-1)
d['annot_id'] <- str_replace_all(d$annot_id, '21502300','2123')
# sort the dataset by time first
time_ordering <- stri_order(d$annot_id,numeric=TRUE)
d_tordered <- d[time_ordering,]

# defien the model  for tfm duration
tfmdurn_model <- lm(tfm_duration~multi_bat, data=d_tordered)
# get the posterior distribution of the estimated coefficients
sim_tfmdurn <- sim(tfmdurn_model, n.sim=2000)
# plot the raw data
plot(tfm_duration~multi_bat, data=d_tordered)
# check the model diagnostics 
plot(tfmdurn_model)

# - is there an influence of the recording hour? 
split.and.join <- function(x){
  splitx <- str_split(x,'_')
  jointx <- paste(splitx[[1]][1],'_',splitx[[1]][2],sep='')
  jointx}
d_tordered$rechour <- sapply(d_tordered$annot_id,split.and.join)

tfmdurn.model.lmm <- lmer(tfm_duration~multi_bat+(1|rechour), data=d_tordered)

# RUN THE SAME MODEL IN RSTARM WITH MCMC -- it's more robust. 
# The results of the lmm show that recording hour does not make any difference. 
# The ACF plot is also something to add to the SI nd also mention it  in the paper. 

# If there is no pattern in the ACF -- you wouldn't expectto see any effect of time - but, it may be 
# interesting to put it into action. You could also just put in the covariate, and report the credible interval 
# to *show* that there is no/low effect. 

# to see if there is any kind of temporal structure
# in the data - run the autocorrelation of the 
# residuals 
# First sort the data according to time and run an autocorrelation
acf(resid(tfmdurn_model))

# Now check the probability that the efect is >0

p_gq_zero = mean(sim_tfmdurn@coef[,2]>0)
# the Fawcett et al. study found a mean difference of 1.8 ms
p_gq_18ms = mean(sim_tfmdurn@coef[,2]>=1.8*10^-3)
p_gq_18ms


# ---# Now, let's include the priors from Fawcett et al. datlm(tfm_duration~multi_bat, data=d_tordered)
library(rstan)
library(rstanarm)

###### The whole audio analysis data set
wholeaudio <- read.csv('obs_nonsilent_measurements_20dBthreshold.csv')
wholeaudio['group.status'] = (wholeaudio$num_bats>1)*1 # 0 is single and 1 is multi
wholeaudio['annot_id'] <- stri_sub(wholeaudio$video_annot_id,8,-1)
wholeaudio['annot_id'] <- str_replace_all(wholeaudio$annot_id, '21502300','2123')

time_ordering <- stri_order(wholeaudio$annot_id,numeric=TRUE)
wholeaudio <- wholeaudio[time_ordering,]
wholeaudio$rechour <- sapply(wholeaudio$annot_id,split.and.join)

# %% Peak-Amplitude 
peak_amp <- subset(wholeaudio, measurement=='rms')
peak_amp$exprms <- log(peak_amp$value)
peak_amp$sqrt <- sqrt(peak_amp$value)


lmm.model <- lmer(exprms ~ group.status + (1|video_annot_id), data=peak_amp)
model2.formula <- formula(exprms ~ group.status + (1|video_annot_id)+(1|rechour))
lmm.model2 <- lmer(model2.formula, data=peak_amp)


plot(resid(lmm.model)~group.status, data=peak_amp)
qqnorm(resid(lmm.model))
qqline(resid(lmm.model))   
boxplot(resid(lmm.model)~group.status, data=peak_amp)

# diagnose for model 2

plot(resid(lmm.model2)~group.status, data=peak_amp)
qqnorm(resid(lmm.model2))
qqline(resid(lmm.model2))   
boxplot(resid(lmm.model2)~group.status, data=peak_amp)


#### There are a few points throwing off the model -- re-run the model including only 80%of all the residuals 
pctile.lim <- 0.1
pctile.80.resid <- quantile(resid(lmm.model2),c(pctile.lim, 1-pctile.lim))
resid80pctile.rows <- resid(lmm.model2)>pctile.80.resid[1] & resid(lmm.model2)<pctile.80.resid[2]

goodpreds.data <- peak_amp[resid80pctile.rows,] # subset only those data that are well predicted
lmm.resid80pct <- lmer(model2.formula, data=goodpreds.data)

qqnorm(resid(lmm.resid80pct))
qqline(resid(lmm.resid80pct))   
boxplot(resid(lmm.resid80pct)~group.status, data=goodp)

## You can trust 


# simulate the bayesian posterior distribution 
lmm.sim <- sim(lmm.model, n.sim=2000)
group.diffs <- quantile(lmm.sim@fixef, c(0.03,0.97)) # report the 97% credible interval
group.diffs.rms <- exp(group.diffs)

# %% Terminal frequency 
fmt<-subset(wholeaudio, measurement=='fm_terminal_freqs')
fmt <- fmt[!is.na(fmt$value),]
fmt$value <- fmt$value*10^-3 # convert it all to kHz
tfm.lmm.mod <- lmer(value ~ group.status + (1|unique_window_id), data=fmt)
tfm.lmm.mod2 <- lmer(value ~ group.status + (1|video_annot_id), data=fmt)
tfm.lmm.mod3 <- lmer(value ~ group.status + (1|video_annot_id/segment_number), data=fmt)

# -- # 
# re-fit the model exluding the points with huuuuge residuals 

# -- inspect model fit
qqnorm(resid(tfm.lmm.mod3))
qqline(resid(tfm.lmm.mod3))   
boxplot(resid(tfm.lmm.mod3)~group.status, ylim=c(-100,20), data=fmt) # a few huge residuals, but residuals are centred around 0. 

# new dataset
fmt2 <- fmt[abs(resid(tfm.lmm.mod3))<20,]
tfm2.lmm.mod3 <- lmer(value ~ group.status + (1|video_annot_id/segment_number), data=fmt2)
# -- inspect model fit after removing big residual points
qqnorm(resid(tfm2.lmm.mod3))
qqline(resid(tfm2.lmm.mod3))   
boxplot(resid(tfm2.lmm.mod3)~group.status,  data=fmt2) # a few huge residuals, but residuals are centred around 0. 


# -- generate posterior distribution
sim.tfm <- sim(tfm2.lmm.mod3, n.sim=2000)
quantile(sim.tfm@fixef[,2],c(0.03,0.5,0.97))

##################### Now do dominant frequency
domfreq <- subset(wholeaudio, measurement=='dominant_frequencies')

domfreq$value <- domfreq$value*10^-3 # convert it all to kHz
domfreq$logfreq <- log(domfreq$value)
domfreq.lmm.mod <- lmer(logfreq ~ group.status + (1|unique_window_id), data=domfreq)
domfreq.lmm.mod2 <- lmer(logfreq ~ group.status + (1|video_annot_id), data=domfreq)
domfreq.lmm.mod3 <- lmer(logfreq ~ group.status + (1|video_annot_id/segment_number), data=domfreq)

# mod3 shows no variation at the segment number level, and so let's just keep mod2

# 1. try using a t-distribution using brms package. 
# 2.  Only include the cenral 80%ile of data and decide if the coefficients change or not. 


#  -- INSPECTION
qqnorm(resid(domfreq.lmm.mod2))
qqline(resid(domfreq.lmm.mod2))   

boxplot(resid(domfreq.lmm.mod2)~group.status,  data=domfreq) # a few huge residuals, but residuals are centred around 0. 

