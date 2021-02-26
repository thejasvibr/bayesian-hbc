# Pondfrog problem 
library(blmeco)
library(rstanarm)
library(arm)

set.seed(252532)

data(pondfrog)

# measurements of the pond Lättloch
newdat <- data.frame(fish=1, 
                     ph=6.4, temp=11.6, height=634, 
                     region=factor("north", levels=levels(pondfrog$region)),
                     waterdepth=1.3, surfacearea=14, vegdensity=1)

str(newdat)   # make sure, all variables are defined in exact the same way (check factor levels!)
str(pondfrog)


# Let's remove the height - because 

eqn.0 <- 'frog~fish+waterdepth+temp+vegdensity+ph'

first.model <- lm(eqn.0,data=pondfrog)
sim.first <- sim(first.model, n.sim=5000)

mean.coefs <- apply(sim.first@coef,2,mean)
lower.pctile.coefs <- apply(sim.first@coef,2,quantile,0.03)
higher.pctile.coefs <- apply(sim.first@coef,2,quantile,0.97)
coefs.range <- rbind(lower.pctile.coefs, mean.coefs, higher.pctile.coefs)

# Run mean prediction for the given model
mean.predn <- coefs.range[2,1] +sum(coefs.range[2,2:6]*newdat[c('fish','waterdepth','temp','vegdensity','ph')])

# Run the prediction interval for the given data across all potentially matching coefficients:
predn.interval <- sim.first@coef[,1] + sim.first@coef[,2:6]*newdat[c('fish','waterdepth','temp','vegdensity','ph')]

coefs.multip <- rep(0,nrow(sim.first@coef))
for (rownum in seq(1,nrow(sim.first@coef)))
{
  coefs.multip[rownum]<-sum(sim.first@coef[rownum,2:6]*newdat[c('fish','waterdepth','temp','vegdensity','ph')])
}
predns.all <- coefs.multip+sim.first@coef[,1]
predn.pctile.range <- quantile(predns.all, c(0.03,0.97))
print(predn.pctile.range)
print(paste('Mean prediction:',mean.predn))

# let's check the actual posterior prediction values:
value.posterior <- rep(0,length(predns.all))
for (i in seq(1,length(predns.all))){
value.posterior[i] = rnorm(1,predns.all[i], sim.first@sigma[i])
}


HPDinterval(as.mcmc(ratioba))
quantile(value.posterior, c(0.03,0.97))

