## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July/Aug 2019
## Purpose: Example script for a simple (no-covariate) geostatistical
##          model fit to point data. 
##          In this example, we use data from Pakistan for the TB 
##
##          Version 1: No covariates, no temporal. Binomial/Gaussian model only
##                     invlink(y_s) = a + Z_s
## ###########################################################



## ###########################################################
## User inputs. 
## TODO: update this to be a config file

# model family
family <- 'binomial'

# number of drawws from the posterior to take. 
# 100 good for testing, 1000 good for final run (but slower and more RAM-intensive)
ndraws <- 100

# outcome variable name
outcome_varname <- 'dpt'

# where your local copy of the generic_mapper repo is (clone from here: https://github.com/rburstein-IDM/generic_mapper)
codepath  <- 'C:/Users/rburstein/OneDrive - IDMOD/code/generic_mapper'


## ###########################################################
## Some Setup: Libraries, paths, loading in data


## Check Libraries and load them
## TODO: do this a docker so this step isnt necessary
libs <- c('data.table','ggplot2','ggridges','lubridate','fasterize','gbm','dismo','rgeos','RCurl',
          'raster','sf','mgcv','gridExtra','gstat','polycor','INLA','matrixStats','RColorBrewer')
for(l in libs){
  if(!require(l,character.only = TRUE, quietly = TRUE)){
    message( sprintf('Did not have the required package << %s >> installed. Downloading now ... ',l))
    install.packages(l) 
  }
  library(l,character.only=TRUE)
}

## locations - relative to user, but in AMUG Dropbox. Make sure you have access to AMUG to use
user      <- Sys.info()[['user']]
root      <- sprintf('C:/Users/%s/Dropbox (IDM)/AMUG/generic_mapper',user)

# Source the utility functions 
source(sprintf('%s/utils.R', codepath))

## Load in data and shapes
# prepped data
d       <- fread(sprintf('%s/data/prepped_input/XXXX.csv',datpath))

# load shapefile
shp <- st_read(dsn = sprintf('%s/WHO_POLIO_GLOBAL_GEODATABASE.gdb',datpath), layer = 'GLOBAL_ADM2')
shp <- subset(shp, as.Date(ENDDATE) >= Sys.Date() & ADM0_VIZ_NAME == "Democratic Republic of the Congo")




## ###########################################################
## Clean up and prep data


# make a raster version of the country shapefile
ext_raster <- trim(fasterize(shp, dpt1))

# rasterize the gadm file
rshp <- fasterize(sf = shp, raster = dpt3[[1]], field = 'id')

# make a prediction frame representing the full raster
predfr <- data.table(idx = 1:ncell(get(rlist[1])), adm2 = as.vector(rshp))

# add lat and lon to prediction frame
predfr[, lon := coordinates(ext_raster)[,1]]
predfr[, lat := coordinates(ext_raster)[,2]]




## ###########################################################
## Basic geostats model using INLA (currently with most settings at default)

# get together all the various objects needed for inla to fit the model
input_obj <- prep_inla_objects(data = d, outcome_varname = outcome_varname, ss_varname = 'N')

  
# fit inla model
res_fit <- inla(input_obj[['formula']],
                data = inla.stack.data(input_obj[['stack']]),
                control.predictor = list(A       = inla.stack.A(input_obj[['stack']]),
                                         link    = 1,
                                         compute = FALSE),
                control.fixed     = list(expand.factor.strategy = 'inla'),
                family            = family,
                num.threads       = 1,
                Ntrials           = input_obj[['N']],
                verbose           = TRUE,
                keep              = FALSE)

summary(res_fit)
  
 

# do prediction from the fitted model
pred <- inla_predict(fitted = res_fit, input = input_obj, ndraws = ndraws, predfr = predfr, ext_raster = ext_raster)


# the outputs of prediction with match ext_raster. In the following, aggregate them to the inputted shapefile
agg_res  <- aggregate_results()





## ###########################################################
## Plot some outputs



# shapefile and plot
shp <- merge(shp, out, by = 'id', all.x = TRUE )

# plots
g1 <- ggplot(shp) + geom_sf(aes(fill = ch_median*100)) + coord_sf(datum = NA) + 
  theme_minimal() + theme(legend.position="bottom") +
  scale_fill_gradientn(colours=brewer.pal(7,"YlGnBu"), 
                       name = 'Seropositive Children (%) \nMedian Estimate') 
g2 <- ggplot(shp) + geom_sf(aes(fill = ch_lower*100)) + coord_sf(datum = NA) + 
  theme_minimal() + #theme(legend.position="bottom")+
  scale_fill_gradientn(limits = c(60,100), colours=brewer.pal(7,"Greys"), name = 'Lower UI') 
g3 <- ggplot(shp) + geom_sf(aes(fill = ch_upper*100)) + coord_sf(datum = NA) + 
  theme_minimal() + #theme(legend.position="bottom")+
  scale_fill_gradientn(limits = c(60,100), colours=brewer.pal(7,"Greys"), name = 'Upper UI') 

g1
grid.arrange(g2,g3,nrow=2)


g1 <- ggplot(shp) + geom_sf(aes(fill = ad_median*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(colours=brewer.pal(7,"YlGnBu"), 
                       name = 'Seropositive Adults (%) \nMedian Estimate') + theme(legend.position="bottom") 
g2 <- ggplot(shp) + geom_sf(aes(fill = ad_lower*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(limits = c(35,82), colours=brewer.pal(7,"Greys"), name = 'Lower UI') 
g3 <- ggplot(shp) + geom_sf(aes(fill = ad_upper*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(limits = c(35,82), colours=brewer.pal(7,"Greys"), name = 'Upper UI') 

g1
grid.arrange(g2,g3,nrow=2)



g1 <- ggplot(shp) + geom_sf(aes(fill = df_median*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(colours=brewer.pal(7,"YlOrRd"), 
                       name = 'Difference in Seroprevalence\nMedian Estimate')  + theme(legend.position="bottom") 
g2 <- ggplot(shp) + geom_sf(aes(fill = df_lower*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(limits = c(-10,60), colours=brewer.pal(7,"Greys"), name = 'Lower UI') 
g3 <- ggplot(shp) + geom_sf(aes(fill = df_upper*100)) + coord_sf(datum = NA) + theme_minimal() +
  scale_fill_gradientn(limits = c(-10,60), colours=brewer.pal(7,"Greys"), name = 'Upper UI') 

g1
grid.arrange(g2,g3,nrow=2)




## Plot Covariate effects ladderplot
tmp <- data.table(summary(ch_inla_mod$model_object)$fixed)
tmp$var <- c('intercept',vars)

ggplot(tmp, aes(x=var,ymin=`0.025quant`,ymax=`0.975quant`,y=mean)) +
  geom_errorbar(width = .1) + theme_minimal() + 
  geom_hline(yintercept =0,color='red') +
  geom_point() +
  xlab('Covariate Name') + ylab('<--- (-)   Effect size (log-odds)   (+) --->')


# univariates
outy <- data.table()
for(v in paste0(tmp$var[-1])){ #},'_cs')){
  form <- as.formula(sprintf('cbind( type2, N-type2 ) ~ %s',v))
  m    <- summary(glm(form,data=dagg[adult==0],family='binomial'))$coefficients[2,1]
  sd   <- summary(glm(form,data=dagg[adult==0],family='binomial'))$coefficients[2,2]
  tmpy <- data.table(var=v,mean=exp(m),lower=exp(m-1.96*sd),upper=exp(m+1.96*sd))
  outy <- rbind(outy,tmpy)
}
  
ggplot(outy, aes(x=var,ymin=lower,ymax=upper,y=mean)) +
  geom_errorbar(width = .1) + theme_minimal() + 
  geom_hline(yintercept =1,color='red') +
  geom_point() +
  xlab('Covariate Name') + ylab('<--- (-)   Odds Ratio   (+) --->')




## Do rank districts











######################################
# how do estimates compare to cluster level data?
res <- data.table(y = dagg[[outcometype]], N = dagg$N,
              raster::extract(pred, cbind(dagg$lon,dagg$lat)))
res[, prev := y/N +rnorm(nrow(res),0,0.01)] # jitter for plotting

ggplot(res, aes(x=prev,y=mean,ymin=lower,ymax=upper)) + theme_minimal() + 
  geom_abline(intercept=0,slope=1,color='red') +
  xlab('Empirical Estimate') + ylab('Model Estimate') +
  geom_errorbar(alpha=0.25) + geom_point(alpha=0.75,size=0.5)




######################################
# population at risk

pred_popatrisk <- inlapredfr*cbind(predfr$u5pop)

# about 1 million children at risk, type 2
sum(pred_popatrisk[,1], na.rm = TRUE)
sum(pred_popatrisk[,2], na.rm = TRUE)
sum(pred_popatrisk[,3], na.rm = TRUE)

plot(insertRaster(ext_raster,cbind(log(pred_popatrisk)[,1])), main='log population at risk')





## TODOS: 
# add a campaign covariate
# quality, residual analysis
# adults map -> joint fit with kids, and with other types















