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
temptrash <- sprintf('%s/temptrash',root)


## Source the utility functions directly from github repo. 
script <- getURL("URL HERE", ssl.verifypeer = FALSE)
eval(parse(text = script))

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


# make a wrapper function to do the model  heavy lifting

fitinlamodel <- function(data = dagg, ndraws = 100, fevars, outcometype){
  
  # make design matrix
  message(' .. data munge')
  narmdagg      <- na.omit(data)
  outcome       <- narmdagg[[outcometype]]
  N             <- narmdagg$N      
  locs          <- cbind(narmdagg$lon, narmdagg$lat)
  narmdagg      <- narmdagg[, fevars, with=FALSE]
  design_matrix <- data.frame(int = 1, narmdagg)
  
  # make a mesh
  mesh <- inla.mesh.2d(loc      = locs,
                       max.edge = c(.5,5),
                       offset   = 1,
                       cutoff   = 0.5)
  plot(mesh);  points(locs, col = 'red')
  
  
  # make spde object
  spde  <- inla.spde2.matern(mesh = mesh,  alpha = 2)
  space <- inla.spde.make.index("space", n.spde = spde$n.spde)
  
  
  # make projector matrix
  A <- inla.spde.make.A(mesh, loc = locs)
  
  
  # make a data stack
  stack.obs <- inla.stack(data    = list(o = outcome),
                          A       = list(A, 1),
                          effects = list(space, design_matrix),
                          tag     = 'est')
  
  
  # inla formula
  formla <- as.formula( sprintf ('o ~ -1 + int + %s + f(space, model = spde)', paste(fevars,collapse='+')) )
  
  
  # fit inla model
  message(' .. model fit')
  res_fit <- inla(formla,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A       = inla.stack.A(stack.obs),
                                           link    = 1,
                                           compute = FALSE),
                  control.compute   = list(dic     = TRUE,
                                           cpo     = TRUE,
                                           config  = TRUE),
                  control.fixed     = list(expand.factor.strategy = 'inla'),
                  family            = 'binomial', #'gaussian',
                  num.threads       = 1,
                  Ntrials           = N,
                  verbose           = FALSE,
                  keep              = FALSE,
                  working.directory = temptrash)
  
  summary(res_fit)
  
  #########
  ##  make a prediction surface
  message(' .. prediction')
  draws <- inla.posterior.sample(ndraws, res_fit)
  
  ## get samples as matrices
  par_names <- rownames(draws[[1]]$latent)
  l_idx <- match(res_fit$names.fixed, par_names) # main effects
  s_idx <- grep('^space.*', par_names) # spatial effects
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  rownames(pred_l) <- res_fit$names.fixed
  
  ## if we fit with a nugget, we also need to take draws of the nugget precision
  #if(length(grep('^IID.ID.*', par_names)) > 0){
  #  pred_n <- sapply(draws, function(x) {
  #    nug.idx <- which(grepl('IID.ID', names(draws[[1]]$hyper)))
  #    x$hyperpar[[nug.idx]]}) ## this gets the precision for the nugget
  #}else{
  #  pred_n <- NULL
  #}
  
  coords <- cbind(predfr$lon,predfr$lat)
  
  # Projector matrix
  A.pred <- inla.spde.make.A(mesh = mesh, loc = coords)
  
  predfr$int <- 1
  vals <- predfr[, c('int',fevars), with=FALSE]
  
  # get ndraws obs at each predframe location
  cell_l <- unname(as.matrix(as(data.matrix(vals), "dgeMatrix") %*% pred_l))
  cell_s <- as.matrix(crossprod(t(A.pred), pred_s))
  
  ## add on nugget effect if applicable
  #if(!is.null(pred_n)){
  #  cell_n <- sapply(pred_n, function(x){
  #    rnorm(n = nrow(cell_l), sd = 1 / sqrt(x), mean = 0)
  #  })
  #} else
  
  # make a cell pred object predframe by ndraws draws
  message(' .. summarize and save predictions')
  cell_pred <- plogis(cell_s+cell_l)
    
  # make summaries for plotting
  pred_q     <- rowQuantiles(cell_pred, probs = c(0.025, 0.975))
  inlapredfr <- data.frame(mean=rowMeans(cell_pred),lower=pred_q[,1],upper=pred_q[,2])
  predsumr   <- insertRaster(ext_raster,inlapredfr) 
  
  
  return(list(model_object = res_fit,
              cell_pred    = cell_pred,
              pred_summary = inlapredfr,
              raster       = predsumr))
  
}






# fit a separate child and adult model:
vars <- c('dpt1', 'dpt_dropout', 'type2_sias','mat_edu', 'imp_sani', 'log_u5pop')
type <- 'type2'

ch_inla_mod <- fitinlamodel(data = dagg[adult==0], fevars = vars, outcometype = type)
ad_inla_mod <- fitinlamodel(data = dagg[adult==1], fevars = vars, outcometype = type)




## plot it
hist(ch_inla_mod$pred_summary[,1]) 
plot(ch_inla_mod$raster) #, zlim = c(0.5, 1.0))
        
hist(ad_inla_mod$inlapredfr[,1]) 
plot(ad_inla_mod$raster) #, zlim = c(0.2, 1.0)) #, main = 'Predicted adult seropositive prevalence')


# children minus adult
plot(ch_inla_mod$raster[[1]]-ad_inla_mod$raster[[1]], main = 'Additional immunity in children over adults')
         
diff_cp    <- ch_inla_mod$cell_pred - ad_inla_mod$cell_pred
diffq      <- rowQuantiles(diff_cp, probs = c(0.025, 0.975))
diffpredfr <- data.frame(mean=rowMeans(diff_cp),lower=diffq[,1],upper=diffq[,2])
diff       <- insertRaster(ext_raster,diffpredfr) 

plot(diff[[2:3]]) #, main = 'Predicted adult seropositive prevalence')


## Aggregate these results at the district level
ch <- data.table(ch_inla_mod$cell_pred,adm2=predfr$adm2,pop=predfr$u5pop)
ad <- data.table(ad_inla_mod$cell_pred,adm2=predfr$adm2,pop=predfr$u5pop)
df <- data.table(diff_cp,adm2=predfr$adm2,pop=predfr$u5pop)

cols <- paste0('V',1:100)
ch <- ch[, lapply(.SD, weighted.mean,w=pop,na.rm=TRUE), .SDcols = cols, by=adm2]
ad <- ad[, lapply(.SD, weighted.mean,w=pop,na.rm=TRUE), .SDcols = cols, by=adm2]
df <- df[, lapply(.SD, weighted.mean,w=pop,na.rm=TRUE), .SDcols = cols, by=adm2]


ch <- ch[,  as.list(quantile(.SD,c(0.025,0.5,0.975),na.rm=TRUE)), by=adm2,.SDcols=cols]
ad <- ad[,  as.list(quantile(.SD,c(0.025,0.5,0.975),na.rm=TRUE)), by=adm2,.SDcols=cols]
df <- df[,  as.list(quantile(.SD,c(0.025,0.5,0.975),na.rm=TRUE)), by=adm2,.SDcols=cols]
names(ch) <- c('id','ch_lower','ch_median','ch_upper')
names(ad) <- c('id','ad_lower','ad_median','ad_upper')
names(df) <- c('id','df_lower','df_median','df_upper')

out <- merge(merge(ch, ad, by = 'id', all = TRUE), df, by = 'id', all = TRUE)


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















