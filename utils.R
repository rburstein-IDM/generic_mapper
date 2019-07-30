## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Utility functions for generic_mapper
## ###########################################################



#' Insert Raster
#' Helper function stolen from seegMBG. Replace values in a raster with those in a column of a matrix.
#' @param raster Raster with target resolution and dimensions
#' @param new_vals New data matrix, each column representing a layer in the resulting raster brick
#'
#' @return a new raster or rasterBrick with the values supplies in new_vals
insertRaster <- function (raster, new_vals) {
  cellIdx <- function (x) 1:ncell(x)
  idx <- cellIdx(raster)
  stopifnot(length(idx) == nrow(new_vals))
  stopifnot(max(idx) <= ncell(raster))
  n <- ncol(new_vals)
  raster_new <- raster::brick(replicate(n,
                                        raster[[1]],
                                        simplify = FALSE))
  names(raster_new) <- colnames(new_vals)

  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  return (raster_new)
}





clean_covariates <- function(){
  
}







## ###########################################################
## Clean up and prep data

# cut down the columns a bit
d <- d2[, c('hh','cluster','adult','province','Sabin_1','Sabin_2','Sabin_3','lat','lon'), with = FALSE]

# na omit for now (note that in Arie's paper they did some raking for missingness)
d <- na.omit(d)
d <- d[lat!=0]

# aggregate to the cluster level (naive to weight for now)
dagg <- d[, .(type1 = sum(Sabin_1>3),
              type2 = sum(Sabin_2>3),
              type3 = sum(Sabin_3>3),
              N     = .N), 
          by = .(cluster,lat,lon,province,adult)]

# make alternative variables
dpt_dropout <- dpt1 - dpt3
u5pop[u5pop==0] <- 0.01
log_u5pop   <- log(u5pop)

# make a raster version of the DRC shape
ext_raster <- trim(fasterize(shp, dpt1))

# clean up rasters and extract values at data locations
rlist <- c('dpt1','dpt3','dpt_dropout','u5pop','log_u5pop','mat_edu','imp_sani')
csfr  <- data.table() # to track center and scaling
for(r in rlist){
  
  message(sprintf(' .. crop it, mask it, extract it, center it, scale it: %s',r))

  # mask and extract
  assign(r, mask(crop(get(r), ext_raster), ext_raster))
  dagg[[r]] <- raster::extract(get(r), SpatialPoints(cbind(dagg$lon,dagg$lat)))
  
  # center and scale
  tmp  <- data.table(var=r,mean=mean(dagg[[r]], na.rm=T),sd=sd(dagg[[r]], na.rm=T))
  csfr <- rbind(csfr, tmp)
  dagg[[paste0(r,'_cs')]] <- ( dagg[[r]]-tmp$mean[1] ) / tmp$sd[1]
  
}

# rasterize the gadm file
rshp <- fasterize(sf = shp, raster = dpt3[[1]], field = 'id')


# make a prediction frame representing the full raster
predfr <- data.table(idx = 1:ncell(get(rlist[1])), adm2 = as.vector(rshp))
for(r in c(rlist)){ 
  predfr[[r]] <- as.vector(get(r))
  predfr[[paste0(r,'_cs')]] <- (predfr[[r]]-csfr[var==r]$mean) / csfr[var==r]$sd
}

# add lat and lon to prediction frame in case we use it
predfr[, lon := coordinates(ext_raster)[,1]]
predfr[, lat := coordinates(ext_raster)[,2]]


## #########################################################
## TODO: Incorporate SIAs as a spatial covariate.

# bring in sia data
sia <- data.table(df_sia_final)

# keep to the 5 years preceding this survey and DRC only
sia <- sia[ADM0_NAME == 'DEMOCRATIC REPUBLIC OF THE CONGO']
sia[, year := year(ymd(start_date))]
sia <- sia[year >= 2009 & year < 2015]

# keep only 0 to 5 and faction >= .5
sia <- sia[agegroup == '0 to 5 years']
sia <- sia[fraction >= 0.5]

# collapse by ADM_2
sia <- sia[, .(type1_sias = sum(vaccinetype %in% c('bOPV','mOPV1','tOPV')),
               type2_sias = sum(vaccinetype %in% c('tOPV')),
               type3_sias = sum(vaccinetype %in% c('bOPV','tOPV'))),
            by = .(ADM0_NAME,ADM1_NAME,ADM2_NAME)]

# add adm2 id info
shpinfo[, ADM0_NAME := toupper(ADM0_VIZ_NAME)]
shpinfo[, ADM1_NAME := toupper(ADM1_VIZ_NAME)]
shpinfo[, ADM2_NAME := toupper(ADM2_VIZ_NAME)]
sia <- merge(sia, shpinfo,  by = c('ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME'), 
             all.x = TRUE)
sia <- sia[, c('id', 'type1_sias', 'type2_sias', 'type3_sias'), with = FALSE]

# make sia rasters and add them to predfr
predfr <- merge(predfr, sia, by.x = 'adm2', by.y = 'id', all.x = TRUE)
predfr <- predfr[order(idx)]

for(r in c('type1_sias', 'type2_sias', 'type3_sias')){
  assign(r, insertRaster(rshp, cbind(as.vector(predfr[[r]]))))
  assign(r, insertRaster(rshp, cbind(as.vector(predfr[[r]]))))
  assign(r, insertRaster(rshp, cbind(as.vector(predfr[[r]]))))
  dagg[[r]] <- raster::extract(get(r), SpatialPoints(cbind(dagg$lon,dagg$lat)))
}

## ###########################################################
## Super simple sanity check model

summary(glm( cbind(type2,N-type2) ~ dpt1,            data = dagg[adult==0], family='binomial'))
summary(glm( cbind(type2,N-type2) ~ type2_sias,      data = dagg[adult==0], family='binomial'))
summary(glm( cbind(type2,N-type2) ~ dpt1+type2_sias, data = dagg[adult==0], family='binomial'))

# simple logistic regression
m1ch <- glm( cbind(type2,N-type2) ~ dpt1+dpt_dropout+type2_sias+mat_edu+imp_sani+log_u5pop, # + lat*lon,
              data = dagg[adult==0], family='binomial')
m1ad <- glm( cbind(type2,N-type2) ~ dpt1+dpt_dropout+type2_sias+mat_edu+imp_sani+log_u5pop, # + lat*lon,
              data = dagg[adult==1], family='binomial')

summary(m1ch)
summary(m1ad)

# logistic gam
m2ch <- gam( cbind(type2,N-type2) ~ s(dpt1) + s(dpt_dropout) + s(type2_sias,k=3) + s(mat_edu) + s(imp_sani) + s(log_u5pop), # + s(lat) + s(lon), 
             data = dagg[adult==0], family='binomial', predict=TRUE )
m2ad <- gam( cbind(type2,N-type2) ~ s(dpt1) + s(dpt_dropout) + s(type2_sias,k=3) +s(mat_edu) + s(imp_sani) + s(log_u5pop), # + s(lat) + s(lon), 
             data = dagg[adult==1], family='binomial', predict=TRUE )
summary(m2ch) # looks like non-linear terms are more predictive.
summary(m2ad)

# how much explanatory power can we get out of a gbm?
incl <- c('type2','dpt1','dpt_dropout','type2_sias','mat_edu','imp_sani','log_u5pop') #,'lat','lon')
run_brt <- function(ddd)
  gbm.step(data            = ddd[,incl,with=F],
           gbm.y           = 1,
           gbm.x           = incl[-1],
           bag.fraction    = 0.5,
           tree.complexity = 4,
           n.trees         = 120,
           learning.rate   = .001,
           offset          = log(ddd$N),
           family          = 'poisson')
m3ch <- run_brt(ddd = dagg[adult==0])
m3ad <- run_brt(ddd = dagg[adult==1]) 

summary(m3ch)
summary(m3ad); plot.gbm(m3ad,4)

# how did the various predictors do?
predtest_ch <- data.table(type2 = dagg[adult==0]$type2/dagg[adult==0]$N,
                          pred1 = plogis(predict(m1ch,newdata=dagg[adult==0])),
                          pred2 = plogis(predict(m2ch,newdata=dagg[adult==0])),
                          pred3 = predict.gbm(m3ch,n.trees=m3ch$gbm.call$best.trees,type='response'))
predtest_ad <- data.table(type2 = dagg[adult==1]$type2/dagg[adult==1]$N,
                          pred1 = plogis(predict(m1ad,newdata=dagg[adult==1])),
                          pred2 = plogis(predict(m2ad,newdata=dagg[adult==1])),
                          pred3 = predict.gbm(m3ad,n.trees=m3ad$gbm.call$best.trees,type='response'))
cor(na.omit(predtest_ch))
cor(na.omit(predtest_ad))
# the covariates are a bit more predictive of child rather than adult immunity.
# the more complex brt model does about the same for both



# map the the 3 predictive covariate only models
pred_ch <- insertRaster(ext_raster,cbind(plogis(predict(m1ch,newdata=predfr)),
                                         plogis(predict.gam(m2ch,newdata=predfr)),
                                         predict(m3ch,n.trees=m3ch$gbm.call$best.trees,
                                                 type='response',newdata=predfr[,incl[-1],with=FALSE])))
pred_ad <- insertRaster(ext_raster,cbind(plogis(predict(m1ad,newdata=predfr)),
                                         plogis(predict.gam(m2ad,newdata=predfr)),
                                         predict(m3ad,n.trees=m3ad$gbm.call$best.trees,
                                                  type='response',newdata=predfr[,incl[-1],with=FALSE])))
plot(pred_ch[[1]])
plot(pred_ad[[1]])  



## ###########################################################
## Basic geostats model - spatial covariates (enter linearly) with a spatial Matern GP spde
##   start with Sabin.2 in under-5s


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















