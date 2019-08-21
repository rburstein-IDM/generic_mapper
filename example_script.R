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

# model family (binomial or gaussian only)
family <- 'binomial'

# number of drawws from the posterior to take. 
# 100 good for testing, 1000 good for final run (but slower and more RAM-intensive)
ndraws <- 100

# outcome variable name
outcome_varname <- 'cookfuel_solid'

# where your local copy of the generic_mapper repo is (clone from here: https://github.com/rburstein-IDM/generic_mapper)
codepath  <- 'C:/Users/rburstein/OneDrive - IDMOD/code/generic_mapper'

# shapefile name to use (make sure its in './Dropbox (IDM)/generic_mapper/data/shp')
shpfname <- 'alhasan-pk-adm-2-with-clusters.shp'

# shapefile name to use (make sure its in './Dropbox (IDM)/generic_mapper/data/pop_rasters')
popfname <- 'wp-adj-pak-2015-5km.tif'

# name of dataset file (make sure its in './Dropbox (IDM)/generic_mapper/data/prepped_input')
# dataset name must follow the standard naming: <cntry>_<svyyr>_<outcome_varname>.csv
# datasets should be prepared at the individual or HH level, with the following variables:
# LATNUM, LONGNUM, CLUSTER, <outcome_varname> (if the outcome is binomial it must be 0, 1 at this level)
cntry <- 'pak'
svyyr <- 2017


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
  library(l, character.only = TRUE, quietly = TRUE)
}

## locations - relative to user, but in AMUG Dropbox. Make sure you have access to AMUG to use
user      <- Sys.info()[['user']]
root      <- sprintf('C:/Users/%s/Dropbox (IDM)/generic_mapper',user)

# Source the utility functions 
source(sprintf('%s/utils.R', codepath))

## Load in data and shapes
# prepped data
d       <- fread(sprintf('%s/data/prepped_input/%s-%s-%s.csv', root, cntry, svyyr, outcome_varname))

# load shapefile
shp <- st_read(sprintf('%s/data/shp/%s', root, shpfname))
shp$id <- 1:nrow(shp)

# load population raster
pop <- raster(sprintf('%s/data/pop_rasters/%s', root, popfname))


## ###########################################################
## Clean up and prep data


# rasterize the gadm file
ext_raster <- fasterize(sf = shp, raster = pop, field = 'id')


# make a prediction frame representing the full raster
predfr <- data.table(idx = 1:ncell(pop), pop = as.vector(pop), collapse_to = as.vector(ext_raster))

# add lat and lon to prediction frame
predfr[, lon := coordinates(ext_raster)[,1]]
predfr[, lat := coordinates(ext_raster)[,2]]


# drop any clusters that say they are at zero zero
d <- d[LATNUM != 0 & LONGNUM != 0]

# collapse the dhs data to the cluster level
if(tolower(family) == 'binomial'){
  dagg <- d[, .(outcome = sum(get(outcome_varname), na.rm = TRUE),
                N       = sum(!is.na(get(outcome_varname)))), 
            by = .(LATNUM, LONGNUM, CLUSTER)]
  dagg[, tomap := outcome/N]
} else if(tolower(family) == 'gaussian'){
  dagg <- d[, .(outcome = mean(get(outcome_varname), na.rm = TRUE),
                N       = sum(!is.na(get(outcome_varname)))), 
            by = .(LATNUM, LONGNUM, CLUSTER)]
  dagg[, tomap := outcome]
} else {
  stop('Family must be binomial or gaussian')
}

# make some quick plots showing data
hist(dagg$tomap)


ggplot(dagg) + geom_sf(data = shp, fill = 'white', colour = 'grey') + theme_minimal() + coord_sf(datum = NA) + 
  geom_point(aes(LONGNUM,LATNUM,size=N,color=tomap),alpha=0.75, stroke = 0, shape = 16) + 
  scale_color_gradientn(values = c(0,0.8,1.0), colours = c("#a6172d","#EFDC05","#4f953b") ) +
  ylab('') + xlab('')


## ###########################################################
## Basic geostats model using INLA (currently with most settings (such as priors) at default)

# get together all the various objects needed for inla to fit the model
input_obj <- prep_inla_objects(data = dagg, outcome_varname = 'outcome', ss_varname = 'N', plot_mesh = TRUE)

# fit inla model
res_fit <- run_inla_model(input = input_obj, model_family = family)
summary(res_fit)
  
# do prediction from the fitted model
pred <- inla_predict(fitted = res_fit, input = input_obj, ndraws = ndraws, predframe = predfr, ext_raster = ext_raster)


# plot the output raster
plot(pred$raster, col= c("#a6172d","#EFDC05","#4f953b"), ncol=3)

