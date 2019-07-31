## ###########################################################
## Author:  Roy Burstein (rburstein@idmod.org)
## Date:    July 2019
## Purpose: Utility functions for generic_mapper
## ###########################################################

message(' <<< LOADING UTILITIES >>> ')

#' Insert Raster
#' Helper function stolen from seegMBG. Replace values in a raster with those in a column of a matrix.
#' @param raster Raster with target resolution and dimensions
#' @param new_vals New data matrix, each column representing a layer in the resulting raster brick
#'
#' @return a new raster or rasterBrick with the values supplies in new_vals
insert_raster <- function (raster, new_vals) {
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





#' Prep all data objects for an INLA geostat model run
#' Note: Version 1: no covariates, no time, just a GP smoothing model
#'
#' @param data data frame where outcome and ss, and coordinates are (later covariates as well)
#'    make sure xy coordinates are called 'LONGNUM' and 'LATNUM'
#' @param outcome_varname String. Name of outcome variable in data 
#' @param ss_varname String. Name of sample size variable in data
#' @param plot_mesh boolean. make a mesh plot, just for show
#'
#' @return list obkect of everything we need to run the inla spatial model
prep_inla_objects <- function(data, outcome_varname, ss_varname, plot_mesh = TRUE){
  
  
  # make a clean design matrix
  data          <- setDT(data)[,c(outcome_varname,ss_varname,'LONGNUM','LATNUM'), with = FALSE]
  message( sprintf('Original dataset has %i rows.', nrow(data)))
  data          <- na.omit(data)
  message( sprintf('NA-removed dataset has %i rows.', nrow(data)))
  outcome       <- data[[outcome_varname]]
  N             <- data[[ss_varname]]     
  locs          <- cbind(data$LONGNUM, data$LATNUM)
  design_matrix <- data.frame(int = rep(1, nrow(data)))
  
  # make a mesh
  # TODO: Allow user to change these args
  mesh <- inla.mesh.2d(loc      = locs,
                       max.edge = c(.5,5),
                       offset   = 1,
                       cutoff   = 0.5)
  
  if(plot_mesh == TRUE){
    plot(mesh)
    points(locs, col = 'red')
  }
  
  
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
  #formla <- as.formula( sprintf ('o ~ -1 + int + %s + f(space, model = spde)', paste(fevars,collapse='+')) )
  formla <- as.formula( 'o ~ -1 + int + f(space, model = spde)' )
  
  # output a list object of everything we need to pass into inla
  return( list( outcome = outcome,
                spde    = spde,
                space   = space,
                N       = N,
                dm      = design_matrix,
                mesh    = mesh,
                stack   = stack.obs,
                formula = formla) )
  
}






#' Wrapper for inla
#'
#' @param input object which is output from prep_inla_objects()
#' @param model_family currently only support for binomial or gaussian
#'
#' @return
#' @export Fitted inla model object
run_inla_model <- function(input, model_family){
  spde  <- input[['spde']]
  space <- input[['space']]
  if(!tolower(model_family) %in% c('binomial','gaussian'))
    stop('Currently only support for binomial or gaussian models')
  
  if(tolower(model_family) == 'binomial'){
    message('Fitting BINOMIAL model')
    fit <- inla(formula           = input[['formula']],
                data              = inla.stack.data(input_obj[['stack']]),
                control.predictor = list(A       = inla.stack.A(input_obj[['stack']]),
                                         link    = 1,
                                         compute = FALSE),
                control.compute   = list(dic     = TRUE,
                                         cpo     = TRUE,
                                         config  = TRUE),
                control.fixed     = list(expand.factor.strategy = 'inla'),
                family            = 'binomial',
                num.threads       = 1,
                Ntrials           = input_obj[['N']],
                verbose           = TRUE,
                keep              = FALSE)
  } else if(tolower(model_family) == 'gaussian'){
    message('Fitting GAUSSIAN model')
    fit <- inla(formula           = input[['formula']],
                data              = inla.stack.data(input_obj[['stack']]),
                control.predictor = list(A       = inla.stack.A(input_obj[['stack']]),
                                         compute = FALSE),
                control.compute   = list(dic     = TRUE,
                                         cpo     = TRUE,
                                         config  = TRUE),
                control.fixed     = list(expand.factor.strategy = 'inla'),
                family            = 'Gaussian',
                num.threads       = 1,
                scale             = input_obj[['N']],
                verbose           = TRUE,
                keep              = FALSE)
  }
  
  return(fit)
}





#' Predict surface out of inla model
#'
#' @param fitted fitted model object
#' @param input data inputs that went into inla mode
#' @param ndraws number of posterior draws to take
#' @param predframe prediction frame
#' @param ext_raster prediction raster matching prediction frame
#'
#' @return list of three output objects. cell_pred: pixel draws, prediction summary, raster formation prediction summary
inla_predict <- function(fitted, input, ndraws, predframe, ext_raster){
  

  ##  get posterior draws for each covariate
  draws <- inla.posterior.sample(ndraws, fitted)
  
  ## get samples as matrices
  par_names <- rownames(draws[[1]]$latent)
  l_idx <- match(fitted$names.fixed, par_names) # main effects
  s_idx <- grep('^space.*', par_names) # spatial effects
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  
  if(length(fitted$names.fixed) == 1) pred_l <- t(pred_l)
  rownames(pred_l) <- fitted$names.fixed
  
  # output coordinates
  coords <- cbind(predframe$lon,predframe$lat)
  
  # Projector matrix from location of REs to prediction frame coordinates
  A.pred <- inla.spde.make.A(mesh = input[['mesh']], loc = coords)
  
  predframe$int <- 1
  vals <- predframe$int #predfr[, c('int',fevars), with=FALSE]
  
  # get ndraws obs at each predframe location
  cell_l <- unname(as.matrix(as(data.matrix(vals), "dgeMatrix") %*% pred_l))
  cell_s <- as.matrix(crossprod(t(A.pred), pred_s))
  
  # make a cell pred object predframe by ndraws draws
  message(' .. summarize and save predictions')
  cell_pred <- plogis(cell_s+cell_l)
  
  # make summaries for plotting
  pred_q     <- rowQuantiles(cell_pred, probs = c(0.025, 0.975))
  inlapredfr <- data.frame(mean=rowMeans(cell_pred),lower=pred_q[,1],upper=pred_q[,2])
  predsumr   <- insert_raster(ext_raster,inlapredfr) 
  predsumr   <- mask(predsumr, ext_raster)
  
  return(list(cell_pred    = cell_pred,
              pred_summary = inlapredfr,
              raster       = predsumr))
  
}











#' Aggregate draws to the adminX level
#'
#' @param pred inla_predict
#' @param predframe prediction frame
#' @param collapse_to String. Name of column in predfr to collapse to
#'
#' @return DT summarizing median, 2.5% and 97.5% quantiles
aggregate_results <- function(pred, predframe){
  
  # make a prediction frame
  tmp <- data.table(pred$cell_pred, adm = predframe[['collapse_to']], pop = predframe$pop)
  
  # define draw columns
  cols <- paste0('V',1:ncol(pred$cell_pred))
  
  # collapse each draw, population weighted, to the admin 2
  tmp  <- tmp[, lapply(.SD, weighted.mean,w=pop,na.rm=TRUE), .SDcols = cols, by=adm]
  
  # summarize across draws
  tmp <- tmp[,  as.list(quantile(.SD,c(0.5,0.025,0.975),na.rm=TRUE)), by=adm, .SDcols=cols]
  
  # rename things
  names(tmp) <- c('id','median','lower','upper')
  
  return(tmp)
  
}











