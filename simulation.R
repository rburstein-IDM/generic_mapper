########################################################################################
# Author:   Roy Burstein
# Date:     August 2019
# Purpose:  For presentation on 21 Aug, simulate and model geostatistical data
########################################################################################


#########  
### SETUP

# where images will get saved
plotdir  <- 'C:/Users/rburstein/Desktop/tmpimgs' 


# simulation paameters
set.seed(12334) 
sims            <- 100   # how many spatial points to simulate
xlim <- ylim    <- c(0,10) # extent of spatial domain  
res             <- 100     # number of pixels per row/column
mean_samplesize <- 100     # mean N for each binomial sample
posterior_draws <- 100    # how many draws to take from the posterior 



# load (and if needed, install) libraries
for(l in c('RandomFields', 'raster', 'data.table', 'ggplot2', 'matrixStats',
           'lhs', 'INLA', 'gridExtra', 'magick','geostatsp')){
  if(!l %in% installed.packages()[,1])
    install.packages(l)
  library(l, character.only = TRUE)
}

# define a helper function to cleanly put table into a raster
insert_raster <- function (raster, new_vals) {
  values(raster)    <- 1
  idx               <- 1:ncell(raster)
  n                 <- ncol(new_vals)
  raster_new        <- raster::brick(replicate(n, raster, simplify = FALSE))
  names(raster_new) <- colnames(new_vals)
  for(i in 1:n) {
    raster_new[[i]][idx] <- new_vals[, i]
  }
  return(raster_new)
}



###############   
#### SIMULATION

# make spatial domain represented as a raster
r <- raster(nrows = res, ncols = res, xmn = xlim[1],ymn = ylim[1], xmx = xlim[2], ymx = ylim[2])

# Simulate an underlying gaussian field
rf    <- RFsimulate(c(var=0.5, range=5, shape=2), x = r, n = 1)
rf_df <- as.data.frame(rf, xy = TRUE) 

# simulate some point locations (s) using latin hypercube sampling
s <- randomLHS(sims, 2) * max(xlim)

# simulate binomial values for those points, get everything into a data.table
d <- data.table(lon = s[,1], lat = s[,2], gp_true = extract(rf, s)) # extract true rf
d[, N      := rpois(.N, mean_samplesize)]  # get a sample size for each point
d[, p_true := plogis( gp_true)]            # true probability surface for comparison
d[, y_obs  := rbinom(.N, N, p_true)]       # simulate the observed value

# plot some of the data
# simulated points
ggplot(d, aes(lon, lat)) + geom_point(aes(size = N), color = 'red') +
  geom_point(aes(size = y_obs)) + theme_minimal()

# simulated points versus true underlying probability
ggplot(d, aes(p_true, y_obs/N)) + geom_point(aes(size = N)) +
  theme_minimal() + geom_abline(intercept=0, slope=1, color='red', lty='dashed')

# true underlying probability surface
ggplot() + geom_raster(data = rf_df, aes(x,y,fill=plogis(sim)))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + labs(fill = 'Truth') +
  geom_point(data= d, aes(x=lon, y=lat, size = y_obs/N), shape = 1) +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank()) 



###############   
#### FIT MODEL

# make mesh and plot it
mesh <- inla.mesh.2d(loc = s, max.edge = c(1,2))
plot(mesh); points(s, col = 'blue', cex = d$N/mean_samplesize)

# make A matrix (aligns mesh and data indices)
A <- inla.spde.make.A(mesh, s)

# make inla spatial model objects
spde  <- inla.spde2.matern(mesh)
space <- inla.spde.make.index("space", n.spde = spde$n.spde)

# make a data stack
stack.obs <- inla.stack(data    = list(Y = d$y_obs),
                        A       = list(A, 1),
                        effects = list(space, as.data.frame(d)))

# fit the model
fit <- inla(formula           = as.formula( 'Y ~ -1 + f(space, model = spde)' ), 
            data              = inla.stack.data(stack.obs),
            control.predictor = list(A       = inla.stack.A(stack.obs)),
            control.compute   = list(config  = TRUE),
            family            = 'binomial',
            Ntrials           = d$N,
            verbose           = FALSE,
            keep              = FALSE)

summary(fit)


#####################
#### MAKE PREDICTIONS

# Use R-INLA function to pull posterior samples
draws     <- inla.posterior.sample(posterior_draws, fit)

# Get draw estimates at mesh vertices
s_idx <- grep('^space.*', rownames(draws[[1]]$latent)) # spatial effects
pred_s <- sapply(draws, function (x) x$latent[s_idx])

# lets look at the mode (replace the first draw with it)
pred_s[,1] <- fit$summary.random$`space`$mode

# Extract output coordinates (prediction frame across the full spatial domain)
coords <- coordinates(r)

# Projector matrix from location of REs to prediction frame coordinates
A.pred <- inla.spde.make.A(mesh = mesh, loc = coords)

# output draws for each pixel
cell_pred <- plogis(as.matrix(A.pred %*% pred_s))

# make summaries for plotting
pred_q     <- rowQuantiles(cell_pred, probs = c(0.025, 0.975))
inlapredfr <- data.frame(mean  = rowMeans(cell_pred),
                         lower = pred_q[,1],
                         upper = pred_q[,2], 
                         x     = coords[,1], 
                         y     = coords[,2])
predsumr   <- insert_raster(r,inlapredfr[,1:3]) 

# get the modal draw another way
ld       <- sapply(draws, function (x) x$logdens$joint )
modedraw <- which(ld == min(ld))

#####################
#### PLOT PREDICTIONS

# plot summary results and truth rasters
g1 <- ggplot() + geom_raster(data = rf_df, aes(x,y,fill=plogis(sim)))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + labs(fill = 'Truth') + 
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank()) 
g2 <- ggplot() + geom_raster(data = inlapredfr, aes(x,y,fill=mean))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) +  ggtitle('Mean') 
g3 <- ggplot() + geom_raster(data = inlapredfr, aes(x,y,fill=lower))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) +  ggtitle('Lower')
g4 <- ggplot() + geom_raster(data = inlapredfr, aes(x,y,fill=upper))  +
  scale_fill_viridis_c(option = "inferno", limits = c(0,1)) +  ggtitle('Upper')

gopt <- theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.background = element_blank()) 

grid.arrange(g1, g2+gopt, g3+gopt, g4+gopt, layout_matrix=rbind(c(1,1,1),c(1,1,1),c(2,3,4)))
grid.arrange(g2+gopt+geom_point(data=d, aes(x=lon,y=lat),color='white',size=.25), g3+gopt, g4+gopt, ncol = 3)

# Plot one vertical slice
out <- data.table()
for(px in 0:(res-1)){
  tmp <- data.table(value = cell_pred[px*res+40,], px = px, draw = 1:posterior_draws)
  tmp$truth <- rf_df[px*res+40,]$sim
  out <- rbind(out,tmp)
}
out <- merge(out, data.table(ld=ld,draw=1:posterior_draws), by = 'draw')
out_mean <- out[, .(value = mean(value)),         by = px]
out_true <- out[, .(value = plogis(mean(truth))), by = px]

ggplot(out, aes(y = value, x = px/res)) + geom_line(aes(group = draw), alpha=0.05) +
  scale_color_viridis_c(option = "inferno", limits = c(0,1)) + ylim(0.25,0.85) +
  theme_bw() + ylab('Posterior Draw Value') + xlab('Lat') +
  geom_line(data=out_mean, lwd = 2, aes( color = value)) + #theme(legend.position = 'none') +
  geom_line(data=out_true, lwd = 4,  color = 'blue', alpha = 0.3) +
  geom_line(data=out[ld%in%min(ld)], color='green') + 
  geom_line(data=out[draw==1],       color='green',lty='dotted')


# plot a few map draws and animate them
drawsr <- data.table(cbind(cell_pred[,1:20]), x= coords[,1], y= coords[,2])
files <- c()
for(i in 1:20){
  f <- sprintf('%s/%s.png',plotdir,i)
  files <- c(files, f)
  png(f)
  tmp <- drawsr[, c(paste0('V',i),'x','y'), with = FALSE]
  names(tmp) <-  c('value','x','y')
  plot(ggplot() + geom_raster(data = tmp, aes(x,y,fill=value))  +
    scale_fill_viridis_c(option = "inferno", limits = c(0,1)) + 
    theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.background = element_blank()))
  dev.off()
}

image_write(image_animate(image_read(files),fps=2), sprintf("%s/animation.gif",plotdir))

g2 + gopt +  geom_point(data=d, aes(x=lat,y=lon))


# posterior for range and variance (using SD which I am fairly sure is what RF uses)
outp.field     <- inla.spde2.result(fit, name = 'space', spde = spde, do.transform = TRUE)
range_postdens <- data.table(inla.tmarginal(function(x) sqrt(8)/x, outp.field$marginals.kappa[[1]]))
var_postdens   <- data.table(inla.tmarginal(function(x) x, outp.field$marginals.variance.nominal[[1]]))

g1 <- ggplot(range_postdens, aes(x,y)) + theme_bw() + geom_line(size=2) + xlab('Matern Range') +
  ylab('Marginal Posterior Density') 
g2 <- ggplot(var_postdens, aes(x,y))   + theme_bw() + geom_line(size=2) + xlab('Nominal Variance') + 
  ylab('Marginal Posterior Density')
grid.arrange(g1,g2)

# joint hyperparameter density
joint <- data.table(fit$joint.hyper)
names(joint) <- c('theta1','theta2','dens')
ggplot() + geom_point(data=joint, aes(x=theta1,y=theta2,color=dens), size = 9)



# lets pretend we have 4 even areas we care to aggregate
adm <- c(rep(rep(c(1,2),each=res/2),res/2),rep(rep(c(3,4),each=res/2),res/2))
plot(insert_raster(r,cbind(adm))) # confirm its 4 squares
cp <- data.table(adm=adm,pop=1,cell_pred)
cols <- paste0('V',1:posterior_draws)
cp  <- cp[, lapply(.SD, weighted.mean,w=pop,na.rm=TRUE), .SDcols = cols, by=adm]
cp <- cp[,  as.list(quantile(.SD,c(0.5,0.025,0.975),na.rm=TRUE)), by=adm, .SDcols=cols]
tru <- data.table(truth=rf_df$sim,adm=adm)[, .(truth = plogis(mean(truth))), by = adm]
cp <- merge(cp, tru, by = 'adm')
names(cp) <- c('admin','median','lower','upper','truth')

ggplot(data=cp, aes(x=factor(admin))) + geom_point(aes(y=truth), size=5, color='blue', shape=1) +
  geom_point(aes(y=median)) + ylim(.35,.75) +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.1) + theme_bw() + xlab('Admin Area') + ylab('')

