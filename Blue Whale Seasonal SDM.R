remove(list = ls())

# Load the required packages
library(biomod2)
library(raster)
library(sf)
library(ggplot2)
library(sp)
library(ggcorrplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

######CURRENT#####
#biological data
blue_seasonal <- st_read("Biological Data/Seasonal Whale Sightings Data/Blue Whale/Blue_Whale_Seasonal_Presence/Blue_Whale_Seasonal_Presence.shp")
blue_seasonal2 <- as.data.frame(blue_seasonal) #turn to dataframe

#envrionmental data
env <- read.csv("Environmental Data/10km Resolution/env_relevant.csv")
env2 <- as.data.frame(env)

# Extract the coordinates of all background data, and the environmental variables of interest
expl.var.coords <- env2[,c(1,2)]
expl.var <- env2[,c(3,5,6,7,18,19)]
#expl.var <- env2[,c(5,7)]

#extract data for model formatting
resp.xy <- blue_seasonal2[,2:3] #extract the long and lat for the sightings
myrespVar <- as.numeric(blue_seasonal2[,7]) #extract the presences for the sightings

#Find the rows with at least one NA in the environmental
rows_with_NAs = which(rowSums(is.na(expl.var)) > 0)

# Remove all rows with NAs, whether or not they have a whale presence
resp.xy = resp.xy[-rows_with_NAs,]
myrespVar = myrespVar[-rows_with_NAs]

#Remove rows which have both a whale presence and at least one NA in the environmental data
expl.var = expl.var[-rows_with_NAs,]
expl.var.coords = expl.var.coords[-rows_with_NAs,]

#set up model parameters
blue_seasonal_options <- BIOMOD_ModelingOptions(GLM = list(type = "quadratic", interaction.level = 1), RF = list(ntree = 1000),
                                       # GAM = list(
                                       #                 algo = "GAM_mgcv",
                                       #                 type = "s_smoother",
                                       #                 k = 25,
                                       #                 interaction.level = 0,
                                       #                 myFormula = NULL,
                                       #                 family = binomial(link = "logit"),
                                       #                 method = "GCV.Cp",
                                       #                 optimizer = c("outer","newton"),
                                       #                 select = FALSE,
                                       #                 knots = NULL,
                                       #                 paraPen = NULL,
                                       #                 control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
                                       #                                , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                       #                                , rank.tol = 1.49011611938477e-08
                                       #                                , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                       #                                , optim = list(factr=1e+07)
                                       #                                , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                       #                                , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE
                                       #                                , efs.lspmax = 15, efs.tol = 0.1, keepData = FALSE, scale.est = "fletcher"
                                       #                                , edge.correct = FALSE)),
                                       MAXENT.Phillips  =  list(
                                         path_to_maxent.jar = "~/Desktop/FOME Github/FOME/Hannah's SDM/maxent/maxent.jar",
                                         memory_allocated =8192,
                                         maximumiterations = 200,
                                         visible = FALSE,
                                         linear = TRUE,
                                         quadratic = TRUE,
                                         product = TRUE,
                                         threshold = TRUE,
                                         hinge = TRUE,
                                         lq2lqptthreshold = 80,
                                         l2lqthreshold = 10,
                                         hingethreshold = 15,
                                         beta_threshold = -1,
                                         beta_categorical = -1,
                                         beta_lqp = -1,
                                         beta_hinge = -1,
                                         defaultprevalence = 0.5))
blue_seasonal_options

myBiomodData <- BIOMOD_FormatingData(resp.var = myrespVar,
                                     expl.var = expl.var,
                                     resp.xy = resp.xy,
                                     resp.name = "Blue_Seasonal", PA.nb.rep = 1,
                                     PA.nb.absences = 10000, PA.strategy = "random")
myBiomodData
plot(myBiomodData)

# ### Crossvalidation ###
# #DataSplitTable.b <- BIOMOD_CrossValidation(bm.format = myBiomodData,
#                                            k = 5,
#                                            nb.rep = 2,
#                                            do.full.models = FALSE)


#running the actual SDM
blue_seasonal_model <- BIOMOD_Modeling(
  bm.format = myBiomodData, models = c("GLM", "RF",
                                       #"GAM"
                                       "MAXENT.Phillips.2"), bm.options = blue_seasonal_options,
                                      nb.rep = 5, data.split.perc = 70, var.import = 3, do.full.models = FALSE,
                                      modeling.id = "Blue Seasonal SDM", metric.eval = c('TSS','ROC'))

blue_seasonal_model

# Get evaluation scores & variables importance
#get_evaluations(blue_seasonal_model)
aa = get_evaluations(blue_seasonal_model)
mean_roc_glm <- mean(aa["ROC", "Testing.data", "GLM", , ])
mean_roc_rf <- mean(aa["ROC", "Testing.data", "RF", , ])
mean_roc_maxent <- mean(aa["ROC", "Testing.data", "MAXENT.Phillips.2", , ])

# Print the mean ROC values
cat("Mean ROC for GLM:", mean_roc_glm, "\n")
cat("Mean ROC for RF:", mean_roc_rf, "\n")
cat("Mean ROC for Maxent:", mean_roc_maxent, "\n")

#TSS
mean_tss_glm <- mean(aa["TSS", "Testing.data", "GLM", , ])
mean_tss_rf <- mean(aa["TSS", "Testing.data", "RF", , ])
mean_tss_maxent <- mean(aa["TSS", "Testing.data", "MAXENT.Phillips.2", , ])

#Print the mean TSS values
cat("Mean TSS for GLM:", mean_tss_glm, "\n")
cat("Mean TSS for RF:", mean_tss_rf, "\n")
cat("Mean TSS for Maxent:", mean_tss_maxent, "\n")

#get variables of importance
get_variables_importance(blue_seasonal_model, as.data.frame = TRUE)

#Plot evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_seasonal_model)
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model, group.by = c('expl.var', 'algo', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model, group.by = c('algo', 'expl.var', 'dataset'))

#Plot response curves
bm_PlotResponseCurves(bm.out = blue_seasonal_model,
                      models.chosen = get_built_models(blue_seasonal_model)[c(1:3, 12:14)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_seasonal_model,
                      models.chosen = get_built_models(blue_seasonal_model)[c(1:3, 12:14)],
                      fixed.var = 'min')
# #bm_PlotResponseCurves(bm.out = blue_seasonal_model,
#                       models.chosen = get_built_models(blue_seasonal_model)[3],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

#project the single models
myBiomodProj <- BIOMOD_Projection(bm.mod = blue_seasonal_model,
                                  proj.name = 'Current',
                                  new.env = expl.var,
                                  new.env.xy = expl.var.coords,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = FALSE)


plot(myBiomodProj)

## Run the ensemble models
# Model ensemble models
blue_model_seasonal_ensemble <- BIOMOD_EnsembleModeling(bm.mod = blue_seasonal_model,
                                               models.chosen = blue_seasonal_model@models.computed,
                                               em.by = 'PA_dataset',
                                               metric.select = c('TSS'),
                                               metric.select.thresh = c(0.7),
                                               var.import = 3,
                                               metric.eval = c('TSS', 'ROC'),
                                               prob.mean = TRUE,
                                               prob.median = TRUE,
                                               prob.cv = FALSE,
                                               prob.ci = TRUE,
                                               prob.ci.alpha = 0.05,
                                               committee.averaging = TRUE,
                                               prob.mean.weight = TRUE,
                                               prob.mean.weight.decay = 'proportional')

blue_model_seasonal_ensemble

# Get evaluation scores & variables importance
bb <- get_evaluations(blue_model_seasonal_ensemble, as.data.frame = TRUE)
ensemble_data <- bb[bb$Model == "EMmeanByTSS", ]
# Calculate the mean ROC value for the ensemble model
mean_roc_ensemble <- mean(ensemble_data$Testing.data[ensemble_data$Eval.metric == "ROC"])
# Print the mean ROC value for the ensemble model
cat("Mean ROC for Ensemble Model:", mean_roc_ensemble, "\n")
# Calculate the mean TSS value for the ensemble model
mean_tss_ensemble <- mean(ensemble_data$Testing.data[ensemble_data$Eval.metric == "TSS"])
# Print the mean TSS value for the ensemble model
cat("Mean TSS for Ensemble Model:", mean_tss_ensemble, "\n")

#variables of importance for ensemble
rankings_blue_seasonal <- get_variables_importance(blue_model_seasonal_ensemble, as.data.frame = TRUE)

# Calculate average MDA values for each expl variable
averageMDA <- tapply(rankings_blue_seasonal$Var.imp, rankings_blue_seasonal$Expl.var, mean)

# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_model_seasonal_ensemble, group.by = 'model')
bm_PlotEvalBoxplot(bm.out = blue_model_seasonal_ensemble, group.by = c('model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_model_seasonal_ensemble, group.by = c('expl.var', 'model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_model_seasonal_ensemble, group.by = c('expl.var', 'model', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_model_seasonal_ensemble, group.by = c('model', 'expl.var', 'dataset'))

# Represent response curves
bm_PlotResponseCurves(bm.out = blue_model_seasonal_ensemble,
                      models.chosen = get_built_models(blue_model_seasonal_ensemble)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_model_seasonal_ensemble,
                      models.chosen = get_built_models(blue_model_seasonal_ensemble)[c(1, 6, 7)],
                      fixed.var = 'min')
# bm_PlotResponseCurves(bm.out = blue_model_seasonal_ensemble,
#                       models.chosen = get_built_models(blue_model_seasonal_ensemble)[7],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

# Project ensemble models (from single projections)
myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = blue_model_seasonal_ensemble,
                                             proj.name = "ensemble",
                                             models.chosen = 'all',
                                             bm.proj = myBiomodProj,
                                             metric.binary = 'all',
                                             metric.filter = 'all')
plot(myBiomodEMProj)


##REPROJECT TO STUDY AREA##
#crop the environmental variables and coordinates to study area to project the ensemble there
lat_min <- 41
lat_max <- 56
lon_min <- -71
lon_max <- -46
cropped_data <- subset(env2, latitude >= lat_min & latitude <= lat_max & longitude >= lon_min & longitude <= lon_max)

# #plotting to see if this makes sense
# world <- ne_countries(scale = "medium", returnclass = "sf") %>%
#   st_transform(crs = 4326)
# ggplot() +
#   #theme_minimal()+ change the theme here
#   geom_point(aes(color = cropped_data$SST, x = cropped_data$longitude, y = cropped_data$latitude)) +
#   geom_sf(data = blue[blue$presence == 1, ], colour = "red", size = 0.7) +
#   geom_sf(data = world, fill = "black", colour = "white", size = 0.2) +
#   coord_sf(xlim = c(-75,-41), ylim = c(39,61)) +
#   xlab("Longitude") + # for the x axis label
#   ylab("Latitude")

cropped.env <- cropped_data[,c(3,5,6,7,18,19)]
#cropped.env <- cropped_data[,c(5,7,19)]
cropped.env.xy <- cropped_data[,c(1,2)]

#remove NAs
rows_with_NAs_Cropped = which(rowSums(is.na(cropped.env)) > 0)

#Remove rows which have both a whale presence and at least one NA in the environmental data
cropped.env = cropped.env[-rows_with_NAs_Cropped,]
cropped.env.xy = cropped.env.xy[-rows_with_NAs_Cropped,]

# Project ensemble models (from single projections)
myBiomodProjEnsemble <- BIOMOD_Projection(bm.mod = blue_seasonal_model,
                                          proj.name = 'Ensemble',
                                          new.env = cropped.env,
                                          new.env.xy = cropped.env.xy,
                                          models.chosen = 'all',
                                          metric.binary = 'all',
                                          metric.filter = 'all',
                                          build.clamping.mask = FALSE)

myBiomodEmReProj <- BIOMOD_EnsembleForecasting(bm.em = blue_model_seasonal_ensemble,
                                               proj.name = "ensemble",
                                               bm.proj = myBiomodProjEnsemble,
                                               models.chosen = 'all',
                                               metric.binary = 'all',
                                               metric.filter = 'all')
plot(myBiomodEmReProj)


#extract habitat suitability values and make data frame
load("~/Desktop/FOME GitHub/FOME/Hannah's SDM/Data/Blue.Seasonal/proj_ensemble/proj_ensemble_Blue.Seasonal_ensemble.RData")
blue_seasonal_hs_results_ensemble = data.frame(cropped.env.xy, HSI = ef.out[,1])
write.csv(blue_seasonal_hs_results_ensemble, file = "blue_seasonal_hs_results_ensemble.csv", row.names = FALSE)
#write.csv(blue_seasonal_hs_results_ensemble, file = "blue_seasonal_hs_results_ensemble_test8.csv", row.names = FALSE)

# #plotting to see if this makes sense
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326)
ggplot() +
  #theme_minimal()+ change the theme here
  geom_point(aes(color = blue_seasonal_hs_results_ensemble$HSI, x = blue_seasonal_hs_results_ensemble$longitude, y = blue_seasonal_hs_results_ensemble$latitude)) +
  scale_color_distiller(palette = "Spectral") +
  geom_sf(data = world, fill = "black", colour = "white", size = 0.2) +
  coord_sf(xlim = c(-71,-46), ylim = c(41,56)) +
  xlab("Longitude") + # for the x axis label
  ylab("Latitude")

write.csv(cropped.env, file = "cropped.env.csv", row.names = FALSE)
#######MID########

#biological data
blue_seasonal <- st_read("Biological Data/Seasonal Whale Sightings Data/Blue Whale/Blue_Whale_Seasonal_Presence/Blue_Whale_Seasonal_Presence.shp")
blue_seasonal2 <- as.data.frame(blue_seasonal)

#environmental mid data
env_mid <- st_read("Environmental Data/10km Resolution/TenKm_AOI_Mid_final/TenKm_AOI_Mid_final.shp")
env_mid2 <- as.data.frame(env_mid)
env_mid2 <- env_mid2 %>% mutate_at(vars(Bathy_Mean, SST, SSS, NPP), as.numeric) #make sure all columns are numeric

# Extract the coordinates of all background data, and the environmental variables of interest
expl.var.coords <- env_mid2[,c(6,7)]
expl.var <- env_mid2[,c(8:13)]

#extract data for model formatting
resp.xy <- blue_seasonal2[,2:3] #extract the long and lat for the sightings
myrespVar <- as.numeric(blue_seasonal2[,7]) #extract the presences for the sightings

#Find the rows with at least one NA in the environmental
rows_with_NAs = which(rowSums(is.na(expl.var)) > 0)

#Remove all rows with NAs, whether or not they have a whale presence
resp.xy = resp.xy[-rows_with_NAs,]
myrespVar = myrespVar[-rows_with_NAs]

#Remove rows which have both a whale presence and at least one NA in the environmental data
expl.var = expl.var[-rows_with_NAs,]
expl.var.coords = expl.var.coords[-rows_with_NAs,]

# #set up the SDM
blue_seasonal_options <- BIOMOD_ModelingOptions(GLM = list(type = "quadratic", interaction.level = 1), RF = list(ntree = 1000),
                                       # GAM = list(
                                       #                 algo = "GAM_mgcv",
                                       #                 type = "s_smoother",
                                       #                 k = 25,
                                       #                 interaction.level = 0,
                                       #                 myFormula = NULL,
                                       #                 family = binomial(link = "logit"),
                                       #                 method = "GCV.Cp",
                                       #                 optimizer = c("outer","newton"),
                                       #                 select = FALSE,
                                       #                 knots = NULL,
                                       #                 paraPen = NULL,
                                       #                 control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
                                       #                                , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                       #                                , rank.tol = 1.49011611938477e-08
                                       #                                , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                       #                                , optim = list(factr=1e+07)
                                       #                                , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                       #                                , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE
                                       #                                , efs.lspmax = 15, efs.tol = 0.1, keepData = FALSE, scale.est = "fletcher"
                                       #                                , edge.correct = FALSE)),
                                       MAXENT.Phillips  =  list(
                                         path_to_maxent.jar = "~/Desktop/FOME Github/FOME/Hannah's SDM/maxent/maxent.jar",
                                         memory_allocated =8192,
                                         maximumiterations = 200,
                                         visible = FALSE,
                                         linear = TRUE,
                                         quadratic = TRUE,
                                         product = TRUE,
                                         threshold = TRUE,
                                         hinge = TRUE,
                                         lq2lqptthreshold = 80,
                                         l2lqthreshold = 10,
                                         hingethreshold = 15,
                                         beta_threshold = -1,
                                         beta_categorical = -1,
                                         beta_lqp = -1,
                                         beta_hinge = -1,
                                         defaultprevalence = 0.5))
blue_seasonal_options

myBiomodData <- BIOMOD_FormatingData(resp.var = myrespVar,
                                     expl.var = expl.var,
                                     resp.xy = resp.xy,
                                     resp.name = "Blue_Seasonal_MID", PA.nb.rep = 1,
                                     PA.nb.absences = 10000, PA.strategy = "random")
myBiomodData
plot(myBiomodData)

### Crossvalidation ###
# DataSplitTable.b <- BIOMOD_CrossValidation(bm.format = myBiomodData,
#                                            k = 5,
#                                            nb.rep = 2,
#                                            do.full.models = FALSE)

#running the actual SDM
blue_seasonal_model_mid <- BIOMOD_Modeling(
  bm.format = myBiomodData, models = c("GLM", "RF",
                                       #"GAM",
                                       "MAXENT.Phillips.2"), bm.options = blue_seasonal_options, nb.rep = 5, data.split.perc = 80,
  var.import = 3, do.full.models = FALSE, modeling.id = "Blue Seasonal MID SDM", metric.eval = c('TSS','ROC'))

blue_seasonal_model_mid

# Get evaluation scores & variables importance
get_evaluations(blue_seasonal_model_mid)
get_variables_importance(blue_seasonal_model_mid, as.data.frame = TRUE)

#Plot evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_seasonal_model_mid)
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_mid, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_mid, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid, group.by = c('expl.var', 'algo', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid, group.by = c('algo', 'expl.var', 'dataset'))

#Plot response curves
bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid,
                      models.chosen = get_built_models(blue_seasonal_model_mid)[c(1:3, 12:14)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid,
                      models.chosen = get_built_models(blue_seasonal_model_mid)[c(1:3, 12:14)],
                      fixed.var = 'min')
# #bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid,
#                       models.chosen = get_built_models(blue_seasonal_model_mid)[3],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

#project the single models
myBiomodProj_Mid <- BIOMOD_Projection(bm.mod = blue_seasonal_model_mid,
                                      proj.name = 'Mid',
                                      new.env = expl.var,
                                      new.env.xy = expl.var.coords,
                                      models.chosen = 'all',
                                      metric.binary = 'all',
                                      metric.filter = 'all',
                                      build.clamping.mask = FALSE)


plot(myBiomodProj_Mid)

## Run the ensemble models
# Model ensemble models
blue_seasonal_model_mid_ensemble <- BIOMOD_EnsembleModeling(bm.mod = blue_seasonal_model_mid,
                                                   models.chosen = blue_seasonal_model_mid@models.computed,
                                                   em.by = 'PA_dataset',
                                                   metric.select = c('TSS'),
                                                   metric.select.thresh = c(0.7),
                                                   var.import = 3,
                                                   metric.eval = c('TSS', 'ROC'),
                                                   prob.mean = TRUE,
                                                   prob.median = TRUE,
                                                   prob.cv = FALSE,
                                                   prob.ci = TRUE,
                                                   prob.ci.alpha = 0.05,
                                                   committee.averaging = TRUE,
                                                   prob.mean.weight = TRUE,
                                                   prob.mean.weight.decay = 'proportional')

blue_seasonal_model_mid_ensemble

# Get evaluation scores & variables importance
get_evaluations(blue_seasonal_model_mid_ensemble, as.data.frame = TRUE)
get_variables_importance(blue_seasonal_model_mid_ensemble, as.data.frame = TRUE)


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_seasonal_model_mid_ensemble, group.by = 'model')
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_mid_ensemble, group.by = c('model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid_ensemble, group.by = c('expl.var', 'model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid_ensemble, group.by = c('expl.var', 'model', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_mid_ensemble, group.by = c('model', 'expl.var', 'dataset'))

# Represent response curves
bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid_ensemble,
                      models.chosen = get_built_models(blue_seasonal_model_mid_ensemble)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid_ensemble,
                      models.chosen = get_built_models(blue_seasonal_model_mid_ensemble)[c(1, 6, 7)],
                      fixed.var = 'min')
# bm_PlotResponseCurves(bm.out = blue_seasonal_model_mid_ensemble,
#                       models.chosen = get_built_models(blue_seasonal_model_mid_ensemble)[7],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

# Project ensemble models (from single projections)
myBiomodEMProj_Mid <- BIOMOD_EnsembleForecasting(bm.em = blue_seasonal_model_mid_ensemble,
                                                 proj.name = "ensemble_mid",
                                                 models.chosen = 'all',
                                                 bm.proj = myBiomodProj_Mid,
                                                 metric.binary = 'all',
                                                 metric.filter = 'all')

##REPROJECT TO STUDY AREA##
#crop the environmental variables and coordinates to study area to project the ensemble there
lat_min <- 41
lat_max <- 56
lon_min <- -71
lon_max <- -46
cropped_data <- subset(env_mid2, latitude >= lat_min & latitude <= lat_max & longitude >= lon_min & longitude <= lon_max)
cropped.env <- cropped_data[,c(8:13)]
cropped.env.xy <- cropped_data[,c(6,7)]

#remove NAs
rows_with_NAs_Cropped = which(rowSums(is.na(cropped.env)) > 0)

#Remove rows which have both a whale presence and at least one NA in the environmental data
cropped.env = cropped.env[-rows_with_NAs_Cropped,]
cropped.env.xy = cropped.env.xy[-rows_with_NAs_Cropped,]

# Project ensemble models (from single projections)
myBiomodProjEnsemble_Mid <- BIOMOD_Projection(bm.mod = blue_seasonal_model_mid,
                                              proj.name = 'Ensemble_MID',
                                              new.env = cropped.env,
                                              new.env.xy = cropped.env.xy,
                                              models.chosen = 'all',
                                              metric.binary = 'all',
                                              metric.filter = 'all',
                                              build.clamping.mask = FALSE)

myBiomodEmReProj_Mid <- BIOMOD_EnsembleForecasting(bm.em = blue_seasonal_model_mid_ensemble,
                                                   proj.name = "ensemble_mid",
                                                   bm.proj = myBiomodProjEnsemble_Mid,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')
plot(myBiomodEmReProj_Mid)

#extract habitat suitability values and make data frame
load("~/Desktop/FOME GitHub/FOME/Hannah's SDM/Data/Blue.Seasonal.MID/proj_ensemble_MID/proj_ensemble_mid_Blue.Seasonal.MID_ensemble.RData")
blue_seasonal_hs_results_ensemble_mid = data.frame(cropped.env.xy, HSI = ef.out[,1])
write.csv(blue_seasonal_hs_results_ensemble_mid, file = "blue_seasonal_hs_results_ensemble_mid.csv", row.names = FALSE)

# #plotting to see if this makes sense
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326)
ggplot() +
  #theme_minimal()+ change the theme here
  geom_point(aes(color = blue_seasonal_hs_results_ensemble_mid$HSI, x = blue_seasonal_hs_results_ensemble_mid$longitude, y = blue_seasonal_hs_results_ensemble_mid$latitude)) +
  scale_color_distiller(palette = "Spectral") +
  geom_sf(data = world, fill = "black", colour = "white", size = 0.2) +
  coord_sf(xlim = c(-71,-46), ylim = c(41,56)) +
  xlab("Longitude") + # for the x axis label
  ylab("Latitude")

#######FUTURE########
#biological data
blue_seasonal <- st_read("Biological Data/Seasonal Whale Sightings Data/Blue Whale/Blue_Whale_Seasonal_Presence/Blue_Whale_Seasonal_Presence.shp")
blue_seasonal2 <- as.data.frame(blue_seasonal)

#environmental mid data
env_fut <- st_read("Environmental Data/10km Resolution/TenKm_AOI_Fut/TenKm_AOI_Fut.shp")
env_fut2 <- as.data.frame(env_fut)
env_fut2 <- env_fut2 %>% mutate_at(vars(Bathy_Mean, SST, SSS, NPP), as.numeric) #make sure all columns are numeric

# Extract the coordinates of all background data, and the environmental variables of interest
expl.var.coords <- env_fut2[,c(6,7)]
expl.var <- env_fut2[,c(8:13)]

#extract data for model formatting
resp.xy <- blue_seasonal2[,2:3] #extract the long and lat for the sightings
myrespVar <- as.numeric(blue_seasonal2[,7]) #extract the presences for the sightings

#Find the rows with at least one NA in the environmental
rows_with_NAs = which(rowSums(is.na(expl.var)) > 0)

#Remove all rows with NAs, whether or not they have a whale presence
resp.xy = resp.xy[-rows_with_NAs,]
myrespVar = myrespVar[-rows_with_NAs]

#Remove rows which have both a whale presence and at least one NA in the environmental data
expl.var = expl.var[-rows_with_NAs,]
expl.var.coords = expl.var.coords[-rows_with_NAs,]

# #set up the SDM
blue_seasonal_options <- BIOMOD_ModelingOptions(GLM = list(type = "quadratic", interaction.level = 1), RF = list(ntree = 1000),
                                       # GAM = list(
                                       #                 algo = "GAM_mgcv",
                                       #                 type = "s_smoother",
                                       #                 k = 25,
                                       #                 interaction.level = 0,
                                       #                 myFormula = NULL,
                                       #                 family = binomial(link = "logit"),
                                       #                 method = "GCV.Cp",
                                       #                 optimizer = c("outer","newton"),
                                       #                 select = FALSE,
                                       #                 knots = NULL,
                                       #                 paraPen = NULL,
                                       #                 control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07
                                       #                                , maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                       #                                , rank.tol = 1.49011611938477e-08
                                       #                                , nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                       #                                , optim = list(factr=1e+07)
                                       #                                , newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0)
                                       #                                , outerPIsteps = 0, idLinksBases = TRUE, scalePenalty = TRUE
                                       #                                , efs.lspmax = 15, efs.tol = 0.1, keepData = FALSE, scale.est = "fletcher"
                                       #                                , edge.correct = FALSE)),
                                       MAXENT.Phillips  =  list(
                                         path_to_maxent.jar = "~/Desktop/FOME Github/FOME/Hannah's SDM/maxent/maxent.jar",
                                         memory_allocated =8192,
                                         maximumiterations = 200,
                                         visible = FALSE,
                                         linear = TRUE,
                                         quadratic = TRUE,
                                         product = TRUE,
                                         threshold = TRUE,
                                         hinge = TRUE,
                                         lq2lqptthreshold = 80,
                                         l2lqthreshold = 10,
                                         hingethreshold = 15,
                                         beta_threshold = -1,
                                         beta_categorical = -1,
                                         beta_lqp = -1,
                                         beta_hinge = -1,
                                         defaultprevalence = 0.5))
blue_seasonal_options

myBiomodData <- BIOMOD_FormatingData(resp.var = myrespVar,
                                     expl.var = expl.var,
                                     resp.xy = resp.xy,
                                     resp.name = "Blue_Seasonal_FUT", PA.nb.rep = 1,
                                     PA.nb.absences = 10000, PA.strategy = "random")
myBiomodData
plot(myBiomodData)

### Crossvalidation ###
# DataSplitTable.b <- BIOMOD_CrossValidation(bm.format = myBiomodData,
#                                            k = 5,
#                                            nb.rep = 2,
#                                            do.full.models = FALSE)

#running the actual SDM
blue_seasonal_model_fut <- BIOMOD_Modeling(
  bm.format = myBiomodData, models = c("GLM", "RF",
                                       #"GAM",
                                       "MAXENT.Phillips.2"), bm.options = blue_seasonal_options, nb.rep = 5, data.split.perc = 80,
  var.import = 3, do.full.models = FALSE, modeling.id = "Blue Seasonal FUT SDM", metric.eval = c('TSS','ROC'))

blue_seasonal_model_fut

# Get evaluation scores & variables importance
get_evaluations(blue_seasonal_model_fut)
get_variables_importance(blue_seasonal_model_fut, as.data.frame = TRUE)

#Plot evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_seasonal_model_fut)
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_fut, group.by = c('algo', 'algo'))
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_fut, group.by = c('algo', 'run'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut, group.by = c('expl.var', 'algo', 'algo'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut, group.by = c('expl.var', 'algo', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut, group.by = c('algo', 'expl.var', 'dataset'))

#Plot response curves
bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut,
                      models.chosen = get_built_models(blue_seasonal_model_fut)[c(1:3, 12:14)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut,
                      models.chosen = get_built_models(blue_seasonal_model_fut)[c(1:3, 12:14)],
                      fixed.var = 'min')
# #bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut,
#                       models.chosen = get_built_models(blue_seasonal_model_fut)[3],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

#project the single models
myBiomodProj_Fut <- BIOMOD_Projection(bm.mod = blue_seasonal_model_fut,
                                      proj.name = 'Fut',
                                      new.env = expl.var,
                                      new.env.xy = expl.var.coords,
                                      models.chosen = 'all',
                                      metric.binary = 'all',
                                      metric.filter = 'all',
                                      build.clamping.mask = FALSE)


plot(myBiomodProj_Fut)

## Run the ensemble models
# Model ensemble models
blue_seasonal_model_fut_ensemble <- BIOMOD_EnsembleModeling(bm.mod = blue_seasonal_model_fut,
                                                   models.chosen = blue_seasonal_model_fut@models.computed,
                                                   em.by = 'PA_dataset',
                                                   metric.select = c('TSS'),
                                                   metric.select.thresh = c(0.7),
                                                   var.import = 3,
                                                   metric.eval = c('TSS', 'ROC'),
                                                   prob.mean = TRUE,
                                                   prob.median = TRUE,
                                                   prob.cv = FALSE,
                                                   prob.ci = TRUE,
                                                   prob.ci.alpha = 0.05,
                                                   committee.averaging = TRUE,
                                                   prob.mean.weight = TRUE,
                                                   prob.mean.weight.decay = 'proportional')

blue_seasonal_model_fut_ensemble

# Get evaluation scores & variables importance
get_evaluations(blue_seasonal_model_fut_ensemble, as.data.frame = TRUE)
get_variables_importance(blue_seasonal_model_fut_ensemble, as.data.frame = TRUE)


# Represent evaluation scores & variables importance
bm_PlotEvalMean(bm.out = blue_seasonal_model_fut_ensemble)
bm_PlotEvalBoxplot(bm.out = blue_seasonal_model_fut_ensemble, group.by = c('model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut_ensemble, group.by = c('expl.var', 'model', 'model'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut_ensemble, group.by = c('expl.var', 'model', 'dataset'))
bm_PlotVarImpBoxplot(bm.out = blue_seasonal_model_fut_ensemble, group.by = c('model', 'expl.var', 'dataset'))

# Represent response curves
bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut_ensemble,
                      models.chosen = get_built_models(blue_seasonal_model_fut_ensemble)[c(1, 6, 7)],
                      fixed.var = 'median')
bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut_ensemble,
                      models.chosen = get_built_models(blue_seasonal_model_fut_ensemble)[c(1, 6, 7)],
                      fixed.var = 'min')
# bm_PlotResponseCurves(bm.out = blue_seasonal_model_fut_ensemble,
#                       models.chosen = get_built_models(blue_seasonal_model_fut_ensemble)[7],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)

# Project ensemble models (from single projections)
myBiomodEMProj_Fut <- BIOMOD_EnsembleForecasting(bm.em = blue_seasonal_model_fut_ensemble,
                                                 proj.name = "ensemble_fut",
                                                 models.chosen = 'all',
                                                 bm.proj = myBiomodProj_Fut,
                                                 metric.binary = 'all',
                                                 metric.filter = 'all')

##REPROJECT TO STUDY AREA##
#crop the environmental variables and coordinates to study area to project the ensemble there
lat_min <- 41
lat_max <- 56
lon_min <- -71
lon_max <- -46
cropped_data <- subset(env_fut2, latitude >= lat_min & latitude <= lat_max & longitude >= lon_min & longitude <= lon_max)
cropped.env <- cropped_data[,c(8:13)]
cropped.env.xy <- cropped_data[,c(6,7)]

#remove NAs
rows_with_NAs_Cropped = which(rowSums(is.na(cropped.env)) > 0)

#Remove rows which have both a whale presence and at least one NA in the environmental data
cropped.env = cropped.env[-rows_with_NAs_Cropped,]
cropped.env.xy = cropped.env.xy[-rows_with_NAs_Cropped,]

# Project ensemble models (from single projections)
myBiomodProjEnsemble_Fut <- BIOMOD_Projection(bm.mod = blue_seasonal_model_fut,
                                              proj.name = 'Ensemble_FUT',
                                              new.env = cropped.env,
                                              new.env.xy = cropped.env.xy,
                                              models.chosen = 'all',
                                              metric.binary = 'all',
                                              metric.filter = 'all',
                                              build.clamping.mask = FALSE)

myBiomodEmReProj_Fut <- BIOMOD_EnsembleForecasting(bm.em = blue_seasonal_model_fut_ensemble,
                                                   proj.name = "ensemble_fut",
                                                   bm.proj = myBiomodProjEnsemble_Fut,
                                                   models.chosen = 'all',
                                                   metric.binary = 'all',
                                                   metric.filter = 'all')
plot(myBiomodEmReProj_Fut)

#extract habitat suitability values and make data frame
load("~/Desktop/FOME GitHub/FOME/Hannah's SDM/Data/Blue.Seasonal.FUT/proj_ensemble_FUT/proj_ensemble_fut_Blue.Seasonal.FUT_ensemble.RData")
blue_seasonal_hs_results_ensemble_fut = data.frame(cropped.env.xy, HSI = ef.out[,1])
write.csv(blue_seasonal_hs_results_ensemble_fut, file = "blue_seasonal_hs_results_ensemble_fut.csv", row.names = FALSE)

# #plotting to see if this makes sense
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs = 4326)
ggplot() +
  #theme_minimal()+ change the theme here
  geom_point(aes(color = blue_seasonal_hs_results_ensemble_fut$HSI, x = blue_seasonal_hs_results_ensemble_fut$longitude, y = blue_seasonal_hs_results_ensemble_fut$latitude)) +
  scale_color_distiller(palette = "Spectral") +
  geom_sf(data = world, fill = "black", colour = "white", size = 0.2) +
  coord_sf(xlim = c(-71,-46), ylim = c(41,56)) +
  xlab("Longitude") + # for the x axis label
  ylab("Latitude")
