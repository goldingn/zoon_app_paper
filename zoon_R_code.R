# This is a simplified version of R code in:
# "The zoon R package for reproducible and shareable species distribution
# modelling"

# You can access the code to reproduce the manuscript in full at:
#  https://github.com/zoonproject/zoon_app_paper
# or in archived form at:
#  https://dx.doi.org/10.5281/zenodo.834875

library(zoon)

mosquito1 <- workflow(occurrence = UKAnophelesPlumbeus,
                      covariate  = UKBioclim,
                      process    = Background(n = 500),
                      model      = MaxEnt,
                      output     = InteractiveMap)

mosquito2 <- workflow(occurrence = UKAnophelesPlumbeus,
                      covariate  = UKBioclim,
                      process    = Chain(Background(n = 500),
                                         StandardiseCov),
                      model      = list(MaxEnt,
                                        GBM,
                                        RandomForest),
                      output     = Appify)


FengPapes <- 
  workflow(occurrence = SpOcc('Dasypus novemcinctus',
                              extent = c(-130, -20, -60, 60)),
           covariate = Bioclim(extent = c(-130, -20, -60, 60),
                               layers = c(1:4, 6, 9, 10, 12, 15)),
           process = Chain(Clean, 
                           MESSMask,  
                           Background(n = 10000,
                                      bias = 200), 
                           Crossvalidate(k = 5)),
           model = MaxEnt,
           output = PrintMap(points = FALSE,
                             threshold = 0.05,
                             thresholdmethod = 'falsenegative',
                             xlim = c(-130, -70),
                             ylim = c(20, 50)))

FengPapesUpdate <- 
  ChangeWorkflow(workflow = FengPapes,
                 output = Chain(ResponseCurve(cov = 1),
                                InteractiveMap)) 

spThin <- function (.data, thin = 50) {
  
  # check these are presence-background data
  stopifnot(all(.data$df$type %in% c('presence', 'background')))
  
  # install & load the package
  zoon::GetPackage('spThin')
  
  # get dataframe & index to presence data
  df <- na.omit(.data$df)
  pres_idx <- which(df$type == 'presence')
  
  # prepare presence data subset and apply thinning
  sub_df <- data.frame(LAT = df$latitude[pres_idx],
                       LONG = df$longitude[pres_idx],
                       SPEC = NA)
  th <- thin(loc.data = sub_df,
             thin.par = thin,
             reps = 1,
             locs.thinned.list.return = TRUE,
             write.files = FALSE,
             write.log.file = FALSE)
  
  # get index to rows in sub_df, update the full dataset and return
  pres_keep_idx <- as.numeric(rownames(th[[1]]))
  .data$df <- rbind(df[pres_idx,][pres_keep_idx, ],
                    df[-pres_idx, ])
  return (.data)
}

BuildModule(object = spThin,
            type = 'process',
            title = 'Spatial thinning of Presence-only Data',
            description = paste('Apply the stochastic spatial thinning',
                                'algorithm implemented in the spThin',
                                'package to presence data in a',
                                'presence-background dataset'),
            details = paste('Full details of the algorithm are available in',
                            'the open-access article by Aiello-Lammens',
                            'et al. (2015): dx.doi.org/10.1111/ecog.01132'),
            author = 'zoon Developers',
            email = 'zoonproject@gmail.com',
            paras = list(thin = paste('Thinning parameter - the required',
                                      'minimum distance (in kilometres) between points',
                                      'after applying the thinning procedure')),
            dataType = 'presence-only',
            check = TRUE)

MaxEntComparison <- workflow(
  occurrence = CarolinaWrenPO,
  covariate  = CarolinaWrenRasters,
  process    = Chain(SubsampleOccurrence(500),
                     Background(n = 10000),
                     CarolinaWrenValidation),
  model      = list(MaxEnt(args = 'threshold=false'),
                    MaxNet,
                    MaxNet(regmult = 0.005),
                    MaxNet(features = 'l'),
                    MaxNet(regmult = 0.005, features = 'l'),
                    LogisticRegression),
  output     = Chain(PrintMap(points = FALSE),
                     PerformanceMeasures))
