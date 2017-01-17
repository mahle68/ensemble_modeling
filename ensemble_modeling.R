
###############################################################################################################################
# title:    "Ensemble modeling for migratory niche of Oriental honey-buzzards (Pernis ptilorhynchus) over the East China Sea" #
# author:   Elham Nourani                                                                                                     #
# location: Nagasaki University, Nagasaki, Japan                                                                              #
# email:    mahle68@gmail.com                                                                                                 #
# date:     17 January 2017                                                                                                   #
# comments:                                                                                                                   #
###############################################################################################################################

#set working directory
setwd("~/ensemble_modeling") 

#create a folder named "input data" within your working directory and save the input presence points and environmental data in it.
#running dir("input data") should return the following:
#[1] "autumn-2km-vars-8.5_present-RCM"         "autumn-2km-vars-future-4.5_46-55-cordex" "autumn-2km-vars-future-4.5_91-00-cordex"
#[4] "autumn-2km-vars-future-8.5_46-55-cordex" "autumn-2km-vars-future-8.5_91-00-cordex" "input_points.csv"  


#load necessary libraries
library(biomod2) # to run niche modleing algorithms and ensemble modeling
library(plyr)    # to convert arrays to dfs
library(raster)  # to open and manipulate spatial data 
library(mgcv)    
library(dismo)   # to generate random background points


################################
## Step1: prepare input files ##
################################

#load species data (presence points)
training<-read.csv("input data/input_points.csv",header=T)

#convert to a spatial file to be used for background selection
sp_training<-training 
coordinates(sp_training)<-~longitude+latitude


#generating background points:

#load masking layer
extent<-raster("input data/autumn-2km-vars-8.5_present-RCM/autumn_ua.asc") # This has to be a raster with no NAs; any of the prepared environmental layers would do

#set layer projection for both the masking layer and presence points

proj4string(extent)<- CRS("+proj=longlat +datum=WGS84")
proj4string(sp_training)<- CRS("+proj=longlat +datum=WGS84")

#generate 5000 random points within the extent, not overlapping with the presence points
background<-data.frame(randomPoints(ext,5000,sp_training))


#define column names for background points
colnames(background)<-c("longitude","latitude")

#put training and background data in a single file.
#add response column to the background (NA) and training points (1) 
background$resp<-rep(NA,nrow(background))
training$resp<-rep(1,nrow(training))

#bind them together
resp<-rbind(training,background)

#convert to a spatial object and define projection
coordinates(resp)<-~longitude+latitude
proj4string(resp)<- CRS("+proj=longlat +datum=WGS84")


#load explanatory variables
Filenames<-list.files("input data/autumn-2km-vars-8.5_present-RCM",pattern=".asc",full.name=T)

#put them in a raster stack and define projection
rstack<-stack(Filenames)
proj4string(rstack)<- CRS("+proj=longlat +datum=WGS84")

#define species name
sp_name<-"pernis_ptilorhynchus"

#initialize the dataset to be used for modeling
myBiomodData<-BIOMOD_FormatingData(resp.var=resp,
                                   expl.var=rstack,
                                   resp.name= sp_name,
                                   PA.nb.rep=3, #repeat pseudrep selection 3 times! This is arbitrary
                                   PA.nb.absences= 4000,
                                   PA.strategy= "random" #randomly picks pseudo-absences from the background points defined as NA in the resp df 
                                   )

############################
## Step 2: model building ##
############################


#configure model options
myBiomodOption <- BIOMOD_ModelingOptions(GLM = list( type = 'quadratic',
                                                     interaction.level = 0,
                                                     myFormula = NULL,
                                                     test = 'AIC',
                                                     family = binomial(link = 'logit'),
                                                     mustart = 0.5,
                                                     control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE) ),
                                         GBM = list( distribution = 'bernoulli',
                                                     n.trees = 2500,
                                                     interaction.depth = 7,
                                                     n.minobsinnode = 5,
                                                     shrinkage = 0.001,
                                                     bag.fraction = 0.5,
                                                     train.fraction = 1,
                                                     cv.folds = 3,
                                                     keep.data = FALSE,
                                                     verbose = FALSE,
                                                     perf.method = 'cv'),
                                         GAM = list( algo = 'GAM_mgcv',
                                                     type = 's_smoother',
                                                     k = -1,
                                                     interaction.level = 0,
                                                     myFormula = NULL,
                                                     family = binomial(link = 'logit'),
                                                     method = 'GCV.Cp',
                                                     optimizer = c('outer','newton'),
                                                     select = FALSE,
                                                     knots = NULL,
                                                     paraPen = NULL,
                                                     control = list(nthreads = 1, irls.reg = 0, epsilon = 1e-07, maxit = 200, trace = FALSE, mgcv.tol = 1e-07, mgcv.half = 15
                                                                    , rank.tol = 1.49011611938477e-08, nlm = list(ndigit=7, gradtol=1e-06, stepmax=2, steptol=1e-04, iterlim=200, check.analyticals=0)
                                                                    , optim = list(factr=1e+07), newton = list(conv.tol=1e-06, maxNstep=5, maxSstep=2, maxHalf=30, use.svd=0), outerPIsteps = 0
                                                                    , idLinksBases = TRUE, scalePenalty = TRUE, keepData = FALSE, scale.est = "pearson") ),
                                         MAXENT.Phillips = list( path_to_maxent.jar = "C:/Users/Elham/Desktop/maxent",
                                                                 memory_allocated = 512,
                                                                 maximumiterations = 1000,
                                                                 visible = FALSE,
                                                                 linear = TRUE,
                                                                 quadratic = TRUE,
                                                                 product = TRUE,
                                                                 threshold = FALSE,
                                                                 hinge = FALSE,
                                                                 lq2lqptthreshold = 80,
                                                                 l2lqthreshold = 10,
                                                                 hingethreshold = NULL,
                                                                 beta_threshold = NULL,
                                                                 beta_categorical = NULL,
                                                                 beta_lqp = -1,
                                                                 beta_hinge = NULL,
                                                                 defaultprevalence = 0.5)
                                        )


# calibrating the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c("GLM","GBM","GAM","MAXENT.Phillips"),
                                    models.options = myBiomodOption,
                                    NbRunEval=10, #number of evaluation runs, i.e. cross-validation folds
                                    DataSplit=80, #% of data to use for calibrating the model; the rest will be used for testing
                                    Prevalence=0.5,#for buidling weighted response weights
                                    VarImport=5, #number of permutations to estimate variable importance
                                    models.eval.meth = c("TSS","ROC"), #evaluation metrics
                                    SaveObj = TRUE, #keep all results and outputs on hard drive
                                    rescal.all.models = TRUE, #all model prediction will be scaled with a binomial GLM
                                    do.full.models = FALSE, # if true, models will be calibrated and evaluated with the whole dataset alone
                                    modeling.id = paste(sp_name,"FirstModeling",sep="") #ID or name of modeling procedure.
                                    )


###############################
### STEP 3: model evaluation ##
###############################

# get all model evaluations
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

#change the array to a dataframe for easier manipulation
Eval.df<-adply(.data=myBiomodModelEval,.margins=c(1:5),.id=c("eval.metric","data.type","model","run","PArep"))

colnames(Eval.df)[6]<-"value"

#create an dataframe including evaluation data
Eval.df<-Eval.df[Eval.df$data.type=="Testing.data",]

#save the evaluation data as a csv file
write.csv(Eval.df,"Eval_df.csv")

#calculate average TSS and ROC for each model and save as csv file
averageEval<-with(Eval.df,tapply(Eval.df$value,list(eval.metric, model),mean))

write.csv(averageEval,"averageEval.csv")


#create separate columns for ROC and TSS
TSS<-Eval.df[Eval.df$eval.metric=="TSS",-1]
colnames(TSS)[5]<-"TSS"
row.names(TSS)<-c(1:120)

ROC<-Eval.df[Eval.df$eval.metric=="ROC",-1]
colnames(ROC)[5]<-"ROC"
row.names(ROC)<-c(1:120)

#put the two data frame into one
Eval.df2<-cbind(ROC,TSS$TSS)
colnames(Eval.df2)[6]<-"TSS"

#extract models with TSS>0.7 and ROC>0.7
GoodModels<-subset(Eval.df2,Eval.df2$TSS>=0.7 & Eval.df2$ROC>=0.8)

#check the summary to make sure everything is good
summary(GoodModels)
table(GoodModels$model) 

#Extract names for the good models
#example name format: "pernis.ptilorhynchus_PA3_RUN10_GLM"
GoodModels$model.name<-paste("pernis.ptilorhynchus",GoodModels$PArep,GoodModels$run,GoodModels$model,sep="_")



#################################
## STEP 4: variable importance ##
#################################

# variable importances 
VarImp<-get_variables_importance(myBiomodModelOut)
dimnames(VarImp)

#convert into a data frame
VarImp.df<-adply(.data=VarImp,.margins=c(1:4),.id=c("variable","model","run","PArep"))
colnames(VarImp.df)[5]<-"value"

#calculate average variable importance across all models for each algorithm and save as a csv file
averageVarImp<-with(VarImp.df,tapply(VarImp.df$value,list(variable, model),mean))
write.csv(averageVarImp,"averageVarImp.csv")

#calculate standar deviations for variable importances across all models for each algorithm and save as a csv file
sdVarImp<-with(VarImp.df,tapply(VarImp.df$value,list(variable, model),sd))
write.csv(sdVarImp,"sdVarImp.csv")



###############################
## STEP 5: ensemble modeling ##
###############################

#ensemble modeling with good models ;)
myBiomodEM <- BIOMOD_EnsembleModeling(
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all', # the way the models will be combined to build the ensemblemodels.
    eval.metric = c("TSS","ROC"), 
    eval.metric.quality.threshold = c(0.7,0.8), #provide a threshhold for all evaluation metrics chosen
    prob.mean = F,
    prob.cv = F,
    prob.ci = F,
    prob.ci.alpha = 0.05,
    prob.median = F,
    committee.averaging = F, #non-weighted average
    prob.mean.weight = T, #weighted mean
    prob.mean.weight.decay = 'proportional', #relative importance of the weights
    VarImport = 5
    )


#########################################
## STEP 6: evaluate the ensemble model ##
#########################################

#get evaluation scores and save as a csv file
myBiomodEMEval<-get_evaluations(myBiomodEM)

write.csv(myBiomodEMEval,"myBiomodEMEval.csv")


#get variable importance
VarImp.Ens<-get_variables_importance(myBiomodEM)

VarImpEns.df<-adply(.data=VarImp.Ens,.margins=c(1:3),.id=c("variable","model","run"))
colnames(VarImpEns.df)[4]<-"value"

#calculate average variable importance and save as a csv file
averageVarImp.Ens<-tapply(VarImpEns.df$value,VarImpEns.df$variable,mean)
write.csv(averageVarImp.Ens,"averageVarImpEns.csv")

#calculate standard deviations for variable importance and save as a csv file
sdVarImp.Ens<-tapply(VarImpEns.df$value,VarImpEns.df$variable,sd)
write.csv(sdVarImp.Ens,"sdVarImpEns.csv")


##################################################
## STEP 7: projecting to the current conditions ##
##################################################


#project only the good models to use for visualizations for each model...
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = rstack, #current conditions
                                  proj.name = 'current',
                                  selected.models = GoodModels$model.name,
                                  compress = FALSE,
                                  build.clamping.mask = FALSE,
                                  output.format = '.grd',
                                  on_0_1000= FALSE)

#separate the model types into raster stacks
GAM.stack<-get_predictions(myBiomodProj,model="GAM")
GLM.stack<-get_predictions(myBiomodProj,model="GLM")
GBM.stack<-get_predictions(myBiomodProj,model="GBM")
MAX.stack<-get_predictions(myBiomodProj,model="MAXENT.Phillips")

#use raster package to avergae predictions for each model type
GAM.average<-overlay(GAM.stack,fun=mean)
GLM.average<-overlay(GLM.stack,fun=mean)
GBM.average<-overlay(GBM.stack,fun=mean)
MAX.average<-overlay(MAX.stack,fun=mean)

windows();plot(GAM.average,main = "GAM")
windows();plot(GLM.average,main = "GLM")
windows();plot(GBM.average,main = "GBM")
windows();plot(MAX.average,main = "MAXENT")

####ensemble model visualization

#project all the individual models on "current" conditions to use later with the ensemble forecast function
myBiomodProj.all <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = rstack, #current conditions
                                      proj.name = 'current',
                                      selected.models = "all",
                                      compress = FALSE,
                                      build.clamping.mask = FALSE,
                                      output.format = '.grd',
                                      on_0_1000= TRUE #0-1 scale
                                      )

#ensemble forecasting for the current conditions
myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj.all
                                         )

#average the forecasted models
myBiomodEF.stack<-get_predictions(myBiomodEF)
myBiomodEF.average<-overlay(myBiomodEF.stack,fun=mean)

#visualize
windows();plot(myBiomodEF.average)

#save as a raster file
writeRaster(myBiomodEF.average,"BiomodEF.average.tif","GTiff")




#################################################
## STEP 8: projecting to the future conditions ##
#################################################


######### PROJECT TO RCP 4.5 ##########


#mid-century conditions (2046-2055):

#Import futre condition files. Filenames must match those of the previous step
Filenames_future_mid_45<-list.files("input data/autumn-2km-vars-future-4.5_46-55-cordex",pattern=".asc",full.name=T)

rstack_future_mid_45<-stack(Filenames_future_mid_45)
proj4string(rstack_future_mid_45)<- CRS("+proj=longlat +datum=WGS84")


myBiomodProjFuture_mid_45 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                               new.env = rstack_future_mid_45,
                                               proj.name = 'future_mid_45',
                                               selected.models = 'all',
                                               binary.meth = 'TSS',
                                               compress = 'xz',
                                               build.clamping.mask = T,
                                               output.format = '.grd'
                                               )

#ensemble forecasting for future conditions
myBiomodEF_Future_mid_45 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                       projection.output = myBiomodProjFuture_mid_45
                                                       )

#plot
myBiomodEF.future.stack_mid_45<-get_predictions(myBiomodEF_Future_mid_45)
windows();plot(myBiomodEF.future.stack_mid_45,main="future_mid_45")

writeRaster(myBiomodEF.future.stack_mid_45,c("future1_mid_45.tiff","future2_mid_45.tiff"),format="GTiff",bylayer=T)

#average the two versions of future together
ave.future_mid_45<-overlay(myBiomodEF.future.stack_mid_45,fun=mean)
writeRaster(ave.future_mid_45,"ave.future_mid.tif","GTiff")



#late century conditions (2091-2100):

#Import futre condition files. Filenames must match those of the previous step
Filenames_future_end_45<-list.files("input data/autumn-2km-vars-future-4.5_91-00-cordex",pattern=".asc",full.name=T)

rstack_future_end_45<-stack(Filenames_future_end_45)
proj4string(rstack_future_end_45)<- CRS("+proj=longlat +datum=WGS84")


myBiomodProjFuture_end_45 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                               new.env = rstack_future_end_45,
                                               proj.name = 'future_end_45',
                                               selected.models = 'all',
                                               binary.meth = 'TSS',
                                               compress = 'xz',
                                               build.clamping.mask = T,
                                               output.format = '.grd'
                                               )

#Ensemble forecasting for future conditions
myBiomodEF_Future_end_45 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                       projection.output = myBiomodProjFuture_end_45
                                                       )

#plot
myBiomodEF.future.stack_end_45<-get_predictions(myBiomodEF_Future_end_45)
windows();plot(myBiomodEF.future.stack_end_45,main="future_end_45")

writeRaster(myBiomodEF.future.stack_end_45,c("future1_end_45.tiff","future2_end_45.tiff"),format="GTiff",bylayer=T)

#average the two versions of future together
ave.future_end_45<-overlay(myBiomodEF.future.stack_end_45,fun=mean)
writeRaster(ave.future_end_45,"ave.future_end_45.tif","GTiff")



######### PROJECT TO RCP 8.5 ##########


#mid century conditions (2046-2055)

#Import futre condition files. Filenames must match those of the previous step
Filenames_future_mid_85<-list.files("input data/autumn-2km-vars-future-8.5_46-55-cordex",pattern=".asc",full.name=T)

rstack_future_mid_85<-stack(Filenames_future_mid_85)
proj4string(rstack_future_mid_85)<- CRS("+proj=longlat +datum=WGS84")


myBiomodProjFuture_mid_85 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                               new.env = rstack_future_mid_85,
                                               proj.name = 'future_mid_85',
                                               selected.models = 'all',
                                               binary.meth = 'TSS',
                                               compress = 'xz',
                                               build.clamping.mask = T,
                                               output.format = '.grd'
                                               )

#Ensemble forecasting for future conditions
myBiomodEF_Future_mid_85 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                       projection.output = myBiomodProjFuture_mid_85
                                                       )

#plot
myBiomodEF.future.stack_mid_85<-get_predictions(myBiomodEF_Future_mid_85)
windows();plot(myBiomodEF.future.stack_mid_85,main="future_mid_85")

writeRaster(myBiomodEF.future.stack_mid_85,c("future1_mid_85.tiff","future2_mid_85.tiff"),format="GTiff",bylayer=T)

#average the two versions of future together

ave.future_mid_85<-overlay(myBiomodEF.future.stack_mid_85,fun=mean)
writeRaster(ave.future_mid_85,"ave.future_mid_85.tif","GTiff")



#late century conditions (2091-2100)

#Import futre condition files. Filenames must match those of the previous step
Filenames_future_end_85<-list.files("input data/autumn-2km-vars-future-8.5_91-00-cordex",pattern=".asc",full.name=T)

rstack_future_end_85<-stack(Filenames_future_end_85)
proj4string(rstack_future_end_85)<- CRS("+proj=longlat +datum=WGS84")


myBiomodProjFuture_end_85 <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                               new.env = rstack_future_end_85,
                                               proj.name = 'future_end_85',
                                               selected.models = 'all',
                                               binary.meth = 'TSS',
                                               compress = 'xz',
                                               build.clamping.mask = T,
                                               output.format = '.grd'
                                               )

#Ensemble forecasting for future conditions
myBiomodEF_Future_end_85 <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                       projection.output = myBiomodProjFuture_end_85
                                                       )

#plot
myBiomodEF.future.stack_end_85<-get_predictions(myBiomodEF_Future_end_85)
windows();plot(myBiomodEF.future.stack_end_85,main="future_end")

writeRaster(myBiomodEF.future.stack_end_85,c("future1_end_85.tiff","future2_end_85.tiff"),format="GTiff",bylayer=T)

#average the two versions of future together
ave.future_end_85<-overlay(myBiomodEF.future.stack_end_85,fun=mean)
writeRaster(ave.future_end_85,"ave.future_end_85.tif","GTiff")

