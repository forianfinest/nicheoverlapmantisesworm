###SCRIPT LIST###

###TRIMMING OCCURENCES AND MODELING###

https://datadryad.org/stash/dataset/doi:10.5061/dryad.zgmsbcc98 (Chang et al. 2021 Dryad folder, reported as Lin 2020 in the manuscript)

https://www.sciencebase.gov/catalog/item/607df429d34e8564d67e3afc (Claunch & Reed 2021, scripts for modeling)

library(sp)
library(raster)
library(rgdal)
library(dismo)
library(usdm)
library(ENMeval)
library(rJava)
library(spThin)

chordo<-read.csv(file = "chordo.csv", header=TRUE) #####FILE WITH THE COORDINATES

chordot<-thin(
chordo,
lat.col = "y",
long.col = "x",
spec.col = "species",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/enm/coordinates", #####CHOOSE OUTPUT FOLDER
out.base = "thinned_data_chordo",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)


pate<-read.csv("patellifera.csv", header = TRUE)

patet<-thin(
pate,
lat.col = "y",
long.col = "x",
spec.col = "species",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/enm/coordinates",
out.base = "thinned_data_pate",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

tiformo <- read.csv("tiformo.csv", header = TRUE)

tiformot<-thin(
tiformo,
lat.col = "y",
long.col = "x",
spec.col = "species",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/enm/coordinates",
out.base = "thinned_data_formo",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

acro <- read.csv("a_japonica.csv", header = TRUE)

acrot<-thin(
acro,
lat.col = "y",
long.col = "x",
spec.col = "species",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/enm/coordinates",
out.base = "thinned_data_acro",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)


###NOTE: WE UPLOADED AS SUP MAT THE ALREADY TRIMMED FILES###

###VIF CALC###

lista_variables2 <-list.files(pattern='*.asc') ####THIS IS FOR MAKING THE RASTERSTACK WITH THE BIOCLIM VARIABLES, YOU FIRST MAKE A LIST WITH THE FILES' NAMES###
a <- stack(lista_variables2)
coordinates(chordot) <- c("x", "y") ####HERE IT'S WITH CHORDODES, BUT IT WORKS SIMILAR WITH THE OTHER SPECIES. ALSO, HERE IT'S X AND Y, BUT YOU CAN CHANGE ACCORDING TO HOW IT'S CALLED ON YOUR DATASET (I.E. LONG AND LAT) ###
sp<-SpatialPoints(chordot)
b<-extract(a,sp)
c<-cbind(chordot,b)
write.csv(c,file="test.csv")
chordovif<-read.csv("test.csv",header=T)
vifstep(chordovif,th=5)

###FOR THE OTHERS###

coordinates(acrot)<-c("longitude", "latitude")
nuacro<-SpatialPoints(acrot)
nuacrob<-extract(a, nuacro, method="bilinear")
acroc<-cbind(acrot, nuacrob)
nuacroc<-cbind(nuthinnedjaponica, nuacrob)
write.csv(nuacroc,file="nuacrovar.csv")
nuacrovar<-read.csv("nuacrovar.csv",header=T)
vifstep(nuacrovar,th=5)

coordinates(patet)<-c("lon", "lat")
nupate<-SpatialPoints(patet)
nupateb<-extract(a, nupate, method="bilinear")
nupatec<-cbind(patet, nupateb)
write.csv(nupatec,file="nupate_variables.csv")
nupatevar<-read.csv("nupate_variables.csv")
vifstep(nupatevar,th=5)

coordinates(tiformot)<-c("x", "y")
nuformo<-SpatialPoints(tiformot)
nuformob<-extract(a, nuformo, method="bilinear")
nuformoc<-cbind(tiformot, nuformob)
write.csv(nuformo2c,file="nuformo_variables_noNAS.csv")
nuformovar<-read.csv("nuformo_variables_noNAS.csv")
vifstep(nuformovar2[3:21],th=5)

###MAKING MODELS###

acro <- read.csv('acro_points.csv', header = TRUE)
coordinates(acro) <- c("x", "y")
datafilesacro <- Sys.glob("*.asc") ###IN THIS CASE WE CREATED OTHER FOLDERS AND PUT THE VIF-SELECTED VARIABLES INSIDE THERE###
acrostack <- stack(datafilesacro)


acrohinge <- ENMevaluate(occ = acro,
envs = acrostack,
bg = NULL,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
n.bg = 2500,
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

pate <- read.csv('pate_points.csv', header = TRUE)
coordinates(pate) <- c("x", "y")
datafilespate <- Sys.glob("*.asc")
patestack <- stack(datafilespate)


patehinge <- ENMevaluate(occ = pate,
envs = patestack,
bg = NULL,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
n.bg = 2500,
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

tiformo <- read.csv('tiformo_points.csv', header = TRUE)
coordinates(tiformo) <- c("x", "y")
datafilestiformo <- Sys.glob("*.asc")
tiformostack <- stack(datafilestiformo)


tiformohinge <- ENMevaluate(occ = tiformo,
envs = tiformostack,
bg = NULL,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
n.bg = 2500,
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

######IF YOU WANT TO USE OUR BACKGROUND DATASET

acrobg<-read.csv('acro_bg.csv', header = TRUE)
coordinates(acrobg) <- c("x", "y")

acrohinge <- ENMevaluate(occ = acro,
envs = acrostack,
bg = acrobg,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

patebg<-read.csv('pate_bg.csv', header = TRUE)
coordinates(patebg) <- c("x", "y")

patehinge <- ENMevaluate(occ = pate,
envs = patestack,
bg = patebg,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

tiformobg<-read.csv('tiformo_bg.csv', header = TRUE)
coordinates(tiformobg) <- c("x", "y")

tiformohinge <- ENMevaluate(occ = tiformo,
envs = tiformostack,
bg = tiformobg,
fc=c('L','LQ','LP','LQP','LQH','H'),
partitions = "block",
RMvalues =c(1,1.5,2,2.5,3,3.5,4,4.5,5),
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

###TO VISUALIZE THE RESULTS AND CHOOSE THE MODEL####

View(acrohinge@results) 
View(patehinge@results) 
View(tiformohinge@results) 

###AFTER THIS, YOU CAN CHECK FOR WHICH MODEL PERFORMS BETTER ACCORDING TO YOUR METRICS. THE, YOU CAN SAVE IT AS A RASTER FILE###
writeRaster(acrohinge@predictions@layers[[38]], 'acro1p5lqh.asc', format='ascii')
writeRaster(patehing@predictions@layers[[53]], 'pate4p5h.asc', format='ascii')
writeRaster(tiformohinge@predictions@layers[[51]], 'tiformo3p5h.asc', format='ascii')


###MODELING FOR CHORDODES FORMOSANUS###

chordot <- read.csv('chordo_points.csv', header = TRUE)
coordinates(chordot) <- c("x", "y")
chordoasc <- Sys.glob("*.asc")
stack_chordo <- stack(chordoasc)

#####GIVEN THE WAY LOWER NUMBER OF OCCURENCES FOR C. FORMOSANUS, WE ALSO MADE A BIAS FILE FOR SAMPLING THE BACKGROUND, FOLLOWING SCRIPTS AND SUGGESTIONS HERE:

https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/

https://github.com/jamiemkass/ENMeval/issues/26

library(dismo)
library(raster)
library(MASS)
library(magrittr)
library(maptools)

test_raster<-raster("wc2_bio_30s_01.asc") ###OR WHICHEVER RASTER HAS THE SAME RESOLUTION CROPPED FOR TAIWAN, THIS IS ONE FROM LIN 2020####
occur.ras <- rasterize(chordot, test_raster, 1)
presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]
dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)))
myraster<-raster(nrows=480, ncols=360, xmn=119.5, xmx=122.5,ymn=21.5,ymx=25.5,vals=dens$z)
deems <- mask(myraster, test_raster)
deems[is.na(deems[])] <- 0
chordobg <- xyFromCell(deems, sample(ncell(deems), 2500, prob=values(deems)))

####WE CAN NOW USE OUR BIASED BG

chordohinge <- ENMevaluate(occ = chordot,
envs = stack_chordo,
bg = chordobg,
fc=c('L','LQ','LP','QP','LQP','LQH','QH','H','QT','LQPT','LQHPT'),
partitions = "jackknife",
RMvalues =c(0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5),
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

####IF YOU WANT TO USE OUR OWN BG

ogchordo_bg <- read.csv("chordo_bg.csv", header = TRUE)
coordinates(ogchordo_bg)<-c("x", "y")

chordohinge <- ENMevaluate(occ = chordot,
envs = stack_chordo,
bg = ogchordo_bg,
fc=c('L','LQ','LP','QP','LQP','LQH','QH','H','QT','LQPT','LQHPT'),
partitions = "jackknife",
RMvalues =c(0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5),
algorithm = 'maxent.jar',
clamp = TRUE, rasterPreds = NULL,
parallel = FALSE, numCores = 1, progbar = TRUE, updateProgress = TRUE)

####CHEK RESULTS

View(chordohinge@results)

writeRaster(chordosus@predictions@layers[[189]],"chordo2p25lqpt.asc",format='ascii')

###PLOTTING MODELS###

library(colorRamps)
plotstack<-stack(chordohinge@predictions@layers[[189]], acrohinge@predictions@layers[[38]],tiformohinge@predictions@layers[[51]], patehing@predictions@layers[[53]])
spplot(plotstack, col.regions=blue2green2red(200))

####BINARY MAPS

https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/

sdm_threshold <- function(sdm, occs, type = "mtp", binary = FALSE){
  occPredVals <- raster::extract(sdm, occs)
  if(type == "mtp"){
    thresh <- min(na.omit(occPredVals))
  } else if(type == "p10"){
    if(length(occPredVals) < 10){
      p10 <- floor(length(occPredVals) * 0.9)
    } else {
      p10 <- ceiling(length(occPredVals) * 0.9)
    }
    thresh <- rev(sort(occPredVals))[p10]
  }
  sdm_thresh <- sdm
  sdm_thresh[sdm_thresh < thresh] <- NA
  if(binary){
    sdm_thresh[sdm_thresh >= thresh] <- 1
  }
  return(sdm_thresh)
}


###SO, BASICALLY YOU HAVE TO PUT THE MODEL AND THE PRESENCE POINTS USED (BELOW FOR C. formosanus "chordot" ARE THE PRESENCE POINTS)

chordo_p10 <- sdm_threshold(chordohinge@predictions@layers[[189]], chordot, "p10", binary=TRUE)
chordo_p10[is.na(chordo_p10[])] <- 0
writeRaster(chordo_p10,"chordodbinp10.asc",format='ascii')

acro_p10 <- sdm_threshold(acrohinge@predictions@layers[[38]], acro, "p10", binary=TRUE)
acro_p10[is.na(acro_p10[])] <- 0
writeRaster(acro_p10,"acrodbinp10.asc",format='ascii')

tiformo_p10 <- sdm_threshold(tiformohinge@predictions@layers[[38]], tiformo, "p10", binary=TRUE)
tiformo_p10[is.na(tiformo_p10[])] <- 0
writeRaster(tiformo_p10,"tiformodbinp10.asc",format='ascii')

pate_p10 <- sdm_threshold(patehinge@predictions@layers[[38]], pate, "p10", binary=TRUE)
pate_p10[is.na(pate_p10[])] <- 0
writeRaster(pate_p10,"patedbinp10.asc",format='ascii')



###ENMTOOLS####

library(ENMTools)

chordo<-raster("chordo2p25lqpt.asc")
pate<-raster("pate4p5h.asc")
tiformo<-raster("formo2p5lqpt.asc")
acro<-raster("tiformo3p5h.asc")

pateformo<-raster.overlap(pate, formo, verbose = FALSE)
write.csv(pateformo, file = "pateformo.csv") ###FOR VISUALIZATION, THE R OUTPUT OMITS DECIMAL NUMBERS FROM THE 8TH ONWARD### 
chordoformo<-raster.overlap(chordo, formo, verbose = FALSE)
write.csv(chordoformo, file = "chordoformo.csv")
acroformo<-raster.overlap(acro, formo, verbose = FALSE)
write.csv(acroformo, file = "acroformo.csv")
acropate<-raster.overlap(acro, pate, verbose = FALSE)
write.csv(acropate, file = "acropate.csv")
acrochordo<-raster.overlap(acro, chordo, verbose = FALSE)
write.csv(acrochordo, file = "acrochordo.csv")
patechordo<-raster.overlap(pate, chordo, verbose = FALSE)
write.csv(patechordo, file = "patechordo.csv")

pearson_chordoacro<-raster.cor(chordo, acro, method = "pearson")
pearson_chordoformo<-raster.cor(chordo, formo, method = "pearson")
pearson_chordopate<-raster.cor(chordo, pate, method = "pearson")
pearson_pateacro<-raster.cor(pate, acro, method = "pearson")
pearson_formoacro<-raster.cor(formo, acro, method = "pearson")
pearson_formopate<-raster.cor(formo, pate, method = "pearson")

