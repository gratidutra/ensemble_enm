library(raster)
library(SSDM)
library(tidyverse)
library(rJava)

# CCSM4

setwd("~/Documents/Biologia/paper_fernando/CCSM4/CCSM4 RCP 8.5/RASTER")
r_8.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/CCSM4/CCSM4 RCP 6.0/RASTER")
r_6.0 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/CCSM4/CCSM4 RCP 4.5/RASTER")
r_4.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

s <- stack(r_6.0, r_4.5, r_8.5)

setwd("~/Documents/Biologia/paper_fernando/CCSM4/CCSM4 RCP 4.5")
occurrence <- read_csv("PONTOS DE OCORRÊNCIA.csv")
occurrence$long <- as.double(occurrence$long)
occurrence$lat <- as.double(occurrence$lat)

CCSM4 <- ensemble_modelling("MAXENT", occurrence, s,
  rep = 100,
  cv = "k-fold",
  cv.param = c(5, 10),
  Xcol = "long", Ycol = "lat",
  ensemble.thresh = c(0.75)
)

# MRICGCM3

setwd("~/Documents/Biologia/paper_fernando/MRICGCM3/MRICGCM3 RCP 8.5/RASTER")
r_8.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/MRICGCM3/MRICGCM3 RCP 6.0/RASTER")
r_6.0 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/MRICGCM3/MRICGCM3 RCP 4.5/RASTER")
r_4.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

s <- stack(r_6.0, r_4.5, r_8.5)

MRICGCM3 <- ensemble_modelling("MAXENT", occurrence, s,
  rep = 100,
  cv = "k-fold",
  cv.param = c(5, 10),
  Xcol = "long", Ycol = "lat",
  ensemble.thresh = c(0.75)
)

# HADGEM2
setwd("~/Documents/Biologia/paper_fernando/HADGEM2-ES/RCP 8.5 HADGEM2-ES/RASTER")
r_8.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/HADGEM2-ES/RCP 6.0 HADGEM2-ES/RASTER")
r_6.0 <- raster("Itapotihyla_langsdorffii_avg.asc")

setwd("~/Documents/Biologia/paper_fernando/HADGEM2-ES/RCP 4.5 HADGEM2-ES/RASTER")
r_4.5 <- raster("Itapotihyla_langsdorffii_avg.asc")

s <- stack(r_6.0, r_4.5, r_8.5)

HADGEM2 <- ensemble_modelling("MAXENT", occurrence, s,
  rep = 100,
  cv = "k-fold",
  cv.param = c(5, 10),
  Xcol = "long", Ycol = "lat",
  ensemble.thresh = c(0.75)
)

# Salvando raster de adequabilidade

setwd("~/Documents/Biologia/paper_fernando/ensemble/suitability")

writeRaster(CCSM4@projection, filename = "ensemble_suitability_ccsm4.grd", bandorder = "BIL", overwrite = TRUE)

writeRaster(MRICGCM3@projection, filename = "ensemble_suitability_mricgm3.grd", bandorder = "BIL", overwrite = TRUE)

writeRaster(HADGEM2@projection, filename = "ensemble_suitability_hadgem2.grd", bandorder = "BIL", overwrite = TRUE)

# Salvando raster binarizado

setwd("~/Documents/Biologia/paper_fernando/ensemble/binary")

writeRaster(CCSM4@binary, filename = "ensemble_binary_ccsm4.grd", bandorder = "BIL", overwrite = TRUE)

writeRaster(MRICGCM3@binary, filename = "ensemble_binary_mricgm3.grd", bandorder = "BIL", overwrite = TRUE)

writeRaster(HADGEM2@binary, filename = "ensemble_binary_hadgem2.grd", bandorder = "BIL", overwrite = TRUE)

# salvando parametros de cada ensemble
setwd("~/Documents/Biologia/paper_fernando/ensemble")

evaluation_ccsm4 <- write_csv(CCSM4@evaluation, 'ccsm4_evaluation.csv')
evaluation_mricgm3 <- write_csv(MRICGCM3@evaluation, 'mricgm3_evaluation.csv')
evaluation_hadgem2 <- write_csv(HADGEM2@evaluation, 'hadgem2_evaluation.csv')

rm(HADGEM2)

#evaluation <- tibble(
#  model = c("CCSM4", "MRICGCM3", "HADGEM2"),
#  treshold = c(CCSM4@evaluation$threshold, MRICGCM3@evaluation$threshold, HADGEM2@evaluation$threshold),
#  AUC = c(CCSM4@evaluation$AUC, MRICGCM3@evaluation$AUC, HADGEM2@evaluation$AUC),
#  omission_rate = c(CCSM4@evaluation$omission.rate, MRICGCM3@evaluation$omission.rate, HADGEM2@evaluation$omission.rate),
#  sensitivity = c(CCSM4@evaluation$sensitivity, MRICGCM3@evaluation$sensitivity, HADGEM2@evaluation$sensitivity),
#  specificity = c(CCSM4@evaluation$specificity, MRICGCM3@evaluation$specificity, HADGEM2@evaluation$specificity),
#  prop_correct = c(CCSM4@evaluation$prop.correct, MRICGCM3@evaluation$prop.correct, HADGEM2@evaluation$prop.correct),
#  Kappa = c(CCSM4@evaluation$Kappa, MRICGCM3@evaluation$Kappa, HADGEM2@evaluation$Kappa)
#)

#write_csv(evaluation, "evaluation.csv")
