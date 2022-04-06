# prodx
# download most recent software version
devtools::install_github("geomorphR/geomorph", ref = "Stable", build_vignettes = TRUE)
devtools::install_github("mlcollyer/RRPP")

# load analysis packages
library(here)
library(StereoMorph)
library(geomorph)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(wesanderson)

# read shape data and define number of sLMs
shapes <- readShapes("shapes")
shapesGM <- readland.shapes(shapes, nCurvePts = c(10,3,5,5,3,10))

# read qualitative data
qdata <- read.csv("qdata.morph.csv",
                  header = TRUE,
                  row.names = 1)

## Generalised Procrustes Analysis

# gpa
Y.gpa <- gpagen(shapesGM, print.progress = FALSE)

## plot
plot(Y.gpa)

# dataframe
gdf <- geomorph.data.frame(shape = Y.gpa$coords,
                           size = Y.gpa$Csize,
                           merged = qdata$merged)

# add centroid size to qdata
qdata$csz <- Y.gpa$Csize

## Boxplot (centroid size)

# attributes
csz <- qdata$csz
merged <- qdata$merged

# palette
pal <- wes_palette("Moonrise2", 6, type = "continuous")

# boxplot of Perdiz arrow points by merged
ggplot(qdata, aes(x = merged, y = csz, color = merged)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  labs(x = 'merged', y = 'Centroid Size')

## Principal Components Analysis

# palette
pal <- wes_palette("Moonrise2", 4, type = "continuous")

# pca
pca <- gm.prcomp(Y.gpa$coords)
summary(pca)

# set plot parameters
pch.gps <- c(1:4)[as.factor(qdata$merged)]
col.gps <- pal[as.factor(qdata$merged)]
col.hull <- c("#798E87", "#C27D38", "#CCC591", "#29211F")

## pca plot
pc.plot <- plot(pca,
                asp = 1,
                pch = pch.gps,
                col = col.gps)
shapeHulls(pc.plot,
           groups = qdata$merged,
           group.cols = col.hull)

### Minima/maxima of PC1/2 with warp grids

# plot x/y maxima/minima
## x - minima
mean.shape <- mshape(Y.gpa$coords)
plotRefToTarget(pca$shapes$shapes.comp1$min, 
                mean.shape)

## x - maxima
plotRefToTarget(pca$shapes$shapes.comp1$max, 
                mean.shape)

## y - minima
plotRefToTarget(pca$shapes$shapes.comp2$min, 
                mean.shape)

## y - maxima
plotRefToTarget(pca$shapes$shapes.comp2$max, 
                mean.shape)

## Procrustes ANOVA: Allometry

# allometry
fit.size <- procD.lm(shape ~ size, 
                     data = gdf, 
                     print.progress = FALSE, 
                     iter = 9999)

# allometry
anova(fit.size)

# unique allometries
fit.unique<-procD.lm(shape ~ size * merged, 
                     data = gdf, 
                     print.progress = FALSE, 
                     iter = 9999)

# unique allometries
anova(fit.unique)

# plot
## PredLine (Adams and Nistri 2010)
plotAllometry(fit.unique, 
              size = gdf$size, 
              logsz = TRUE, 
              method = "PredLine", 
              pch = pch.gps, 
              col = col.gps)

## Procrustes ANOVA: Shape and size

# shape
fit.sh.reg <- procD.lm(shape ~ merged,
                       data = gdf,
                       print.progress = FALSE,
                       iter = 9999)

# difference in shape by size class?
anova(fit.sh.reg)

# pairwise comparison of LS means = which differ?
pw.sh.reg <- pairwise(fit.sh.reg,
                      groups = qdata$merged)
summary(pw.sh.reg, 
        confidence = 0.95, 
        test.type = "dist")

# size
fit.sz.reg <- procD.lm(size ~ merged,
                       data = gdf,
                       print.progress = FALSE,
                       iter = 9999)

# size
anova(fit.sz.reg)

# pairwise comparison of LS means = which differ?
pw.sz.reg <- pairwise(fit.sz.reg,
                      groups = qdata$merged)
summary(pw.sz.reg, 
        confidence = 0.95, 
        test.type = "dist")

## Trajectory analysis

# trajectory analysis::shape
TA <- trajectory.analysis(fit.sh.reg,
                          groups = qdata$region, 
                          traj.pts = qdata$size, 
                          print.progress = FALSE)

# magnitude difference
summary(TA,
        attribute = "MD",
        show.trajectories = TRUE)

# plot
TP <- plot(TA, 
           pch = as.numeric(qdata$region), 
           bg = as.numeric(qdata$size))
add.trajectories(TP, 
                 traj.pch = c(15, 17),
                 traj.cex = 1)


# subset landmark coordinates to produce mean shapes
new.coords <- coords.subset(A = Y.gpa$coords,
                            group = qdata$merged)
names(new.coords)

# group shape means
mean <- lapply(new.coords, mshape)

## plot mean shape
plot(mean$north_L)
plot(mean$north_S)
plot(mean$south_L)
plot(mean$south_S)

## comparison plot
plotRefToTarget(mean$south_S,
                mean$south_L,
                method = "TPS",
                mag = 1,
                useRefPts = TRUE)
