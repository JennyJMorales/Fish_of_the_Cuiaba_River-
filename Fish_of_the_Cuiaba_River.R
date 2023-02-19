# Ecomorphology as a predictor of fish trophic guilds in streams of a Brazilian tropical savanna

# Jenny J. Morales1, Lucia Mateus2, & Jerry Penha2

# 1 Programa de Pós-graduação em Ecologia e Conservação da Biodiversidade, Instituto de Biociências, Universidade Federal de Mato Grosso, Cuiabá, Brasil.
# 2 Laboratório de Ecologia e Manejo de Recursos Pesqueiros, Instituto de Biociências, Universidade Federal de Mato Grosso, Cuiabá, Brasil.
#----------------------------------------------------------------------------

library(vegan)
library(adespatial)  
library(car)
library(ggplot2)
library(MASS)
library (rrcov)
library (ellipse)


##-------REPOSITOTY: Fish_of_the_Cuiabá_River

#####----RESOURCES:

#Alimentary Importance Index values (%) for 55 fish species:
trophica<-read.table("trophic.txt",header=TRUE) 
trophic<-(trophica[-1])
rownames(trophic) <- 1:nrow(trophic)

#Mean values for 12 ecomorphological traits with corrected body size (standard length) for 55 fish species:
morphological<-read.table("morpho_with_size_correction.txt",header=TRUE)
morphol<-(morphological[-1])
rownames(morphol) <- 1:nrow(morphol)

#Phylogenetic distance matrix between 55 fish species:
nudos.filogenetica<-read.table("filo100.txt",header=TRUE)


####--------Redundancy Analysis (RDA): the mean values for each of the 12 ecomorphological traits 
#form all 55 species formed a predictor variables matrix of original (non-standardised) descriptors, 
#while IAi % values provide the matrix with response variables.------------------------------------------

RDA<-rda(trophic~.,morphol)             
summary(RDA)
anova(RDA,permutations=how(nperm=999))                           #to identify the significance level (p-value = 0.05) of the model. 
(R2<- RsquareAdj(RDA)$r.squared)               
(r2adj <- RsquareAdj(RDA)$adj.r.squared)                         #ddjusted coefficient of determination r2_adj, i.e. % of variance explained by the model
coef(RDA)
vif.cca(RDA)                                                     #to assess whether there is collinearity (VIF>10) between the predictors
RDA_for<-forward.sel(trophic,morphol,adjR2thresh=R2, nperm=9999) #selects the variables that best explain dietary variation 

plot(RDA, scaling = 1,                                           #to plot the model
     cex = 2,
     xlab="RDA 1 = 29.95 %", 
     ylab = "RDA 2 = 11.98 %",
     xlim=c(-25,25), ylim=c(-15,15),
     font = 1,bty = "L",type = "p", family = "serif",
     cex.lab=1.3, cex.axis=1.3, cex.sub=1.3)


#Partial Mantel Tests: were performed to correlate the diet distance matrix (response)
# with the ecomorphological distance matrices and taxonomic distances

dist.trophica<-vegdist(trophic,method = "horn")    # to obtain a matrix of diet distance consumed among the 55 fish species
dist.morpho<-vegdist(morphol,method = "euclidean") # to obtain a distance matrix of the ecomorphological traits between 55 fish species
filo100<-(nudos.filogenetica[-1])      # the phylogeny between species is already a distance matrix 



#to test if morphology of species explains the diet consumed by species, removing the influence of phylogeny:
mantel.partial(dist.trophica,dist.morpho,filo100,method = "pearson",permutations = 1000)

# to test if relatedness between species explains the diet they consume, eliminating the influence of morphology:
mantel.partial(dist.trophica,filo100,dist.morpho,method = "pearson",permutations = 1000)







