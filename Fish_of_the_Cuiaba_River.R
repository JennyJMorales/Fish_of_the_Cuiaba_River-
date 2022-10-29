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
nudos.filogenetica<-read.table("filogenetica.txt",header=TRUE)
filo<-(nudos.filogenetica[-1])


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


#Linear Discriminant Analysis (LDA): to evaluate the ecomorphological traits 
#that were the indicators of whether body shape could be useful for predicting 
#the trophic guilds to which fish species belong in tropical savanna streams 

###----
plot.lda = function(lda.out, groups, colour.vec=NULL, plot.sites=1, plot.centroids=0, xax=1, yax=2, plot.env=TRUE, plot.ell=TRUE, title="LDA predicted classes", mul.coef=2, pos.names=NULL, col.env="black", xlim=NULL, ylim=NULL) 
  ### 
  # lda.out : Output file of function lda() in {MASS}.
  # groups  : Vector listing the group number for each object (factor or numeric).
  # colour.vec : personal vector with colour names, to be used for the groups in
  #    the plot instead of the standard colours. Example vector with 7 colours:
  #    col.vec = c("gray60","bisque","brown4","red","blue","darkgreen","orange4")
  #    Function colors() makes 657 colours available to R users.
  # plot.sites: 0 = Do not plot the sites.
  #             1 = plot symbols for the sites
  #             2 = print the site names
  # plot.centroids: 0 = Do not plot the group centroids.
#                 1 = Plot the group centroids (triangles) with assigned 
#                     group colours.
# xax : The axis to be used for abscissa of the plot.
# yax : The axis to be used for ordinate of the plot.
# plot.env=TRUE: plot the explanatory variables to the plot. FALSE: Do not 
#                plot them.
# plot.ell=TRUE: plot the 95% coverage ellipses of the groups. FALSE: Do not
#                plot them.
# title : Allows user to customize the title that will print above the plot.
# mul.coef : Multiplication factor for the length of the variable arrows. Some
#            trial and error is needed here.
# pos.names : Offset the names of the binary variables: 1=bottom, 2=left, 
#             3=top, 4=right.
# col.env="black" : Colour for the environmental variable arrows and names.
# xlim, ylim : Vectors describing the minimum and maximum values of the plotted region.
#
# License: GPL-2
# Authors: Daniel Borcard and Pierre Legendre

{
  if(class(lda.out) != "lda") stop("File lda.out was not produced by function lda of MASS")
  if(min(summary(as.factor(groups))) < 2) stop("There is at least one group with less than 2 observations")
  library(ellipse)
  coef = lda.out$scaling   # Standardized discriminant function coefficients
  k <- ncol(coef)  # Number of canonical axes
  lev <- length(levels(as.factor(groups)))  # Number of groups
  # print(c(k,lev))
  Fp <- predict(lda.out)$x
  centre <- matrix(NA,lev,k)
  for(i in 1:lev) { # Compute the group centroids in LDA space
    #	centre[i,] <- apply(Fp[groups==i,], 2, mean) }  
    centre[i,] <- apply(Fp[groups==levels(as.factor(gr))[i],], 2, mean) }  
  # print(centre)
  class.num <- as.numeric(predict(lda.out)$class) # Assignment of sites to classes
  # print(class.num)
  if(xax > k) stop("Their are not enough canonical axes; change the xax value")
  if(yax > k) stop("Their are not enough canonical axes; change the yax value")
  xlab=paste("LDA axis",xax," ")
  ylab=paste("LDA axis",yax," ")
  plot(Fp[,xax], Fp[,yax], type="n", main=title, xlab=xlab, ylab=ylab, xlim, ylim, asp=1)
  if(plot.sites==1) {   # Plot symbols for the sites
    if(is.null(colour.vec)) {
      points(Fp[,xax], Fp[,yax], pch=21, bg=class.num+1)
    } else {
      colour.sel <- colour.vec[class.num]
      points(Fp[,xax], Fp[,yax], pch=21, bg=colour.sel)
    }
  } else if(plot.sites==2) {
    if(is.null(colour.vec)) {
      text(Fp[,xax], Fp[,yax], row.names(Fp), col=class.num+1)
    } else {
      colour.sel <- colour.vec[class.num]
      text(Fp[,xax], Fp[,yax], row.names(Fp), col=colour.sel)
    }
  }	
  if(plot.centroids) {
    if(is.null(colour.vec)) {
      points(centre[,xax], centre[,yax], pch=24, bg=(1:lev)+1, cex=3)
    } else {
      colour.sel <- colour.vec[1:lev]
      # print(colour.sel)
      points(centre[,xax], centre[,yax], pch=24, bg=colour.sel, cex=3)
    }
  }
  abline(v=0, lty="dotted")
  abline(h=0, lty="dotted")
  # Draw 95% ellipses around the groups
  if(plot.ell) {
    for(i in 1:length(levels(as.factor(groups)))) { 
      #		cov <- cov(Fp[groups==i,c(xax,yax)])
      cov <- cov(Fp[groups==levels(as.factor(gr))[i],
                    c(xax,yax)])
      #		centre <- apply(Fp[groups==i,c(xax,yax)], 2, mean)
      centre <- apply(Fp[groups==levels(as.factor(gr))[i],
                         c(xax,yax)], 2, mean)
      lines(ellipse(cov, centre=centre, level=0.95))
    }
  }
  if(plot.env) { 
    arrows(x0=0, y0=0, x1=coef[,xax]*mul.coef, y1=coef[,yax]*mul.coef, 
           col=col.env, code=2, lty=1, length=0.1, lwd=1) 	
    if(!is.null(rownames(coef))) {
      text(x=coef[,xax]*mul.coef, y=coef[,yax]*mul.coef, 
           rownames(coef), col=col.env, cex=1, pos=pos.names) }
  }
}                                        #function to plot the LDA model
######---------------------------

trophic.euc<-decostand (trophic,"standardize")   #to standardise the dietary values for the 55 species
#They represent the trophic groups previously classified given the food consumed by the fish species:
gr <- c(1,2,2,2,2,2,5,2,1,1,3,1,2,2,2,2,2,2,2,2,2,2,2,3,2,4,4,2,5,4,2,2,2,1,2,2,2,2,2,4,2,4,2,2,2,2,2,2,2,4,2,2,3,4,1) 
# 1 = omnivore
# 2 = invertivore
# 3 = detritivore
# 4 = algivore/periphytivore
# 5 = piscivore

morp <- decostand (morphol,"normalize")             # to standardise the 12 ecomorphological traits
Wilks.test(morphol, gr)                             # to test with wilks' lambda that predictor variables have different means
                                                    # p-value = 0.001; i.e. the means are different

morpho_<-as.data.frame(morphol)
morp.sc <- as.data.frame (scale(morpho_))
morf.lda <-lda(gr~., data = morp.sc)                # to obtain discrimination functions: based on the standardised variables
summary(morf.lda)
morf.lda$means                                      # mean values of the 5 trophic groups for the ecomorphological traits
(c<-morf.lda$scaling)                               # to explaer classification functions
morf.lda$svd^2                                      # to calculate the canonical eigenvalues 
(Fp <- predict(morf.lda)$x)                         # to position the objects in the space of the canonical variants
(morf.class <-predict(morf.lda)$class)              # classification of objects
(morf.post <- round(predict(morf.lda)$posterior,2)) # posterior probability that objects belong to groups
(morf.table <-table(gr,morf.class))                 # contingency table of previous versus predicted classifications 
diag(prop.table(morf.table, 1))                     # proportion of correct classification (classification success)


plot.lda(lda.out = morf.lda,
         groups = gr,
         plot.sites = 1,
         plot.centroids = 0,
         xax=1, yax=2, 
         plot.env=TRUE, 
         plot.ell=TRUE, 
         title="Predicted groups by LDA model", 
         mul.coef=4, 
         pos.names= 4, 
         col.env="black", 
         xlim=c(-10,10), ylim= c(-10,10))


#Partial Mantel Tests: were performed to correlate the diet distance matrix (response)
# with the ecomorphological distance matrices and taxonomic distances

dist.trophica<-vegdist(trophic,method = "horn")    # to obtain a matrix of diet distance consumed among the 55 fish species
dist.morpho<-vegdist(morphol,method = "euclidean") # to obtain a distance matrix of the ecomorphological traits between 55 fish species
filo<-(nudos.filogenetica[-1])                     # the phylogeny between species is already a distance matrix 

#to test if morphology of species explains the diet consumed by species, removing the influence of phylogeny:
mantel.partial(dist.trophica,dist.morpho,filo,method = "pearson",permutations = 1000)

# to test if relatedness between species explains the diet they consume, eliminating the influence of morphology:
mantel.partial(dist.trophica,filo,dist.morpho,method = "pearson",permutations = 1000)







