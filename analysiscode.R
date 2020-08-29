#####
## Code to analyze environmental and clonal varation and the link between them in white top in N Colorado
## Authors: Natalie West and Amanda Stahlke
## Last updated 2020-08-25

#### Libraries ####
library(BiodiversityR)
library(vegan)
library(goeveg)
library(lme4)
library(MuMIn)
library(rpart)
library(party)
library(subselect)
library(MASS)


library(dplyr)
library(reshape2)
library(ggplot2)
library(ggstatsplot)

set.seed(123)

#### Import data and subset ####
gen_env<-na.omit(read.csv("geno_enviro.csv"))
names(gen_env)
#gen_env$BL<-gen_env$bare+gen_env$litter #composite background variable
gen_env$Site<-as.factor(gen_env$Site)
head(gen_env)

site_GN_rank <- gen_env %>%
  group_by(Site) %>%
  summarise(unique(GN_Rank))
site_GN_rank 


clone_rate <- gen_env %>%
  group_by(Site) %>%
  summarise(sum(Clone)) %>%
  mutate(clone_rate=1-`sum(Clone)`/15) %>%
  as.data.frame()
merge(clone_rate, site_GN_rank)

stemcts <- gen_env %>%
  # group_by(Site) %>%
  summarise(avg=mean(stem), sd=sd(stem)) %>%
  mutate(avg/(100))

arrange(gen_env, stem) %>%
  dplyr::select(site.plot.subplot, stem)

#### Veg summary ####
## Aggregate site variables for the average and standard deviation across a site. ##

veg_mat<-aggregate(list("grass"=gen_env$grass,"bare"=gen_env$bare,"draba"=gen_env$Ldraba.cvr,"litter"=gen_env$litter, "forb"=gen_env$forb),
                   by=list("Site"=gen_env$Site),FUN=mean)

veg_mat.sd<-aggregate(list("grass"=gen_env$grass,"bare"=gen_env$bare,"draba"=gen_env$Ldraba.cvr,"litter"=gen_env$litter),
                      by=list("Site"=gen_env$Site),FUN=sd)

GN<-read.csv("G_N_bysite.csv") ### this clonal richness is different from what we ultimately report, but rank is the same
veg_mat$R_GN<-GN[,"R_GN"]
## R_GN is the rank of G/N across sites.
## From least to most clonal: Site 1, 6, 3, 5, 2, 4. 

veg_mat$GN<-GN$site.G.N
veg_mat

plot(veg_mat$R_GN,veg_mat$draba,xaxt="n",ylim=c(0,60),ylab="% Draba",xlab="Prevalence of Clones")
plot(veg_mat$R_GN,veg_mat$grass,xaxt="n",ylim=c(0,60),ylab="% Grass", xlab="Prevalence of Clones")
plot(veg_mat$R_GN,veg_mat$bare,xaxt="n",ylim=c(0,60),ylab="% Bare",xlab="Prevalence of Clones")


#### Are there general site-level differences among measured variables we think are probably important? ####
ggplot(data = gen_env, aes(x=canopy)) +
  geom_histogram()

ggplot(data = gen_env, aes(x=Ldraba.cvr)) +
  geom_histogram(bins=100)

ggplot(data = gen_env, aes(x=stem)) +
  geom_histogram(bins=100)

ggplot(data = gen_env, aes(x=Site, y=pH)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)

ggplot(data = gen_env, aes(x=Site, y=canopy)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggsave("canopy_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=woody)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggsave("woody_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=K)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)

ggplot(data = gen_env, aes(x=Site, y=litter)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggplot(data = gen_env, aes(x=Site, y=grass)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggsave("grass_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=forb)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggsave("forb_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=bare)) +
  geom_boxplot() + 
  geom_point(alpha=0.3)
ggsave("bare_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=stem)) +
  geom_boxplot() + 
  geom_point(alpha=0.3) 
ggsave("stem_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=Ldraba.cvr)) +
  geom_boxplot() + 
  geom_point(alpha=0.3) 
ggsave("Ldraba.cvr_site.jpg")

ggplot(data = gen_env, aes(x=Site, y=Sand)) +
  geom_boxplot() + 
  geom_point(alpha=0.3) 
ggsave("Sand_site.jpg")

#### Correlation among vriables ####
cormat <- round(cor(gen_env[,c(6:9,11:dim(gen_env)[2])],  method = "spearman"),2) # drop texture and site

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE)

ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   hjust = 1))+
  coord_fixed()

melted_cormat %>%
  filter(abs(value) >.5 & value <1) %>%
  arrange(desc(value))

melted_cormat %>%
  filter(Var1 == "Clone" |  Var2 == "Clone") %>%
  arrange(desc(value))

melted_cormat %>%
  filter(Var1 == "stem" |  Var2 == "stem") %>%
  arrange(desc(abs(value)))

melted_cormat %>%
  filter(Var1 == "bare" |  Var2 == "bare") %>%
  arrange(desc(abs(value)))

melted_cormat %>%
  filter(Var1 == "grass" |  Var2 == "grass") %>%
  arrange(desc(abs(value)))

melted_cormat %>%
  filter(Var1 == "forb" |  Var2 == "forb") %>%
  arrange(desc(abs(value)))

melted_cormat %>%
  filter(Var1 == "Sand" |  Var2 == "Sand") %>%
  arrange(desc(abs(value)))

melted_cormat %>%
  filter(Var1 == "canopy" |  Var2 == "canopy") %>%
  arrange(desc(abs(value)))

### set a threshold - remove variables that have greater that 0.5
melted_cormat %>%
  filter(abs(value) >.5 & value <1) %>%
  arrange(desc(abs(value)))

#### NMDS of subplots ####
all_env_fact<-gen_env[,11:dim(gen_env)[2]]
names(all_env_fact)

env.mds<-metaMDS(all_env_fact, distance = "gower", k=4)
stressplot(env.mds)

env.mds2<-metaMDS(all_env_fact, distance = "gower", k=4, previous.best = env.mds)
stressplot(env.mds2)
env.mds2$species

dim.result <- dimcheckMDS(all_env_fact, distance="gower", k=8, trymax = 20)
data.frame(dimension = seq_along(dim.result), stress=dim.result)

### 2 dims is very close to rule of thumb stress below 0.206

#### nMMDS viz ####
env.plot<-ordiplot(env.mds,display = "sites",type="n", choices = c(1,2))
ordisymbol(env.plot, gen_env, "Site", legend=TRUE,cex=2, lwd=2, choices = c(1,2))
env.factors<-envfit(env.mds, all_env_fact, permutations = 999)
plot(env.factors, p.max=0.01)
### this is pretty crowded; make this in ggplot with  the factors we won't use later in grey
### https://jkzorz.github.io/2020/04/04/NMDS-extras.html
### Keep OM, Zn, pH, Fe, Sand ###
fact_keep <- c("pH", "OM", "Zn", "Fe", "P", "Sand", "grass", "forb", "canopy")
fact_drops  <- setdiff(names(all_env_fact), fact_keep)
env_fact<-data.frame(all_env_fact[,-which(names(all_env_fact)%in%fact_drops)])
head(env_fact)

#### create a composite variable of soil variables via NMDS? ####
fact_drops  <- c("canopy", "grass", "forb", "woody", "bare", "Sand", "litter", "Clay", "Silt")
soilz <- data.frame(all_env_fact[,-which(names(all_env_fact)%in%fact_drops)])
head(soilz)

soilz.mds<-metaMDS(soilz,distance = "gower", k=4)
soilz.mds<-metaMDS(soilz,distance = "gower", k=4,
                   previous.best = soilz.mds)
dim.soil.result <- dimcheckMDS(soilz,distance="gower",k=4,trymax = 20)
data.frame(dimension = seq_along(dim.soil.result), stress=dim.soil.result)

nut.plot<-ordiplot(soilz.mds,display = "sites",type="n")
ordisymbol(nut.plot, gen_env, "Site", legend=TRUE,cex=2,lwd=2)
nut.factors<-envfit(soilz.mds,soilz)
plot(nut.factors)

soil.nut1<-soilz.mds$points[,"MDS1"]
soil.nut2<-soilz.mds$points[,"MDS2"]

##### Examine multivariate dispersion among sites ####
## Betadispersion - sites. Create a dissimilarity matrix for betdispersion analysis; Looks at multivariate dispersion (variance) among SITES
commax<-vegdist(all_env_fact, method = "gower")
fum<-with(gen_env,betadisper(commax,Site))
plot(fum)
fum
boxplot(fum) # variation among sampled points
permutest(fum) #-- Difference in multivariate dispersion between sites
plot(TukeyHSD(fum))


##### Do we need to account for site? ####
env2 <- as.matrix(all_env_fact)
colnames(env2)

anosim <- anosim(env2, grouping = gen_env$Site, distance = "gower", permutations = 9999) 
anosim
summary(anosim) ### yep

adonis <- adonis2(env2 ~ gen_env$Site)
adonis ### yep

##### Can we model a linear relationship btw environment and subplot clonality? ####

### Use uncorrelated OM, Zn, pH, Fe, Sand
### I think I can actually just skip this because of the above anosim/adonis results

### remove? 
# cl_global_glm <- glm(Clone ~
#                        I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
#                        OM + Zn + Fe + pH  + P,
#                      data=gen_env, family="binomial",
#                      na.action = "na.fail",
#                      weights = rep(15,dim(gen_env)[1]))
# summary(cl_global_glm)
# ## Sand and OM not significant anymore? :(
# 
# ## cooks distance to identify overly influential points
# cooksd<-cooks.distance(cl_global_glm)
# plot(cooksd)
# abline(h = 4*mean(cooksd, na.rm=T), col="red")
# text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
# influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])
# 
# cbind(gen_env[influential,])
# 
# gen_env_pruned <- gen_env[-influential,]
# 
# cl_global_glm <- glm(Clone ~ 
#                        I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
#                        OM + Zn + Fe + pH  + P,
#                      data=gen_env_pruned, family="binomial",
#                      na.action = "na.fail",
#                      weights = rep(15-length(influential),dim(gen_env_pruned)[1]))
# summary(cl_global_glm) 
### after removing influential subplots, pH is not significant

# cl_global_glm_pruned <- glm(Clone ~ 
#                        I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
#                        OM + Zn + Fe  + P,
#                      data=gen_env_pruned, family="binomial",
#                      na.action = "na.fail",
#                      weights = rep(15-length(influential),dim(gen_env_pruned)[1]))
# summary(cl_global_glm_pruned)

#### GLMER ####
### what should we do about residulas? 

cl_global_glmer <- glmer(Clone ~ 
                       I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                       OM + Zn  + Fe  + pH + 
                         (1|gen_env$Site),
                     data=gen_env, family="binomial",
                     na.action = "na.fail",
                     control = lme4::glmerControl(
                       optimizer = "Nelder_Mead",
                       calc.derivs = FALSE,
                       boundary.tol = 1e-7),
                     weights = rep(15,dim(gen_env)[1]))
summary(cl_global_glmer) ## everything but canopy

### 1. Test for overly influential points
cl_global_glm <- glm(Clone ~
                        I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                        OM + Zn + Fe + pH, 
                      data=gen_env, family="binomial",
                      na.action = "na.fail",
                      weights = rep(15,dim(gen_env)[1]))
summary(cl_global_glm) ## residual deviance >>>> df :(
## everything but Sand, OM, Zn

## cooks distance to identify overly influential points
cooksd<-cooks.distance(cl_global_glm)
plot(cooksd)
abline(h = 4*mean(cooksd, na.rm=T), col="red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])

cbind(gen_env[influential,]) ## high pH subplots that are all clonal

gen_env_pruned <- gen_env[-influential,]

cl_global_glm_pruned <- glm(Clone ~
                       I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                       OM + Zn + Fe + pH,
                     data=gen_env_pruned, family="binomial",
                     na.action = "na.fail",
                     weights = rep(13,dim(gen_env_pruned)[1]))
summary(cl_global_glm_pruned) # +/- everything

cl_global_glmer_pruned <- glmer(Clone ~ 
                           I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                           OM + Zn  + Fe  + pH + 
                           (1|gen_env_pruned$Site),
                         data=gen_env_pruned, family="binomial",
                         na.action = "na.fail",
                         control = lme4::glmerControl(
                           optimizer = "Nelder_Mead",
                           calc.derivs = FALSE,
                           boundary.tol = 1e-7),
                         weights = rep(13,dim(gen_env_pruned)[1]))
summary(cl_global_glmer_pruned) ## everything but OM

## re-run cooks distance to identify overly influential points
cooksd<-cooks.distance(cl_redux1_glm_pruned)
plot(cooksd)
abline(h = 4*mean(cooksd, na.rm=T), col="red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
## cooksd not SO huge

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])
cbind(gen_env_pruned[influential,]) ## now some clonal, some not; all from less clonal sites 3, 5, 6
gen_env_pruned2 <- gen_env_pruned[-influential,]

cl_redux1_glm_pruned2 <- glm(Clone ~
                              I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                              OM + Zn + Fe + pH,
                            data=gen_env_pruned2, family="binomial",
                            na.action = "na.fail",
                            weights = rep(13,dim(gen_env_pruned2)[1]))
summary(cl_redux1_glm_pruned2) # everuthing but Zn

glm_cl <- dredge(cl_redux1_glm_pruned2)
importance(subset(glm_cl, delta<=2))

cl_global_glmer_pruned2 <- glmer(Clone ~ 
                           I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                           OM + Zn  + Fe +  pH + 
                           (1|gen_env_pruned2$Site),
                         data=gen_env_pruned2, family="binomial",
                         na.action = "na.fail", 
                         control = lme4::glmerControl(
                           optimizer = "Nelder_Mead",
                           calc.derivs = FALSE,
                           boundary.tol = 1e-7),
                         weights = rep(13, dim(gen_env_pruned2)[1]))
summary(cl_global_glmer_pruned2) ## everything but OM
summary(cl_global_glmer_pruned)
## pruned and pruned2 yield qualitatively the same results for glm
## but very different for glmer
## carry on without dropping those 'outliers'???

### controlling for site as a mixed effect has a strong effect on predictors:
### relative effect sizes are the same, but many more variables are significant 
### without controlling for site effects. Hmm....

### Start back at the global glmer, pruned, without OM

cl_redux1_glmer <- glmer(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                            Zn + pH + 
                           (1|gen_env_pruned$Site),
                         data=gen_env_pruned, family="binomial",
                         na.action = "na.fail",
                         control = lme4::glmerControl(
                           optimizer = "Nelder_Mead",
                           calc.derivs = FALSE,
                           boundary.tol = 1e-7),
                         weights = rep((15-rm),dim(gen_env_pruned)[1]))
summary(cl_redux1_glmer)

d_pruned <- dredge(cl_redux1_glmer) 
importance(subset(d_pruned, delta<=4)) ## canopy, grass, P, Zn

cl_redux1_glmer <- glmer(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                           Zn + pH + 
                           (1|gen_env_pruned$Site),
                         data=gen_env_pruned, family="binomial",
                         na.action = "na.fail",
                         control = lme4::glmerControl(
                           optimizer = "Nelder_Mead",
                           calc.derivs = FALSE,
                           boundary.tol = 1e-7),
                         weights = rep((15-rm),dim(gen_env_pruned)[1]))
summary(cl_redux1_glmer)


#### viz regression models ####
ggstatsplot::ggcoefstats(
  x = cl_redux1_glmer,
  point.args = list(color = "red", size = 3, shape = 15),
  vline.args = list(size = 1, color = "#CC79A7", linetype = "dotdash"),
  # stats.label.color = c("#0072B2", "#D55E00", "darkgreen"),
  ggstatsplot.layer = FALSE
) + # note the order in which the labels are entered
  ggplot2::labs(x = "regression coefficient", y = NULL)
ggsave("subplot_clonality_glmer.jpg", height = 4.5, width = 6, units = "in", dpi = 300)
dev.off()

#### Do the same predictors of clonality predict Stem density? ####
stem_global_glm <-glm(stem ~ I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                            OM + Zn  + Fe + pH,
                          data=gen_env_pruned, 
                          family = "poisson",
                          na.action = "na.fail")
summary(stem_global_glm)

stem_global_glmer <-glmer(stem ~ I(Sand/100) + I(grass/100) + 
                            I(forb/100)+ I(canopy/100) +
                            OM + Zn  + Fe + pH +
                            (1|gen_env_pruned$Site),
                          data=gen_env_pruned, 
                          control = lme4::glmerControl(
                            optimizer = "Nelder_Mead",
                            calc.derivs = FALSE,
                            boundary.tol = 1e-7),
                          family = "poisson",
                          na.action = "na.fail")
summary(stem_global_glmer)

d_stem_pruned <- dredge(stem_global_glmer) ### 4 convergence warnings
importance(subset(d_stem_pruned, delta<=2)) ## Fe, grass, OM, pH, Sand, Zn

stem_redux_glmer <-glmer(stem ~ I(Sand/100) + I(grass/100) + 
                            OM + Fe + pH +
                            (1|gen_env_pruned$Site),
                          data=gen_env_pruned, 
                          control = lme4::glmerControl(
                            optimizer = "Nelder_Mead",
                            calc.derivs = FALSE,
                            boundary.tol = 1e-7),
                          family = "poisson",
                          na.action = "na.fail")
summary(stem_redux_glmer)

ggstatsplot::ggcoefstats(
  x = stem_redux_glmer,
  point.args = list(color = "red", size = 3, shape = 15),
  vline.args = list(size = 1, color = "#CC79A7", linetype = "dotdash"),
  # stats.label.color = c("#0072B2", "#D55E00", "darkgreen"),
  ggstatsplot.layer = FALSE
) + # note the order in which the labels are entered
  ggplot2::labs(x = "regression coefficient", y = NULL)
ggsave("subplot_stem_glmer.jpg", height = 4.5, width = 6, units = "in", dpi = 300)
dev.off()

ggstatsplot::ggcoefstats(
  x = cl_global_mds_glmer_pruned,
  point.args = list(color = "red", size = 3, shape = 15),
  vline.args = list(size = 1, color = "#CC79A7", linetype = "dotdash"),
  # stats.label.color = c("#0072B2", "#D55E00", "darkgreen"),
  ggstatsplot.layer = FALSE
) + # note the order in which the labels are entered
  ggplot2::labs(x = "regression coefficient", y = NULL)
ggsave("subplot_clonality_nmds12_glmer.jpg", height = 4.5, width = 6, units = "in", dpi = 300)
dev.off()


#### use nMDS instead in regression? ####
## http://r-sig-ecology.471788.n2.nabble.com/NMDS-axes-scores-td7579648.html
## axis1 + axis2 + axis1:axis2 ???
gen_env$soil.nut1 <- soil.nut1
gen_env$soil.nut2 <- soil.nut2

cl_global_mds_glmer <- glmer(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                               soil.nut1 + soil.nut2 + soil.nut1:soil.nut2 +
                               (1|gen_env$Site),
                             data=gen_env, family="binomial",
                             na.action = "na.fail",
                             control = lme4::glmerControl(
                               optimizer = "Nelder_Mead",
                               calc.derivs = FALSE,
                               boundary.tol = 1e-7),
                             weights = rep((15-rm),dim(gen_env)[1]))
summary(cl_global_mds_glmer)

## glm
cl_global_mds_glm <- glm(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                           soil.nut1 + soil.nut2 + soil.nut1:soil.nut2,
                         data=gen_env, family="binomial",
                         na.action = "na.fail",
                         weights = rep(15, dim(gen_env)[1]))
summary(cl_global_mds_glm)

## remove outsized
cooksd<-cooks.distance(cl_global_mds_glm)
plot(cooksd)
abline(h = 4*mean(cooksd, na.rm=T), col="red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])
cbind(gen_env[influential,]) ## mostly clonal; all from less clonal sites 1, 3, 5 (same as before?)
length(influential)

gen_env_pruned <- gen_env[-influential,]
soil.nut1_pruned <- soilz.mds$points[-influential,"MDS1"]
soil.nut2_pruned <- soilz.mds$points[-influential,"MDS2"]
gen_env_pruned$soil.nut1_pruned <- soil.nut1_pruned
gen_env_pruned$soil.nut2_pruned <- soil.nut2_pruned

## rerun glm
cl_global_mds_glm_pruned <- glm(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                                  soil.nut1_pruned + soil.nut2_pruned + soil.nut1_pruned:soil.nut2_pruned,
                                data=gen_env_pruned, family="binomial",
                                na.action = "na.fail",
                                weights = rep(13, dim(gen_env_pruned)[1]))
summary(cl_global_mds_glm_pruned) ## qualitatively the same but interaction term not sig

cl_global_mds_glmer_pruned <- glmer(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                                      soil.nut1_pruned + soil.nut2_pruned + soil.nut1_pruned:soil.nut2_pruned +
                                      (1|gen_env_pruned$Site),
                                    data=gen_env_pruned, family="binomial",
                                    na.action = "na.fail",
                                    control = lme4::glmerControl(
                                      optimizer = "Nelder_Mead",
                                      calc.derivs = FALSE,
                                      boundary.tol = 1e-7),
                                    weights = rep(13,dim(gen_env_pruned)[1]))
summary(cl_global_mds_glmer_pruned)

## nmds and stem density ####
cl_global_mds_glmer <- glmer(Clone ~  I(grass/100) + I(forb/100)+ I(canopy/100) +
                               soil.nut1 + soil.nut2 + soil.nut1:soil.nut2 +
                               (1|gen_env$Site),
                             data=gen_env, family="binomial",
                             na.action = "na.fail",
                             control = lme4::glmerControl(
                               optimizer = "Nelder_Mead",
                               calc.derivs = FALSE,
                               boundary.tol = 1e-7),
                             weights = rep((15-rm),dim(gen_env)[1]))
summary(cl_global_mds_glmer)

## glm
stem_global_mds_glm <-glm(stem ~ I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                        soil.nut1 + soil.nut2 + soil.nut1:soil.nut2,
                      data=gen_env, 
                      family = "poisson",
                      na.action = "na.fail")
summary(stem_global_mds_glm)

## remove outsized
cooksd<-cooks.distance(stem_global_mds_glm)
plot(cooksd)
abline(h = 4*mean(cooksd, na.rm=T), col="red")
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])
cbind(gen_env[influential,]) ## mostly clonal; all from less clonal sites 1, 3, 5 (same as before?)
length(influential)

gen_env_pruned <- gen_env[-influential,]
soil.nut1_pruned <- soilz.mds$points[-influential,"MDS1"]
soil.nut2_pruned <- soilz.mds$points[-influential,"MDS2"]
gen_env_pruned$soil.nut1_pruned <- soil.nut1_pruned
gen_env_pruned$soil.nut2_pruned <- soil.nut2_pruned

## rerun glm
stem_global_mds_glm_pruned <-glm(stem ~ I(Sand/100) + I(grass/100) + I(forb/100)+ I(canopy/100) +
                            soil.nut1_pruned + soil.nut2_pruned + soil.nut1_pruned:soil.nut2_pruned,
                          data=gen_env_pruned, 
                          family = "poisson",
                          na.action = "na.fail")
summary(stem_global_mds_glm_pruned)

stem_global_glmer_pruned <-glmer(stem ~ I(Sand/100) + I(grass/100) + 
                            I(forb/100)+ I(canopy/100) +
                              soil.nut1_pruned + soil.nut2_pruned + soil.nut1_pruned:soil.nut2_pruned +
                            (1|gen_env_pruned$Site),
                          data=gen_env_pruned, 
                          control = lme4::glmerControl(
                            optimizer = "Nelder_Mead",
                            calc.derivs = FALSE,
                            boundary.tol = 1e-7),
                          family = "poisson",
                          na.action = "na.fail")
summary(stem_global_glmer_pruned)

