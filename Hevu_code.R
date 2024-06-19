## Install packages


#install.pacakges("virtualspecies")
#install.packages("devtools")
#devtools::install_github("marlonecobos/kuenm")
#install.pacakges("spThin")
#install.packages("hierfstat")
#install.packages("vegan")
#install.packages("raster")
#install.packages("related", repo="https://R-Forge.R-project.org")
#install.packages("dplyr")
#if(!requireNamespace("corMLPE", quietly = TRUE)) 
#  remotes::install_github("nspope/corMLPE")
#install.packages"adegenet")
#install.packages("ggplot2")
#install.packages("leaflet")
#install.packages("ape")
#install.packages("ggdendro")
#install.packages("tidyverse")


## Libraries


library(virtualspecies)
library(kuenm)
library(spThin)
library(hierfstat)
library(vegan)
library(raster)
library(related)
library(dplyr)
library(corMLPE)
library(adegenet)
library(ggplot2)
library(leaflet)
library(ape)
library(ggdendro)
library(tidyverse)

## Spatial thinning - occurrence points


occs <- read.csv("full.csv")


occs <- occs[,-c(1)] # remove species column
colnames(occs) # double check

cts <- thin.algorithm(rec.df.orig = occs, thin.par = 1, reps = 5)

View(cts[[1]])

l <- cts[[1]]

write.csv(l, "joint.csv")

train_bin <- sample_frac(cts[[1]], size = 0.70, replace = FALSE) # randomly sample 70% of data points from the thinned data set.
test_bin  <- dplyr::anti_join(cts[[1]], train_bin) # assign the other 30%.

write.csv(train_bin, "train.csv")
write.csv(test_bin, "test.csv")


## Correlation - environmental layers


files <- list.files(pattern = "asc")

# stack the rasters into one object

bio_cor <- stack(files)

removeCollinearity(
  bio_cor,
  multicollinearity.cutoff = 0.7,
  select.variables = TRUE,
  method = "pearson")


## SDM construction

# kuenm - model calibration


dir <- getwd()

# The next chunk of code is for preparing the arguments for using the function following the modularity principle. 
# These variables can be changed according to each case.


occ_joint <- "joint.csv" # all spatial coordinates
occ_tra <- "train.csv" # training data set created from spThin
M_var_dir <- "M_variables" # directory where buffered variables are
batch_cal <- "Candidate_models" # batch script to be run 
out_dir <- "Candidate_Models" # output directory
reg_mult <- c(1,2,3,4,5,6,7,8,9,10) # regularization multipliers from 1-10
f_clas <- "all"
args <- NULL
maxent_path <- dir # where the maxent jar file is located
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          maxent.path = maxent_path, wait = wait, run = run)


# kuenm - model evaluation


# Evaluation is a crucial step in model calibration. This step centers on selecting candidate models and their 
# associated parameters to identify the best models for the purposes of the study. The kuenm_ceval function evaluates
# candidate models based on three distinct criteria: statistical significance (based on partial ROC analyses), 
# prediction ability (we use omission rates, but other metrics, such as overall correct classification rate, 
# can also be used), and model complexity (here evaluated using AICc). The following code chunk calls the function 
# help window.


occ_test <- "test.csv" # test data set created from spThin
out_eval <- "Calibration_results" # output directory
threshold <- 5
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform MOP calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for the previous function
if(utils::packageVersion("dplyr") > "0.50"){
  
} else {
  
}

# This code also allows evaluating candidate models that were created previously, selecting those with best 
# performance based on the three criteria.

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, 
                        batch = batch_cal, out.eval = out_eval, threshold = threshold, 
                        rand.percent = rand_percent, iterations = iterations, kept = kept, 
                        selection = selection, parallel.proc = paral_proc)


# kuenm - model creation


# After selecting parameters that produce best models, the next step is to create the final models,
# and if needed transfer them to other environmental data sets (e.g., to other time periods or other geographic
# regions). The function help is called via this code


batch_fin <- "Final_models" # batch script
mod_dir <- "Final_Models" # output directory
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "logistic"
project <- TRUE
G_var_dir <- "G_variables" # directory with subdirectories for projected scenarios and extents
ext_type <- "all"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL
# Again, some of the variables used here as arguments were already created for the previous functions

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)


# kuenm - final model evalutation


occ_ind <- "ind.csv" # independent spatial coordinates not used in the joint, test, or train data set
replicates <- TRUE
out_feval <- "Final_Models" # output directory
mod_dir <- "Final_Models" # batch script
# Most of the variables used here as arguments were already created for previous functions


fin_eval <- kuenm_feval(path = mod_dir, occ.joint = occ_joint, occ.ind = occ_ind, replicates = replicates,
                        out.eval = out_feval, threshold = threshold, rand.percent = rand_percent,
                        iterations = iterations, parallel.proc = paral_proc)



## SDM thresholding 


#HTML provided is the website where this code was adapted from. 

#https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/#:~:text=10th%20percentile%20training%20presence&text=It%20assumes%20that%20the%2010,greater%20region%20than%20the%20MTP.

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

sdm <- raster("sdm.asc")

occs <- read.csv("train.csv")

plot(sdm)
points(occs[,2:3], pch = 19, cex = 0.5)

occPredVals <- raster::extract(sdm, occs[,2:3])

p10 <- sdm_threshold(sdm, occs[,2:3], "p10")

plot(p10)


writeRaster(p10, "p10.asc")

## SNPfiltR

h <- read.vcfR("file.vcf")

h

popmap<-read.csv("popmap.txt", sep = "\t")

popmap <- popmap[-c(43,60:78),] # remove Americana, Zion, and Pinaleno individuals

h <- hard_filter(h, gq = 30, depth = 5)

h

h <- filter_allele_balance(h)

h

max_depth(h)

h <- max_depth(h, maxdepth = 50)

h

missing_by_sample(vcfR=h, popmap = popmap)

#run function to drop samples above the threshold we want from the vcf
h<-missing_by_sample(vcfR=h, cutoff = .9)

#subset popmap to only include retained individuals
popmap<-popmap[popmap$id %in% colnames(h@gt),]

#remove invariant sites generated by dropping individuals
h<-min_mac(h, min.mac = 4)

h

vcfR::write.vcf(h, "./filtered.vcf")

## Relatedness


# This code was modified from J.D.M.

x <- readgenotypedata("file.txt") # modified text file of Structure output from populations module


species <- "vulnerata"

# input colony number (e.g., 1, 2, or 3)
colony_number <- 1

# define number of simulations
nsims <- 100

# simulations
# first, obtain overall allele freqs for each locus
# allow 5% variation of allele frequencies (i.e., a little randomness)
n_loci <- (ncol(x) - 1) / 2
alleles <- list()
allele_freqs <- list()

# add all real and simulated data into a single data frame for relatedness estimates
total <- x
#offspring_sims2 <- do.call("rbind", offspring_sims)
#half_sibs2 <- do.call("rbind", half_sibs)
#unrelated2 <- do.call("rbind", unrelated)
#total <- rbind(total, offspring_sims2, half_sibs2, unrelated2)

# clean up objects no longer needed and clean up garbage
#rm(offspring_sims2, half_sibs2, unrelated2, offspring_sims, half_sibs, unrelated)

# measure relatedness of all real and simulated individuals (this may take up to 30-60 minutes)
total_relate <- coancestry(x$gdata, lynchli=2, wang=2)
View(total_relate)


## Heterozygosity

gen <- adegenet::read.genepop("populations.snps.gen") # from populations module
bs <- basic.stats(gen, diploid = T) # this will generate several, total summary statistics
# it can be modified to have individual populations to run basic stats on for He and Ho estimates
# I split up populations by opening up the gen file in a text editor and then iteratively importing them into R and using basic.stats

bs$overall

## Dendrogram

k <- read.csv("file.csv", header = TRUE, row.names = "Individuals")

View(k)

k <- k[,-1]

mtc <- dist(scale(k), method = "euclidean")

hc <- hclust(mtc, method = "ward.D2")

plot(hc)

plot(as.phylo(hc),  cex = 1,
     edge.color = "black", edge.width = 2,
     tip.color = "black")

## FST


gen <- adegenet::read.genepop("populations.snps.gen") # from populations module
gen.87 <- genet.dist(gen, method = "Nei87") 
gen.87
gen.87 <- as.matrix(gen.87)
write.table(gen.87, "gen.txt")


## PCA

# outputs from SNPfiltR and Plink here

pca <- read_table2("hevu_13June2024.eigenvec", col_names = FALSE)
eigenval <- scan("hevu_13June2024.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot pca
b <- ggplot(pca, aes(X3, X4, col = X2)) + geom_point(size = 3) + scale_color_manual(values = c("#4D99FF", "#FF9999", "#FF9999", "#008080", "#FF9999", "#7326E5"))
b <- b + coord_equal() +stat_ellipse() + theme()
b + xlab(bquote("PC1 (28.56%)")) + ylab(bquote("PC2 (12.10%)")) + guides(color = guide_legend(title = "Locality"))
b

# outputs from Stacks here

hevu <- read.table("populations.structure", sep = "\t", header = TRUE) # from populations module
View(hevu)
hevu <- hevu %>% distinct(Individuals, .keep_all = TRUE) # in the case there are repeat individuals

hevu_c <- hevu[,-c(1:2)] # removes the columns with individual names and pop IDs
View(hevu_c)
pca_h <- prcomp(hevu_c)
summary(pca_h)
pci <- data.frame(pca_h$x, Locality = hevu$Locale) # create a data frame with the pca output, appending locality assignments for ggplot
h <- ggplot2::ggplot(pci, aes(x = PC1, y = PC2, color = hevu$Locale)) + 
  geom_point(size = 2) +
  scale_color_manual(values = c("#4D99FF", "#FF9999", "#FF9999", "#008080", "#FF9999", "#7326E5")) +
  stat_ellipse() +
  xlab(bquote("PC1 (7.38%)")) + 
  ylab(bquote("PC2 (7.12%)")) +
  guides(color = guide_legend(title = "Locality")) +
  theme()


## IBD


Dgeo <- read.table("geo.txt") # generate from decimal degrees between each population
Dgen <- read.table("gen.txt") # generated from FST between each population from above code
Dgeo <- as.matrix(Dgeo)
Dgen <- as.matrix(Dgen)

ibd.v <- vegan::mantel(Dgen, Dgeo, method = "pearson", permutations = 9999)


## MLPE


points_to_line <- function(data, long, lat, id_field = NULL, sort_field = NULL) {
  
  # Convert to SpatialPointsDataFrame
  coordinates(data) <- c(long, lat)
  
  # If there is a sort field...
  if (!is.null(sort_field)) {
    if (!is.null(id_field)) {
      data <- data[order(data[[id_field]], data[[sort_field]]), ]
    } else {
      data <- data[order(data[[sort_field]]), ]
    }
  }
  
  # If there is only one path...
  if (is.null(id_field)) {
    
    lines <- SpatialLines(list(Lines(list(Line(data)), "id")))
    
    return(lines)
    
    # Now, if we have multiple lines...
  } else if (!is.null(id_field)) {  
    
    # Split into a list by ID field
    paths <- sp::split(data, data[[id_field]])
    
    sp_lines <- SpatialLines(list(Lines(list(Line(paths[[1]])), "line1")))
    
    # I like for loops, what can I say...
    for (p in 2:length(paths)) {
      id <- paste0("line", as.character(p))
      l <- SpatialLines(list(Lines(list(Line(paths[[p]])), id)))
      sp_lines <- spRbind(sp_lines, l)
    }
    
    return(sp_lines)
  }
}

# This code was modified from https://bookdown.org/hhwagner1/LandGenCourse_book/WE_12.html as part of a landscape genetics course

dat <- list.files(pattern = "csv") # import all csv files with decimal degrees between populations to create vectors between each one. There should be a csv file per pair. 


v_lines <- points_to_line(data = dat, # create the vector lines
                          long = "long", 
                          lat = "lat",
                          sort_field = "sort")

# visualize the lines
leaflet(data = v_lines) %>%
  addTiles() %>%
  addPolylines()

# create a buffer (in meters) around each line
v_buff <- buffer(v_lines, width = 2000)

# once again, visualize the buffered lines
leaflet(data = v_buff) %>%
  addTiles() %>%
  addPolylines()

files <- list.files(pattern = "asc")

# stack them for quicker processing
env <- stack(files)

# take the values of each environmental variables along the buffer
vals <- extract(env, v_buff, df = T)


# model1 --> grinnellian niche 
# model2 --> bioclim niche

MLPE <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = TRUE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

mod1 <- MLPE(log.fst ~ bio2 + bio8 + bio9 + bio15 +
               direction + stream + dem + tcc +
               (1|pop1), data)

mod2 <- MLPE(log.fst ~ bio2 + bio8 + bio9 + 
               bio17 + bio18 + bio19 + (1|pop1), data)

MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}

mod1noREML <- MLPEnoREML(log.fst ~ bio2 + bio8 + bio9 + bio15 +
                           direction + stream + dem + tcc + (1|pop1), data)

mod2noREML <- MLPEnoREML(log.fst ~ bio2 + bio8 + bio9 + 
                           bio17 + bio18 + bio19 + (1|pop1), data)


Models <- list(Grinnellian = mod1noREML, Bioclim = mod2noREML)
IC <- data.frame(AIC = sapply(Models, AIC),
                 BIC = sapply(Models, BIC)) 
IC

IC <- data.frame(IC, k = sapply(Models, function(ls) attr(logLik(ls), "df")))
IC

N = nrow(data)
IC$AICc <- IC$AIC + 2*IC$k*(IC$k+1)/(48-IC$k-1)
IC$AICc <- IC$AIC + 2*IC$k*(IC$k+1)/(N-IC$k-1)
IC

# Calculate model weights for AICc

AICcmin <- min(IC$AICc)
RL <- exp(-0.5*(IC$AICc - AICcmin))
sumRL <- sum(RL)
IC$AICcmin <- RL/sumRL

# Calculate model weights for BIC

BICmin <- min(IC$BIC)
RL.B <- exp(-0.5*(IC$BIC - BICmin))
sumRL.B <- sum(RL.B)
IC$BICew <- RL.B/sumRL.B

# Check it out

round(IC,3)

# 2.f Confidence intervals for predictors

ModelsREML <- list(Grinnellian=mod1, Bioclim=mod2)

# 95% CI for these estimates
confint(ModelsREML$Grinnellian, level = 0.95, method = "Wald")

confint(ModelsREML$Bioclim, level = 0.95, method = "Wald")



