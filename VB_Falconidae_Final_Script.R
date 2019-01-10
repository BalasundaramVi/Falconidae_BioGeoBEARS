# Falconidae C174 Final Project Script
# Created by Vignesh Balasundaram (504495938)
# 03/20/2018

# Set Working Directory
setwd("~/Documents/CS/R/Falconidae_BioGeoBears")

# Import Libraries
library(ape)
library(ggtree)
library(EBImage)
library(BAMMtools)
library(coda)
library(phytools)
library(geiger)
library(optimx)
library(FD)
library(snow)
library(parallel)
library(BioGeoBEARS)

##################################################
#####     PART I: CREATING THE PHYLOGENY     #####
##################################################

# Load Phylogeny of all Birds
load("bird.phylo")
zoom(bird.phylo, grep("Falconidae", bird.phylo$Family), cex = 0.5)

# Find node of Falconidae from bird.phylo
MRCA <- getMRCA(bird.phylo, tip = grep("Falconidae", bird.phylo$Family))
MRCA # Node is 16322
dev.off()

# Extract and plot Falconidae clade
Falconidae <- extract.clade(bird.phylo, node = 16322)
plot(Falconidae, cex = 0.4)
add.scale.bar()

# Label Falconidae species found in South America
CladeA <- getMRCA(Falconidae, tip = c("Phalcoboenus_australis", "Caracara_cheriway"))
CladeA # to find appropriate node
dev.off()

# Supplemental Figure B1: Falcanidae Clade with Highlighted Caracaras
PT <- ggtree(Falconidae)
PT <- PT %>% phylopic("6cebf754-cb71-448d-a5bb-947157205264", color ="black", alpha = .2) +
  geom_hilight(112, fill = "indianred2", alpha = 0.4) +
  geom_cladelabel(node = 112, label = "Caracaras", align = T, angle = 270, hjust = 'center')
print(PT)

# Create a .tre File with only the necessary tree in question (Falconidae)
write.tree(Falconidae, "FalconidaeTree.tre")

# Test to see if .tre file has been made correctly
# Supplemental Figure B2: Creation of Plot of Labeled South American Species
test_tree <- read.tree("FalconidaeTree.tre")
test_tree <- ladderize(test_tree)
SL <- read.csv("SpeciesList.csv")
test_tree$tip.label <- SL$?..Number
plot(test_tree, cex = 0.5, show.tip.label = TRUE, tip.col = c(rep("blue", 45), rep("red", 19)))
title("Plot of South American Species to Species in Other Continents")
legend(x = 0, y = 20, legend = c("South American", "Other"), 
       fill = c(rep("red", 1), rep("blue", 2)), cex = 1)
axisPhylo()
dev.off()

# Load Falconidae Tree
tree <- read.tree("FalconidaeTree.tre")
tree <- ladderize(tree)

## Finding Nodes for Specific Subclades: ##
# O_MRCA subclades of "Other" - Falconidae species inhabiting multiple regions
# A_MRCA subclade of Falconidae species inhabiting Asia (with one outlier in Africa)
# S1 + S2 MRCA of two South American subclades
O_MRCA <- getMRCA(tree, tip = c("Falco_biarmicus", "Falco_alopex"))
A_MRCA <- getMRCA(tree, tip = c("Microhierax_caerulescens", "Polihierax_semitorquatus"))
S1_MRCA <- getMRCA(tree, tip = c("Micrastur_gilvicollis", "Herpetotheres_cachinnans"))
S2_MRCA <- getMRCA(tree, tip = c("Phalcoboenus_australis", "Spiziapteryx_circumcincta"))

# Creation of Heatmap of Falconidae Phylogeny to Geographic Area
geodata <- read.delim("FalconGeo2.txt")
Africa <- geodata$Africa
Asia <- geodata$Asia
Europe <- geodata$Europe
SouthAmerica <- geodata$South.America
NorthAmerica <- geodata$North.Central.America
Australia <- geodata$Australia

names(Africa) <- GeoTree$tip.label # naming each row with correct species
g <- cbind(Africa, Asia, Europe, SouthAmerica, NorthAmerica, Australia)

# Figure 1: Visual Representation of Geographic Ranges of Falconidae Species
p <- ggtree(GeoTree) + geom_tiplab(size=3.5, align=TRUE, linesize=.5) + theme_tree2()
pp <- (p + scale_y_continuous(expand=c(0, 0.3))) %>%
  gheatmap(g, offset=8, width=0.6, colnames=F, high = "steelblue", low = "lightgrey") %>%
  scale_x_ggtree()
pp + theme(legend.position="right")


##################################################
#####          PART II: BAMM SCRIPT          #####
##################################################

# Import tree
dev.off()
tree <- read.tree("FalconidaeTree.tre")
tree <- ladderize(tree)
plot(tree, cex = 0.4)
axisPhylo()

# Preparing BAMM control file for Falconidae Tree
# My Windows XPS 15 Pro has 4 cores
is.ultrametric(tree)
is.binary(tree)
min(tree$edge.length)
# estimate priors and create control files
priors <- setBAMMpriors(tree, outfile = NULL)
generateControlFile(file = "Falconidaecontrolfile.txt", params = list(
  treefile = "FalconidaeTree.tre",
  globalSamplingFraction = "1", # 100% sampling
  seed = sample(1:1000000, 1),
  overwrite = "0",
  expectedNumberOfShifts = "1",
  lambdaInitPrior = as.numeric(priors["lambdaInitPrior"]),
  lambdaShiftPrior = as.numeric(priors["lambdaShiftPrior"]),
  muInitPrior = as.numeric(priors["muInitPrior"]),
  numberOfGenerations = "50000000", # if you have fewer cores, make this 5 mil instead of 10 mil
  mcmcWriteFreq = "1000",
  eventDataWriteFreq = "1000",
  printFreq = "1000",
  acceptanceResetFreq = "1000",
  outName = "Falconidae",
  numberOfChains = "4")) # set to number of CPUs

# Run all the BAMM commands on Command Prompt at this point

#####    COMMAND PROMPT CODE     #####
##
## cd \Users\Vignesh B\Documents\RSTUDIO\FinalProject\C174FinalReport
## bamm -c Falconidaecontrolfile.txt
##
##### END OF COMMAND PROMPT CODE #####

# Import BAMM results into R
edata <- getEventData(tree, eventdata = "Falconidae_event_data.txt", burnin = 0.1)
summary(edata)

# Check quality of BAMM results
# Assess MCMC convergence
# Supplemental Figure B3: MCMC Convergence
mcmcout <- read.csv("Falconidae_mcmc_out.txt")
plot(mcmcout$logLik ~ mcmcout$generation)

# Test for convergence of the MCMC chains with the coda package for R
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout),]
# Effective Size Post Burn of N_shifts and LogLik
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)

##### Analyze BAMM results with BAMMtools #####

# Analysis of rate shifts
computeBayesFactors("Falconidae_mcmc_out.txt", expectedNumberOfShifts = 1, burnin = 0.1)
# above output diagonal values (1~1, 2~2,..., 6~6) all equal 1.00

# Plot histograms of the prior and posterior probability distributions
# Figure 2: Plot of Posterior and Prior Probability Distributions
plotPrior("Falconidae_mcmc_out.txt", expectedNumberOfShifts = 1)

# Analysis of Rates
# Supplemental Figure B5: Mean Speciation Rate of Falconidae Tree
s <- plot.bammdata(edata, spex = "s", labels = F, font = 3, cex = 0.3)
title(main = "Mean speciation rate", sub = "time before present")
addBAMMlegend(s, location = "topleft", nTicks = 1)
axisPhylo()

# plot rate shift configurations and compare frequency of 
# alternative configurations in posterior distribution
css <- credibleShiftSet(edata, expectedNumberOfShifts = 1, threshold = 5, set.limit = 0.95)
css$number.distinct
summary(css)

sss <- plot.credibleshiftset(css, border = F)

dev.off() # reset the plotting dimensions specified by previous commands
plot.new() # double check the plotting dimensions and clear previous data
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts = 1)

# Figure 3: Best Shift Configuration of Falconidae from BAMM
ss <- plot.bammdata(best, labels = T, font = 3, cex = 0.4)
title(main = "Best shift configuration", sub = "time before present")
addBAMMlegend(ss, location = "left", nTicks = 1)
addBAMMshifts(best, cex = 3, pch = 1)
axisPhylo()

# plotting probability of branches containing rate shift
# the longer the branch, the higher the probability of a rate shift
dev.off()
par(font = 1)
marg_probs <- marginalShiftProbsTree(edata)
# Figure 6: Marginal Shift Probability of Falconidae separated into three major regions:
# 1. South American subclades in Red
# 2. Asia in Green (and one African species in black)
# 3. Species inhabiting multiple continental regions in Blue
plot.phylo(marg_probs, cex = 0.5, tip.color = 
             c(rep("steelblue", 38), rep("darkgreen",6), rep("black",1), rep("indianred2", 19)))
legend(x = 0, y = 23, legend = c("South America","Africa", "Asia", "Other Continents"),
       fill = c(rep("indianred2", 1),rep("black",1),rep("darkgreen", 1),rep("steelblue", 1)),cex = 1)
title(main = "Marginal shift probability")
add.scale.bar("bottomleft", font = 0.7)

# Falconidae global clade speciation and extinction rates (mean and variance)
global_rates <- getCladeRates(edata)
mean(global_rates$lambda) # Average rate of speciation of 0.166
mean(global_rates$mu) # Average rate of extinction of 0.0790
var(global_rates$lambda) # speciation variance of 0.00136
var(global_rates$mu) # extinction variance of 
# the speciation rate estimate is 0.17 new species per million years
quantile(global_rates$lambda, c(0.05, 0.95))
# there is a 90% probability that the speciation rate of this clade is between 0.117 and 0.235

# Individual Subclade Analysis of Rates:
# O = Falconidae clade with species inhabiting multiple regions, A = Asia/Africa, S = South America
O_rates <- getCladeRates(edata, node = O_MRCA)
A_rates <- getCladeRates(edata, node = A_MRCA)
S1_rates <- getCladeRates(edata, node = S1_MRCA)
S2_rates <- getCladeRates(edata, node = S2_MRCA)
# average of the two South American subclades
Sl_rates <- (S1_rates$lambda + S2_rates$lambda)/2
Sm_rates <- (S1_rates$mu + S2_rates$mu)/2

# Table 2: Results of Speciation Rates of Subclades by Geographic Region
## SPECIATION RATES ##
mean(O_rates$lambda) # Average rate of speciation of 0.168
mean(A_rates$lambda) # Average rate of speciation of 0.166
mean(Sl_rates) # Average rate of speciation of 0.165
var(O_rates$lambda) # speciation variance of 0.00134
var(A_rates$lambda) # speciation variance of 0.00141
var(Sl_rates) # speciation variance of 0.00147

## EXTINCTION RATES ##
mean(O_rates$mu) # Average rate of extinction of 0.0787
mean(A_rates$mu) # Average rate of extinction of 0.0790
mean(Sm_rates) # Average rate of extinction of 0.0795
var(O_rates$mu) # extinction variance of 0.00286
var(A_rates$mu) # extinction variance of 0.00287
var(Sm_rates) # extinction variance of 0.00289

mean(O_rates$mu) # Average rate of death of 0.0787
mean(A_rates$mu) # Average rate of death of 0.0790
mean(Sm_rates) # Average rate of death of 0.0795


# Rate-through-time analyses
dev.off()
par(font = 1)
# Figure 4: Rate-through-time analysis of entire Falconidae Tree
plotRateThroughTime(edata, ratetype = "speciation", avgCol = "black", intervalCol = "gray80", 
                    intervals = c(0.05, 0.95), opacity = 1)
# Figure 5: Side-by-side Rate-through-time plots of speciation of entire Falconidae clade (black)
# and speciation of Falcon species inhabiting multiple regions (blue)
plotRateThroughTime(edata, ratetype = "speciation", avgCol = "black", intervalCol = "gray40", 
                    node = O_MRCA, nodetype = 'exclude', intervals = c(0.05, 0.95), opacity = 0.3,
                    start.time = 15, xlim = c(15,0), ylim = c(0.155,0.165)) +
  plotRateThroughTime(edata, ratetype = "speciation", avgCol = "steelblue", intervalCol = "gray40",
                      intervals = c(0.05, 0.95), opacity = 0.3, node = O_MRCA, add = TRUE,
                      start.time = 15, xlim = c(15,0), ylim = c(0.155,0.165))

##################################################
#####        PART III: MCCR ANALYSIS         #####
##################################################

# Double check tree is still correct
tree <- read.tree("FalconidaeTree.tre")
fgamma <- ltt(tree, plot = FALSE)$gamma

# Find the gamma statistic for project tree
ltt(tree)
temp <- ltt(tree, log = TRUE)
lines(c(0, max(nodeHeights(tree))), c(log(2), log(length(tree$tip.label))), 
      lty = "dashed", lwd = 2, col = "red")

# Transparency function for log through time plot
makeTransparent<-function(someColor,alpha=10){
  newColor<-col2rgb(someColor)
  apply(newColor,2,function(curcoldata){
    rgb(red=curcoldata[1],green=curcoldata[2],blue=curcoldata[3],
        alpha=alpha,maxColorValue=255)
  })
}

# Supplemental Figure B6: Plot of Falconidae Tree and Lineage Through Time
obj<-ltt(tree)
plotTree(tree,color=makeTransparent("blue",alpha=50),ftype="off",add=TRUE, mar=par()$mar)
# Finding the gamma value of Falconidae
fgamma <- ltt(tree)$gamma # Falconidae gamma is 1.240405

# Find the age of the clade
age <- branching.times(tree)[1] # age at origin branching node (65) is 38.56723
# specify the total richness of the clade (including species not included in tree)
richness <- 64
# Estimate speciation rate for all species in clade
falconbirth = (log(richness) - log(2))/age
falconbirth # Estimated speciation rate for the clade (all species in tree)

# Monte Carlo Constant Rates (MCCR) Test
# simulate gamma values when trees are undersampled
num_simulations<-100 #number of simulations
g1_null<-numeric(num_simulations) #g1_null will hold the simulated gamma values for pruned trees
for(i in 1:num_simulations) {
  sim.bdtree(falconbirth, d=0, stop = "taxa", n=64)->sim_tree 
  drop.random(sim_tree, 0)->prune # Here we drop 0 species RANDOMLY from the tree
  gammaStat(prune)->g1_null[i] # this stores the gamma values from the pruned trees
}
dev.off()

# Figure 7: Histogram of null distribution of gamma statistics 
# wtih Falconidae gamma marked with red arrow

# create a histogram of the null distribution
hist(g1_null, xlim = c(-3.5, 3.5))
#arrow indicates where the observed gamma falls in the null just generated
arrows(fgamma, 40, fgamma, 0, col="red", lwd=2) 

# Table 1: Results of teh MCCR Analysis of Falconidae
fgamma <- ltt(tree)$gamma # Falconidae gamma is 1.240405
mean(g1_null) # mean of the null distribution g1 is -0.1385018
# Calculating the p-value through count of all values falling under Falconidae gamma
smallerNull<-g1_null<=fgamma
count<-sum(smallerNull)
mccr_pval<-(count+1)/(num_simulations+1)
mccr_pval # p-value of 0.9306931

##################################################
##### PART IV: Analysis Through BioGeoBEARS  #####
##################################################
dev.off()
plot.new()

# Biogeographical Models of Falconidae Ancestral Biogeographic States
# Figure 8: BAYAREALIKE +J Model
source("FalconidaeBioGeoBEARS.R")
# Ttable 3: AIC Results from BioGeoBEARS
read.table("restable_AIC_rellike_formatted.txt")

