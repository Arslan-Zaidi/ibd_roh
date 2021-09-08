
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(c("-a", "--sibs"), type = "double", help = "sibs"),
  make_option(c("-b", "--cousin1"), type = "double", help = "1st cousin prop"),
  make_option(c("-c", "--cousin2"), type = "double", help = "2nd cousin prop"),
  make_option(c("-d", "--cousin3"), type = "double", help = "3rd cousin prop"),

  make_option(c("-t", "--texact"), type = "double", default = 10, help = "max_T_exact_model"),
  make_option(c("-s", "--kinship_scale"), type = "double", default = 1, help = "scaling factor for kinship"),
  make_option(c("-k", "--okinship"), type = "double", default = 0, help = "old kinship (before t)"),
  make_option(c("-p", "--pop"), type = "integer", default = 1, help = "which population to use (1 = pathan, 2 = jatt)"),
  make_option(c("-n", "--ne"), type = "integer", default = 1, help = "ne trajectory (1,2,3) or constant ne(>3)"),
  make_option(c("-o", "--output", type = "character", help = "path to output directory", default = "output/"))
)

args <- parse_args(OptionParser(option_list = option_list))

#rm(list=ls())
#setwd("/lustre/scratch115/projects/autozyg/BiB/code_for_popgen_paper")
min.cm=5
max.cm=30

###load("Ne_trajectories_from_IBDNe_spliced_onto_SMCpp.self_declared_groups.before_brute_for_filtering_on_IBDseq.RData")
### load a table of consanguinity rates estimated from self-reported data in the Bradford Pakistani mothers, split by homogeneous sub-group defined in FineStructure --> this will be used to get a naive estimate of kinship
#load("consanguinity_rates_estimated_naively_from_self_reported_information.homogenous_groups_from_fineSTRUCTURE.RData")
#rownames(consang.by.bir2.robust)=gsub("Jatt/Choudhary","Jatt/Choudhry",rownames(consang.by.bir2.robust))

### load IBDNe results for these homogeneous sub-groups --> these will be used as a model for the number of mating pairs over the last 50 generations
#ibdne.new=read.delim('IBDNe_HomogeneousGroups_from_fineSTRUCTURE.with_brute_force_plus_HLA_plus_centromere_filtering.txt',header=T, as.is=T)
ibdne.new = read.delim('/Users/Azaidi/endogamy/data/BiB_CoreExome_IBDNe_different_FROH_filters_Dec20.txt', header = T, as.is = T)

#args <- commandArgs(trailingOnly = TRUE)

#rates vector (sibs, 1st cousins, 2nd cousins, 3rd cousins)
#rates = as.numeric(args[ c(2 : length(args) ) ])
rates = c(args$sibs, args$cousin1, args$cousin2, args$cousin3)
rates_name = paste(rates, collapse="_")

# Specify the number of generations during which we run the exact model
#max_T_exact_model = as.numeric(args[1])
#max_T_exact_model = 10
max_T_exact_model = args$texact

#specify the kinship prior to max_T_exact_model
#old.kinship = as.numeric(args[2])
#old.kinship = 0
old.kinship = args$okinship

#specify which model of recent kinship you want to use
#which.recent.consang=as.numeric(args[3])
#which.recent.consang = 1 #scaling factor to increase kinship
which.recent.consang = args$kinship_scale

#define which Ne trajectory you will use
#which.ne.trajectory=as.numeric(args[4])
which.ne.trajectory = args$ne

#define which population you are going to model
#which.pop=as.numeric(args[1])
#which.pop = 1 #Jatt
which.pop = args$pop

# This is the maximal generation to track when computing the distribution of the TMRCA
max_T_coal = 45055 # Just need to be large. We set this to be the number of generations back in time that the SMC++ estimates of population size in the HGDP paper go.

### Specify which population is being modelled and pull out the IBDNe estimates for that population
if(which.pop==1){
pop="PATHAN"
pop.label="PATHAN_in_Cluster_8"
ibdne.pop=ibdne.new[ibdne.new$pop=="Pathan" & ibdne.new$GEN<50,]
ibdne.pop$t=ibdne.pop$GEN+1
} else {
pop="Jatt/Choudhry"
pop.label="JattChoudhry_in_cluster_10"
ibdne.pop=ibdne.new[ibdne.new$pop=="Choudhry_Jatt" & ibdne.new$GEN<50,]
ibdne.pop$t=ibdne.pop$GEN+1
}

### We will define a data frame to store the historical Ne back in time
historic.ne=list()
historic.ne[[pop]]=data.frame(t=1:max_T_coal,Ne=NA)

### Specify the Ne trajectory to be used
if(which.ne.trajectory==1){
    #use the IBDNe estimate for last 50 generations then constant
    ne.label = "Ne_trajectory_from_IBDNe_for_50_gens_then_constant"
    historic.ne[[pop]][1:50,2] = ibdne.pop[1:50,"NE"]
    historic.ne[[pop]][50:nrow(historic.ne[[pop]]),2] = ibdne.pop[50,"NE"]
}
if(which.ne.trajectory==2){
 #use the lower bound of the 95% CI of the IBDNe estimate for last 50 generations then constant
  ne.label = "Ne_trajectory_from_IBDNe_LWR_CI_for_50_gens_then_constant"
  historic.ne[[pop]][1:50,2] = ibdne.pop[1:50,"LWR.95.CI"]
  historic.ne[[pop]][50:nrow(historic.ne[[pop]]),2] =	ibdne.pop[50,"LWR.95.CI"]
}
if(which.ne.trajectory==3){
 #use the upper bound of the 95% CI of the IBDNe estimate for last 50 generations then constant
  ne.label = "Ne_trajectory_from_IBDNe_UPR_CI_for_50_gens_then_constant"
  historic.ne[[pop]][1:50,2] = ibdne.pop[1:50,"UPR.95.CI"]
  historic.ne[[pop]][50:nrow(historic.ne[[pop]]),2] = ibdne.pop[50,"UPR.95.CI"]
}
if(which.ne.trajectory>3){
### use a constant value for t>0
    #constant.ne.values=c(1e3,5e3,10e3,20e3,30e3,40e3,50e3,100e3,6e3,7e3,8e3,9e3) ### can change this to include any constant Ne value you want
    #my.constant.ne=constant.ne.values[which.ne.trajectory-3]
    #historic.ne[[pop]][,2]=my.constant.ne
    historic.ne[[pop]][,2] = which.ne.trajectory
    ne.label=paste("constant_Ne_",which.ne.trajectory/1e3,"k",sep="")
}

### set the average kinship between spouses to be used for t<=max_T_exact_model
# "rates" is a vector of probabilities of the relationship between a husband and a wife (a mating pair). The first entry in the vector is for sibling matings - we will assume 0.
# which.recent.consang is just a scaling factor we will multiply by the empirical estimates of consanguinity rates from the self-reported data, to increase the kinship value.
#rates = c(0,consang.by.bir2.robust[pop,] * which.recent.consang)
rates = rates * which.recent.consang
n = length(rates)
kinship = sum(rates * 4^(-(1:n)))

### set the average kinship between spouses for t>max_T_exact_model
# This may have been specified on the commandline but if it was given as -1, this step just sets it to the same value to be used for  t<=max_T_exact_model
if(old.kinship<0){
    old.kinship=kinship
}

### Make a label for the output file
if(which.recent.consang==1){
    label=paste(pop.label,".",ne.label,".empirical_consang_rates_for_",max_T_exact_model,"_gens",sep="")
} else {
    label=paste(pop.label,".",ne.label,".mean_kinship_",kinship,"_for_",max_T_exact_model,"_gens",sep="")
}
label=paste(label,".then_mean_kinship_",old.kinship,sep="")

cat("\nmax_T_exact_model = ",max_T_exact_model,"; recent kinship = ",kinship,"; old.kinship=",old.kinship,"\n")

# function specifying the average kinship between spouses over time
get_q = function(t,old.kinship){
if (t<=max_T_exact_model)
  return(kinship)
else
  return(old.kinship)
}


# function specifying the population size over time --> this just extracts it from the data frame we created above as historic.ne[[pop]]
get_N = function(t,historic.ne.pop){
  return(historic.ne.pop[historic.ne.pop$t==t,"Ne"])
}

### Run code to calculate transition matrix, run Markov chain to get TMRC distribution for the last max_T_exact_model generations,
#then use exponential distribution to get TMRCA distribution prior to that, splice these two parts back together,
#then calculate the probability of an IBD or ROH segment being in a given length bin

# First read the code containing the calc_ell_dist() function.
source("LengthDistribution.flexible.R")

#specify the bins for length in cM over which you want to compute the statistics
ell_seq = seq(5,30,by=1)

#Run the code
if(which.ne.trajectory>3){
save.output=calc_ell_dist(ell_seq,get_N,rates,get_q,max_T_coal,max_T_exact_model,num_runs=100000,old.kinship,ibdne=FALSE)
} else {
save.output=calc_ell_dist(ell_seq,get_N,rates,get_q,max_T_coal,max_T_exact_model,num_runs=100000,old.kinship,ibdne=TRUE)
}

# Add some other variables to the saved output in case you want to use them later.
save.output[["old.kinship"]]=old.kinship
save.output[["ells"]]=ell_seq
save.output[["N_func"]]=get_N
save.output[["max_T_coal"]]=max_T_coal
save.output[["max_T_exact_model"]]= max_T_exact_model

# Save the output
save(save.output,file=paste(args$output,"/output.",label,"_r",rates_name,".RData",sep=""))
