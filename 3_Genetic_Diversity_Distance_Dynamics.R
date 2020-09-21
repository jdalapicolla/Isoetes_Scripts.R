###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

###########   STEP 03: GENETIC DIVERSITY, DISTANCE, AND TAJIMA D   ############

#Script prepared by Jeronymo Dalapicolla and Rodolfo Jaffé for the manuscript:
# Dalapicolla et al. - Conservation implications genetic structure in the narrowest endemic quillwort from the Eastern Amazon


########## PRE-ANALYSIS AND INFORMATION ----

##A. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

##B. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP.
source("functions_LanGen.R")


##C. INSTALL THE PACKAGES
if("remotes" %in% rownames(installed.packages()) == FALSE){install.packages("remotes")
} else {print (paste0("'remotes' has already been installed in library"))}
if("BiocManager" %in% rownames(installed.packages()) == FALSE){install.packages("BiocManager")
} else {print (paste0("'BiocManager' has already been installed in library"))}
if("pacman" %in% rownames(installed.packages()) == FALSE){install.packages("pacman")
} else {print (paste0("'pacman' has already been installed in library"))}
if("devtools" %in% rownames(installed.packages()) == FALSE){install.packages("devtools")
} else {print (paste0("'devtools' has already been installed in library"))}

if("r2vcftools" %in% rownames(installed.packages()) == FALSE){remotes::install_github("nspope/r2vcftools")
} else {print (paste0("'r2vcftools' has already been installed in library"))}
if("LEA" %in% rownames(installed.packages()) == FALSE){BiocManager::install("LEA")
} else {print (paste0("'LEA' has already been installed in library"))}
if("vcfR" %in% rownames(installed.packages()) == FALSE){install.packages("vcfR")
} else {print (paste0("'vcfR' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("reshape2" %in% rownames(installed.packages()) == FALSE){install.packages("reshape2")
} else {print (paste0("'reshape2' has already been installed in library"))}
if("vegan" %in% rownames(installed.packages()) == FALSE){install.packages("vegan")
} else {print (paste0("'vegan' has already been installed in library"))}
if("mmod" %in% rownames(installed.packages()) == FALSE){install.packages("mmod")
} else {print (paste0("'mmod' has already been installed in library"))}
if("poppr" %in% rownames(installed.packages()) == FALSE){install.packages("poppr")
} else {print (paste0("'poppr' has already been installed in library"))}
if("ecodist" %in% rownames(installed.packages()) == FALSE){install.packages("ecodist")
} else {print (paste0("'ecodist' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("strataG" %in% rownames(installed.packages()) == FALSE){devtools::install_github('ericarcher/strataG', build_vignettes = TRUE, force = T)
} else {print (paste0("'strataG' has already been installed in library"))}

#C. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vegan, ecodist, vcfR, adegenet, poppr, mmod, reshape2, ggplot2, dartR, strataG)

##D. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Diversity", "./Results_Distance", "./Results_NeEstimator"))



###### 1. Loading Files ------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "isoetes"

###1.2. LOAD VCF FILES: 
#A. Load neutral .vcf file with geographical information and genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps_neutral)
names(snps_neutral@meta) #verify col names in metafile

  #B. Number of cluster and method:
optimal_K = 1
method = "sNMF"

#B. Position of samples by populion by sNMF approach. Choose one method and change it on script:
for (i in 1:length(unique(snps_neutral@meta$PopID_snmf))){
  pop = which(snps_neutral@meta$PopID_snmf == i)
  assign(paste0("pop_sNMF_", i), pop)
}


###1.3. VERIFY THE MEAN COVERAGE DEPTH: 
#A. Mean coverage in all dataset
site.depth = Query(snps_neutral, type="site-mean-depth")
summary(site.depth$MEAN_DEPTH) #Mean = 83.81 / Median = 71.49
hist(site.depth$MEAN_DEPTH, breaks=30)

#B. Mean coverage by individual
coverage_ind = c()
for (p in 1:length(snps_neutral@sample_id)){
beta = Query(Subset(snps_neutral, samples = p), type="site-mean-depth")
coverage_ind[p] = mean(beta$MEAN_DEPTH, na.rm = T)}
#verify
coverage_ind
#save as metafile
snps_neutral@meta$coverage = coverage_ind

write.csv(coverage_ind, "Results_Metafiles/coverage.csv")






###### 2. Genetic Diversity ---------------
###2.1. ESTIMATE GENETIC DIVERSITY:
#A. By species/taxon (all samples together)
Overall = GenDiv(snps_neutral)
write.csv(Overall, file=paste0("Results_Diversity/Diversity_Overall_neutral_", project_name, ".csv"))

#B. By individual
ind = GenDiv_IND(snps_neutral)
write.csv(ind, file=paste0("Results_Diversity/Diversity_Individual_neutral_", project_name , ".csv"))








###### 3. Genetic Distance by Individual ------
###3.1. RELATEDNESS IN ALL SAMPLES
#Yang's Relatedness: distance values among all pairs are calculated once (below diagonal) with distance within individuals
#A. Calculate Relatedness
REL_YANG = Relatedness(snps_neutral, type = "yang",verbose = TRUE)
head(REL_YANG)
nrow(REL_YANG) 
colnames(REL_YANG)<-c("INDV1","INDV2","RELATEDNESS_AJK_Yang")

#B. Save results:
write.csv(REL_YANG, file=paste0("Results_Distance/Yang_Reletedness_IND_neutral_", project_name, ".csv"))


###3.2. Dps DISTANCE
#A. Convert VCF to genind
snps =  read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_LEA_DAPC_TESS.vcf"), verbose = T)

#B. Convert data into genind
geninddata = vcfR2genind(snps)
geninddata@tab[1:10,1:10]

#C. Calculate proportions of shared alleles between pairs of individuals
dps_shared = propShared(geninddata)
head(dps_shared)

#D. Save results:
write.csv(dps_shared, file=paste0("Results_Distance/Dps_Shared_IND_neutral_", project_name, ".csv"))


###3.3. DISTANCE HISTOGRAM
#A. Values for Dps:
x = dps_shared[lower.tri(dps_shared, diag = FALSE)]
length(x)

#B. Exclude Relatedness between same individual and create a copy of REL data frame
REL2 = REL_YANG[REL_YANG$INDV1 != REL_YANG$INDV2,]
#verify
head(REL2)
length(REL2[,1]) #1485 comparisons
length(REL_YANG[,1])
length(REL_YANG[,1]) - length(REL2[,1])
#define a dataframe 
df = REL2[,-1]
df[,1] = x
dim(df)
colnames(df)<-c("Dps","Rel")
df2 = melt(df)
#verify
summary(df)
sd(df$Rel)

#C. Theme for graphs
theme_2020 = theme(legend.text = element_text(face = "italic",
                                              colour="black",
                                              family = "Helvetica",
                                              size = rel(1)), 
                   axis.title = element_text(colour="black",
                                             face = "bold",
                                             family = "Helvetica",
                                             size = rel(1.2)), 
                   axis.text = element_text(family = "Helvetica",
                                            colour = "black",
                                            angle = 0,
                                            hjust = 0.5,
                                            vjust = 0.5,
                                            size = rel(1)), 
                   axis.line = element_line(size = 1,colour = "black"), 
                   axis.ticks = element_line(colour="black",size = rel(1)),
                   
                   panel.grid.minor = element_blank(), 
                   panel.background = element_rect(fill = "whitesmoke"), 
                   panel.grid.major = element_line(colour="black",size = rel(0.2), linetype = "dotted"),
                   legend.key = element_blank(), 
                   legend.title = element_text(colour = "black",
                                               face = "bold",
                                               size = rel(1.5),
                                               hjust = 0.5,
                                               vjust = 0.5,
                                               family = "Helvetica"), 
                   plot.title = element_text(colour = "black",
                                             face = "bold",
                                             hjust = 0.5, #alingment
                                             size = rel(1.7),
                                             family = "Helvetica"))

#D. Plot the graph
graph2 = ggplot(data=df2, mapping = aes(x= value, fill = variable)) +
  geom_histogram(binwidth = 0.0005) +
  facet_wrap(~variable, nrow = 2, scales = "free") +
  theme_2020 +
  theme(legend.position="none")

#E. Save as pdf:
pdf("histogram_related_dps.pdf", onefile = F)
plot.new()
graph2
dev.off()








####### 5. Converting VCF to Genepop to run NeEstimator -----
###5.1. LOAD VCF
snpsR = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral.vcf"), verbose = T)


###5.2. DEFINING POPULATIONS USING THE VCF METAFILE:
snps = vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID=T)
VCFsummary(snps)
population = as.factor(snps@meta$PopID_snmf)


##5.3. CONVERTING FILES:
#A. VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#B. Adding strata (pops) into Genind
snps_genind@pop = population

#C. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

#D. Converting Gtypes to GENEPOP to run in NeEstimator and save it:
genepopWrite(snps_gtypes, "isoetes_genepop_regular.txt")

#-----------------------------------------------

##5.4. CONVERTING FILES WITHOUT MISSING DATA:
#A. VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#B. Adding strata (pops) into Genind
snps_genind@pop = population

#C. Removing Missing Data
snps_genind = missingno(snps_genind, type = "loci", cutoff = 0)

#D. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(snps_genind)
class(snps_gtypes)

#E. Converting Gtypes to GENEPOP to run in NeEstimator and save it:
genepopWrite(snps_gtypes, "isoetes_genepop_NOMISSING.txt")

#-----------------------------------------------

##5.4. CONVERTING FILES WITHOUT MISSING DATA RESAMPLING BY SNPS NUMBER:
#A. VCF to Genind
snps_genind = vcfR2genind(snpsR)
class(snps_genind)

#B. Adding strata (pops) into Genind
snps_genind@pop = population

#C. Removing Missing Data
snps_genind = missingno(snps_genind, type = "loci", cutoff = 0)

#D. Genind to Genlight
gl = gi2gl(snps_genind)

#E. Create a subset of snps or individuals
#by loci number. 30,000 random SNPs
set.seed(123)
index.loc = sample(nLoc(gl), 30000, replace = F)
index.loc

#F. Subsetting
glsub = gl[,index.loc] # only snps

#G. Convert to GENEIND
genind = gl2gi(glsub)
genind

#H. Converting Genind to Gtypes
snps_gtypes = genind2gtypes(genind)
class(snps_gtypes)

#I. Converting Gtypes to GENEPOP to run in NeEstimator and save it:
genepopWrite(snps_gtypes, "isoetes_genepop_30000.txt")

#J. Load this file in NeEstimator to run Ne analyses.

##END