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
if("gghighlight" %in% rownames(installed.packages()) == FALSE){install.packages("gghighlight")
} else {print (paste0("'gghighlight' has already been installed in library"))}
if("related" %in% rownames(installed.packages()) == FALSE){install.packages("related", repos="http://R-Forge.R-project.org")
} else {print (paste0("'related' has already been installed in library"))}


#C. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vegan, ecodist, vcfR, adegenet, poppr, mmod, reshape2, ggplot2, dartR, strataG, gghighlight, related)

##D. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_Diversity", "./Results_Distance", "./Results_NeEstimator"))



###### 1. Loading Files ------
###1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "isoetes"

###1.2. LOAD VCF FILES: 
#A. Load neutral .vcf file with geographical information and genetic clusters ID:
snps_neutral = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf"), overwriteID=T)
VCFsummary(snps_neutral) #55 individuals and 35638 SNPs.
names(snps_neutral@meta) #verify col names in metafile


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


#C. Graphs
# add results from individual genetic metrics
ind_dataset = cbind(snps_neutral@meta, ind) 
names(ind_dataset)

#create a df for graphs
df = melt(ind_dataset[, c(3,11,12,13)],id.vars=c("local_ID"))
df
facets_ind = c("He_O","He_E", "FIS")

#graphs
he0 = ggplot(df[df$variable %in% facets_ind[1],], aes(value, fill = local_ID))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ local_ID, nrow = 1) +
  scale_fill_manual (values = c("grey", "red", "blue", "yellow", "green"),
                     labels = c("Center", "East", "North", "West", "South")) +
  theme_bw()

he = ggplot(df[df$variable %in% facets_ind[2],], aes(value, fill = local_ID))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ local_ID, nrow = 1) +
  scale_fill_manual (values = c("grey", "red", "blue", "yellow", "green"),
                     labels = c("Center", "East", "North", "West", "South")) +
  theme_bw()

fis = ggplot(df[df$variable %in% facets_ind[3],], aes(value, fill = local_ID))+
  geom_histogram(bins = 30, color="black") +
  gghighlight() +
  facet_wrap(~ local_ID, nrow = 1) +
  scale_fill_manual (values = c("grey", "red", "blue", "yellow", "green"),
                     labels = c("Center", "East", "North", "West", "South")) +
  theme_bw()

#verify
he0
he
fis

#save
pdf("./Histogram_IND_He0.pdf", onefile =F, width = 7, height =2 )
he0
dev.off()

pdf("./Histogram_IND_He.pdf", onefile =F, width = 7, height =2)
he
dev.off()

pdf("./Histogram_IND_FIS.pdf", onefile =F, width = 7, height =2)
fis
dev.off()


###2.2. COMPARING GENETIC DIVERSITY AMONG AREAS:
#A. Number of areas:
optimal_K = 5
areas = unique(snps_neutral@meta$local_ID)

#B. Position of samples by area:
for (i in areas){
  pop = which(snps_neutral@meta$local_ID == i)
  assign(paste0("areas_", i), pop)
}

#C. By cluster. ONLY IF YOUR K >= 2
#define a list of positions for each cluster
clusters = list(areas_centro, areas_L, areas_N, areas_O, areas_S)
#subset
for (i in 1:optimal_K){
  UNIND = snps_neutral@sample_id[clusters[[i]]]
  pop = Subset(snps_neutral, samples=UNIND)
  assign(paste0("pop_", i), pop)
}


iso_diver1 = GenDiv_IND(pop_1)
iso_diver1$SPE = c("center")
iso_diver2 = GenDiv_IND(pop_2)
iso_diver2$SPE = c("east")
iso_diver3 = GenDiv_IND(pop_3)
iso_diver3$SPE = c("north")
iso_diver4 = GenDiv_IND(pop_4)
iso_diver4$SPE = c("west")
iso_diver5 = GenDiv_IND(pop_5)
iso_diver5$SPE = c("south")



table_diversity = rbind(iso_diver1, iso_diver2, iso_diver3, iso_diver4, iso_diver5)
table_diversity

#create a data frame to analysis
variaveis = table_diversity[, 1:4]
species = as.factor(table_diversity[ ,5])
plan_anova = data.frame(species,variaveis)
classes = as.factor(plan_anova$species)

#calculate ANOVA
sink("ANOVA_TUKEYS.doc")
for(i in 2:ncol(plan_anova)){
  column<- names (plan_anova[i])
  result_anova<- aov(plan_anova[,i]~classes, data= plan_anova)
  result_anova2<- summary(aov(plan_anova[,i]~classes, data= plan_anova))
  tk<- TukeyHSD (result_anova)
  df <- tibble::rownames_to_column(as.data.frame(tk$classes), "Areas")
  graph = ggplot(data=df, aes(y=Areas, x=diff, xmin=lwr, xmax=upr))+
    geom_point() +
    geom_errorbarh(height=.2) +
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=1)+
    theme_bw()
  pdf(paste0("Tukey_diver_", names(plan_anova[i]) ,".pdf"))
  plot(graph)
  dev.off()
  }
sink()

#calculate P-value
sink("ANOVA_p-values.doc")
for(i in 2:ncol(plan_anova)){
  column<- names (plan_anova[i])
  result_anova<- summary(aov(plan_anova[,i]~classes, data= plan_anova))[[1]][["Pr(>F)"]]
  print(column)
  print(result_anova) 
}
sink()




###### 3. Genetic Distance by Individual ------
###3.1. RELATEDNESS IN ALL SAMPLES
#Yang's Relatedness:
#A. Calculate Relatedness
REL_YANG = Relatedness(snps_neutral, type = "yang",verbose = TRUE)
head(REL_YANG)
nrow(REL_YANG) #1540 comparisons
colnames(REL_YANG)<-c("INDV1","INDV2","RELATEDNESS_AJK_Yang")

#B. Exclude Relatedness between same individual
REL_YANG = REL_YANG[REL_YANG$INDV1 != REL_YANG$INDV2,]
#verify
head(REL_YANG)
nrow(REL_YANG) #1485 comparisons

#C. Save results:
write.csv(REL_YANG, file=paste0("Results_Distance/Yang_Reletedness_IND_neutral_", project_name, ".csv"))



#D. Plot the graph
graph2 = ggplot(data=REL_YANG, mapping = aes(x= RELATEDNESS_AJK_Yang)) +
  geom_histogram(binwidth = 0.0005) +
  theme_bw()

graph2

#E. Save as pdf:
pdf("histogram_related.pdf", onefile = F)
graph2
dev.off()



####### 5. Converting VCF to Genepop to run NeEstimator -----
###5.1. LOAD VCF
snpsR = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf"), verbose = T)


###5.2. DEFINING POPULATIONS USING THE VCF METAFILE:
snps = vcfLink(paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf"), overwriteID=T)
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
set.seed(3455)
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
set.seed(3455)
index.loc = sample(nLoc(gl), 500, replace = F)
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
genepopWrite(snps_gtypes, "isoetes_genepop_500.txt")
