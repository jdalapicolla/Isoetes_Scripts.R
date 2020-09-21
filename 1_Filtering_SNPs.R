###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

########################## STEP 01: FILTERING SNPS ############################

#Script prepared by Jeronymo Dalapicolla and Rodolfo Jaffé for the manuscript:
# Dalapicolla et al. - Conservation implications genetic structure in the narrowest endemic quillwort from the Eastern Amazon


########## PRE-ANALYSIS AND INFORMATION ----

##A. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

##B. DOWNLOAD .vcf FILES USING THE LINK IN THE MANUSCRIPT.

##C. TABLE S1 IS THE "coords" FILE USED HERE.

##D. INSTALL AND LOAD THE PACKAGES
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
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")
} else {print (paste0("'tidyverse' has already been installed in library"))}
if("vcfR" %in% rownames(installed.packages()) == FALSE){install.packages("vcfR")
} else {print (paste0("'vcfR' has already been installed in library"))}
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("poppr" %in% rownames(installed.packages()) == FALSE){install.packages("poppr")
} else {print (paste0("'poppr' has already been installed in library"))}
if("qvalue" %in% rownames(installed.packages()) == FALSE){BiocManager::install("qvalue")
} else {print (paste0("'qvalue' has already been installed in library"))}
if("psych" %in% rownames(installed.packages()) == FALSE){install.packages("psych")
} else {print (paste0("'psych' has already been installed in library"))}

#E. Load multiple packages using the package 'pacman'. If the package is missing "p_load" will download it from CRAN. "" in packages names is not mandatory.
pacman::p_load(r2vcftools, tidyverse, vcfR, dartR, poppr, adegenet, LEA, qvalue, psych)

##F. LOAD AUXILIARY FUNCTIONS. DOWNLOAD IT FROM THE GITHUB INFORMED IN MANUSCRIPT
source("functions_LanGen.R")

##G. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./adapt_var_mapping/Filtering", "./Results_Metafiles", "./Results_Filters/sNMF_FST"))






###### 1. Load Files and Verify Data Quality ----
###1.1. CHOOSE A NAME FOR THE PROJECT AND THE PATH FOR VCF/GVCF AND COORDS FILES:
#A. vcf/gvcf file path:
vcf_file = "vcf/isoetes_original_30k.vcf"

#B. Project name:
project_name = "isoetes"

#C. Coordenates file path and load it:
coord_file = "coords/isoetes_coords_certo.csv"
#load the file
coords = read.csv(coord_file)
#verify if it loaded correctly
head(coords)
tail(coords)

#D. Define the columns positions for localities ID (local), sample ID (sample_ID), longitude, and latitude. Add more columns if you need.
local = 7
sample_ID = 3
longitude = 1
latitude = 2

###1.2. LOAD THE VCF FILE
#A. VCF
snps_raw = vcfLink(vcf_file, overwriteID=T)
VCFsummary(snps_raw)
raw_SZ = capture.output(VCFsummary(snps_raw)) 








###### 2. Filtering Neutral Loci - By Quality -----
###2.1. FILTERING BY QUALITY
#A. Filter for neutral SNPs
snps = Filter(snps_raw, filterOptions(max.alleles=2, min.alleles=2, minQ=30, max.missing = 0.9, maf=0.05, min.meanDP=50, max.meanDP=200, hwe=0.0001), indels="remove")
VCFsummary(snps)
filtered_SZ = capture.output(VCFsummary(snps))

#B. Define R² value:
snps_fil_hwe = snps
r2 = 0.4

###2.2. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) WITHIN CONTIGS:
## If you have more than one population, you can use this to identify SNPs deviating from HW equilibrium within each population, and then removes those SNPs that are in desequilibrium in all populations. You just need to subset your samples:

#A. remove snps with R² value
ld_within <- Linkage(snps_fil_hwe, type="geno-r2", linkageOptions(min.r2=r2))
#ld_within <- read.csv(paste0("adapt_var_mapping/Filtering/ld_within_",r2, "_hwe_test2.csv")) #load this file if it was saved before
head(ld_within)
hist(ld_within$R.2)
write.csv(ld_within, file=paste0("adapt_var_mapping/Filtering/ld_within_",r2, "_hwe_test2.csv"))

#B. Select one set of the correlated snps (ID1 or ID2)
ld_snps <- ld_within$ID1
nold_snps <- snps_fil_hwe@site_id[!(snps_fil_hwe@site_id %in% ld_snps)] 
snps_fil_ld <- Subset(snps_fil_hwe, sites=nold_snps) # Keep snps that are not in LD.
neutralLDWithin_SZ = capture.output(VCFsummary(snps_fil_ld))
neutralLDWithin_SZ ##277 individuals and 11104 SNPs.

###2.3. FILTER DATASET BY LINKAGE DISEQUILIBRIUM (LD) BETWEEN CONTIGS:
#A. remove snps with R²
ld_between <- Linkage(snps_fil_ld, type="interchrom-geno-r2", linkageOptions(min.r2=r2)) 
#ld_between <- read.csv(paste0("adapt_var_mapping/Filtering/ld_between_", r2, "_hwe_test2.csv")) #load this file if it was saved before
head(ld_between)
hist(ld_between$R.2)
write.csv(ld_between, file= paste0("adapt_var_mapping/Filtering/ld_between_", r2, "_hwe_test2.csv"))

#B. Select one set of the correlated snps (ID1 or ID2)
ld2_snps <- ld_between$ID1
nold2_snps <- snps_fil_ld@site_id[!(snps_fil_ld@site_id %in% ld2_snps)]
snps_fil_ldF <- Subset(snps_fil_ld, sites=nold2_snps) # Keep snps that are not in LD.
neutralLDBetween_SZ = capture.output(VCFsummary(snps_fil_ldF)) 
neutralLDBetween_SZ ##277 individuals and 6602 SNPs.

###2.4. SAVE THE .VCF FILE WITH ONLY NEUTRAL SNPS:
Save(snps_fil_ldF, paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"))







###### 3. Adding Geographical and other Information in vcf's Metafile -----
###3.1.  ESTIMATING MISSING DATA
#A. Missing per locus
Missing <- apply(GenotypeMatrix(snps_fil_ldF), 2, function(x) sum(x < 0)/length(x)*100)
summary(Missing) ## Max missing = 9.09%
hist(Missing)

#B. Missing per individual
Missing_ind <- apply(GenotypeMatrix(snps_fil_ldF),1, function(x) sum(x<0)/length(x)*100)
summary(Missing_ind) ## Max missing = 0.402%
hist(Missing_ind)


###3.2. ADD GEOGRAPHICAL INFORMATION AND MISSING DATA BY INDIVIDUALS TO A VCF FILE:
#A. Save missind data information as data frame
Missingind.df = as.data.frame(Missing_ind)
Missingind.df$ID = row.names(Missingind.df)
Missingind.df

#B. Remove geographical information of individuals that you excluded in step #3 
missing_coord = coords[coords[,sample_ID] %in% Missingind.df$ID,]
head(missing_coord)
tail(missing_coord)
length(missing_coord[,1])

#C. Sort the geographical information file according to the sample order of the missing data file:
missing_coord_sorted = missing_coord[order(match(missing_coord[,sample_ID], Missingind.df$ID)),]
head(missing_coord_sorted)

#D. Check if there are some difference between files. If they are corrected, function will return "character (0)"
setdiff(Missingind.df$ID, missing_coord_sorted[,sample_ID])

#E. Check if identical samples were selected in the same order. If they are corrected, function will return "TRUE"
identical(as.character(Missingind.df$ID), as.character(missing_coord_sorted[,sample_ID]))

#F. Create a data frame merging missing and geographical information
missing_local = cbind(missing_coord_sorted[,c(local, longitude, latitude)], Missingind.df)
head(missing_local)
tail(missing_local)

#G. Add coords and missing information to vcf file. Be careful to the order!
#verify order. If they are right, function will return "TRUE"
identical(as.character(missing_local$ID), as.character(snps_fil_ldF@meta$sample_name))
#add the information. You can add more columns if you need, like ecological traits or other grouping classification 
snps_fil_ldF@meta$local_ID = missing_local[,1] #local ID
snps_fil_ldF@meta$longitude = missing_local[,2] #longitude
snps_fil_ldF@meta$latitude = missing_local[,3] #latitude
snps_fil_ldF@meta$missing_PC = missing_local[,4] #% of missing data
snps_fil_ldF@meta$ind_ID = missing_local[,5] #sample ID without "_sorted"
#verify the file
head(snps_fil_ldF@meta) #verify
tail(snps_fil_ldF@meta) #verify

#H. Save metafile with missing data and geographical information:
write_csv(snps_fil_ldF@meta, paste0("./Results_Metafiles/missing_data_coords_", project_name, ".csv"))


###3.3. SAVE THE .VCF FILE WITH ONLY NEUTRAL SNPS:
Save(snps_fil_ldF, paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"))








####### 4. Filtering Neutral Loci - By FST Outliers -----
###4.1. CREATE A geno OBJECT:
#A. Load the .VCF file with only neutral SNPs:
snps = vcfLink(paste0("vcf/", project_name,"_filtered_neutral_partial.vcf"), overwriteID=T)
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#B. Convert VCF to geno object. You need to specify the output file. It will automatically subset the vcf file and assign it as a new object.
snps = Geno(snps, output.file = paste0("./Results_Filters/sNMF_FST/", project_name, "_filtered_neutral_partial.geno"))
VCFsummary(snps) ##277 individuals and 6602 SNPs.

#C. Choose the 6 values of alpha. Alpha is the value of the regularization parameter (by default: 10). The results can depend on the value of this parameter, especially for small data sets. Less than 10,000 SNPs you can use values from 1000 to 10000. More than 10,000 SNP between 1 to 2,000. You can test different values.
alpha_values = c(10, 100, 500, 1000, 2000, 4000)

#D. Create folders for alpha values and copy .geno object in each folder:
for (i in alpha_values){
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i)
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("./Results_Filters/sNMF_FST/", project_name, "_filtered_neutral_partial.geno"), path )
}

#E. Set parameters to run SNMF (LEA) using different alpha values.
K = c(1:10) # set the number of K to be tested
replications = 10 # number of replication in each K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores for run in parallel

###4.2. RUN sNMF (LEA)
#A.Run a loop for all alpha values 
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop), pro_snmf)
}

#B. To load the SNMF projects. This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop), pro_snmf)
}

#C. Summary of the project
summary(project1)
summary(project2)
summary(project3)
summary(project4)
summary(project5)
summary(project6)

#D. View cross-Entropy plot
PlotK(project1) #5
PlotK(project2) #5
PlotK(project3) #5
PlotK(project4) #4
PlotK(project5) #4
PlotK(project6) #5
#another way to view the results
#plot(project1, lwd = 5, col = "red", pch=1)
#plot(project2, lwd = 5, col = "red", pch=1)
#plot(project3, lwd = 5, col = "red", pch=1)
#plot(project4, lwd = 5, col = "red", pch=1)
#plot(project5, lwd = 5, col = "red", pch=1)
#plot(project6, lwd = 5, col = "red", pch=1)

#E. Save Cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results_Filters/sNMF_FST/Cross_Entropy_sNMF_Alpha_",  i, ".pdf"), onefile = F)
  path = paste0("./Results_Filters/sNMF_FST/Alpha", i,"/", project_name, "_filtered_neutral_partial.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
}

#F. Select optimal K value
optimal_K = 2

#G. ATTENTION, if your dataset is K = 1 force a K = 2 to be able to filter SNPs with FST outliers.

#H. Select best run (lowest cross-entropy)  
best_run = Best.run(nrep=replications, optimalK=optimal_K, p1=project1, p2=project2, p3=project3, p4=project4, p5=project5, p6=project6)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("./Results_Filters/sNMF_FST/Alpha", alpha_values[as.numeric(best_run_split[6])],"/", project_name, "_filtered_neutral_partial.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])


###4.3. IDENTIFY AND REMOVE FST OUTLIERS
#A. Compute the FST statistics using best run
FST = fst(project, run, optimal_K) # you need at least 2 populations for a population-based test, so K>1.

#B. Compute the GIF - genomic inflation factor
lambda <- GIF(project, run, optimal_K, fst.values=FST)
lambda #1.222415

#C. Compute adjusted p-values from the combined z-scores and plot histogram of p-values
n = dim(Q(project, run, optimal_K))[1]
z.scores = sqrt(FST*(n-optimal_K)/(1-FST))
adj.p.values = pchisq(z.scores^2/lambda, df = optimal_K-1, lower = FALSE)
hist(adj.p.values, col = "red")

#D. Test different lambda values and plot histogram of p-values
adj.p.values = pchisq(z.scores^2/0.6, df = optimal_K-1, lower = FALSE) ## it is the best, but is still strange #try best value until close to one, uniform distribution with a peak at 0, and max 10% of snps 
hist(adj.p.values, col = "green")
C_fst <- candidates(alpha=0.05, adj.p.values)
ManPlot(adj.p.values, C_fst,"Fst")

#E. after you choose one, save the result
pdf("./Results_Filters/sNMF_FST/sNMF_FST_Outliers_adj_p_values.pdf", onefile = T)
hist(adj.p.values, col = "green")
dev.off()

#F. Candidate loci for FDR control: Benjamini-Hochberg at level q
C_fst <- candidates(alpha=0.05, adj.p.values)

#G. save the Manhatan plot
pdf("./Results_Filters/sNMF_FST/sNMF_FST_Outliers_Manhatan_plot.pdf", onefile = T) #add the number of snps removed in the name of pdf
ManPlot(adj.p.values, C_fst, paste0("Fst - Removing ",  length(C_fst), " SNPs"))
dev.off()

#H. Exclude candidate FST outlier
snps_fil_ldF_candidate <- Subset(snps, sites=C_fst)
snps_fil_ldF_candidate@site_id ## These are all the candidate SNPs
B_snmf = snps_fil_ldF_candidate@site_id
Chrom(snps_fil_ldF_candidate)

candidates_fst <- snps_fil_ldF_candidate@site_id
All_snp <- snps@site_id
N_snp <- All_snp[!(All_snp %in% candidates_fst)] ###Exclude all candidate loci

snps_neutral <- Subset(snps, sites=N_snp)
VCFsummary(snps_neutral) #55 individuals and 35638 SNPs.
length(N_snp)
length(snps@site_id)-length(C_fst) 
length(snps_neutral@site_id)

neutral_after_fst_sNMF = capture.output(VCFsummary(snps_neutral))
neutral_after_fst_sNMF #277 individuals and 5277 SNPs.

#I. Save the vcf
Save(snps_neutral, paste0("vcf/", project_name, "_filtered_neutral.vcf"))

##END