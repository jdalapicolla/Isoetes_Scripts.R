###############################################################################
####################### VALE INSTITUTE OF TECHNOLOGY ##########################
###############################################################################

#######################   STEP 02: GENETIC STRUCTURE   ########################

#Script prepared by Jeronymo Dalapicolla and Rodolfo Jaffé for the manuscript:
# Dalapicolla et al. - Conservation implications genetic structure in the narrowest endemic quillwort from the Eastern Amazon


########## PRE-ANALYSIS AND INFORMATION ----

##A. THIS PIPELINE WAS DESIGNED IN UBUNTU 18.04 LTS, USING RSTUDIO 1.2.1335 AND R 3.6.3

##B. LOAD THE FILE "functions_LanGen.R" WITH FUNCTIONS TO BE USED ON THIS STEP.
source("functions_LanGen.R")


##C. INSTALL PACKAGES
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
if("dartR" %in% rownames(installed.packages()) == FALSE){install.packages("dartR")
} else {print (paste0("'dartR' has already been installed in library"))}
if("adegenet" %in% rownames(installed.packages()) == FALSE){install.packages("adegenet")
} else {print (paste0("'adegenet' has already been installed in library"))}
if("tidyverse" %in% rownames(installed.packages()) == FALSE){install.packages("tidyverse")
} else {print (paste0("'tidyverse' has already been installed in library"))}
if("ggplot2" %in% rownames(installed.packages()) == FALSE){install.packages("ggplot2")
} else {print (paste0("'ggplot2' has already been installed in library"))}


#D. load packages multiple packages use the package: 'pacman'. If the package is missing "p_load" will download it from CRAN. Using "" in the name of packages isn't mandatory.
pacman::p_load(r2vcftools, LEA, vcfR, dartR, adegenet, tidyverse, ggplot2)

##E. CREATE FOLDERS AND DIRECTORIES TO SAVE THE RESULTS:
create_dir(c("./Results_snmf", "./Results_DAPC", "./Results_PCA"))






###### 1. Genetic Structure using sNMF -----
##1.1. CHOOSE A NAME FOR THE PROJECT. MUST BE THE SAME ONE THAN FILTERING STEP:
#A. Project name:
project_name = "isoetes"

###1.2. POPULATION ASSIGNMENT ANALYSIS USING sNMF WITHOUT FST OUTLIERS:
#Now we will carry out population assignment again, but using only neutral loci dataset.
#A1. Load neutral .vcf file after remove the outlier SNPs
snps_fil_ldF_neutral <- vcfLink(paste0("vcf/", project_name, "_filtered_neutral.vcf"), overwriteID = T) 
VCFsummary(snps_fil_ldF_neutral) #55 individuals and 35638 SNPs.

##B. Convert to geno object.You need to specify the output file. It will automatically subset the vcf file and assign it as a new object
snps_fil_ldF_neutral <- Geno(snps_fil_ldF_neutral, output.file = paste0("vcf/", project_name, "_filtered_neutral.geno"))
VCFsummary(snps_fil_ldF_neutral) #55 individuals and 35638 SNPs.

#C. Create folders for alpha values and copy .geno object in each folder:
alpha_values = c(10, 100, 500, 1000, 2000, 4000)
for (i in alpha_values){
  path = paste0("./Results_snmf/Alpha", i, "n")
  if (dir.exists(file.path(getwd(), path)) == FALSE)
  {dir.create(path, recursive = T, showWarnings = F)} else (print (paste0(path, " has already been created. Be careful with overwritting")))
  file.copy(paste0("vcf/", project_name, "_filtered_neutral.geno"), path )
}

#D. Set parameters to run SNMF (LEA) using different alpha.
K = c(1:10) # K to be tested
replications = 10 # numeber of replication by K
ploidy = 2 # species ploidy
CPU = 4 #Number of cores

#E. Run sNMF (LEA) using different alpha.
set.seed(123)
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.geno")
  pro_snmf = snmf(path, K = K, rep = replications, alpha = i, entropy = T, ploidy = ploidy , project = "new", CPU= CPU)
  assign(paste0("project_snmf", loop, "n"), pro_snmf)
}

#F. To load the SNMF projects in a new R session (after quitting R), use:  project = load.snmfProject("Alpha1000//pilocarpus_filtered_ld_hw.snmfProject") ##This allows you to save time because you do not need to run SNMF again!
loop = 0 #set ALWAYS as 0.
for (i in alpha_values){
  loop = loop +1
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  pro_snmf = load.snmfProject(path)
  assign(paste0("project", loop, "n"), pro_snmf)
}

#G. summary of the project
summary(project1n)
summary(project2n)
summary(project3n)
summary(project4n)
summary(project5n)
summary(project6n)

#H. View Cross-Entropy plots
PlotK(project1n) #4
PlotK(project2n) #4
PlotK(project3n) #4
PlotK(project4n) #4
PlotK(project5n) #3
PlotK(project6n) #4
#another way to view the results
#plot(project1n, lwd = 5, col = "red", pch=1)
#plot(project2n, lwd = 5, col = "red", pch=1)
#plot(project3n, lwd = 5, col = "red", pch=1)
#plot(project4n, lwd = 5, col = "red", pch=1)
#plot(project5n, lwd = 5, col = "red", pch=1)
#plot(project6n, lwd = 5, col = "red", pch=1)

#I. Save graphs of cross-Entropy plot with standard deviation error bars
for (i in alpha_values){
  pdf(paste0("./Results_snmf/Cross_Entropy_sNMF_Alpha_",  i, "n.pdf"), onefile = F)
  path = paste0("Results_snmf/Alpha", i,"n/", project_name, "_filtered_neutral.snmfProject")
  print(PlotK(load.snmfProject(path)))
  dev.off()
}

#J. Select optimal K value
optimal_K = 1

#K. Select best run (lowest cross-entropy)
best_run = Best.run(nrep=10, optimalK=optimal_K, p1=project1n, p2=project2n, p3=project3n, p4=project4n, p5=project5n, p6=project6n)
#load the best project
best_run_split = scan(text = best_run, what = "")
path_best_run = paste0("Results_snmf/Alpha", alpha_values[as.numeric(best_run_split[6])],"n/", project_name, "_filtered_neutral.snmfProject")
#set the values
project = load.snmfProject(path_best_run)
run=as.numeric(best_run_split[9])

#L. Barplot replace best run information
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
order = barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
write.table(order[[1]], "./Results_Metafiles/snmf_bestK_pipegraph_order_neutral.txt")

pdf("./Results_snmf/sNMF_Pipegraph_neutral_version2.pdf", onefile = F)
plot.new()
barplotK(Qfile=project, Pop = optimal_K, Run_B = run)
dev.off()

pdf("./Results_snmf/sNMF_Pipegraph_neutral.pdf", onefile = F)
my.colors <- rainbow(optimal_K)
LEA::barchart(project, K = optimal_K, run = run, border = NA, space = 0, col = my.colors, xlab = "Individuals", ylab = "Ancestry proportions", main = "Ancestry matrix") -> bp
axis(1, at = 1:length(bp$order),labels = bp$order, las=1, cex.axis = .3)
dev.off()

#M. Add admixture coefficient and replace the population ID to vcf file
Qmat <- as.data.frame(Q(project, run=run, K=optimal_K))
head(Qmat)

columns = c() 

for (i in 1:ncol(Qmat)){
  columns[i] = paste0("Adx_Coeff_", i)
}

colnames(Qmat) = columns
head(Qmat)
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

for (i in 1:optimal_K){
  j = ncol(snps_fil_ldF_neutral@meta)+1
  snps_fil_ldF_neutral@meta[,j] = Qmat[i]
}
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#N. Verify individuals and maximum Adx_Coeff and add to metafile:
popIds = apply(Qmat, 1, which.max)
snps_fil_ldF_neutral@meta$PopID_snmf <- popIds
head(snps_fil_ldF_neutral@meta)
tail(snps_fil_ldF_neutral@meta)

#O. Save new vcf file with sNMF results and pop ID
Save(snps_fil_ldF_neutral, paste0("vcf/", project_name, "_filtered_neutral_LEA.vcf")) #save vcf
write.csv(snps_fil_ldF_neutral@meta, paste0("Results_Metafiles/ancestry_coef_LEA_", project_name, ".csv"), quote = F) #save result as table
write.csv(as.data.frame(snps_fil_ldF_neutral@meta$PopID_snmf), file= paste0("Results_Metafiles/", project_name, "_neutral_popsID_LEA.csv")) #save only pop ID from sNMF

VCFsummary(snps_fil_ldF_neutral) #277 individuals and 5268 SNPs.







###### 2. Genetic Structure using DAPC -----
###2.1. POPULATION ASSIGNMENT ANALYSIS USING DAPC:
#A. Load neutral .vcf file after removing outlier SNPs:
vcf = read.vcfR(paste0("vcf/", project_name, "_filtered_neutral.vcf"), verbose = FALSE)

#B. Convert "VCF" to "GENIND"
input = vcfR2genind(vcf)
input

#C. Perform a PCA to choose the number of PC in the DAPC:
input_scaled = scaleGen (input, center = TRUE, scale = TRUE, NA.method = "mean")
pca_input = dudi.pca(input_scaled, center = TRUE, scannf = FALSE)

#D. % of PC variation
pc_pca = as.data.frame(pca_input$eig)
pc_pca[,2] = (pc_pca/sum(pc_pca))*100

#D1. Rule of 100% of variance:
index_100 = length(rownames(input@tab))-1
index_100 # number of PC to reach to 100% of explained variance

#D2. Rule of at least 95% of variance:
index_95 = length(which(cumsum(pc_pca[,2]) <= 95))
index_95

#D3. Rule of at least 70% of variance:
index_70 = length(which(cumsum(pc_pca[,2]) <= 70))
index_70 

#D4. Rule of minimum variance per PCs:
variance_pc = 100/(nrow(input@tab)-1)
variance_pc #PCs that increase the explained variance bellow this threshould will be removed
#calculate number of PCs
index_min = length(pc_pca[,2][pc_pca[,2] >= variance_pc])
index_min

#E. Identification of the clusters (We specify that we want to evaluate up to k = 10 groups (max.n.clust=40)
index_100 #54 PCs - k = 1 #For me, 100% of variation is more coherent.
index_95  #51 PCs - k = 1
index_70  #36 PCs - k = 1
index_min #25 PCs - k = 1
#If you see a plateau on the graph, you can choose a number of PCs in the begining of this plateau. In Policarpus there's no plateau.
set.seed(13) #set a seed
grp = find.clusters(input, max.n.clust=10, scale = TRUE) # center=T by default.
#Digit the number of PCs to retain.
#verify Group (DPCA) of the first 10 individuals and the Size of Groups
head(grp$grp, 10)
grp$size


#F. Select the ’best’ BIC is often indicated by an below in the curve of BIC values as a function of k
#data frame with results 
df = as.data.frame(grp$Kstat)
df$K = c(1:10)
colnames(df) = c("BIC", "K")
df

#graph
dapc_k_graph = ggplot(df, aes(x=K, y=BIC)) + 
  geom_line() + 
  geom_point(size=4, shape=21, fill="red", color="darkred") + xlab("Number of Genetic Clusters") + ylab("BIC")+
  theme_pca +
  theme(axis.text.x=element_text(angle=0, vjust=0.5), axis.text = element_text(size=12, color = "black"), axis.title.x = element_text(size=12, color="black"),axis.title.y = element_text(size=12, color="black")) +
  scale_x_continuous(breaks = seq(0,10, by=2))

#save best k graph
pdf("./Results_DAPC/Bestk_DAPC.pdf", onefile = F)
dapc_k_graph
dev.off()







####### 3. Genetic Structure using PCA -----
###3.1. RUN THE PCA:
#A. You ran the PCA before:
pca_input

#B. % of contribution for three first PCs
PC1 = paste0("PC1 (",round(pca_input$eig[1]/sum(pca_input$eig)*100,2),"%)")
PC2 = paste0("PC2 (",round(pca_input$eig[2]/sum(pca_input$eig)*100,2),"%)")
PC3 = paste0("PC3 (",round(pca_input$eig[3]/sum(pca_input$eig)*100,2),"%)")

#C. Min and max values for axes
pc1_range = c(min(pca_input$li[,1]), max(pca_input$li[,1]))
pc2_range = c(min(pca_input$li[,2]), max(pca_input$li[,2]))
pc3_range = c(min(pca_input$li[,3]), max(pca_input$li[,3]))

#D. save results as .CSV
#coordenates in the PCA by axis
write.csv(pca_input$li, file = "Results_PCA/Axis_coord.csv")
#% of PC variation
perc_pca = as.data.frame(pca_input$eig)
soma = sum (perc_pca)
perc_pca[,2] = (perc_pca/soma)*100
colnames(perc_pca) = c("Eigenvalues", "Contribution (%)")
write.csv(perc_pca, file = "Results_PCA/contribution_pc_eig.csv")


###3.2. PCA GRAPHS:
#data frame:
df2 = pca_input$li
df2$site = snps_fil_ldF_neutral@meta$local_ID
df2

#graph
graph_pca = ggplot(data=df2) +
  geom_point(mapping = aes(x=Axis1, y = Axis2, fill=site), size=4, colour="black", pch=21) +
  theme_pca +
  labs(x = PC1, y = PC2) +
  scale_fill_manual (values = c("grey", "red", "blue", "yellow", "green"),
                     labels = c("Center", "East", "North", "West", "South")) +
  guides(fill=guide_legend(title="Area"))

#save graph
pdf(paste0("./Results_PCA/PCA_pca_", project_name,".pdf"))
graph_pca
dev.off()

##END