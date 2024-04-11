main_dir <- "/Users/siddhantkalra/Desktop/TGP/"
setwd(main_dir)

#Libraries----
library(pacman)
p_load(Biobase,ape,clusterProfiler,ggnewscale,ggupset,
       corrplot,data.table,dplyr,GEOquery,ggplot2,heatmap3,
       magrittr,PCAtools,pheatmap,qpcR,RColorBrewer,sjmisc,
       stringr,scmamp,tidyr,tidyverse,reshape2,ggplot2,
       org.Dm.eg.db,topGO,Rgraphviz,ggplot2,dplyr,hrbrthemes,
       viridis,enrichplot)
library(clusterProfiler)

#Part 1 (Merging files)----
#Combining all the files into one that
#will be used for analysis.
#I am using a file "1_Combining_files.csv" that contains path of the sub folder
#where multiple .txt files are stored.

source("Custom_functions.R")

#Looping over the "1_Combining_files.csv" (containing Path and Prefix info)
path_file <- read.csv("1_Combining_files.csv")
for (j in 1:nrow(path_file)){
  temp_path <- path_file[j,1]
  file <- list.files(path = temp_path, pattern=".txt")
  file
  variable_name <- paste0("df",j)
  setwd(temp_path)
  assign(variable_name,merging(file, path_file))
  setwd(main_dir)
}

#Reading 
Ortholog_info <- read.csv("2_orthologs.txt", sep ="\t")
f <- ortholog(Ortholog_info) #This will return the combined file
dim(f) #7836 genes and 97 columns (Experimental condition + gene name)
write.table(f, file="Final_combined_FPKM.csv", 
            sep = ",", 
            quote = FALSE,
            row.names = FALSE)
#End of Part1
#Part 2 Principal Component Analysis and Correlation graphs-----
setwd(main_dir)
rm(list=setdiff(ls(), "main_dir"))
source("Custom_functions.R")

dir.create("1_PCA")
f <-read.csv("Final_combined_FPKM.csv",check.names = FALSE,
                   row.names = 1)
#Generating Metadata file
metadata <- as.data.frame(matrix(ncol=5, nrow = ncol(f)))
names(metadata) <- c("Condition", "Stage",
                     "Species", "Previous_gen", "Current_gen")
metadata$Condition <- names(f)
metadata$Previous_gen[grepl("CC",metadata$Condition)] <- "C"
metadata$Current_gen[grepl("CC",metadata$Condition)] <- "C"

metadata$Previous_gen[grepl("CO",metadata$Condition)] <- "C"
metadata$Current_gen[grepl("CO",metadata$Condition)] <- "O"

metadata$Previous_gen[grepl("OC",metadata$Condition)] <- "O"
metadata$Current_gen[grepl("OC",metadata$Condition)] <- "C"

metadata$Previous_gen[grepl("OO",metadata$Condition)] <- "O"
metadata$Current_gen[grepl("OO",metadata$Condition)] <- "O"

metadata$Stage <- ifelse(grepl("Adult", metadata$Condition,
                               ignore.case = TRUE),
                         "Adult", "Larvae")
metadata$Species <- ifelse(grepl("sim", metadata$Condition,
                                 ignore.case = TRUE),
                           "D.Simulans", "D.Sechellia")
write.table(metadata, file="Metadata_generated.csv",
            sep=",", quote= FALSE, row.names = FALSE)
rownames(metadata) <- metadata$Condition
metadata$Condition <- NULL
#For combined
dir.create("1_PCA/Combined")
pdf(file="1_PCA/Combined/Plots.pdf")
my_pca(f, metadata)
dev.off()

#For Hybrid
dir.create("1_PCA/Hybrid")
f <-read.csv("Final_combined_FPKM.csv",check.names = FALSE,
                  row.names = 1)
f <- f %>%
  dplyr::select(contains("Hybrid"))
metadata <- read.csv("Metadata_generated.csv")
metadata <- subset(metadata, grepl("Hybrid", metadata$Condition))
rownames(metadata) <- metadata$Condition
metadata$Condition <- NULL
pdf(file="1_PCA/Hybrid/Plots.pdf")
my_pca(f, metadata)
dev.off()

#For Parent
dir.create("1_PCA/Parent")
f <-read.csv("Final_combined_FPKM.csv",check.names = FALSE,
             row.names = 1)
f <- f %>%
  dplyr::select(!contains("Hybrid"))
metadata <- read.csv("Metadata_generated.csv")
metadata <- subset(metadata, !grepl("Hybrid", metadata$Condition))
rownames(metadata) <- metadata$Condition
metadata$Condition <- NULL
pdf(file="1_PCA/Parent/Plots.pdf")
my_pca(f, metadata)
dev.off()
#end of part2

#Part 3 Creating Factors and data file----
setwd(main_dir)
rm(list=setdiff(ls(), "main_dir"))
source("Custom_functions.R")
f <-read.csv("Final_combined_FPKM.csv",check.names = FALSE,
             row.names = 1)
my_factors_and_data(f)
#end of part 3
#Part 4 Linear Model----
setwd(main_dir)
rm(list=setdiff(ls(), "main_dir"))
source("Custom_functions.R")
dir.create("2_LinearModel")

#For Combined
f <-read.csv("Factors_and_data.csv",check.names = FALSE)
dir.create("2_LinearModel/Combined")
setwd("2_LinearModel/Combined")
pdf(file="BarGraphs.pdf")
my_linear_model(f,7)
dev.off()
setwd(main_dir)

#For Parent
f <- read.csv("Factors_and_data.csv",check.names = FALSE,
              stringsAsFactors = FALSE)
dir.create("2_LinearModel/Parent")
setwd("2_LinearModel/Parent")
f$Generation <- 0
f <- f[!grepl("Hybrid", f$categories),]
pdf(file="BarGraphs.pdf")
my_linear_model(f,7)
dev.off()
setwd(main_dir)

#For Hybrid
f <- read.csv("Factors_and_data.csv",check.names = FALSE,
              stringsAsFactors = FALSE)
dir.create("2_LinearModel/Hybrid")
setwd("2_LinearModel/Hybrid/")
f$Generation <- 0
f <- f[grepl("Hybrid", f$categories),]
pdf(file="BarGraphs.pdf")
my_linear_model(f,7)
dev.off()
setwd(main_dir)
#End of part 4
#Part 5 Cis/Trans Interactions----
setwd(main_dir)
rm(list=setdiff(ls(), "main_dir"))
source("Custom_functions.R")
dir.create("3_Allele_Interactions")
f <-read.csv("Final_combined_FPKM.csv",check.names = FALSE,
             row.names = 1)

f1<-read.csv("2_LinearModel/Parent/7_QValues.csv")
a<-as.data.frame(f1$gene_name)
b<-f1$Species
Parent_species<-cbind(a,b)
Parent_species<-as.data.frame(Parent_species)
colnames(Parent_species)<-c("gene_name","P_Species_Qvalues")

f1<-read.csv("2_LinearModel/Hybrid/7_QValues.csv")
a<-as.data.frame(f1$gene_name)
b<-f1$Species
Hybrid_species<-cbind(a,b)
Hybrid_species<-as.data.frame(Hybrid_species)
colnames(Hybrid_species)<-c("gene_name","H_Species_Qvalues")

f1<-read.csv("2_LinearModel/Combined/7_QValues.csv")
a<-as.data.frame(f1$gene_name)
b<-f1$Species.Generation
All_species<-cbind(a,b)
All_species<-as.data.frame(All_species)
colnames(All_species)<-c("gene_name","A_Species_Gen_Qvalues")

data<-Reduce(function(x,y) merge(x,y,by="gene_name",all=TRUE) ,list(Parent_species,Hybrid_species,All_species))
rownames(data) <- data$gene_name
data$gene_name <- NULL


#For Adults
setwd("3_Allele_Interactions/")
#setwd("3_Allele_updated_2/")
f1 <- f %>%
  dplyr::select(contains("Adult"))
pdf(file="Adult.pdf")
my_allele_interaction(f1,"Adult",data) 
dev.off()

#For Larvae
f1 <- f %>%
  dplyr::select(!contains("Adult"))
pdf(file="Larvae.pdf")
my_allele_interaction(f1,"Larvae",data) 
dev.off()

#Bar graphs
pdf(file="Graphs.pdf")
f <- list.files(pattern = ".csv")
message(f)
df <- as.data.frame(matrix(ncol=2,nrow=0))
colnames(df) <- c("Check","Category")
for (i in f){
  f1 <- read.csv(i)
  Check <- f1$Check
  Check <- as.data.frame(Check)
  Category <- f1$Category
  Category <- as.data.frame(Category)
  temp <- cbind(Check,Category)
  df <- rbind(df,temp)
  write.table(df, file="Bar_graph1.csv", sep=",",
              quote= FALSE, row.names = FALSE)
}
data <- as.data.frame(table(df))

for (i in unique(data$Category)){
  temp <- subset(data, data$Category == i)
  my_bar <- barplot(height=temp$Freq,
                    names=temp$Check,
                    col="#69b3a2",
                    main=i)
  text(my_bar, temp$Check,paste(temp$Freq))
}
dev.off()

#Part 6 Fisher Exact test Odd ratios----
main_dir <- "/Users/siddhantkalra/Desktop/TGP/"
setwd(main_dir)
pattern <- "(Adult|Larvae).*\\.csv$"
file_list <- list.files(paste0(main_dir,"3_Allele_Interactions"),
                        full.names = TRUE,
                        pattern = pattern)
# file_list <- list.files(paste0(main_dir,"3_Allele_updated_2"),
#                         full.names = TRUE,
#                         pattern = pattern)
file_list
dir.create("4_OddRatios")

#For Parents
f1 <- read.csv("2_LinearModel/Parent/11_QValue_gene_names.csv")
dir.create("4_OddRatios/Parents")
setwd("4_OddRatios/Parents")
# setwd("4_OddRatios_updated/")
pdf(file="Plots.pdf")
my_oddRatios(f1,file_list)
#my_oddRatios_updated(f1,file_list)
dev.off()
setwd(main_dir)

#For Hybrids
f1 <- read.csv("2_LinearModel/Hybrid/11_QValue_gene_names.csv")
dir.create("4_OddRatios/Hybrids")
setwd("4_OddRatios/Hybrids")
pdf(file="Plots.pdf")
my_oddRatios(f1,file_list)
dev.off()
setwd(main_dir)

#For Combined
f1 <- read.csv("2_LinearModel/Combined/11_QValue_gene_names.csv")
dir.create("4_OddRatios/Combine")
setwd("4_OddRatios/Combine")
pdf(file="Plots.pdf")
my_oddRatios(f1,file_list)
dev.off()
setwd(main_dir)

#Part 7 Represntative Heatmap----
setwd(main_dir)
source("Custom_functions.R")
pattern <- "Odd_ratios"
file_list <- list.files(paste0(main_dir,"4_OddRatios/Parents"),
                        full.names = TRUE,
                        pattern = pattern)
file_list
pdf("Main_represntative_HeatMap.pdf")
representative_map(file_list, "All 8 conditions overlapped")
pattern <- "Odd_ratios_Larvae"
file_list <- list.files(paste0(main_dir,"4_OddRatios/Parents"),
                        full.names = TRUE,
                        pattern = pattern)
file_list
representative_map(file_list, "All Larvae conditions overlapped")
pattern <- "Odd_ratios_Adult"
file_list <- list.files(paste0(main_dir,"4_OddRatios/Parents"),
                        full.names = TRUE,
                        pattern = pattern)
file_list
representative_map(file_list, "All Adult conditions overlapped")
dev.off()
