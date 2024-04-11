setwd("/Users/siddhantkalra/Desktop/TGP/")
dir.create("Updated_counts")
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
setwd("Updated_counts/")
f1 <- f %>%
  dplyr::select(contains("Adult"))
keyword <- "Adult"
pdf(file="Adult.pdf")

# #For Larvae
# setwd("Updated_counts/")
# f1 <- f %>%
#   dplyr::select(!contains("Adult"))
# keyword <- "Larvae"
# pdf(file="Larvae.pdf")

f <- f1
combinations <- c("CC", "CO", "OC", "OO")
for (val in combinations){
  f1 <- f %>%
    dplyr::select(contains(val))
  parent <- f1 %>%
    dplyr::select(!contains("Hybrid"))
  hybrid <- f1 %>%
    dplyr::select(contains("Hybrid"))
  
  #P1->Adult Sim
  P1 <- as.data.frame(
    rowMeans(parent[,1:3], 
             na.rm = FALSE, dims = 1)
  )  
  colnames(P1) <- "P1"
  
  #P2->Adult Sec
  P2 <- as.data.frame(
    rowMeans(parent[,4:6], 
             na.rm = FALSE, dims = 1)
  )  
  colnames(P2) <- "P2"
  
  #A1-> Hybrid Sim
  A1 <- as.data.frame(
    rowMeans(hybrid[,1:3], 
             na.rm = FALSE, dims = 1)
  )  
  colnames(A1) <- "A1"
  
  #A2-> Hybrid Sec
  A2 <- as.data.frame(
    rowMeans(hybrid[,4:6], 
             na.rm = FALSE, dims = 1)
  )  
  colnames(A2) <- "A2"
  df <- cbind(P1,P2,A1,A2)
  df$'P1/P2' <- (df$P1+0.000001)/(df$P2+0.000001)
  df$'log2(P1/P2)' <- log2(df$`P1/P2`)
  df$'A1/A2' <- (df$A1+0.000001)/(df$A2+0.000001)
  df$'log2(A1/A2)' <- log2(df$`A1/A2`)
  df$Check <- paste0(keyword,"_",val)
  df <- cbind(df,data)
  
  df$Category[(df$P_Species_Qvalues>0.05) & (df$H_Species_Qvalues>0.05) & (df$A_Species_Gen_Qvalues>0.05)] <- "Conserved"
  df$Category[(df$P_Species_Qvalues<0.05) & (df$H_Species_Qvalues<0.05) & (df$A_Species_Gen_Qvalues>0.05)] <- "Only_Cis"
  df$Category[(df$P_Species_Qvalues<0.05) & (df$H_Species_Qvalues>0.05) & (df$A_Species_Gen_Qvalues<0.05)] <- "Only_Trans"
  df$Category[(df$P_Species_Qvalues<0.05) & (df$H_Species_Qvalues<0.05) & (df$A_Species_Gen_Qvalues<0.05) & ((df$`log2(P1/P2)`/df$`log2(A1/A2)`)>1)] <- "Cis_plus_Trans"
  df$Category[(df$P_Species_Qvalues<0.05) & (df$H_Species_Qvalues<0.05) & (df$A_Species_Gen_Qvalues<0.05) & ((df$`log2(P1/P2)`/df$`log2(A1/A2)`)<1)] <- "Cis_by_Trans"
  df$Category[(df$P_Species_Qvalues>0.05) & (df$H_Species_Qvalues<0.05) & (df$A_Species_Gen_Qvalues<0.05)] <- "Cis_equal_Trans"
  df$Category[(df$P_Species_Qvalues<0.05) & (df$H_Species_Qvalues>0.05) & (df$A_Species_Gen_Qvalues>0.05)] <- "Ambiguous"
  df$Category[(df$P_Species_Qvalues>0.05) & (df$H_Species_Qvalues<0.05) & (df$A_Species_Gen_Qvalues>0.05)] <- "Ambiguous"
  df$Category[(df$P_Species_Qvalues>0.05) & (df$H_Species_Qvalues>0.05) & (df$A_Species_Gen_Qvalues<0.05)] <- "Ambiguous"
  
  df <- df[which(df$P1>0 & df$P2 > 0),]
  df <- df[which(df$P1>1 | df$P2 > 1),]
  df <- df[which(df$A1>0 & df$A2 > 0),]
  df <- df[which(df$A1>1 | df$A2 > 1),]
  
  a<-as.data.frame(table(df$Category))
  a <- subset(a, a$Var1 != "Ambiguous") #Comment this line if needed to include ambiguous
  lbls <- paste(a$Var1, a$Freq)
  lbls
  pie(a$Freq,labels = lbls, main=paste0(keyword,"_",val))
  
  #Conserved = Black
  Conserved<-df[which(df$Category=="Conserved"),]
  plot(Conserved$`log2(P1/P2)`, Conserved$`log2(A1/A2)`, pch = 19, col = "#000000", 
       xlim=c(-10,10), ylim=c(-10,10), main=paste0(keyword,"_",val), xlab = "log2(P1/P2)", ylab = "log2(A1/A2)")
  
  # Ambiguous<-df1[which(df1$Category=="Ambiguous"),]
  # points(Ambiguous$`log2(P1/P2)`,Ambiguous$`log2(A1/A2)`,pch = 19, col = "#8BE064", xlim=c(-10,10), ylim=c(-10,10))
  
  #Compensating = Green
  Compensating<-df[which(df$Category=="Cis_equal_Trans"),]
  points(Compensating$`log2(P1/P2)`,Compensating$`log2(A1/A2)`,pch = 19, col = "#8BE064", xlim=c(-10,10), ylim=c(-10,10))
  
  #Only Cis = Skin
  Only_Cis<-df[which(df$Category=="Only_Cis"),]
  points(Only_Cis$`log2(P1/P2)`,Only_Cis$`log2(A1/A2)`,pch = 19, col = "#FDC086", xlim=c(-10,10), ylim=c(-10,10))
  
  #Only Trans = Red
  Only_Trans<-df[which(df$Category=="Only_Trans"),]
  points(Only_Trans$`log2(P1/P2)`,Only_Trans$`log2(A1/A2)`,pch = 19, col = "#E31A1C", xlim=c(-10,10), ylim=c(-10,10))
  
  #Cis_plus_Trans = Blue
  Cis_plus_Trans<-df[which(df$Category=="Cis_plus_Trans"),]
  points(Cis_plus_Trans$`log2(P1/P2)`,Cis_plus_Trans$`log2(A1/A2)`,pch = 19, col = "#0096FF", xlim=c(-10,10), ylim=c(-10,10))
  
  Cis_by_Trans<-df[which(df$Category=="Cis_by_Trans"),]
  points(Cis_by_Trans$`log2(P1/P2)`,Cis_by_Trans$`log2(A1/A2)`,pch = 19, col = "pink", xlim=c(-10,10), ylim=c(-10,10))
  
  abline(v=0)
  abline(h=0)
  abline(coef = c(0,1))
  abline(coef = c(0,-1))
  abline(v=0)
  abline(h=0)
  abline(coef = c(0,1))
  abline(coef = c(0,-1))

  df <- setDT(df, keep.rownames = "Gene")[]
  write.table(df, file=paste0(keyword,"_",val,".csv"),
              sep=",", quote = FALSE,
              row.names = FALSE)
  
}
dev.off()

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
}
  write.table(df, file="Bar_graph1.csv", sep=",",
            quote= FALSE, row.names = FALSE)
  
df_for_bar_graph <- as.data.frame(matrix(ncol = 3, nrow = 0))
names(df_for_bar_graph) <- c("Experimental_Condition", "Category", "Count")
for(i in unique(df$Check)){
  temp <- subset(df, df$Check==i)
  ambi <- nrow(subset(temp, temp$Category =="Ambiguous"))
  con <- nrow(subset(temp, temp$Category =="Conserved"))
  all_cis <- nrow(subset(temp, grepl("Cis", temp$Category)))
  all_trans <- nrow(subset(temp, grepl("Trans", temp$Category)))
  Category <- c("Ambiguous","Conserved", "All_Cis", "All_Trans")
  Count <- c(ambi,con, all_cis,all_trans)
  Experimental_Condition <- rep(temp[1,1],4)
  temp1 <- cbind(Experimental_Condition,Category,Count)
  df_for_bar_graph <- rbind(df_for_bar_graph,temp1)
}

write.table(df_for_bar_graph, file="BarGraphs_2.csv",
            sep=",", quote = FALSE, row.names = FALSE)
df_for_bar_graph <- read.csv("BarGraphs_2.csv")

for (i in unique(df_for_bar_graph$Category)){
  temp <- subset(df_for_bar_graph, df_for_bar_graph$Category == i)
  my_bar <- barplot(height=temp$Count,
                    names=temp$Experimental_Condition,
                    col="#69b3a2",
                    main=i)
  text(my_bar, temp$Count,paste(temp$Count))
}  
dev.off()

#Odd ratios and heatmap
pdf(file="OddRatios_Heatmap.pdf")
pattern <- "(Adult|Larvae).*\\.csv$"
file_list <- list.files(pattern = pattern)
file_list

#For parents
f1 <- read.csv("/Users/siddhantkalra/Desktop/TGP/2_LinearModel/Parent/11_QValue_gene_names.csv")
for (i in file_list){
  f <- read.csv(i)
  f <- as.data.frame(cbind(f$Gene, f$Check, f$Category))
  names(f) <- c("Gene", "Check", "Category")
  df <- data.frame(matrix(ncol = ncol(f1), nrow = 4))
  names(df) <- names(f1)
  row.names(df) <- c("Size", "Conserved", "All_Cis",
                     "All_Trans")
  colnames(df)[1] <- "Genome"
  col_names <- colnames(f1)
  col_names
  for (j in 1:length(col_names)){
    x <- f1[,j]
    x <- na.omit(x)
    x <- as.data.frame(x)
    colnames(x) <- col_names[j]
    df[1,j] <- nrow(x)
    mer <- merge(f, x, by.x="Gene", by.y= names(x))
    df[2,j] <- nrow(subset(mer,grepl("Conserved",mer$Category)))
    df[3,j] <- nrow(subset(mer,grepl("Cis",mer$Category)))
    df[4,j] <- nrow(subset(mer,grepl("Trans",mer$Category)))
  }
  df_pvalue <- df
  df_odd_ratios <- df
  df_pvalue <- df_pvalue[-1,-1]
  df_odd_ratios <- df_odd_ratios[-1,-1]
  
  for (m in 2:nrow(df)){
    for (n in 2:ncol(df)){
      mini <- df[c(1,m), c(1,n)]
      test_results <- fisher.test(mini)
      df_pvalue[m-1,n-1] <- test_results[["p.value"]]
      df_odd_ratios[m-1,n-1] <- test_results[["estimate"]]
    }
  }
  
  #Log2 transforming odd ratios
  for (m in 1:nrow(df_odd_ratios)){
    for (n in 1:ncol(df_odd_ratios)){
      df_odd_ratios[m,n] <-log2(df_odd_ratios[m,n] + 0.000000001)
    }
  }
  
  #Converting Non-significant values to zero
  for (m in 1:nrow(df_odd_ratios)){
    for (n in 1:ncol(df_odd_ratios)){
      if (df_pvalue[m,n] > 0.05){
        df_odd_ratios[m,n] <- 0
      }
    }
  }
  paletteLength <- 101
  df_odd_ratios_2 <- df_odd_ratios
  
  df_odd_ratios_2 <- ifelse(df_odd_ratios > 0, 1, ifelse(df_odd_ratios < 0, -1, 0))
  
  #df_odd_ratios_2 <- ifelse(df_odd_ratios > 0, 1, -1)

  #color_palette <- colorRampPalette(c("blue", "white", "red"))(n = 101)  # 101 colors for a smooth gradient
  #my_color <- colorRampPalette(c("yellow", "white", "blue"))(paletteLength)
  my_color <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  
  # myBreaks <- c(seq(min(df_odd_ratios), 0,
  #                   length.out = ceiling(paletteLength/2)+1),
  #               seq(max(df_odd_ratios)/paletteLength, max(df_odd_ratios),
  #                   length.out = floor(paletteLength/2)+50))
  
  #1 original
  pheatmap(as.matrix(df_odd_ratios_2), main = f[1,2]) #, breaks = myBreaks)
  #2 original without clustering
  pheatmap(as.matrix(df_odd_ratios_2), main = f[1,2],
           cluster_rows = FALSE, cluster_cols = FALSE) #, breaks = myBreaks)
  
    # df[2,j] <- nrow(mer[which(mer$Category == "Conserved"),])
    # df[3,j] <- nrow(mer[which(mer$Category == "Only_Cis"),])
    # df[4,j] <- nrow(mer[which(mer$Category == "Only_Trans"),])
    # df[5,j] <- nrow(mer[which(mer$Category == "Enhancing"),])
    # df[6,j] <- nrow(mer[which(mer$Category == "Compensating"),])
  
}
dev.off()
