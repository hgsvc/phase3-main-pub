require("dbscan")
require("fpc")
library("factoextra")
library("SNPRelate")
library("data.table")
library("raster")


# 1kg info file
hgsvc_samples <- read.delim("1KG_3202_samples.ped",header=T)

#hgsvc_samples_2504 <- read.table("all_2504_1KG_samples.txt",header=F)
#colnames(hgsvc_samples_2504) <- c("SampleID")
#cbind(rep("2504_unrelated",2504),hgsvc_samples_2504$SampleID)



##### IGNORE THIS BLOCK -- RELEVANT TO HPRC AND HGSVC SAMPLES ONLY
# Subsets of samples of interest (NOTE THAT I HAVE 2 VERSIONS OF FILES WITH THESE SAMPLES)
year_1_samples <- read.table("../SNV_genotypes/YEAR_1_LongRead_samples_v2.txt",header=T)
year_2_samples <- read.delim("../SNV_genotypes/YEAR_2_LongRead_samples_v2.txt",header=T)
colnames(year_1_samples) <- "Sample"
colnames(year_2_samples) <- "Sample"

hgsvc_lr_YR1and2_samples <- rbind(year_1_samples,year_2_samples)
hgsvc_lr_YR1and2_samples <- unique(hgsvc_lr_YR1and2_samples)
hprc_samples <- read.delim("../SNV_genotypes/HPRC_combined_120samples.txt",header=T)
colnames(hprc_samples)[2] <- "Sample"

# Create a merged set for HGSVC and HPRC sample PEDs
merged_1 <- merge(hgsvc_samples,hgsvc_samples_2504,by="SampleID",all=T)
write.table(merged_1,"merged_1.tab",quote=F,sep="\t",row.names=F)

# Combine all HPRC and HGSVC samples into one array
hgsvc_hprs_sampleIDs <- c(year_1_samples$Sample,year_2_samples$Sample,hprc_samples$Sample)
hgsvc_hprs_sampleIDs <- unique(sort(hgsvc_hprs_sampleIDs))


# PCA ANALYSIS STARTS HERE
vcf.fn <- "pangenie_merged_bi_nosnvs.vcf.gz"

#snpgdsClose(genofile)
# Convert VCF to GDS
snpgdsVCF2GDS(vcf.fn, "hgsvc_sv.gds", method="copy.num.of.ref",verbose = T)
genofile <- openfn.gds("hgsvc_sv.gds")

# Start PCA 
partial_pca <- snpgdsPCA(genofile,num.thread = 4,algorithm = "randomized")
## EIGENVECT containcs PCS you need to plot
partial_eigenvect <- partial_pca$eigenvect

# Annotate with population labels (1KG)
out_mx <- data.frame(cbind(partial_pca$sample.id,partial_eigenvect))
colnames(out_mx) <- c("SampleID",paste0("PC",1:16))
merged_out <- merge(out_mx,hgsvc_samples,by="SampleID")
write.table(merged_out,"HGSVC_3202_Autosomes_PCA.tab",quote=F,sep='\t',row.names=F)


###### LOAD SAVED PCA MATRIX ######
# When rerunning, read the tab output file:
merged_out <- read.delim("HGSVC_3202_Autosomes_PCA.tab",header=T)
merged_out <- read.delim("../SNV_genotypes/HGSVC_2504_chr22GATK_PCA.tab",header=T)

merged_out$PC1 <- as.numeric(as.character(merged_out$PC1))
merged_out$PC2 <- as.numeric(as.character(merged_out$PC2))
merged_out$PC3 <- as.numeric(as.character(merged_out$PC3))
merged_out$PC4 <- as.numeric(as.character(merged_out$PC4))
merged_out$PC5 <- as.numeric(as.character(merged_out$PC5))
merged_out$PC6 <- as.numeric(as.character(merged_out$PC6))
merged_out$PC7 <- as.numeric(as.character(merged_out$PC7))

#merged_out$Superpopulation_colors <- rep("black",3202)
merged_out$Superpopulation_colors <- rep("black",2504)

merged_out[merged_out$Superpopulation=="EUR",]$Superpopulation_colors <- rep("#e41a1c",length(which(merged_out$Superpopulation=="EUR")))
merged_out[merged_out$Superpopulation=="EAS",]$Superpopulation_colors <- rep("#377eb8",length(which(merged_out$Superpopulation=="EAS")))
merged_out[merged_out$Superpopulation=="SAS",]$Superpopulation_colors <- rep("#4daf4a",length(which(merged_out$Superpopulation=="SAS")))
merged_out[merged_out$Superpopulation=="AFR",]$Superpopulation_colors <- rep("#984ea3",length(which(merged_out$Superpopulation=="AFR")))
merged_out[merged_out$Superpopulation=="AMR",]$Superpopulation_colors <- rep("#ff7f00",length(which(merged_out$Superpopulation=="AMR")))



# PCA panel (1x2) for the Pangenie genotypes
par(mfrow=c(1,1))
plot(merged_out$PC1,merged_out$PC2,pch=19,col=merged_out$Superpopulation_colors,xlab="PC1",ylab="PC2")
plot(merged_out$PC2,merged_out$PC3,pch=19,col=merged_out$Superpopulation_colors,xlab="PC2",ylab="PC3",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065))



# PCA panel (2x2) for the GATK genotypes
par(mfrow=c(2,2))
plot(merged_out$PC1,merged_out$PC2,pch=19,col=merged_out$Superpopulation_colors,xlab="PC1 (1.29%)",ylab="PC2 (0.42%)",ylim=c(-0.05,0.05),xlim=c(-0.02,0.05))
plot(merged_out$PC2,merged_out$PC3,pch=19,col=merged_out$Superpopulation_colors,xlab="PC2 (0.42%)",ylab="PC3 (0.21%)",xlim=c(-0.05,0.05),ylim=c(-0.05,0.065))
plot(merged_out$PC3,merged_out$PC4,pch=19,col=merged_out$Superpopulation_colors,xlab="PC3 (0.21%)",ylab="PC4 (0.16%)",xlim=c(-0.05,0.065),ylim=c(-0.04,0.12))
plot(merged_out$PC4,merged_out$PC5,pch=19,col=merged_out$Superpopulation_colors,xlab="PC4 (0.16%)",ylab="PC5 (0.14%)",xlim=c(-0.04,0.12),ylim=c(-0.25,0.052))
plot(merged_out$PC5,merged_out$PC6,pch=19,col=merged_out$Superpopulation_colors,xlab="PC5",ylab="PC6")
plot(merged_out$PC6,merged_out$PC7,pch=19,col=merged_out$Superpopulation_colors,xlab="PC6",ylab="PC7")



########## LABEL SPECIFICALLY THE SUBSET OF 34 YEAR 1 SAMPLES and 24 YEAR 2 samples
merged_out_yr1 <- merged_out[merged_out$SampleID %in% year_1_samples$YR_1_SAMPLE,]
merged_out_yr2 <- merged_out[merged_out$SampleID %in% year_2_samples$Sample,]

# Combined YR1 and YR2 samples
merged_out_yr1_and_2 <- rbind(merged_out_yr1,merged_out_yr2)


# HPRC subset
merged_out_hprc <- merged_out[merged_out$SampleID %in% hprc_samples$ChildID,]


# Combine both HPRC and HGSVC
combined_hprc_hgsvc <- rbind(merged_out_hprc,merged_out_yr1_and_2)




## Part 1/3: HPRC samples filled in black diamond
par(new=T); plot(merged_out_hprc$PC1,merged_out_hprc$PC2,pch=23,col=merged_out_hprc$Superpopulation_colors,xlab="",ylab="",bg="darkgrey",ylim=c(-0.05,0.03),xlim=c(-0.02,0.05),cex=1.2)

par(new=T); plot(merged_out_hprc$PC2,merged_out_hprc$PC3,pch=23,col=merged_out_hprc$Superpopulation_colors,xlab="",ylab="",bg="grey",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065),cex=1.2)



## Part 2/3: YEAR 1 and 2 Long read sequencing samples annotated as triangles filled in black color
par(new=T); plot(merged_out_yr1_and_2$PC1,merged_out_yr1_and_2$PC2,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",bg="black",ylim=c(-0.05,0.03),xlim=c(-0.02,0.05),cex=1.2)

par(new=T); plot(merged_out_yr1_and_2$PC2,merged_out_yr1_and_2$PC3,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",bg="black",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065),cex=1.2)



## Part 3/3: Centroid approach samples

merged_out_centroidA <- merged_out[merged_out$SampleID%in%centroid_merged_optionA$CENTROID_SAMPLES_MIN,]
merged_out_centroidB <- merged_out[merged_out$SampleID%in%centroid_merged_optionB$CENTROID_SAMPLES_MAX,]

### MIN from centroid samples are CIRCLES
par(new=T); plot(merged_out_centroidA$PC1,merged_out_centroidA$PC2,pch=21,col=merged_out_centroidA$Superpopulation_colors,xlab="",ylab="",ylim=c(-0.05,0.03),xlim=c(-0.02,0.05),cex=1.2,bg="black")

### MAX from centroid samples are squares
par(new=T); plot(merged_out_centroidB$PC1,merged_out_centroidB$PC2,pch=22,col=merged_out_centroidB$Superpopulation_colors,xlab="",ylab="",ylim=c(-0.05,0.03),xlim=c(-0.02,0.05),cex=1.2,bg="black")


### MIN from centroid samples are CIRCLES
par(new=T); plot(merged_out_centroidA$PC2,merged_out_centroidA$PC3,pch=21,col=merged_out_centroidA$Superpopulation_colors,xlab="",ylab="",bg="black",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065),cex=1.2)


### MAX from centroid samples are squares
par(new=T); plot(merged_out_centroidB$PC2,merged_out_centroidB$PC3,pch=22,col=merged_out_centroidB$Superpopulation_colors,xlab="",ylab="",bg="black",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065),cex=1.2)



# Combined HPRC and HGSVC projections

# Take 2 on YR 1 and 2 sample annotation (this time using GATK's genotypes, not pangenies)
par(new=T);plot(merged_out_yr1_and_2$PC1,merged_out_yr1_and_2$PC2,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",ylim=c(-0.05,0.05),xlim=c(-0.02,0.05),bg="black")
par(new=T);plot(merged_out_yr1_and_2$PC2,merged_out_yr1_and_2$PC3,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",xlim=c(-0.05,0.05),ylim=c(-0.05,0.065),bg="black")
par(new=T);plot(merged_out_yr1_and_2$PC3,merged_out_yr1_and_2$PC4,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",xlim=c(-0.05,0.065),ylim=c(-0.04,0.12),bg="black")
par(new=T);plot(merged_out_yr1_and_2$PC4,merged_out_yr1_and_2$PC5,pch=25,col=merged_out_yr1_and_2$Superpopulation_colors,xlab="",ylab="",xlim=c(-0.04,0.12),ylim=c(-0.25,0.052),bg="black")



# YEAR 1 Long read sequencing samples annotated as triangled filled in black color
par(new=T); plot(merged_out_yr1$PC1,merged_out_yr1$PC2,pch=25,col=merged_out_yr1$Superpopulation_colors,xlab="",ylab="",bg="black")
par(new=T); plot(merged_out_yr1$PC2,merged_out_yr1$PC3,pch=25,col=merged_out_yr1$Superpopulation_colors,xlab="",ylab="",bg="black",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065))



# YEAR 2 Long read sequencing samples annotated as triangled filled in black color
par(new=T); plot(merged_out_yr2$PC1,merged_out_yr2$PC2,pch=25,col=merged_out_yr2$Superpopulation_colors,xlab="",ylab="",bg="black",ylim=c(-0.05,0.03),xlim=c(-0.02,0.05))
par(new=T); plot(merged_out_yr2$PC2,merged_out_yr2$PC3,pch=25,col=merged_out_yr2$Superpopulation_colors,xlab="",ylab="",bg="black",xlim=c(-0.05,0.03),ylim=c(-0.05,0.065))





###### CLUSTERING APPROACHES ######

# Start simple: kmeans

# Sutract each superpopulation separately

# Determine nr of clusters for each superpop
merged_kmeans <- kmeans(merged_out[which(merged_out$Superpopulation=="EUR"),2:3],centers = 1)


# Calculate optimal number of centers (k values)
factoextra::fviz_nbclust(x = merged_out[which(merged_out$Superpopulation=="AMR"),2:3], kmeans, method = "gap_stat")


# Plot each kmeans cluster (manually enter optimal number of clusters, i.e., centers)
km.opt.eur <- kmeans(merged_out[which(merged_out$Superpopulation=="EUR"),2:3],centers = 3,nstart=25)
km.opt.afr <- kmeans(merged_out[which(merged_out$Superpopulation=="AFR"),2:3],centers = 7,nstart=25)
km.opt.eas <- kmeans(merged_out[which(merged_out$Superpopulation=="EAS"),2:3],centers = 1,nstart=25)
km.opt.sas <- kmeans(merged_out[which(merged_out$Superpopulation=="SAS"),2:3],centers = 4,nstart=25)
km.opt.amr <- kmeans(merged_out[which(merged_out$Superpopulation=="AMR"),2:3],centers = 6,nstart=25)


fviz_cluster(km.opt.amr, data = merged_out[which(merged_out$Superpopulation=="AMR"),2:3],
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())



#Start of scripts for measuring distance of each sample to its respective centroid

merged_out_distances_from_centroids <- merged_out[,c(1,18:24)]
merged_out_distances_from_centroids$DIST_from_CENTROID <- rep(99999,nrow(merged_out))
merged_out_distances_from_centroids$CENTROID_ID <- rep(0,nrow(merged_out))


# Iterate through the array of bootstrap-optimized kmeans objects
kmeans_array <- c("km.opt.afr","km.opt.eur","km.opt.eas","km.opt.sas","km.opt.amr")


for(i in 1:length(kmeans_array)){
  
  current_kmean_obj <- eval(parse(text=kmeans_array[i]))
  current_clusters <- eval(parse(text=paste0(kmeans_array[i],"$cluster")))

  # Go into each cluster to calculate distance from centroid
  for(j in sort(unique(current_clusters))){
    
  # Create slice for AMR and k-cluster j=1
  slice_idx <- as.numeric(as.character(names(which(current_clusters==j))))
    
  # For each cluster samples in each population, measure distance from the center
  clust_spec_center_coords <- current_kmean_obj$centers[j,]
  
  #sample_coord <- (merged_out[slice_idx,2:3],merged_out[slice_idx,3]),nrow = 2))
  clust_spec_center_coords <- matrix(c(rep(clust_spec_center_coords[1],length(slice_idx)),rep(clust_spec_center_coords[2],length(slice_idx))),ncol=2)
  
  # calculate distance from cluster-center for each appropriate sample
  clust_cpecific_distances <- pointDistance(merged_out[slice_idx,2:3],clust_spec_center_coords,lonlat = F)
  
  # Write to the specific rows and the last two columns the centroid id and distance to it
  merged_out_distances_from_centroids[slice_idx,]$DIST_from_CENTROID <- clust_cpecific_distances
  merged_out_distances_from_centroids[slice_idx,]$CENTROID_ID <- rep(j,length(slice_idx))
  #
  }

}


write.csv(merged_out_distances_from_centroids,"merged_out_distances_from_centroids_21.csv",row.names=F,quote=F)

merged_out_distances_from_centroids <- read.csv("merged_out_distances_from_centroids_21.csv",header=T)
merged_out_distances_from_centroids_removedSelectedSamples <- merged_out_distances_from_centroids[-which(merged_out_distances_from_centroids$SampleID %in% hgsvc_hprs_sampleIDs),]
merged_out_distances_from_centroids <- merged_out_distances_from_centroids_removedSelectedSamples

# Plot chosen samples in the big PCA (visual validation)
# Go back to the plotting section

centroid_merged_optionA <- c(NA)
centroid_merged_optionB <- c(NA)


for (superpop in c("AFR","AMR","EUR","EAS","SAS")){
  
  slice <- merged_out_distances_from_centroids[merged_out_distances_from_centroids$Superpopulation==superpop,]
  for (myc in sort(unique(slice$CENTROID_ID))){
    slice_c <- slice[slice$CENTROID_ID==as.numeric(as.character(myc)),]
    sample_a <- slice_c[which.min(slice_c$DIST_from_CENTROID),]$SampleID
    sample_b <- slice_c[which.max(slice_c$DIST_from_CENTROID),]$SampleID
    
    centroid_merged_optionA <- c(centroid_merged_optionA,sample_a)
    centroid_merged_optionB <- c(centroid_merged_optionB,sample_b)
  }
  
}


centroid_merged_optionA <- data.frame(CENTROID_SAMPLES_MIN=na.omit(as.character(centroid_merged_optionA)))
centroid_merged_optionB <- data.frame(CENTROID_SAMPLES_MAX=na.omit(as.character(centroid_merged_optionB)))

