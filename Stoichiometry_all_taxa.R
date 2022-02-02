#Script for analyzing TPM data

#Name working directory
working_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/Polyploidy_project/Stoicheometry/TPM_analysis/"

#Set working directory
setwd(working_dir)

#devtools::install_github("karthik/wesanderson")

##Load packages
package_list<-c("dplyr", "plyr", "seqinr", "ggplot2", "gplots", "RColorBrewer", "ape", "insect", 
                "phytools", "ggplot2", "grid", "treeio", "phangorn", "Biostrings", "reshape", 
                "wesanderson", "matrixStats", "ggridges", "scales", "gridExtra", "reshape2")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#If installing and librarying the wesanderson color palette didnt work, do the following and then try again
#devtools::install_github("karthik/wesanderson")

#Set the colors for nuc/cp/mt
mt_col<-"#CC6677"
nuc_col<-"#DDCC77"
cp_col<-"#44AA99"

#We will loop through the genera below to analyze each
#Set a string to identify each genus
study_taxa<-c("Arabidopsis", "Arachis", "Chenopodium", "Gossypium")

#Indicate whether the polyploid genome was originally split into two different files (on for each subgenome)
split_subgenomes<-c("TRUE", "FALSE", "FALSE", "TRUE")

###Some random notes about input files
##Arachis:
#only has 4 replicates for the polyploid
#The 'target_id' field  is different because some seq ids are concatinations of two seq id's. I need to figure out that the deal is with that.

##Gossypium:
#AD1=Gossypium hirsutum, AD2=Gossypium barbadense

##Arabidopsis:
#Arabidopsis Aare had 7 replicates. Removed reps 6 and 7 from the tpm file to standardize things across species.

#All genera:
#Diploid 1 = maternal; Diploid 2 = paternal

#Loop through each genus
for(t in 1:length(study_taxa)){
  #Set the directory for the species
  taxon_dir<-paste("Files_", study_taxa[t], sep = "")
  
  #Read in the one-file (this file contains information about quartets)
  onefile<-read.table(file = paste(taxon_dir, "/Onefile.txt", sep = ""), header = TRUE, sep = "\t")
  
  #Read in the tpm file (this files contains the normalized expression for each gene)
  tpm_file<-read.table(file = paste(taxon_dir, "/Tpm.tbl", sep = ""), header = TRUE, sep = "\t")
  
  #Standardize the names in the tpm file
  names(tpm_file)<-c("target_id","Diploid1.1","Diploid1.2","Diploid1.3","Diploid1.4","Diploid1.5",
                     "Polyploid1.1","Polyploid1.2","Polyploid1.3","Polyploid1.4","Polyploid1.5",
                     "Diploid2.1","Diploid2.2","Diploid2.3","Diploid2.4","Diploid2.5")
  
  #Set columns as numeric
  i<-2:16
  tpm_file[,i] <- apply(tpm_file[,i], 2,function(x) as.numeric(as.character(x)))
  
  #Remove the concatenated gene names to only retain the first one
  tpm_file$target_id<-sapply(strsplit(as.character(tpm_file$target_id), "__"), `[`, 1)
  
  ##Get the outlier genes (genes with expression >1000 TPM that are not cp/mt encoded genes)
  #We looked through these manually and removed (see below)
  outliers<-paste(tpm_file$target_id[which(rowMeans(tpm_file[,2:16], na.rm = TRUE)>1000)][-
  grep("cp_|mt_",
  tpm_file$target_id[which(rowMeans(tpm_file[,2:16], na.rm = TRUE)>1000)]
  )])
  
  ## Pull out the sequences for these outliers
  #read in the fasta file (this is a cds file for every gene model in the polyploid)
  cds_file<-read.fasta(paste(taxon_dir, "/polyploid_cds_catfile.fasta", sep = ""))
  #get the seqs that match outliers
  outlier_seqs<-cds_file[grep(paste(outliers, collapse = "\\b|"), names(cds_file))]
  #write fasta of outliers
  write.fasta(sequences = outlier_seqs, names = names(outlier_seqs), file.out = paste(taxon_dir,"/all_outlier_seqs.fasta", sep = ""))

  ### NOTE: At this point we did some manual cleaning to remove the outliers.
  #In brief: we blasted each of the 'outlier' sequences against NCBI. 
  #If the top hits had some sort of non-coding RNA annotation (i.e. tRNA, ncRNA, rRNA), we put the name of that sequence in the file: outliers_to_remove.csv (inside of the particular genus folder)
  #This manual step needs to be done in order for the below code to run
  #A quick and dirty way to do that is run the loop once to generate the all_outlier_seqs.fasta file, manually make the outliers_to_remove.csv, and then run the whole loop again.
  
  #read in csv file for bad outliers
  bad_genes<-paste(read.csv(file = paste(taxon_dir, "/outliers_to_remove.csv", sep = ""), header = TRUE)$Bad_outlier_genes)
  #Remove bad outliers
  tpm_file_outlier<-tpm_file[-grep(paste(bad_genes, collapse = "\\b|"), tpm_file$target_id),]
  #remove rows with NA (note that I am just looking at one column because some species have NA cols)
  tpm_file_outlier<-tpm_file_outlier[which(complete.cases(tpm_file_outlier$Diploid1.1)),]
  #Rescale the file back to TPMs (TPM means columns should always add up to 1 million)
  for(s in 2:ncol(tpm_file_outlier)){
    tpm_file_outlier[,s]<-tpm_file_outlier[,s]/(sum(tpm_file_outlier[,s])/1000000)
  }
  
  #Generate melted file for proportion expression fig
  cp_temp_df<-tpm_file_outlier[grep("cp_", tpm_file_outlier$target_id),]
  mt_temp_df<-tpm_file_outlier[grep("mt_", tpm_file_outlier$target_id),]
  nuc_temp_df<-tpm_file_outlier[-grep("cp_|mt_", tpm_file_outlier$target_id),]
  
  tpm_file_outlier_temp<-rbind(
    cbind(data.frame(Genome=rep("Nuc", nrow(nuc_temp_df))), nuc_temp_df),
    cbind(data.frame(Genome=rep("Mt", nrow(mt_temp_df))), mt_temp_df),
    cbind(data.frame(Genome=rep("Cp", nrow(cp_temp_df))), cp_temp_df)
  )
  
  #Add the averages accross replicates
  tpm_file_outlier_full<-cbind(tpm_file_outlier_temp,
        data.frame(Diploid1_av=rowMeans(tpm_file_outlier_temp[,3:7], na.rm = TRUE), 
                    Polyploid_av=rowMeans(tpm_file_outlier_temp[,8:12], na.rm = TRUE),
                    Diploid2_av=rowMeans(tpm_file_outlier_temp[,13:17], na.rm = TRUE))
  )
  
  #Write rescaled tpm files
  write.csv(tpm_file_outlier_full, file = paste0("Rescaled_tpm_outliersGone_", study_taxa[t], ".csv"), quote = FALSE, row.names = FALSE)
  
  ## Generate date that do not include plastid genes
  #Remove plastid genes
  tpm_file_noCp<-subset(tpm_file_outlier_full, tpm_file_outlier_full$Genome!="Cp")
  
  #Rescale the file back to TPMs
  for(s in 3:ncol(tpm_file_noCp)){
    tpm_file_noCp[,s]<-tpm_file_noCp[,s]/(sum(tpm_file_noCp[,s])/1000000)
  }
            
  #initialize matrix
  props_mat<-matrix(,nrow = (ncol(tpm_file_outlier)-1), ncol = 4)

  #Get the proportions for each sample
  for(x in 2:(ncol(tpm_file_outlier))){
    props_mat[(x-1),]<-c(names(tpm_file_outlier)[x],sum(nuc_temp_df[,x]), sum(cp_temp_df[,x]), sum(mt_temp_df[,x]))
  }

  Sum_genomes_df<-data.frame(Species=sapply(strsplit(as.character(props_mat[,1]), "\\."), `[`, 1),
                             Replicate=sapply(strsplit(as.character(props_mat[,1]), "\\."), `[`, 2),
                             Nuc_sum=as.numeric(paste(props_mat[,2])),
                             Cp_sum=as.numeric(paste(props_mat[,3])),
                             Mt_sum=as.numeric(paste(props_mat[,4]))
  )
  
  #Write genome proportions file
  write.csv(Sum_genomes_df, file = paste0("Genome_proportions_", study_taxa[t], ".csv"), quote = FALSE, row.names = FALSE)
  
  #Melt the summed dataframe
  Sum_genomes_melt<-melt(Sum_genomes_df, id=c("Species", "Replicate"))

  #The order of the levels in the Sum_genomes_melt$variable determine the order of the stacking in the stacked bar plot
  Sum_genomes_melt$variable<-relevel(Sum_genomes_melt$variable, ref = "Nuc_sum")
  Sum_genomes_melt$variable<-relevel(Sum_genomes_melt$variable, ref = "Mt_sum")
  
  ##incorperate the quartet info
  #Note that for our purpuses we want to sum each pair of homeologs together so that we can compare diploid to polyploids
  #Since the mapping for all species was done using the polyploid genome, we need to do this summing for diploids and polyploids
  
  #Ask whether the polyploid subgenomes are split or not (the naming conventions are slightly different depending)
  if(as.logical(split_subgenomes[t])){
  ### Get quartet info
  #Filter 1-file to only contain quartets
  quart_onefile<-onefile[
    which(onefile$numSeqs_paternalDiploid==1 & onefile$numSeqs_paternalTetraploid==1 & onefile$numSeqs_maternalTetraploid==1 & onefile$numSeqs_maternalDiploid==1),
    ]
  #Find the homeolog pairs
  quart_tpm_mat<-matrix(, nrow = nrow(quart_onefile), ncol = 33)
  #loop through quartets
  for(q in 1:nrow(quart_onefile)){
    patPoly_temp<-paste(quart_onefile$paternalTetraploid[q])
    matPoly_temp<-paste(quart_onefile$maternalTetraploid[q])
    #generate matrix
    quart_tpm_mat[q,]<-c(paste(quart_onefile$OGID[q]),
                         patPoly_temp,
                         matPoly_temp,
                         as.numeric(paste(tpm_file_outlier[which(tpm_file_outlier$target_id==patPoly_temp),2:16])),
                         as.numeric(paste(tpm_file_outlier[which(tpm_file_outlier$target_id==matPoly_temp),2:16]))
                         )}
  }else{
    #Filter 1-file to only contain quartets
    quart_onefile<-onefile[
      which(onefile$numSeqs_paternalDiploid==1 & onefile$numSeqs_tetraploid==2 & onefile$numSeqs_maternalDiploid==1),
      ]
    #Find the homeolog pairs
    quart_tpm_mat<-matrix(, nrow = nrow(quart_onefile), ncol = 33)
    #loop through quartets
    #q<-17
    for(q in 1:nrow(quart_onefile)){
      firstPoly_temp<-unlist(strsplit(paste(quart_onefile$tetraploid[q]), ","))[1]
      secondPoly_temp<-unlist(strsplit(paste(quart_onefile$tetraploid[q]), ","))[2]
      #generate matrix
      quart_tpm_mat[q,]<-c(paste(quart_onefile$OGID[q]),
                           firstPoly_temp,
                           secondPoly_temp,
                           as.numeric(paste(tpm_file_outlier[which(tpm_file_outlier$target_id==firstPoly_temp),2:16])),
                           as.numeric(paste(tpm_file_outlier[which(tpm_file_outlier$target_id==secondPoly_temp),2:16])))}
  }#End subgenome split if/else

  #Convert to dataframe
  quart_tpm_df<-as.data.frame(quart_tpm_mat)
  #add names
  names(quart_tpm_df)<-c(
    "OGID", "patPolyploid", "matPolyploid",
    "pat_Diploid1.1","pat_Diploid1.2","pat_Diploid1.3","pat_Diploid1.4","pat_Diploid1.5",
    "pat_Polyploid1.1","pat_Polyploid1.2","pat_Polyploid1.3","pat_Polyploid1.4","pat_Polyploid1.5",
    "pat_Diploid2.1","pat_Diploid2.2","pat_Diploid2.3","pat_Diploid2.4","pat_Diploid2.5",
    "mat_Diploid1.1","mat_Diploid1.2","mat_Diploid1.3","mat_Diploid1.4","mat_Diploid1.5",
    "mat_Polyploid1.1","mat_Polyploid1.2","mat_Polyploid1.3","mat_Polyploid1.4","mat_Polyploid1.5",
    "mat_Diploid2.1","mat_Diploid2.2","mat_Diploid2.3","mat_Diploid2.4","mat_Diploid2.5")

  #Change the numeric cols to class=numeric
  i<-4:33
  quart_tpm_df[,i] <- apply(quart_tpm_df[,i], 2,function(x) as.numeric(as.character(x)))

  #Sum the tpm values mapped to pat and mat homeologs (see rationale above)
  tpm_homeo_sums<-data.frame(
    OGID=paste(quart_tpm_df[,1]),
    patPoly=paste(quart_tpm_df[,2]),
    matPoly=paste(quart_tpm_df[,3]),
    Diploid1.1=quart_tpm_df[,4] + quart_tpm_df[,19],
    Diploid1.2=quart_tpm_df[,5] + quart_tpm_df[,20],
    Diploid1.3=quart_tpm_df[,6] + quart_tpm_df[,21],
    Diploid1.4=quart_tpm_df[,7] + quart_tpm_df[,22],
    Diploid1.5=quart_tpm_df[,8] + quart_tpm_df[,23],
    Polyploid1.1=quart_tpm_df[,9] + quart_tpm_df[,24],
    Polyploid1.2=quart_tpm_df[,10] + quart_tpm_df[,25],
    Polyploid1.3=quart_tpm_df[,11] + quart_tpm_df[,26],
    Polyploid1.4=quart_tpm_df[,12] + quart_tpm_df[,27],
    Polyploid1.5=quart_tpm_df[,13] + quart_tpm_df[,28],
    Diploid2.1=quart_tpm_df[,14] + quart_tpm_df[,29],
    Diploid2.2=quart_tpm_df[,15] + quart_tpm_df[,30],
    Diploid2.3=quart_tpm_df[,16] + quart_tpm_df[,31],
    Diploid2.4=quart_tpm_df[,17] + quart_tpm_df[,32],
    Diploid2.5=quart_tpm_df[,18] + quart_tpm_df[,33]
  )
  #remove rows with NA (note that I am just looking at one column for NAs)
  tpm_homeo_sums<-tpm_homeo_sums[which(complete.cases(tpm_homeo_sums$Diploid1.1)),]

  #Add cymira data (this is in the 'one-file')
  tpm_file_quartets_prescale<-left_join(tpm_homeo_sums, quart_onefile)

  ## RESCALE (no orgaanelles)
  #rename object
  tpm_file_quartets_temp<-tpm_file_quartets_prescale
  #Rescale quartet file without organelle genes
  for(s in 4:18){
    tpm_file_quartets_temp[,s]<-tpm_file_quartets_temp[,s]/(sum(tpm_file_quartets_temp[,s])/1000000)
  }
  
  #Add averages
  tpm_file_quartets<-cbind(tpm_file_quartets_temp,
                               data.frame(Diploid1_av=rowMeans(tpm_file_quartets_temp[,4:8], na.rm = TRUE), 
                                          Polyploid_av=rowMeans(tpm_file_quartets_temp[,9:13], na.rm = TRUE),
                                          Diploid2_av=rowMeans(tpm_file_quartets_temp[,14:18], na.rm = TRUE),
                                          data.frame(Targeting_short=sapply(strsplit(as.character(tpm_file_quartets_temp$Targeting), "_"), `[`, 1))
                                          )
  )
  
  #Reorder the levels of the Targeting (this is for aesthetics when plotting)
  tpm_file_quartets$Targeting_short<-relevel(tpm_file_quartets$Targeting_short, ref = "Dual-targeted")
  tpm_file_quartets$Targeting_short<-relevel(tpm_file_quartets$Targeting_short, ref = "Mitochondria-targeted")
  tpm_file_quartets$Targeting_short<-relevel(tpm_file_quartets$Targeting_short, ref = "Plastid-targeted")
  
  #Combine the quartet and organelle data and rescale (note that tpm_homeo_sums and tpm_file_outlier are both rescaled after outliers so they're still compatible with eachother)
  org_temp<-tpm_file_outlier[grep("cp_|mt_", tpm_file_outlier$target_id),]
  tpm_file_quartets_organells<-rbind(tpm_homeo_sums,
  cbind(data.frame(OGID=org_temp$target_id, patPoly=org_temp$target_id, matPoly=org_temp$target_id),org_temp[,2:16])
  )

  #Get the organelle tpm files
  cp_tpm<-tpm_file_quartets_organells[grep("cp_", tpm_file_quartets_organells$OGID),]
  mt_tpm<-tpm_file_quartets_organells[grep("mt_", tpm_file_quartets_organells$OGID),]

  ## RESCALE (with orgaanelles)
  for(s in 4:18){
    tpm_file_quartets_organells[,s]<-tpm_file_quartets_organells[,s]/(sum(tpm_file_quartets_organells[,s])/1000000)
  }

  #Add cymira data to the df
  tpm_file_quartets_organells_cymira_temp<-left_join(tpm_file_quartets_organells, quart_onefile)
  
  #Add averages
  tpm_file_quartets_organells_cymira<-cbind(tpm_file_quartets_organells_cymira_temp,
                           data.frame(Diploid1_av=rowMeans(tpm_file_quartets_organells_cymira_temp[,4:8], na.rm = TRUE), 
                                      Polyploid_av=rowMeans(tpm_file_quartets_organells_cymira_temp[,9:13], na.rm = TRUE),
                                      Diploid2_av=rowMeans(tpm_file_quartets_organells_cymira_temp[,14:18], na.rm = TRUE),
                                      data.frame(Targeting_short=sapply(strsplit(as.character(tpm_file_quartets_organells_cymira_temp$Targeting), "_"), `[`, 1))
                                      )
  )
  
  # #Reorder the levels of the Targeting
  # tpm_file_quartets_organells_cymira$Targeting_short<-relevel(tpm_file_quartets_organells_cymira$Targeting_short, ref = "Dual-targeted")
  # tpm_file_quartets_organells_cymira$Targeting_short<-relevel(tpm_file_quartets_organells_cymira$Targeting_short, ref = "Mitochondria-targeted")
  # tpm_file_quartets_organells_cymira$Targeting_short<-relevel(tpm_file_quartets_organells_cymira$Targeting_short, ref = "Plastid-targeted")
  # 

  ## Pull out the cytonuclear interactions
  #Read in cymira associations file
  #This file was generated manually. 
  #Cymira only covers nuclear genes, so we needed a reference to describe what function each organelle gene is involved in (for most it's obvious)
  #Read file
  int_df<-read.csv(file = "Interactions.csv", header = TRUE)

  #Make empty matrix
  int_res_mat<-matrix(,ncol = 31, nrow = length(levels(int_df$Cymira_association)))
  
  #Kind of a wonky loop (but it works!) to get the average expression for each nuclear and organllar gene associated with each functional category
  for(i in 1:length(levels(int_df$Cymira_association))){
    keeper_og_df<-filter(tpm_file_quartets_organells_cymira,grepl(paste(int_df$Organelle_gene[which(int_df$Cymira_association==levels(int_df$Cymira_association)[i])], collapse="|"), tpm_file_quartets_organells_cymira$OGID))
    keeper_nuc_df<-tpm_file_quartets_organells_cymira[which(tpm_file_quartets_organells_cymira$CyMIRA_interaction==paste(levels(int_df$Cymira_association)[i])),]

    int_res_mat[i,]<-c(paste(levels(int_df$Cymira_association)[i]),
                       colMeans(keeper_nuc_df[,4:18]),
                       colMeans(keeper_og_df[,4:18])
    )}
  
  #Convert to df
  int_res_df<-as.data.frame(int_res_mat)

  #Rename columns
  names(int_res_df)<-c("Cymira_category",
                       "nuc_Diploid1.1","nuc_Diploid1.2","nuc_Diploid1.3","nuc_Diploid1.4","nuc_Diploid1.5",
                       "nuc_Polyploid1.1","nuc_Polyploid1.2","nuc_Polyploid1.3","nuc_Polyploid1.4","nuc_Polyploid1.5",
                       "nuc_Diploid2.1","nuc_Diploid2.2","nuc_Diploid2.3","nuc_Diploid2.4","nuc_Diploid2.5",
                       "org_Diploid1.1","org_Diploid1.2","org_Diploid1.3","org_Diploid1.4","org_Diploid1.5",
                       "org_Polyploid1.1","org_Polyploid1.2","org_Polyploid1.3","org_Polyploid1.4","org_Polyploid1.5",
                       "org_Diploid2.1","org_Diploid2.2","org_Diploid2.3","org_Diploid2.4","org_Diploid2.5")
  
  #Shorten the long name (for aesthetics with plotting)
  int_res_df$Cymira_category<-paste(int_res_df$Cymira_category)
  int_res_df$Cymira_category[which(int_res_df$Cymira_category=="mt-Transcription_and_Transcript_Maturation")]<-"mt-Transcription"
  int_res_df$Cymira_category[which(int_res_df$Cymira_category=="pt-Transcription_and_Transcript_Maturation")]<-"pt-Transcription"
  
  #Change numeric cols to class=numeric
  i<-2:31
  int_res_df[,i] <- apply(int_res_df[,i], 2,function(x) as.numeric(as.character(x)))
  
  #Write the interactions dataset (used for figure emulating Havird and Sloan 2016 [Fig. 3])
  write.csv(int_res_df, file = paste0("Cytonuclear_interaction_averages_", study_taxa[t], ".csv"), quote = FALSE, row.names = FALSE)
  
  cytonuc_ratios_df<-data.frame(Cymira_category=paste(int_res_df[,1]),
  Diploid1.1=int_res_df[,17]/int_res_df[,2],
  Diploid1.2=int_res_df[,18]/int_res_df[,3],
  Diploid1.3=int_res_df[,19]/int_res_df[,4],
  Diploid1.4=int_res_df[,20]/int_res_df[,5],
  Diploid1.5=int_res_df[,21]/int_res_df[,6],
  Polyploid1.1=int_res_df[,22]/int_res_df[,7],
  Polyploid1.2=int_res_df[,23]/int_res_df[,8],
  Polyploid1.3=int_res_df[,24]/int_res_df[,9],
  Polyploid1.4=int_res_df[,25]/int_res_df[,10],
  Polyploid1.5=int_res_df[,26]/int_res_df[,11],
  Diploid2.1=int_res_df[,27]/int_res_df[,12],
  Diploid2.2=int_res_df[,28]/int_res_df[,13],
  Diploid2.3=int_res_df[,29]/int_res_df[,14],
  Diploid2.4=int_res_df[,30]/int_res_df[,15],
  Diploid2.5=int_res_df[,31]/int_res_df[,16])
  
  #rearrange this df
  cytonuc_ratios_melt<-melt(cytonuc_ratios_df, id.vars = "Cymira_category")
  
  #Add a species column
  cytonuc_ratios_melt_plus<-cbind(cytonuc_ratios_melt, data.frame(Species=sapply(strsplit(as.character(cytonuc_ratios_melt[,2]), "\\."), `[`, 1)))

  #Rearrange to do the Havird and Sloan-like figure
  melt_int<-melt(int_res_df, id.vars = c("Cymira_category"))
  
  #Add Genome and Species cols 
  havird_df<-cbind(melt_int,
                   data.frame(Genome=sapply(strsplit(as.character(melt_int[,2]), "_"), `[`, 1),
                              Species=sapply(strsplit(sapply(strsplit(as.character(melt_int[,2]), "\\."), `[`, 1), "_"),`[`, 2)))

  #Change level order (for aesthetic)
  havird_df$Genome<-relevel(havird_df$Genome, ref = "org") #save this object

  #Prepare data for a cytonuclear correlation plot
  cytonuc_correlate_temp<-cbind(subset(havird_df, havird_df$Genome=="org"), subset(havird_df, havird_df$Genome=="nuc"))
  cytonuc_correlate<-cytonuc_correlate_temp[,c(1, 2, 3, 5, 8)]
  names(cytonuc_correlate)<-c("Cymira_category", "Sample_name_temp", "Org_average_tpm", "Species", "Nuc_average_tpm")
  cytonuc_correlate<-cbind(cytonuc_correlate, data.frame(Sample_name=sapply(strsplit(as.character(cytonuc_correlate$Sample_name_temp), "_"), `[`, 2),Replicate=sapply(strsplit(as.character(cytonuc_correlate$Sample_name_temp), "\\."), `[`, 2)))[,-2]

  #Get the averages accross replicates
  cytonuc_correlate_avs<-data.frame(
    Cymira_category=paste(int_res_df[,1]),
    nuc_Diploid1_av=rowMeans(int_res_df[,2:6], na.rm = TRUE),
    nuc_Polyploid_av=rowMeans(int_res_df[,7:11], na.rm = TRUE),
    nuc_Diploid2_av=rowMeans(int_res_df[,12:16], na.rm = TRUE),
    org_Diploid1_av=rowMeans(int_res_df[,17:21], na.rm = TRUE),
    org_Polyploid_av=rowMeans(int_res_df[,22:26], na.rm = TRUE),
    org_Diploid2_av=rowMeans(int_res_df[,27:31], na.rm = TRUE),
    Organelle=sapply(strsplit(as.character(int_res_df$Cymira_category), "-"), `[`, 1)
  )

  #Cytonuclear ratios averaged
  cytonuc_correlate_avs_ratios<-data.frame(
    Cymira_category=paste(cytonuc_correlate_avs[,1]),
    Diploid1_ratio=cytonuc_correlate_avs[,5]/cytonuc_correlate_avs[,2],
    Polyploid_ratio=cytonuc_correlate_avs[,6]/cytonuc_correlate_avs[,3],
    Diploid2_ratio=cytonuc_correlate_avs[,7]/cytonuc_correlate_avs[,4],
    Organelle=cytonuc_correlate_avs[,8]
  )
  
  #polyploid/diploid ratios (nuclear vs organellar)
  polyploid_ratios_cytonuc<-data.frame(
    Cymira_category=paste(cytonuc_correlate_avs[,1]),
    polyploid_diploid1_nuc=cytonuc_correlate_avs[,3]/cytonuc_correlate_avs[,2],
    polyploid_diploid2_nuc=cytonuc_correlate_avs[,3]/cytonuc_correlate_avs[,4],
    polyploid_diploid1_org=cytonuc_correlate_avs[,6]/cytonuc_correlate_avs[,5],
    polyploid_diploid2_org=cytonuc_correlate_avs[,6]/cytonuc_correlate_avs[,7],
    diploid2_diploid1_nuc=cytonuc_correlate_avs[,4]/cytonuc_correlate_avs[,2],
    diploid2_diploid1_org=cytonuc_correlate_avs[,7]/cytonuc_correlate_avs[,5],
    Organelle=cytonuc_correlate_avs[,8]
  )
  
  ###Assign all of the output DFs as new objects (with names specifying the genus)
  assign(paste0("tpm_file_outlier_", study_taxa[t]), tpm_file_outlier_full)
  assign(paste0("tpm_file_noCp_", study_taxa[t]), tpm_file_noCp)
  assign(paste0("sum_genomes_melt_", study_taxa[t]), Sum_genomes_melt)
  assign(paste0("tpm_file_quartets_", study_taxa[t]), tpm_file_quartets)
  assign(paste0("tpm_file_quartets_organells_cymira_", study_taxa[t]), tpm_file_quartets_organells_cymira)
  assign(paste0("cp_tpm_", study_taxa[t]), cp_tpm)
  assign(paste0("mt_tpm_", study_taxa[t]), mt_tpm)
  assign(paste0("int_res_df_", study_taxa[t]), int_res_df)
  assign(paste0("cytonuc_ratios_melt_plus_", study_taxa[t]), cytonuc_ratios_melt_plus)
  assign(paste0("havird_df_", study_taxa[t]), havird_df)
  assign(paste0("cytonuc_correlate_", study_taxa[t]), cytonuc_correlate)
  assign(paste0("cytonuc_correlate_avs_", study_taxa[t]), cytonuc_correlate_avs)
  assign(paste0("cytonuc_correlate_avs_ratios_", study_taxa[t]), cytonuc_correlate_avs_ratios)
  assign(paste0("polyploid_ratios_cytonuc_", study_taxa[t]), polyploid_ratios_cytonuc)

  }#End study taxon loop
  

########################################################
########################################################
########################################################
########### Make plots with all species ################
########################################################
########################################################
########################################################

#### Plot the proportions from each genome
#Store as objects and plot as multi-panel plot (using grid.arrange() below)
g1<-ggplot(sum_genomes_melt_Arabidopsis, aes(fill=variable, y=value, x=Replicate))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("Arabidopsis")+
  scale_fill_manual(values = c(mt_col, nuc_col, cp_col))+
  theme_classic()+
  facet_wrap(Species~.)

g2<-ggplot(sum_genomes_melt_Arachis, aes(fill=variable, y=value, x=Replicate))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("Arachis")+
  scale_fill_manual(values = c(mt_col, nuc_col, cp_col))+
  theme_classic()+
  facet_wrap(Species~.)

g3<-ggplot(sum_genomes_melt_Chenopodium, aes(fill=variable, y=value, x=Replicate))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("Chenopodium")+
  scale_fill_manual(values = c(mt_col, nuc_col, cp_col))+
  theme_classic()+
  facet_wrap(Species~.)

g4<-ggplot(sum_genomes_melt_Gossypium, aes(fill=variable, y=value, x=Replicate))+
  geom_bar(position="fill", stat="identity")+
  ggtitle("Gossypium")+
  scale_fill_manual(values = c(mt_col, nuc_col, cp_col))+
  theme_classic()+
  facet_wrap(Species~.)

#Make the multi-panel plot
grid.arrange(g1, g2, g3, g4, ncol=1)

#Make distribution plots (Note, we omitted nuclear genes for aesthetics)
p1<-ggplot(tpm_file_outlier_Arabidopsis, aes(x=Diploid2_av, y=Genome, fill=Genome))+
  geom_density_ridges(rel_min_height = 0.01, alpha = 0.7)+
  geom_density_ridges(data =subset(tpm_file_outlier_Arabidopsis, Genome!="Nuc"), rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2, alpha = 0.7)+
  scale_fill_manual(values = c(nuc_col, mt_col, cp_col))+
  scale_x_continuous(trans="pseudo_log", breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(-5,500000))+
  ggtitle("Arabidopsis")+ xlab("expression/gene (tpm)")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 8, angle = 90), axis.title = element_text(size = 10))+
  theme_classic()

p2<-ggplot(tpm_file_outlier_Arachis, aes(x=Diploid2_av, y=Genome, fill=Genome))+
  geom_density_ridges(rel_min_height = 0.01, alpha = 0.7)+
  geom_density_ridges(data =subset(tpm_file_outlier_Arachis, Genome!="Nuc"), rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2, alpha = 0.7)+
  scale_fill_manual(values = c(nuc_col, mt_col, cp_col))+
  scale_x_continuous(trans="pseudo_log", breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(-5,500000))+
  ggtitle("Arachis")+ xlab("expression/gene (tpm)")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 8, angle = 90), axis.title = element_text(size = 10))+
  theme_classic()

p3<-ggplot(tpm_file_outlier_Chenopodium, aes(x=Diploid2_av, y=Genome, fill=Genome))+
  geom_density_ridges(rel_min_height = 0.01, alpha = 0.7)+
  geom_density_ridges(data =subset(tpm_file_outlier_Chenopodium, Genome!="Nuc"), rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2, alpha = 0.7)+
  scale_fill_manual(values = c(nuc_col, mt_col, cp_col))+
  scale_x_continuous(trans="pseudo_log", breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(-5,500000))+
  ggtitle("Chenopodium")+ xlab("expression/gene (tpm)")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 8, angle = 90), axis.title = element_text(size = 10))+
  theme_classic()

p4<-ggplot(tpm_file_outlier_Gossypium, aes(x=Diploid2_av, y=Genome, fill=Genome))+
  geom_density_ridges(rel_min_height = 0.01, alpha = 0.7)+
  geom_density_ridges(data =subset(tpm_file_outlier_Gossypium, Genome!="Nuc"), rel_min_height = 0.01, jittered_points = TRUE,
                      position = position_points_jitter(width = 0.5, height = 0),
                      point_shape = "|", point_size = 2, alpha = 0.7)+
  scale_fill_manual(values = c(nuc_col, mt_col, cp_col))+
  scale_x_continuous(trans="pseudo_log", breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(-5,500000))+
  ggtitle("Gossypium")+ xlab("expression/gene (tpm)")+
  theme(legend.title = element_blank(), axis.text = element_text(size = 8, angle = 90), axis.title = element_text(size = 10))+
  theme_classic()

#Make the multi-panel plot
grid.arrange(p1, p2, p3, p4, ncol=1)


###Combined cp and mt plots###

#Arabidopsis
cp_subset_Arabidopsis<-subset(tpm_file_outlier_Arabidopsis, tpm_file_outlier_Arabidopsis$Genome=="Cp")
mt_subset_Arabidopsis<-subset(tpm_file_noCp_Arabidopsis, tpm_file_noCp_Arabidopsis$Genome=="Mt")
tpm_file_combo_Arabidopsis<-rbind(cp_subset_Arabidopsis, mt_subset_Arabidopsis)

cm01<-ggplot(tpm_file_combo_Arabidopsis, aes(y=Polyploid_av, x=Diploid1_av))+
  #ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm02<-ggplot(tpm_file_combo_Arabidopsis, aes(y=Polyploid_av, x=Diploid2_av))+
  #ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Paternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm03<-ggplot(tpm_file_combo_Arabidopsis, aes(y=Diploid2_av, x=Diploid1_av))+
  #ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arabidopsis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Paternal diploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

#Arachis
cp_subset_Arachis<-subset(tpm_file_outlier_Arachis, tpm_file_outlier_Arachis$Genome=="Cp")
mt_subset_Arachis<-subset(tpm_file_noCp_Arachis, tpm_file_noCp_Arachis$Genome=="Mt")
tpm_file_combo_Arachis<-rbind(cp_subset_Arachis, mt_subset_Arachis)

cm04<-ggplot(tpm_file_combo_Arachis, aes(y=Polyploid_av, x=Diploid1_av))+
  #ggtitle("Arachis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm05<-ggplot(tpm_file_combo_Arachis, aes(y=Polyploid_av, x=Diploid2_av))+
  #ggtitle("Arachis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Paternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm06<-ggplot(tpm_file_combo_Arachis, aes(y=Diploid2_av, x=Diploid1_av))+
  #ggtitle("Arachis")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Arachis, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Paternal diploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

#Chenopodium
cp_subset_Chenopodium<-subset(tpm_file_outlier_Chenopodium, tpm_file_outlier_Chenopodium$Genome=="Cp")
mt_subset_Chenopodium<-subset(tpm_file_noCp_Chenopodium, tpm_file_noCp_Chenopodium$Genome=="Mt")
tpm_file_combo_Chenopodium<-rbind(cp_subset_Chenopodium, mt_subset_Chenopodium)

cm07<-ggplot(tpm_file_combo_Chenopodium, aes(y=Polyploid_av, x=Diploid1_av))+
  #ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm08<-ggplot(tpm_file_combo_Chenopodium, aes(y=Polyploid_av, x=Diploid2_av))+
  #ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Paternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm09<-ggplot(tpm_file_combo_Chenopodium, aes(y=Diploid2_av, x=Diploid1_av))+
  #ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Chenopodium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Paternal diploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

#Gossypium
cp_subset_Gossypium<-subset(tpm_file_outlier_Gossypium, tpm_file_outlier_Gossypium$Genome=="Cp")
mt_subset_Gossypium<-subset(tpm_file_noCp_Gossypium, tpm_file_noCp_Gossypium$Genome=="Mt")
tpm_file_combo_Gossypium<-rbind(cp_subset_Gossypium, mt_subset_Gossypium)

cm10<-ggplot(tpm_file_combo_Gossypium, aes(y=Polyploid_av, x=Diploid1_av))+
  #ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm11<-ggplot(tpm_file_combo_Gossypium, aes(y=Polyploid_av, x=Diploid2_av))+
  #ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Polyploid")+xlab("Paternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

cm12<-ggplot(tpm_file_combo_Gossypium, aes(y=Diploid2_av, x=Diploid1_av))+
  #ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Genome)), size=1, alpha=0.7)+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(10,1000000), breaks = c(1, 10, 100, 1000, 10000, 100000,1000000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Cp"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_combo_Gossypium, Genome== "Mt"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  ylab("Paternal diploid")+xlab("Maternal diploid")+
  theme(legend.position = "none", panel.background = element_blank())

#Make the multi-panel plot
grid.arrange(cm01, cm02,cm03,  
             cm04, cm05, cm06,
             cm07, cm08, cm09,
             cm10, cm11, cm12,
             ncol=3)

### Plot nuclear gene expression (quartets) in diploid vs polyploid
#Arabidopsis
tpm_file_quartets_Arabidopsis<-rbind(
  subset(tpm_file_quartets_Arabidopsis, Targeting_short=="Not-organelle-targeted"),
  subset(tpm_file_quartets_Arabidopsis, Targeting_short=="Dual-targeted"),
  subset(tpm_file_quartets_Arabidopsis, Targeting_short=="Mitochondria-targeted"),
  subset(tpm_file_quartets_Arabidopsis, Targeting_short=="Plastid-targeted")
         )

n1<-ggplot(tpm_file_quartets_Arabidopsis, aes(y=Polyploid_av, x=Diploid1_av))+
  ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n2<-ggplot(tpm_file_quartets_Arabidopsis, aes(y=Polyploid_av, x=Diploid2_av))+
  ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Paternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n3<-ggplot(tpm_file_quartets_Arabidopsis, aes(y=Diploid2_av, x=Diploid1_av))+
  ggtitle("Arabidopsis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Paternal diploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arabidopsis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

#Arachis
tpm_file_quartets_Arachis<-rbind(
  subset(tpm_file_quartets_Arachis, Targeting_short=="Not-organelle-targeted"),
  subset(tpm_file_quartets_Arachis, Targeting_short=="Dual-targeted"),
  subset(tpm_file_quartets_Arachis, Targeting_short=="Mitochondria-targeted"),
  subset(tpm_file_quartets_Arachis, Targeting_short=="Plastid-targeted")
)

n4<-ggplot(tpm_file_quartets_Arachis, aes(y=Polyploid_av, x=Diploid1_av))+
  ggtitle("Arachis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n5<-ggplot(tpm_file_quartets_Arachis, aes(y=Polyploid_av, x=Diploid2_av))+
  ggtitle("Arachis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Paternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n6<-ggplot(tpm_file_quartets_Arachis, aes(y=Diploid2_av, x=Diploid1_av))+
  ggtitle("Arachis")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Paternal diploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Arachis, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

#Chenopodium
tpm_file_quartets_Chenopodium<-rbind(
  subset(tpm_file_quartets_Chenopodium, Targeting_short=="Not-organelle-targeted"),
  subset(tpm_file_quartets_Chenopodium, Targeting_short=="Dual-targeted"),
  subset(tpm_file_quartets_Chenopodium, Targeting_short=="Mitochondria-targeted"),
  subset(tpm_file_quartets_Chenopodium, Targeting_short=="Plastid-targeted")
)

n7<-ggplot(tpm_file_quartets_Chenopodium, aes(y=Polyploid_av, x=Diploid1_av))+
  ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n8<-ggplot(tpm_file_quartets_Chenopodium, aes(y=Polyploid_av, x=Diploid2_av))+
  ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Paternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n9<-ggplot(tpm_file_quartets_Chenopodium, aes(y=Diploid2_av, x=Diploid1_av))+
  ggtitle("Chenopodium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Paternal diploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Chenopodium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

#Gossypium
tpm_file_quartets_Gossypium<-rbind(
  subset(tpm_file_quartets_Gossypium, Targeting_short=="Not-organelle-targeted"),
  subset(tpm_file_quartets_Gossypium, Targeting_short=="Dual-targeted"),
  subset(tpm_file_quartets_Gossypium, Targeting_short=="Mitochondria-targeted"),
  subset(tpm_file_quartets_Gossypium, Targeting_short=="Plastid-targeted")
)

n10<-ggplot(tpm_file_quartets_Gossypium, aes(y=Polyploid_av, x=Diploid1_av))+
  ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n11<-ggplot(tpm_file_quartets_Gossypium, aes(y=Polyploid_av, x=Diploid2_av))+
  ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Paternal diploid")+ylab("Polyploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())

n12<-ggplot(tpm_file_quartets_Gossypium, aes(y=Diploid2_av, x=Diploid1_av))+
  ggtitle("Gossypium")+
  geom_point(aes(colour = factor(Targeting_short)), size=0.5, alpha=0.4)+
  xlab("Maternal diploid")+ylab("Paternal diploid")+
  scale_color_manual(values = c(cp_col, mt_col, "brown", "light gray"))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Plastid-targeted"), method='lm', se=FALSE, colour=cp_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Mitochondria-targeted"), method='lm', se=FALSE, colour=mt_col, size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Dual-targeted"), method='lm', se=FALSE, colour="brown", size = 0.5)+
  geom_smooth(data=subset(tpm_file_quartets_Gossypium, Targeting_short== "Not-organelle-targeted"), method='lm', se=FALSE, colour="light gray", size = 0.5)+
  theme(legend.position = "none", panel.background = element_blank())


#Print big grid of figures
grid.arrange(n1, n2 ,n3,  
             n4, n5 , n6,
             n7, n8 , n9,
             n10, n11 , n12,
             ncol=3)


### Make cytonuclear interaction figure similar to Havird et al

#Arabidopsis
h1<-ggplot(havird_df_Arabidopsis, aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size=1.5, alpha=0.6, stroke=0) +
  ggtitle("Arabidopsis")+
  theme_bw() + scale_color_manual(values = c("black", nuc_col)) + 
  ylab("TPM") + 
  stat_summary(havird_df_Arabidopsis, geom = "point", fun = "mean", size = 6, shape = 95, position=position_dodge(width = 0.8), mapping=aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  theme(legend.position="right", axis.text = element_text(size = 6, angle = 90), axis.title = element_text(size = 8), legend.text = element_text(size = 8, face="italic"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(Species ~ .)

#Arachis
h2<-ggplot(havird_df_Arachis, aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size=1.5, alpha=0.6, stroke=0) +
  ggtitle("Arachis")+
  theme_bw() + scale_color_manual(values = c("black", nuc_col)) + 
  ylab("TPM") + 
  stat_summary(havird_df_Arachis, geom = "point", fun = "mean", size = 6, shape = 95, position=position_dodge(width = 0.8), mapping=aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  theme(legend.position="right", axis.text = element_text(size = 6, angle = 90), axis.title = element_text(size = 8), legend.text = element_text(size = 8, face="italic"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(Species ~ .)

#Chenopodium
h3<-ggplot(havird_df_Chenopodium, aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size=1.5, alpha=0.6, stroke=0) +
  ggtitle("Chenopodium")+
  theme_bw() + scale_color_manual(values = c("black", nuc_col)) + 
  ylab("TPM") + 
  stat_summary(havird_df_Chenopodium, geom = "point", fun = "mean", size = 6, shape = 95, position=position_dodge(width = 0.8), mapping=aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  theme(legend.position="right", axis.text = element_text(size = 6, angle = 90), axis.title = element_text(size = 8), legend.text = element_text(size = 8, face="italic"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(Species ~ .)

#Gossypium
h4<-ggplot(havird_df_Gossypium, aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size=1.5, alpha=0.6, stroke=0) +
  ggtitle("Gossypium")+
  theme_bw() + scale_color_manual(values = c("black", nuc_col)) + 
  ylab("TPM") + 
  stat_summary(havird_df_Gossypium, geom = "point", fun = "mean", size = 6, shape = 95, position=position_dodge(width = 0.8), mapping=aes(y=value, x=Cymira_category, group=Genome, color=Genome)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  theme(legend.position="right", axis.text = element_text(size = 6, angle = 90), axis.title = element_text(size = 8), legend.text = element_text(size = 8, face="italic"), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  facet_wrap(Species ~ .)

#Print big grid of figures
grid.arrange(h1, h2 ,h3, h4, ncol=1)

###Correlation plot of org vs nuc
###

#Arabidopsis Dip1
i1<-ggplot(cytonuc_correlate_avs_Arabidopsis, aes(x=nuc_Diploid1_av, y=org_Diploid1_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arabidopsis; Diploid1")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid1_av, y=org_Diploid1_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()
  
#Arabidopsis Poly
i2<-ggplot(cytonuc_correlate_avs_Arabidopsis, aes(x=nuc_Polyploid_av, y=org_Polyploid_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arabidopsis; Polyploid")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Polyploid_av, y=org_Polyploid_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Arabidopsis Dip2
i3<-ggplot(cytonuc_correlate_avs_Arabidopsis, aes(x=nuc_Diploid2_av, y=org_Diploid2_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arabidopsis; Diploid2")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid2_av, y=org_Diploid2_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()


### ARACHIS
#Arachis Dip1
i4<-ggplot(cytonuc_correlate_avs_Arachis, aes(x=nuc_Diploid1_av, y=org_Diploid1_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arachis; Diploid1")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid1_av, y=org_Diploid1_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Arachis Poly
i5<-ggplot(cytonuc_correlate_avs_Arachis, aes(x=nuc_Polyploid_av, y=org_Polyploid_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arachis; Polyploid")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Polyploid_av, y=org_Polyploid_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Arachis Dip2
i6<-ggplot(cytonuc_correlate_avs_Arachis, aes(x=nuc_Diploid2_av, y=org_Diploid2_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Arachis; Diploid2")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid2_av, y=org_Diploid2_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

###CHENOPODIUM
#Chenopodium Dip1
i7<-ggplot(cytonuc_correlate_avs_Chenopodium, aes(x=nuc_Diploid1_av, y=org_Diploid1_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Chenopodium; Diploid1")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid1_av, y=org_Diploid1_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Chenopodium Poly
i8<-ggplot(cytonuc_correlate_avs_Chenopodium, aes(x=nuc_Polyploid_av, y=org_Polyploid_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Chenopodium; Polyploid")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Polyploid_av, y=org_Polyploid_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Chenopodium Dip2
i9<-ggplot(cytonuc_correlate_avs_Chenopodium, aes(x=nuc_Diploid2_av, y=org_Diploid2_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Chenopodium; Diploid2")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid2_av, y=org_Diploid2_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()


###GOSSYPIUM
#Gossypium Dip1
i10<-ggplot(cytonuc_correlate_avs_Gossypium, aes(x=nuc_Diploid1_av, y=org_Diploid1_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Gossypium; Diploid1")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid1_av, y=org_Diploid1_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Gossypium Poly
i11<-ggplot(cytonuc_correlate_avs_Gossypium, aes(x=nuc_Polyploid_av, y=org_Polyploid_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Gossypium; Polyploid")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Polyploid_av, y=org_Polyploid_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

#Gossypium Dip2
i12<-ggplot(cytonuc_correlate_avs_Gossypium, aes(x=nuc_Diploid2_av, y=org_Diploid2_av, color=Organelle, shape=Cymira_category))+
  ggtitle("Gossypium; Diploid2")+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,100000), breaks = c(1, 10, 100, 1000, 10000, 100000))+
  xlab("Nuclear expression (TPM)") + 
  ylab("Organellar expression (TPM)")+ 
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(x=nuc_Diploid2_av, y=org_Diploid2_av), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  theme_classic()

###

#Print large grid of figures
grid.arrange(i1, i2 ,i3,  
             i4, i5 , i6,
             i7, i8 , i9,
             i10, i11 , i12,
             ncol=3)

#Some alternative ways of plotting
#grid.arrange(i2, i5, i8, i11, ncol=4)

# grid.arrange(i1, i4, i7, i10, 
#              i3, i6, i9, i12, ncol=4)

##### 'Ratio of ratios' plots #####
##organellar/nuclear
#Arabidopsis
r1<-ggplot(cytonuc_correlate_avs_ratios_Arabidopsis, aes(y=Polyploid_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 50))+
  scale_y_continuous(limits = c(10, 50))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r2<-ggplot(cytonuc_correlate_avs_ratios_Arabidopsis, aes(y=Polyploid_ratio, x=Diploid2_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Paternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(0, 150))+
  scale_y_continuous(limits = c(0, 150))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid2_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r3<-ggplot(cytonuc_correlate_avs_ratios_Arabidopsis, aes(y=Diploid2_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Paternal diploid")+
  scale_x_continuous(limits = c(0, 150))+
  scale_y_continuous(limits = c(0, 150))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Diploid2_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,1])/
                                                                         summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arabidopsis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                       )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Arachis
r4<-ggplot(cytonuc_correlate_avs_ratios_Arachis, aes(y=Polyploid_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 150))+
  scale_y_continuous(limits = c(10, 150))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r5<-ggplot(cytonuc_correlate_avs_ratios_Arachis, aes(y=Polyploid_ratio, x=Diploid2_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Paternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 150))+
  scale_y_continuous(limits = c(10, 150))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid2_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r6<-ggplot(cytonuc_correlate_avs_ratios_Arachis, aes(y=Diploid2_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Paternal diploid")+
  scale_x_continuous(limits = c(10, 150))+
  scale_y_continuous(limits = c(10, 150))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Diploid2_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,1])/
                                                                      summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Arachis))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Chenopodium
r7<-ggplot(cytonuc_correlate_avs_ratios_Chenopodium, aes(y=Polyploid_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 100))+
  scale_y_continuous(limits = c(10, 100))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r8<-ggplot(cytonuc_correlate_avs_ratios_Chenopodium, aes(y=Polyploid_ratio, x=Diploid2_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Paternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 100))+
  scale_y_continuous(limits = c(10, 100))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid2_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r9<-ggplot(cytonuc_correlate_avs_ratios_Chenopodium, aes(y=Diploid2_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Paternal diploid")+
  scale_x_continuous(limits = c(10, 100))+
  scale_y_continuous(limits = c(10, 100))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Diploid2_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Chenopodium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Gossypium
r10<-ggplot(cytonuc_correlate_avs_ratios_Gossypium, aes(y=Polyploid_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 300))+
  scale_y_continuous(limits = c(10, 300))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r11<-ggplot(cytonuc_correlate_avs_ratios_Gossypium, aes(y=Polyploid_ratio, x=Diploid2_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Paternal diploid")+
  ylab("Polyploid")+
  scale_x_continuous(limits = c(10, 175))+
  scale_y_continuous(limits = c(10, 175))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Polyploid_ratio, x=Diploid2_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

r12<-ggplot(cytonuc_correlate_avs_ratios_Gossypium, aes(y=Diploid2_ratio, x=Diploid1_ratio, color=Organelle, shape=Cymira_category))+
  geom_point(size = 2)+
  scale_shape_manual(values=c(0,1,2,0,1,2,5,6))+
  scale_color_manual(values = c(mt_col, cp_col))+
  xlab("Maternal diploid")+
  ylab("Paternal diploid")+
  scale_x_continuous(limits = c(10, 300))+
  scale_y_continuous(limits = c(10, 300))+
  geom_abline(intercept = 0, slope = 1)+
  geom_smooth(aes(y=Diploid2_ratio, x=Diploid1_ratio), method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, inherit.aes = FALSE)+
  geom_text(size=1, col="black", aes(x = -Inf, y = Inf, hjust = -0.5, vjust = 1,
                                     label = paste0(
                                       "R-squared: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$r.squared,
                                       "\nSlope: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2],
                                       "\nSlope=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,4],
                                       "\nIntercept=0 P-val: ", summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[1,4],
                                       "\nSlope=1 P-val: ", 2*pt(q=
                                                                   (abs(1-summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,1])/
                                                                      summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = cytonuc_correlate_avs_ratios_Gossypium))$coefficients[2,2]),
                                                                 df=6, lower.tail=FALSE)
                                     )))+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend


#Print big grid of figures
grid.arrange(r1, r2, r3,
             r4, r5, r6,
             r7, r8, r9,
             r10, r11, r12,
             ncol=3)


### Polyploid ratios
#Arabidopsis
pp1<-ggplot(polyploid_ratios_cytonuc_Arabidopsis, aes(x=polyploid_diploid1_nuc, y=polyploid_diploid1_org, color=Cymira_category))+
  #ggtitle("Arabidopsis")+
  xlab("Polyploid nuclear / Diploid1 nuclear")+
  ylab("Polyploid organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp2<-ggplot(polyploid_ratios_cytonuc_Arabidopsis, aes(x=polyploid_diploid2_nuc, y=polyploid_diploid2_org, color=Cymira_category))+
  #ggtitle("Arabidopsis")+
  xlab("Polyploid nuclear / Diploid2 nuclear")+
  ylab("Polyploid organellar / Diploid2 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp3<-ggplot(polyploid_ratios_cytonuc_Arabidopsis, aes(x=diploid2_diploid1_nuc, y=diploid2_diploid1_org, color=Cymira_category))+
  #ggtitle("Arabidopsis")+
  xlab("Diploid2 nuclear / Diploid1 nuclear")+
  ylab("Diploid2 organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Arachis
pp4<-ggplot(polyploid_ratios_cytonuc_Arachis, aes(x=polyploid_diploid1_nuc, y=polyploid_diploid1_org, color=Cymira_category))+
  #ggtitle("Arachis")+
  xlab("Polyploid nuclear / Diploid1 nuclear")+
  ylab("Polyploid organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp5<-ggplot(polyploid_ratios_cytonuc_Arachis, aes(x=polyploid_diploid2_nuc, y=polyploid_diploid2_org, color=Cymira_category))+
  #ggtitle("Arachis")+
  xlab("Polyploid nuclear / Diploid2 nuclear")+
  ylab("Polyploid organellar / Diploid2 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp6<-ggplot(polyploid_ratios_cytonuc_Arachis, aes(x=diploid2_diploid1_nuc, y=diploid2_diploid1_org, color=Cymira_category))+
  #ggtitle("Arachis")+
  xlab("Diploid2 nuclear / Diploid1 nuclear")+
  ylab("Diploid2 organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Chenopodium
pp7<-ggplot(polyploid_ratios_cytonuc_Chenopodium, aes(x=polyploid_diploid1_nuc, y=polyploid_diploid1_org, color=Cymira_category))+
  #ggtitle("Chenopodium")+
  xlab("Polyploid nuclear / Diploid1 nuclear")+
  ylab("Polyploid organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp8<-ggplot(polyploid_ratios_cytonuc_Chenopodium, aes(x=polyploid_diploid2_nuc, y=polyploid_diploid2_org, color=Cymira_category))+
  #ggtitle("Chenopodium")+
  xlab("Polyploid nuclear / Diploid2 nuclear")+
  ylab("Polyploid organellar / Diploid2 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp9<-ggplot(polyploid_ratios_cytonuc_Chenopodium, aes(x=diploid2_diploid1_nuc, y=diploid2_diploid1_org, color=Cymira_category))+
  #ggtitle("Chenopodium")+
  xlab("Diploid2 nuclear / Diploid1 nuclear")+
  ylab("Diploid2 organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Gossypium
pp10<-ggplot(polyploid_ratios_cytonuc_Gossypium, aes(x=polyploid_diploid1_nuc, y=polyploid_diploid1_org, color=Cymira_category))+
  #ggtitle("Gossypium")+
  xlab("Polyploid nuclear / Diploid1 nuclear")+
  ylab("Polyploid organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp11<-ggplot(polyploid_ratios_cytonuc_Gossypium, aes(x=polyploid_diploid2_nuc, y=polyploid_diploid2_org, color=Cymira_category))+
  #ggtitle("Gossypium")+
  xlab("Polyploid nuclear / Diploid2 nuclear")+
  ylab("Polyploid organellar / Diploid2 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

pp12<-ggplot(polyploid_ratios_cytonuc_Gossypium, aes(x=diploid2_diploid1_nuc, y=diploid2_diploid1_org, color=Cymira_category))+
  #ggtitle("Gossypium")+
  xlab("Diploid2 nuclear / Diploid1 nuclear")+
  ylab("Diploid2 organellar / Diploid1 organellar")+
  scale_x_continuous(limits = c(0, 2))+
  scale_y_continuous(limits = c(0, 2))+
  geom_point(size = 1)+
  #geom_text(aes(label=Cymira_category),hjust=1, vjust=1, size=2)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 1)+
  geom_smooth(method='lm', formula= y~x, colour="Black", size=0.5, se = FALSE, linetype="dashed")+
  theme(panel.background = element_blank(), axis.title = element_text(size = 8), axis.text = element_text(size = 8), axis.line = element_line(colour = "black"), legend.position = "none") #Remove legend.position argument to show legend

#Print big grid of figures
grid.arrange(pp1, pp2, pp3,
             pp4, pp5, pp6,
             pp7, pp8, pp9,
             pp10, pp11, pp12,
             ncol=3)

############################################################
#################        DATA TABLES       #################
############################################################

## Generate the statistics/data tables for various figures
#Figure 5
for(t in 1:length(study_taxa)){

mat_temp<-matrix(, ncol = 8, nrow = 3)

#Get LM
lm_pd1<-summary(lm(formula =  Polyploid_ratio~Diploid1_ratio, data = eval(as.name(paste0("cytonuc_correlate_avs_ratios_", study_taxa[t])))))
lm_pd2<-summary(lm(formula =  Polyploid_ratio~Diploid2_ratio, data = eval(as.name(paste0("cytonuc_correlate_avs_ratios_", study_taxa[t])))))
lm_d1d2<-summary(lm(formula =  Diploid2_ratio~Diploid1_ratio, data = eval(as.name(paste0("cytonuc_correlate_avs_ratios_", study_taxa[t])))))

#Genus
mat_temp[,1]<-c(study_taxa[t], study_taxa[t], study_taxa[t])

#y variable
mat_temp[,2]<-c("Polyploid", "Polyploid", "Paternal_diploid")

#x variable
mat_temp[,3]<-c("Maternal_diploid", "Paternal_diploid", "Maternal_diploid")

#Rsquared
mat_temp[,4]<-c(lm_pd1$r.squared, lm_pd2$r.squared, lm_d1d2$r.squared)

#Slope
mat_temp[,5]<-c(lm_pd1$coefficients[2], lm_pd2$coefficients[2], lm_d1d2$coefficients[2])

#P-value (slope = 0)
mat_temp[,6]<-c(lm_pd1$coefficients[2,4], lm_pd2$coefficients[2,4], lm_d1d2$coefficients[2,4])

#P-value (slope = 1)
mat_temp[,7]<-c(2*pt(q=(abs(1-lm_pd1$coefficients[2,1])/lm_pd1$coefficients[2,2]), df=6, lower.tail=FALSE),
  2*pt(q=(abs(1-lm_pd2$coefficients[2,1])/lm_pd2$coefficients[2,2]), df=6, lower.tail=FALSE),
  2*pt(q=(abs(1-lm_d1d2$coefficients[2,1])/lm_d1d2$coefficients[2,2]), df=6, lower.tail=FALSE))

#P-value (intercept = 0)
mat_temp[,8]<-c(lm_pd1$coefficients[1,4], lm_pd2$coefficients[1,4], lm_d1d2$coefficients[1,4])

df_temp<-as.data.frame(mat_temp)

names(df_temp)<-c("Genus", "y-variable", "x-variable", "R2", "Slope", "P-value (slope = 0)", "P-value (slope = 1)", "P-value (intercept = 0)")

#Write it
write.table(df_temp, file = "Fig_5_LM_stats.csv", quote = FALSE, sep = ",", row.names = FALSE, append = TRUE)
}


#Figure 3b
for(t in 1:length(study_taxa)){
  
  temp_df<-eval(as.name(paste0("cytonuc_correlate_avs_", study_taxa[t])))
  
  temp_df[,2:7]<-log10(temp_df[,2:7])
  
  mat_temp<-matrix(, ncol = 7, nrow = 3)
  
  #Get LM
  lm_pol<-summary(lm(formula =  org_Polyploid_av~nuc_Polyploid_av, data = temp_df))
  lm_mat<-summary(lm(formula =  org_Diploid1_av~nuc_Diploid1_av, data = temp_df))
  lm_pat<-summary(lm(formula =  org_Diploid2_av~nuc_Diploid2_av, data = temp_df))
  
  #Genus
  mat_temp[,1]<-c(study_taxa[t], study_taxa[t], study_taxa[t])
  
  #Species
  mat_temp[,2]<-c("Polyploid", "Maternal_diploid", "Paternal_diploid")
  
  #Rsquared
  mat_temp[,3]<-c(lm_pol$r.squared, lm_mat$r.squared, lm_pat$r.squared)
  
  #Slope
  mat_temp[,4]<-c(lm_pol$coefficients[2], lm_mat$coefficients[2], lm_pat$coefficients[2])
  
  #Intercept
  mat_temp[,5]<-c(lm_pol$coefficients[1,1], lm_mat$coefficients[1,1], lm_pat$coefficients[1,1])
  
  #P-value (slope = 0)
  mat_temp[,6]<-c(lm_pol$coefficients[2,4], lm_mat$coefficients[2,4], lm_pat$coefficients[2,4])
  
  #P-value (slope = 1)
  mat_temp[,7]<-c(2*pt(q=(abs(1-lm_pol$coefficients[2,1])/lm_pol$coefficients[2,2]), df=6, lower.tail=FALSE),
                  2*pt(q=(abs(1-lm_mat$coefficients[2,1])/lm_mat$coefficients[2,2]), df=6, lower.tail=FALSE),
                  2*pt(q=(abs(1-lm_pat$coefficients[2,1])/lm_pat$coefficients[2,2]), df=6, lower.tail=FALSE))
  
  # #P-value (intercept = 0)
  # mat_temp[,7]<-c(lm_pol$coefficients[1,4], lm_mat$coefficients[1,4], lm_pat$coefficients[1,4])
  # 
  
  df_temp<-as.data.frame(mat_temp)
  
  names(df_temp)<-c("Genus", "Species", "R2", "Slope", "Intercept", "P-value (slope = 0)", "P-value (slope = 1)")
  
  #Write it
  write.table(df_temp, file = "Fig_3b_LM_stats.csv", quote = FALSE, sep = ",", row.names = FALSE, append = TRUE)
}

#Figure 4
#tpm_file_combo_
mat_temp2<-matrix(, ncol = 7, nrow = 4)

for(t in 1:length(study_taxa)){
  
  combo_df<-eval(as.name(paste0("tpm_file_combo_", study_taxa[t])))[,c(1,2,18,19,20)]
  combo_df[,3:5]<-log10(combo_df[,3:5])
  cp_temp<-subset(combo_df, Genome=="Cp")
  mt_temp<-subset(combo_df, Genome=="Mt")
  
  mat_temp2[t,]<-c(study_taxa[t],
  summary(lm(formula =  Polyploid_av~Diploid1_av, data = cp_temp))$r.squared,
  summary(lm(formula =  Polyploid_av~Diploid2_av, data = cp_temp))$r.squared,
  summary(lm(formula =  Diploid1_av~Diploid2_av, data = cp_temp))$r.squared,
  summary(lm(formula =  Polyploid_av~Diploid1_av, data = mt_temp))$r.squared,
  summary(lm(formula =  Polyploid_av~Diploid2_av, data = mt_temp))$r.squared,
  summary(lm(formula =  Diploid1_av~Diploid2_av, data = mt_temp))$r.squared)

}
res_df2<-as.data.frame(mat_temp2)

names(res_df2)<-c("Genus", "CP_Poly_vs_MatDip","CP_Poly_vs_PatDip", "CP_PatDip_vs_MatDip",
                    "MT_Poly_vs_MatDip","MT_Poly_vs_PatDip", "MT_PatDip_vs_MatDip")

#Write it
write.table(res_df2, file = "Fig_4_R2_values.csv", quote = FALSE, sep = ",", row.names = FALSE, append = FALSE)

#### Figure S5
#tpm_file_quartets_

mat_temp3<-matrix(, ncol = 16, nrow = 4)

for(t in 1:length(study_taxa)){
  
  df3<-eval(as.name(paste0("tpm_file_quartets_", study_taxa[t])))
  
  df4<-df3[,c(which(names(df3)=="OGID"),
    which(names(df3)=="Diploid1_av"),
    which(names(df3)=="Polyploid_av"),
    which(names(df3)=="Diploid2_av"),
    which(names(df3)=="Targeting_short")
    )]
  
  df4[,2:4]<-log10(df4[,2:4])
  
  df4$Diploid1_av[which(df4$Diploid1_av==-Inf)]<-0
  df4$Diploid2_av[which(df4$Diploid2_av==-Inf)]<-0
  df4$Polyploid_av[which(df4$Polyploid_av==-Inf)]<-0
  
  Ncp_temp<-subset(df4, Targeting_short=="Plastid-targeted")
  Nmt_temp<-subset(df4, Targeting_short=="Mitochondria-targeted")
  Ndual_temp<-subset(df4, Targeting_short=="Dual-targeted")
  Nnot_temp<-subset(df4, Targeting_short=="Not-organelle-targeted")
  
  
  mat_temp3[t,]<-c(study_taxa[t],
                   summary(lm(formula =  Polyploid_av~Diploid1_av, data = df4))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid1_av, data = Ncp_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid1_av, data = Nmt_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid1_av, data = Ndual_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid1_av, data = Nnot_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid2_av, data = df4))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid2_av, data = Ncp_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid2_av, data = Nmt_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid2_av, data = Ndual_temp))$r.squared,
                   summary(lm(formula =  Polyploid_av~Diploid2_av, data = Nnot_temp))$r.squared,
                   summary(lm(formula =  Diploid2_av~Diploid1_av, data = df4))$r.squared,
                   summary(lm(formula =  Diploid2_av~Diploid1_av, data = Ncp_temp))$r.squared,
                   summary(lm(formula =  Diploid2_av~Diploid1_av, data = Nmt_temp))$r.squared,
                   summary(lm(formula =  Diploid2_av~Diploid1_av, data = Ndual_temp))$r.squared,
                   summary(lm(formula =  Diploid2_av~Diploid1_av, data = Nnot_temp))$r.squared)
}
  
resdf4<-as.data.frame(mat_temp3)      

names(resdf4)<-c("Genus", "AllNuc_Poly_vs_MatDip", "nCP_Poly_vs_MatDip", "nMT_Poly_vs_MatDip", "nDual_Poly_vs_MatDip", "Not_Poly_vs_MatDip",
                 "AllNuc_Poly_vs_PatDip", "nCP_Poly_vs_PatDip", "nMT_Poly_vs_PatDip", "nDual_Poly_vs_PatDip", "Not_Poly_vs_PatDip",
                 "AllNuc_PatDip_vs_MatDip", "nCP_PatDip_vs_MatDip", "nMT_PatDip_vs_MatDip", "nDual_PatDip_vs_MatDip", "Not_PatDip_vs_MatDip")
                 
#Write it
write.table(resdf4, file = "Fig_S5_R2_values.csv", quote = FALSE, sep = ",", row.names = FALSE, append = FALSE)


