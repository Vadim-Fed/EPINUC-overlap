#..............Read to PEAK OVERLAP ANALYSIS............................................
library(dplyr)
library(readr)

Chromosome_length <- as.data.frame(read_csv("Insert path to chromosome lenght file.csv",col_names = F)) #Load this file as well
names(Chromosome_length)<-c("chr","chr_length")
a<-list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","Y")

####.....Organizing Datasets....####

Ab_Plasma_reads<-read.delim("Path to EPINuc antibody reads.BED",header = F)[1:3]
names(Ab_Plasma_reads)<-c("chr","Start","End")
Ab_Plasma_reads<-arrange(Ab_Plasma_reads,chr)

#Colon
Colon1<-read.delim("Path to 1st ENCODE Colon chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3]
Colon2<-read.delim("Path to 2nd ENCODE Colon chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Colon1)<-c("chr","Start","End")
names(Colon2)<-c("chr","Start","End")
Colon<-rbind(Colon1,Colon2)
Colon<-arrange(Colon,chr)

Overlappp<-function(y){
  Plasma<-subset(Colon1,chr==paste0("chr",y))
  Tissue<-subset(Colon2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Colon_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Lung
Lung1<-read.delim("Path to 1st ENCODE Lung chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
Lung2<-read.delim("Path to 2nd ENCODE Lung chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Lung1)<-c("chr","Start","End")
names(Lung2)<-c("chr","Start","End")
Lung<-rbind(Lung1,Lung2)
Lung<-arrange(Lung,chr)

Overlappp<-function(y){
  Plasma<-subset(Lung1,chr==paste0("chr",y))
  Tissue<-subset(Lung2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Lung_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Heart
Heart1<-read.delim("Path to 1st ENCODE Heart chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
Heart2<-read.delim("Path to 2nd ENCODE Heart chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Heart1)<-c("chr","Start","End")
names(Heart2)<-c("chr","Start","End")
Heart<-rbind(Heart1,Heart2)
Heart<-arrange(Heart,chr)

Overlappp<-function(y){
  Plasma<-subset(Heart1,chr==paste0("chr",y))
  Tissue<-subset(Heart2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Heart_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Brain
Brain1<-read.delim("Path to 1st ENCODE Brain chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
Brain2<-read.delim("Path to 2nd ENCODE Brain chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Brain1)<-c("chr","Start","End")
names(Brain2)<-c("chr","Start","End")
Brain<-rbind(Brain1,Brain2)
Brain<-arrange(Brain,chr)

Overlappp<-function(y){
  Plasma<-subset(Brain1,chr==paste0("chr",y))
  Tissue<-subset(Brain2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Brain_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#BM
BM1<-read.delim("Path to 1st ENCODE Bone-marrow chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
BM2<-read.delim("Path to 2nd ENCODE Bone-marrow chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(BM1)<-c("chr","Start","End")
names(BM2)<-c("chr","Start","End")
BM<-rbind(BM1,BM2)
BM<-arrange(BM,chr)

Overlappp<-function(y){
  Plasma<-subset(BM1,chr==paste0("chr",y))
  Tissue<-subset(BM2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
BM_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Breast
Breast1<-read.delim("Path to 1st ENCODE Breast chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
Breast2<-read.delim("Path to 2nd ENCODE Breast chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Breast1)<-c("chr","Start","End")
names(Breast2)<-c("chr","Start","End")
Breast<-rbind(Breast1,Breast2)
Breast<-arrange(Breast,chr)

Overlappp<-function(y){
  Plasma<-subset(Breast1,chr==paste0("chr",y))
  Tissue<-subset(Breast2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Breast_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#B-cell
BC1<-read.delim("Path to 1st ENCODE B-cell chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
BC2<-read.delim("Path to 2nd ENCODE B-cell chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(BC1)<-c("chr","Start","End")
names(BC2)<-c("chr","Start","End")
BC<-rbind(BC1,BC2)
BC<-arrange(BC,chr)

Overlappp<-function(y){
  Plasma<-subset(BC1,chr==paste0("chr",y))
  Tissue<-subset(BC2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
BC_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Liver
Liver1<-read.delim("Path to 1st ENCODE Liver chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
Liver2<-read.delim("Path to 2nd ENCODE Liver chip-seq dataset.BED compatible to EPINuc antibody",header = F)[1:3] 
names(Liver1)<-c("chr","Start","End")
names(Liver2)<-c("chr","Start","End")
Liver<-rbind(Liver1,Liver2)
Liver<-arrange(Liver,chr)

Overlappp<-function(y){
  Plasma<-subset(Liver1,chr==paste0("chr",y))
  Tissue<-subset(Liver2,chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
Liver_overlaping_Peaks<-do.call(rbind,lapply(a, Overlappp))

#Aggregating data
Tissues<-list(Colon,Lung,Heart,Brain,BM,Breast,BC,Liver)
names(Tissues)<-c("Colon","Lung","Heart","Brain","BM","Breast","BC","Liver")

Overlaping_Peaks<-list(Colon_overlaping_Peaks,Lung_overlaping_Peaks,Heart_overlaping_Peaks,Brain_overlaping_Peaks,
                       BM_overlaping_Peaks,Breast_overlaping_Peaks,BC_overlaping_Peaks,Liver_overlaping_Peaks)
names(Overlaping_Peaks)<-names(Tissues)

###..............Generating Tissue unique profile......####

w<-list()
Unique_profiles<-function(k){
  Tissue_interest_unique<-Tissues[[names(Tissues)[k]]]
  Tissue_no_interest<-Tissues[names(Tissues)[-k]]
  
  for (z in 1:7) {
    Unique<-function(y){
      v<-list()
      Original<-subset(as.data.frame(Tissue_interest_unique),chr==paste0("chr",y))
      Subtract<-subset(as.data.frame(Tissue_no_interest[[z]]),chr==paste0("chr",y))
      
      Start_Overlap<-lapply(1:length(Original$Start), function(x) sum(Subtract$Start<Original$Start[x]&Original$Start[x]<Subtract$End))
      Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
      
      Middle_Overlap<-lapply(1:length(Original$Start), function(x) sum(Original$Start[x]>Subtract$Start&Original$End[x]<Subtract$End))
      Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
      
      Middle_Overlap2<-lapply(1:length(Original$Start), function(x) sum(Original$Start[x]<Subtract$Start&Original$End[x]>Subtract$End))
      Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
      
      End_Overlap<-lapply(1:length(Original$End), function(x) sum(Subtract$Start<Original$End[x]&Original$End[x]<Subtract$End))
      End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
      
      
      v[[y]]<-Original[-(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
    }
    Tissue_interest_unique<-do.call(rbind,lapply(a, Unique))
    #print(dim(Tissue_interest_unique))
  } 
  w[[k]]<-Tissue_interest_unique
}

All_Tissue_unique_profiles<-lapply(1:8,Unique_profiles)
names(All_Tissue_unique_profiles)<-names(Tissues)

###.......Plasma reads intersect with unique profiles.....###
J<-list()
Intersect<-function(z){
Overlappp<-function(y){
  v<-list()
  Plasma<-subset(Ab_Plasma_reads,chr==paste0("chr",y))
  Tissue<-subset(All_Tissue_unique_profiles[[z]],chr==paste0("chr",y))
  
  Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
  Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
  Plasma[which(Start_Overlap!=0),]
  
  Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
  Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
  Plasma[which(Middle_Overlap!=0),]
  
  Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
  Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
  
  End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
  End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
  Plasma[which(End_Overlap!=0),]
  
  v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
}
J[[z]]<-do.call(rbind,lapply(a, Overlappp))
}
Plasma_Intersect<-lapply(1:8, Intersect)
names(Plasma_Intersect)<-names(Tissues)

###...QC Step - Testing and removing duplicated intersection of Plasma reads with chip-seq peaks
#Overlap between chip-seq shared peaks and plasma intersected reads
j<-list()
Intersect<-function(z){
  Overlappp<-function(y){
    Plasma<-subset(Overlaping_Peaks[[z]],chr==paste0("chr",y))
    Tissue<-subset(Plasma_Intersect[[z]],chr==paste0("chr",y))
    
    Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
    Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
    Plasma[which(Start_Overlap!=0),]
    
    Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
    Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
    Plasma[which(Middle_Overlap!=0),]
    
    Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
    Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
    
    End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
    End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
    Plasma[which(End_Overlap!=0),]
    
    v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
  }
  j[[z]]<-do.call(rbind,lapply(a, Overlappp))
}
Plasma_Chip_dual_overlap<-lapply(1:8, Intersect)
names(Plasma_Chip_dual_overlap)<-names(Tissues)

#Overlap between plasma reads intersected shared peaks and unique profiles 
j<-list()
Intersect<-function(z){
  Overlappp<-function(y){
    Plasma<-subset(Plasma_Chip_dual_overlap[[z]],chr==paste0("chr",y))
    Tissue<-subset(All_Tissue_unique_profiles[[z]],chr==paste0("chr",y))
    
    Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
    Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
    Plasma[which(Start_Overlap!=0),]
    
    Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
    Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
    Plasma[which(Middle_Overlap!=0),]
    
    Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
    Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
    
    End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
    End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
    Plasma[which(End_Overlap!=0),]
    
    v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
  }
  j[[z]]<-do.call(rbind,lapply(a, Overlappp))
}
Plasma_Chip_overlap_Unique<-lapply(1:8, Intersect)
names(Plasma_Chip_overlap_Unique)<-names(Tissues)

###Plasma intersect values corrected
Plasma_Intersect_values<-unlist(lapply(Plasma_Intersect, nrow))
double_intersect<-unlist(lapply(Plasma_Chip_overlap_Unique, nrow))
Plasma_Intersect_values_Corrected<-Plasma_Intersect_values-double_intersect

#######..... Assessing overlap significance....######
iter=10000 #number of iterations ##############Change to 10K
h<-list()
By_chance<-function(k){
Bootstrap<-function(y){
  
  W<-list()
  Sampling<-function(y){
    DF<-data.frame(chr=rep(paste0("chr",y),unname(table(Ab_Plasma_reads$chr)[names(table(Ab_Plasma_reads$chr))==paste0("chr",y)])),
                   Start=sample(1:Chromosome_length[Chromosome_length$chr==y,2],unname(table(Ab_Plasma_reads$chr)[names(table(Ab_Plasma_reads$chr))==paste0("chr",y)])))
    DF<-arrange(DF,Start)           
    DF$End<-DF$Start+220
    W[[y]]<-DF
  }
  
  Sample.DF<-do.call(rbind,lapply(a, Sampling))
  
  v<-list()
  Overlappp<-function(y){
    Plasma<-subset(Sample.DF,chr==paste0("chr",y))
    Tissue<-subset(All_Tissue_unique_profiles[[k]],chr==paste0("chr",y))
    
    Start_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Tissue$Start<Plasma$Start[x]&Plasma$Start[x]<Tissue$End))
    Start_Overlap<-as.data.frame(do.call(rbind,Start_Overlap))
    Plasma[which(Start_Overlap!=0),]
    
    Middle_Overlap<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]>Tissue$Start&Plasma$End[x]<Tissue$End))
    Middle_Overlap<-as.data.frame(do.call(rbind,Middle_Overlap))
    Plasma[which(Middle_Overlap!=0),]
    
    Middle_Overlap2<-lapply(1:length(Plasma$Start), function(x) sum(Plasma$Start[x]<Tissue$Start&Plasma$End[x]>Tissue$End))
    Middle_Overlap2<-as.data.frame(do.call(rbind,Middle_Overlap2))
    
    End_Overlap<-lapply(1:length(Plasma$End), function(x) sum(Tissue$Start<Plasma$End[x]&Plasma$End[x]<Tissue$End))
    End_Overlap<-as.data.frame(do.call(rbind,End_Overlap))
    Plasma[which(End_Overlap!=0),]
    
    v[[y]]<-Plasma[unique(c(which(Start_Overlap!=0),which(End_Overlap!=0),which(Middle_Overlap!=0),which(Middle_Overlap2!=0))),]
  }
  
  h[[y]]<-do.call(rbind,lapply(a, Overlappp))
}
P<-lapply(1:iter, Bootstrap)
P1<-as.data.frame(do.call(rbind,lapply(1:length(P), function (x) dim(P[[x]])[1])))[,1]
mean=mean(P1)
Sd=sd(P1)
n=iter
xbar=unname(Plasma_Intersect_values_Corrected[k])
z=(xbar-mean)/(Sd)
Pvalue<-2*pnorm(-abs(z)) #two tail
hist(P1,main=paste(names(Tissues)[k],"random overlap","","p=",round(Pvalue,15)),
     xlab="Overlap",xlim = c(min(P1)-IQR(P1)*1.5,max(P1)+IQR(P1)*1.5))
abline(v=xbar, col="red", lwd=3, lty=2)
}

lapply(1:8, By_chance)


