###################################################
### This program is to check overlap between UW dataset and PHSKC dataset
### The goal is to take two separate data sources of HIV gene sequences linked with some demographic
### variables and make a merged dataset that includes all overlaps and all unique individuals.
### For right now we are only going to compare UW and PHSKC, although we also have the WA state data.
##################################################

# Set your working directory and/or path the data files

# One way to do it
cur.dir <-"I:/groups/SOM/ICRC/DataTeam/Genomics/SeattleHIVSources_Josh/"
setwd(cur.dir)

# Another way to do it
# this is my own path, you will need to edit this
path <- "~/Dropbox/R01_Seattle/RecordLinkage"
# And then you load the files using that path to the directory
load(paste0(path, "/alignments and trees.rdata"))

##################################
## UW data
## There are N individuals that have a ProRT seq in this file.
####################################

# wherever you have stored the file locally  
uw.rtseq <- read.csv(file=" ... /ProtRt_Sequence_2018.csv", header=T)
  
uw.rt <- uw.rtseq[, c(1,4)]
names(uw.rt) <- c("ptid", "PrtRT.seq")

##assign a seqid, the dataset does not have it.
uw.rt$SeqID <- paste(uw.rt$ptid, rownames(uw.rt), sep="-")
dim(uw.rt)
#[1] 4745    3  # these are old results

# Again, wherever you have stored the file locally 
uw.person = read.csv(file="... Demographic_2018.csv", header=T)
dim(uw.person) 
#[1] 2795   11 
uw.pt <- uw.person[, c("?..PtID","BirthYear","PresentSex")]
names(uw.pt) <- c("ptid","birth.year","sex")

##merge seq data with demo info
uw.rt1 <- merge(uw.rt, uw.pt, by="ptid",all.x=T)
dim(uw.rt1)
#[1] 4745    5

uw.rt1 <- uw.rt1[, c("ptid","PrtRT.seq","birth.year","sex","SeqID")]


###########################################
## PHSKC data 
###########################################

# wherever you have stored the file locally 
# Note that this file is made by the other script, which takes the FASTA sequence file
# and makes it into a .csv file. I don't think this is the best way to do this... I'd rather
# have the UW .csv file of sequences be made into a FASTA file, actually.
ph.rtseq <- read.csv(".... PHSKC_PR_RTSeq.csv",header=T)

# wherever you have stored the file locally  
meta = readRDS(file="... sequences_meta.rds")
meta <- meta[, c("seqID", "newnum")]
names(meta) <- c("SeqID","ptid")
ph.rtseq <- merge(ph.rtseq, meta, by="SeqID", all.x=T)
dim(ph.rtseq)
#[1] 16227     3

# wherever you have stored the file locally  
person <- readRDS(file=".... person.rds")
ph.pt <- person[, c("newnum", "b_yr","sex")]
names(ph.pt) <- c("ptid","birth.year","sex")

##merge seq with demo info
ph.rt1 <- merge(ph.rtseq, ph.pt, by="ptid",all.x=T)
dim(ph.rt1)
#[1] 16227     5
  
length(unique(ph.rt1$ptid))
#[1] 9445


##ph.rt1$SeqID <- NULL
ph.rt1 <- ph.rt1[, c("ptid","PrtRT.seq","birth.year","sex","SeqID")]

#############################
########Start to compare uw.rt1 & ph.rt1 by (birth year, sex, seq)
#############################
library(RecordLinkage)

dim(uw.rt1)
#[1] 4745    4
length(unique(uw.rt1$ptid))
#[1] 2737

dim(ph.rt1)
#[1] 16227     5
length(unique(ph.rt1$ptid))
#[1] 9445


#compare uw rt seq with phskc rt seq by (birth year, sex and PR/RT seq)
res <- compare.linkage (uw.rt1, ph.rt1, blockfld=c(2,3,4),phonetic=c(3,4))
matchones <- res$pairs
dim(matchones)
#[1] 3868    8

table(matchones$sex, useNA='always')
# 1 <NA> 
#3868    0 
table(matchones$birth.year, useNA='always')
#   1 <NA> 
#3868    0 
table(matchones$PrtRT.seq, useNA='always')
#   1 <NA> 
#3868    0 

#########start to retrieve overlap ones
uw.rt1$row.indx <- as.numeric(rownames(uw.rt1))
names(uw.rt1) <- c("uw.ptid","uw.rtseq","uw.byear","uw.sex","uw.seqid","uw.row.indx")

matched1 <- merge(matchones, uw.rt1, by.x="id1",by.y="uw.row.indx",all.x=T)

###########
ph.rt1$row.indx <- as.numeric(rownames(ph.rt1))
names(ph.rt1) <- c("ph.ptid","ph.rtseq","ph.byear","ph.sex","ph.seqid","ph.row.indx")

matched2 <- merge(matched1, ph.rt1, by.x="id2",by.y="ph.row.indx",all.x=T)

##found 3868 matched seq, 2368 individuals between UW RT seq and PHSKC RT seq
table(matched2$PrtRT.seq, useNA='always')
#   1 <NA> 
#3868    0 

#> length(unique(out$ph.ptid))
#[1] 2368
#> length(unique(out$uw.ptid))
#[1] 2368

out <- matched2[, c("uw.seqid","uw.ptid","uw.byear","uw.sex","uw.rtseq","ph.seqid","ph.ptid","ph.byear","ph.sex","ph.rtseq")]
#> out[as.character(out$uw.rtseq) != as.character(out$ph.rtseq),]
# [1] uw.seqid uw.ptid  uw.byear uw.sex   uw.rtseq ph.seqid ph.ptid  ph.byear ph.sex   ph.rtseq
#<0 rows> (or 0-length row.names)
#> out[out$uw.sex != out$ph.sex,]
# [1] uw.seqid uw.ptid  uw.byear uw.sex   uw.rtseq ph.seqid ph.ptid  ph.byear ph.sex   ph.rtseq
#<0 rows> (or 0-length row.names)
#> out[out$uw.byear != out$ph.byear,]
# [1] uw.seqid uw.ptid  uw.byear uw.sex   uw.rtseq ph.seqid ph.ptid  ph.byear ph.sex   ph.rtseq
#<0 rows> (or 0-length row.names)

##save
write.csv(out , file="./DataExport/OverlapPtids_UW_PHSKC.csv",quote=F,row.names=F)

# This is the overlap

# We are also interested in all the unique sequences found in either data set, so that we can
# make a single merged data set without duplicates.

#############################################
#############################################End