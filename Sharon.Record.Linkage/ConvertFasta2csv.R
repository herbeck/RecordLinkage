###################################################
###
### PHSKC has fasta file, this program is to convert PR/RT seq to csv file
##################################################


cur.dir <-"I:/groups/SOM/ICRC/DataTeam/Genomics/SeattleHIVSources_Josh/"
setwd(cur.dir)


####################################
input <- readLines("./DataOriginal/2019_02_20 - for Josh/Sequences.fas")
output <- file("./DataIntermediate/PHSKC_PR_RTSeq.csv","w")

currentSeq <- 0
newLine <- 0
delimiter = ","
writeLines(paste("SeqID",delimiter,"PrtRT.seq","\n"), output, sep="")

for(i in 1:length(input)) {
  ##only process PR/RT seq, ignore INT seq
  if(strtrim(input[i], 1) == ">" & strtrim(input[i], 6) == ">PR/RT") {
  
    seqid = substring(input[i],2)
    if(currentSeq == 0) {
      writeLines(paste(seqid,delimiter), output, sep="")
      currentSeq <- currentSeq + 1
    } else {
      writeLines(paste("\n",seqid,delimiter, sep =""), output, sep="")
    }
  } else if (strtrim(input[i], 1) != ">") {
    writeLines(paste(input[i]), output, sep="")  ##write seq
  }
}

close(output)






#############################################
##############################################End