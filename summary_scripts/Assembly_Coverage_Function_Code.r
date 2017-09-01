## R CODE FOR ANALYSIS OF RESULTING ASSEMBLY COVERAGE FILES ##

library(hexbin)
library(data.table)
library(doBy)
library(lattice)
library(gridExtra)
#load appropriate libraries. 

temp_data <- function(file) {
  temp_dataset <- read.table(file, header=F, sep=" ")
  #Imports the data
  temp_dataset <- temp_dataset[ which(temp_dataset$V1!='genome'), ]
  #remove anything in the file where 'genome' is in column V1
  temp_dataset$loc = as.character(lapply(strsplit(as.character(temp_dataset$V1), split="_"), "[", 1))
  #split each line using '_' as the split variable. assign the second part of the split to loc (locus)
  temp_dataset$A <- basename(file)
  temp_dataset$taxon = as.character(lapply(strsplit(as.character(temp_dataset$A), split="-"), "[", 1)) 
  #split each line using '_' as the split variable. assign the second part of the split to taxon
  return(temp_dataset)
  #return the dataset
}
#The above is a function that opens each file sent to it and cleans up the data a bit.

make_dataset <- function(working_dir) {
  #a function to generate the dataset. 
  setwd(working_dir)
  files <- list.files(pattern = '.per.contig.summary')
  #set working directory and get files.

  for (file in files){
  #simple for loop traversing all files the loop picked up earlier. 
    if (exists("dataset")){
      #if the dataset already exists, this will start the append process. 
      Sys.sleep(0.1)
      print(file)
      flush.console() 
      #provides a way to see what files are being worked on.
      
      temp_dataset <- temp_data(file)
      #run the function setup above on the file.
   
      dataset<-rbind(dataset, temp_dataset)
      #join the new data to the old dataset. 
      
      rm(temp_dataset)
      #as a precaution, remove the temp_dataset from memory. 
  
    } else {
      #if dataset does not exist, we need to make it. 
      Sys.sleep(0.1)
      print(file)
      flush.console()
      #provides a way to see what files are being worked on.
      
      dataset <- temp_data(file)
      #this runs the function on the file and generates 'dataset' 
    }
  }
new_data <- dataset
rm(dataset)
#this is a precaution against if you run the function twice. It will remove 'dataset' so that if the loop is run again, a new dataset will be made. 
return(new_data)
}
#The above is a function loops over files in the directory provided to the function.

new_data <- make_dataset("/directory/path/to/your/files")
#make the dataset using the functions 'make_dataset' and 'temp_data'. 
#to get this statement to work properly, you will need to execute the code for the functions into active memory (similar to loading a pacakge)

#########
# now that you have the data, the statements below, summarize that data in various ways. 

summary_avg <- summaryBy(V2 ~ taxon, data = new_data, FUN = function(x) { c(n = length(x), m = mean(x), sd = sd(x), min=min(x), max=max(x)) } )

