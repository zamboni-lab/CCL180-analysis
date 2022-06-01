# Basic Function to fix data format

#list to character
fixlc<-function(obj){as.character(unlist(obj))}
#list to numeric
fixln<-function(obj){as.numeric(as.character(unlist(obj)))} # should encode text as factors firstf
#factor to numeric
fixlf<-function(obj){as.numeric(as.factor(unlist(obj)))} # should encode text as factors first