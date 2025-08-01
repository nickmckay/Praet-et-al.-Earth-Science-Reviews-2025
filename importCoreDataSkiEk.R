#import files and build core list for varveModeling
getCoreData <- function(inDir="~/Dropbox/VarveModel/SkilakEklutna/Input/"){
setwd(inDir)
#load thicknesses
thickIn = read.csv("thicknesses.csv") #csv file with years in first column, and different cores in each subsequent column



#load OC and UC priors
countingPriors = read.csv("countPriors.csv")


#load marker layer information
markLayers =read.csv("markerLayers.csv")

#need to create a list, that has entries for each core (numbered), and contains six variables:
#years (BP), vector
#thicknesses , vector
#markerLayers, sparse vector with marker Layer names at marker layers
#UCprior
#OCprior
#MLprior,  sparse vector with marker Layer probabilities at marker layers

#how many cores?
nCores = ncol(thickIn)-1
coreNames=names(thickIn[,2:ncol(thickIn)])


cores = vector(mode = "list",length = nCores)
for(c in 1:nCores){
  #grab years and thicknesses
  year = thickIn[,1]
  thicks = thickIn[,c+1]
  yrIn = which(!is.na(thicks) & !is.na(year))
  year = year[yrIn]
  thicks = thicks[yrIn]
  
  #create marker layers
  mls = markLayers[,c+1]
  mlsNames = as.character(markLayers[,1])
  mlsYrin=c()
  for(m in 1:length(mls)){
    if(!is.na(mls[m])){
      mlsYrin[m] = which(year==mls[m])
    }
  }
  goodmls=which(!is.na(mlsYrin))
  mlsYrin = mlsYrin[goodmls]
  mlsNames = mlsNames[goodmls]
  
  #setup empty markerLayers same size as thicknesses
  markerLayers = matrix(nrow=length(thicks))
  markerLayers[mlsYrin]=mlsNames
  
  #repeat with markerLayer prioris
  mlPriors = markLayers[,nCores+2]/100
  if(max(mlPriors)>1 | min(mlPriors)<0){
    stop("problem with your marker layer priors")
  }
  mlPriors = mlPriors[goodmls]
  mlPrior = matrix(nrow=length(thicks))
  mlPrior[mlsYrin]=mlPriors
  
  #ok, now overcounting and undercounting
  ocPrior = countingPriors[,c+1][yrIn]
  ucPrior = countingPriors[,c+1][yrIn]
  
  
  
  coreData = data.frame(year, thickness = thicks,markerLayers,mlPrior,ocPrior, ucPrior)  
  
  cores[[c]] = coreData
  
}


save(cores,file="SkilakEklutnaInput.Rdata")
return(cores)
}
