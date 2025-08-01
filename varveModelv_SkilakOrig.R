#.libPaths("/home/npm4/R/x86_64-redhat-linux-gnu-library/3.4/")
library(here)
library(tidyverse)
library(matrixStats)

#load in list of core dataframes

####DIVIDE BY two TO MAKE IT CLOSER TO REALITY - we artifically inflated earlier.... 
# for(c in 1:length(cores)){
#   cores[[c]]$ocPrior = cores[[c]]$ocPrior/2
#   cores[[c]]$ucPrior = cores[[c]]$ucPrior/2
#   
# }



inputFolders <- list.dirs(here("SkilakEklutna","Input"))[str_detect(list.dirs(here("SkilakEklutna","Input")),pattern = "\\d$")]

#inputFolders <- inputFolders[-c(1:11)]

source(here("importCoreDataSkiEk.R"))
source(here("syntheticVarves.R"))
source(here("simulateOverAndUndercounting.R"))
nSim =1000
nMLit=1000


#FOR NOW, CORES have to have a hard top! Top marker layer must have a probability of 1
#this can be fixed, but needs some work. 
for(inputFolder in inputFolders){
  print(glue::glue("Now working on {inputFolder}"))
  cores <- getCoreData(inputFolder)
  nCores = length(cores)
  
  #need a list of all marker layers in order
  allMarkerLayers = read.csv(file.path(inputFolder,"MarkerLayers.csv"))[,1]
  
  #setup Post of markerLayers
  MLpostAge = matrix(NA,nrow=length(allMarkerLayers),ncol=nSim)
  
  #get some info about the thicknesses
  startYears = sapply(cores,FUN = function(x)min(x$year))
  endYears = sapply(cores,FUN = function(x)max(x$year))
  
  #make a matrix
  overlapMat = sapply(cores,FUN = function(x)x$thickness[which(x$year==max(startYears)):which(x$year==min(endYears))])
  
  #any NAs
  overlapMat = overlapMat[which(!apply(is.na(overlapMat),MARGIN = 1,any)),]
  
  
  meanThick=colMeans(overlapMat)
  sdThick=colSds(overlapMat)
  
  totYears=sapply(cores,FUN = function(x)length(x$thickness))
  
  #Setup composites matrix
  compMat = matrix(NA,nrow = round(2*max(totYears)),ncol = nSim)
  
  #Update priors based on years between marker layers
  MLsim=list()
  #how many iterations to update priors from marker layers?
  
  
  #pull out marker layer order
  # allML = lapply(cores,FUN = function(x)x$markerLayers[which(!is.na(x$markerLayers))])
  # uniqueML = unlist(unique(allML))
  # nML = sapply(allML,length)
  # meanML = c()
  # for(nm in 1:length(uniqueML)){
  #   meanML[nm] = mean(sapply(cores,FUN = function(x)which(x$markerLayers==uniqueML[nm])),na.rm=TRUE)
  # }
  # names(meanML)=uniqueML
  # meanML = meanML[order(meanML)]
  
  #Setup storage for new prior information
  for(c in 1:nCores){
    cores[[c]]$bigOCA=matrix(data = NA,nrow = length(cores[[c]]$ocPrior),ncol=nMLit)
    cores[[c]]$bigUCA=matrix(data = NA,nrow = length(cores[[c]]$ucPrior),ncol = nMLit)
  }
  
  print("Simulate Marker Laters")
  #runs simulation
  for(it in 1:nMLit){
    if(it%%50 == 0){
      cat(round(it/nMLit * 100),"% \r")
    }
    for(c in 1:nCores){
      isML=which(!is.na(cores[[c]]$markerLayers))
      realTest=runif(length(isML))<cores[[c]]$mlPrior[isML]
      MLsim[[c]]=cores[[c]]$markerLayers
      MLsim[[c]][isML[which(!realTest)]]=NA
    }
    
    #identify which tie points are found in all cores in this simulation
    allML = lapply(cores,FUN = function(x)x$markerLayers[which(!is.na(x$markerLayers))])
    lastML = sapply(cores,FUN = function(x)x$markerLayers[max(which(!is.na(x$markerLayers)))])
    firstML = sapply(cores,FUN = function(x)x$markerLayers[min(which(!is.na(x$markerLayers)))])
    
    #create matrix of where we have marker layers
    simML = matrix(FALSE,ncol = nCores,nrow= length(allMarkerLayers))
    for(a in 1:length(allMarkerLayers)){
      for(c in 1:nCores){
        if(length(which(allML[[c]]==allMarkerLayers[a])) > 0 ){
          simML[a,c]=TRUE
        }
      }
    }
    #now fill in FALSEs that are surrounded by TRUES
    for(c in 1:nCores){
      #which are false
      isF=which(!simML[,c])
      for(f in isF){
        if(f<max(which(simML[,c])) &  f>min(which(simML[,c]))){
          simML[f,c]=TRUE
        }
      }
    }
    nMLsim = cbind(allMarkerLayers,rowSums(simML))
    
    
    TPlist = lapply(MLsim,unique)
    allTP =unique(unlist(TPlist))
    allTP=allTP[which(!is.na(allTP))]
    inAll=c()
    inTwo=c()
    for (t in 1:length(allTP)){
      inAll[t]=all(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))
      inTwo[t]=sum(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))==nMLsim[which(allTP[t]==nMLsim[,1]),2]
      
    }
    
    #which Marker Layers are we simulating here
    ML2sim=as.vector(allTP[which(inTwo)])
    ML2sim=na.exclude(ML2sim)
    MLyears = matrix(NA,nrow = length(ML2sim),ncol = nCores)
    for(m in 1:nrow(MLyears)){
      allMLList = (sapply(cores,FUN = function(x)which(x$markerLayers==ML2sim[m])))
      for(c in 1:nCores){
        thisML = allMLList[[c]]
        #remove MLs that didn't pass 
        if(length(thisML==1)){
          MLyears[m,c] = thisML
        }
      }
    }
    
    #make sure that the MLyears is in the correct order!!!
    ody = apply(MLyears,2,diff)
    if(any(ody[is.finite(ody)]<0)){
      stop("There are reversals in your marker layers. Check your input files.")
    }
    #print(MLyears)
    
    
    #WARNING: AS CURRENTLY CONSTRUCTED, THIS ONLY MAKES SENSE FOR SYMMETRIC OC UC PRIORS!
    #calculate diffs
    dML=diff(MLyears)
    
    #ASSUME A NORMAL DISTRIBUTION FOR THE DURATION of the YEARS and DRAW FROM IT
    #THIS MAY NOT BE IDEAL!!!!!!!!!
    diffMeans = rowMeans(diff(MLyears),na.rm=TRUE)
    diffSds = rowSds(diff(MLyears),na.rm=TRUE)
    simTrueDiff = rnorm(n = length(dML),mean = diffMeans,sd = diffSds)
    OCUC = dML/simTrueDiff
    #rowSds(diff(MLyears)) 
    OCadd = OCUC-1
    OCadd[OCadd<=0]=0
    UCadd = -1*(OCUC-1)
    UCadd[UCadd<=0]=0
    
    
    for (c in 1:nCores){
      for (i in 1:(nrow(MLyears)-1)){
        mlstart=which(ML2sim[i]==cores[[c]]$markerLayers )
        mlstop=which(ML2sim[i+1]==cores[[c]]$markerLayers )
        if(length(mlstart) ==1 & length(mlstop)==1){
          #find index over which to update priors
          pri=which(ML2sim[i]==cores[[c]]$markerLayers ):which(ML2sim[i+1]==cores[[c]]$markerLayers )
          cores[[c]]$bigOCA[pri,it]=OCadd[i,c]
          cores[[c]]$bigUCA[pri,it]=UCadd[i,c]
        }
      }
    }
    
  }
  
  
  #use simulations on marker layers to get distributions on yearly priors.
  for (c in 1:nCores){
    #this way makes single estimates based on the mean. A distribuiton is better.
    # cores[[c]]$ocPrior=rowSums(cbind(cores[[c]]$ocPrior,rowMeans(cores[[c]]$bigOCA,na.rm=TRUE)),na.rm=TRUE) 
    # cores[[c]]$ucPrior=rowSums(cbind(cores[[c]]$ucPrior,rowMeans(cores[[c]]$bigUCA,na.rm=TRUE)),na.rm=TRUE)
    
    #create a matrix of priors
    cores[[c]]$ocPriorMat = matrix(rep(cores[[c]]$ocPrior,times = nMLit),ncol=nMLit)+cores[[c]]$bigOCA
    cores[[c]]$ucPriorMat = matrix(rep(cores[[c]]$ucPrior,times = nMLit),ncol=nMLit)+cores[[c]]$bigUCA
    
    
    
  }
  
  
  
  #setup posterior matrix storage.
  corePost = vector("list",length = nCores)
  for(c in 1:nCores){
    corePost[[c]]$mlPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
    corePost[[c]]$ucPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
    corePost[[c]]$ocPost=matrix(0,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
    corePost[[c]]$meanCorr=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
    corePost[[c]]$thicks=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
    corePost[[c]]$thick2=c()
    corePost[[c]]$mlCum=matrix(NA,nrow = round(length(cores[[c]]$year)*1.5),ncol=nSim) 
  }
  
  

# Simulate counts ---------------------------------------------------------

  
  #Simulate segments, update priors and create composites.
  print("Simulate counts")
  
  for(i in 1:nSim){
    #empty out this year and last year
    lastMLyear = 0
    thisMLyear = 0
    
    MLsim=list()
    #simulate whether or not tie points are "TRUE"
    for(c in 1:nCores){
      isML=which(!is.na(cores[[c]]$markerLayers))
      realTest=runif(length(isML))<cores[[c]]$mlPrior[isML]
      MLsim[[c]]=cores[[c]]$markerLayers
      MLsim[[c]][isML[which(!realTest)]]=NA
    }
    
    #identify which tie points are found in at least 2?
    TPlist = lapply(MLsim,unique)
    allTP =unique(unlist(TPlist))
    allTP=allTP[which(!is.na(allTP))]
    inAll=c()
    inTwo=c()
    for (t in 1:length(allTP)){
      inAll[t]=all(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))
      inTwo[t]=sum(sapply(TPlist,FUN = function(x)any(is.element(x,allTP[t]))))==nMLsim[which(allTP[t]==nMLsim[,1]),2]
      
    }
    #which Marker Layers are we simulating here
    ML2sim=as.vector(allTP[which(inTwo)])
    ML2sim=na.exclude(ML2sim)
    MLyears = matrix(NA,nrow = length(ML2sim),ncol = nCores)
    for(m in 1:nrow(MLyears)){
      allMLList = (sapply(cores,FUN = function(x)which(x$markerLayers==ML2sim[m])))
      for(c in 1:nCores){
        thisML = allMLList[[c]]
        #remove MLs that didn't pass 
        if(length(thisML==1)){
          MLyears[m,c] = thisML
        }
      }
    }
    #Do not reorder! THIS MEANS THAT Marker Layers MUST be in stratigraphic order
    # newOrder = order(rowMeans(MLyears,na.rm=TRUE))
    # MLyears = MLyears[newOrder,]
    # ML2sim = ML2sim[newOrder]
    
    #make sure that there are no reversals!!!
    MLyears
    ody = apply(MLyears,2,diff)
    
    #store Marker Layer posterior data
    for(c in 1:nCores){
      for(ml in 1:length(ML2sim)){
        wmlp=which(ML2sim[ml]==cores[[c]]$markerLayers)
        corePost[[c]]$mlPost[wmlp,i]=ML2sim[ml]
      }
    }
    
    #create a vector that gives instructions about matching seg lengths
    mustMatchLengths = matrix(TRUE,nrow=(nrow(MLyears)-1))
    
    #add to ML years if need be.
    #does the first segment start at 1? (is the surface a marker layer?)
    ML1test = MLyears[1,!is.na(MLyears[1,])]==1
    if (!all(ML1test) & any(ML1test)){
      stop("You cant have a marker layer that starts at 1 for some cores and not all of them")
    }
    if (!all(ML1test) ){
      #if no top seg, add one
      MLyears=rbind(matrix(1,ncol=nCores),MLyears)
      mustMatchLengths=rbind(TRUE,mustMatchLengths)#ALERT!!!!! THIS MEANS THAT WE'VE GIVEN A prob of 1 to the top layer!
    }
    
    lastMLtest = MLyears[nrow(MLyears),!is.na(MLyears[nrow(MLyears),])]
    totYearTest = totYears[!is.na(MLyears[nrow(MLyears),])]
    if(any((totYearTest - lastMLtest)==0) & !all((totYearTest - lastMLtest)==0)){
      stop("Your last marker layer must either be the oldest years, or at least 2 years a way from your oldest years")
    }
    if(!any((totYearTest - lastMLtest)==0)){
      #add in the bottom MLyears layer
      MLyears=rbind(MLyears,totYears)
      ML2sim=append(ML2sim,"bottom")
      
      mustMatchLengths=rbind(mustMatchLengths,FALSE)
      
    }
    
    #set up segments
    seg=vector("list",nrow(MLyears)-1)
    
    
    for(s in 1:(nrow(MLyears)-1)){
      for(c in 1:nCores){
        if(s==(nrow(MLyears)-1) | !is.na(MLyears[s,c]) & !is.na(MLyears[s+1,c])){
          if(s==1){
            seg[[s]][[c]]=cores[[c]][(MLyears[s,c]):MLyears[s+1,c],]#include the first marker layer the first year
            
          }else if(s< (nrow(MLyears)-1)){
            seg[[s]][[c]]=cores[[c]][(MLyears[s,c]+1):MLyears[s+1,c],]
          }else{#IF ITS THE LAST ROW, make the segment from the last, non-NA marker layer
            lastRealMarker = max(MLyears[-nrow(MLyears),c],na.rm=T)
            seg[[s]][[c]]=cores[[c]][(lastRealMarker+1):MLyears[s+1,c],]
          }
        }
      }
    }
    segTotalThicks = matrix(NA,nrow = 5,ncol=length(seg))
    segTotalYears = matrix(NA,nrow = 5,ncol=length(seg))
    
    ###CHECK THAT SEGMENT THICKNESS ARE CORRECT
    for(c in 1:nCores){
      for(s in 1:(nrow(MLyears)-1)){
        if(length(seg[[s]])>=c){
          segTotalThicks[c,s] = sum(seg[[s]][[c]]$thickness)
          segTotalYears[c,s] = length(seg[[s]][[c]]$year)
        }
      }
    }
    
    
    #   if(!all(sapply(seg[[1]],FUN = function(x)min(x$year)) == 1)){
    #     #make a new marker layer that goes to the top, that doesn't have the same length requirement.
    #     topSeg[[c]] = cores[[c]][MLyears[s,c]:MLyears[s+1,c],]
    #   }
    #segment  
    
    newSeg=vector("list",nrow(MLyears)-1)
    #simulate between segments 
    tSim=vector("list",nCores)
    
    composite = c()
    for(s in 1:length(seg)){#for each segment...
      #which cores are we simulating?
      core2sim=which(!sapply(seg[[s]],is.null))
      orig.core2sim = core2sim
      
      sameLength = FALSE
      t=0
      while(!sameLength){#
        t=t+1
        #prior picker - pick a number to randomly sample from the OC/UC priors
        ppi = sample.int(nMLit,size = 1)
        
        #build in failsafe
        if(t>1e4){
          if(length(core2sim>2)){
           # print(paste("segment",as.character(s),"hit count limit.", as.character(length(core2sim)),"cores. Reducing number of cores in simulation"))
            core2sim = core2sim[-sample(core2sim)[1]]#remove one of the cores.
            t=0
          }else{        
            stop("Too many iterations..no more cores to remove.. stopped.")
          }
        }
        for(c in core2sim){
          #simulate new segments given the UC and OC priors.
          tSim[[c]]= simulateOverAndUndercounting(thicks = seg[[s]][[c]]$thickness,OCP = seg[[s]][[c]]$ocPriorMat[,ppi],UCP = seg[[s]][[c]]$ucPriorMat[,ppi])
        }
        
        newSegLengths = sapply(tSim,FUN = function(x)length(x$newThicks))[core2sim]
        sameLength= (all(newSegLengths==newSegLengths[1]) | !mustMatchLengths[s])
        
        
      }
      
      #see if any didn't run
      didnt.sim = orig.core2sim[!orig.core2sim %in% core2sim]
      if(length(didnt.sim)>0){#if any didn;t run...
        for(dd in didnt.sim){#one by one, simulate the length for any that got excluded
          
          ddt=0
          while(TRUE){#and simulate until you match the length
            ppi = sample.int(nMLit,size = 1)
            ddt=ddt+1
            tSim[[dd]]= simulateOverAndUndercounting(thicks = seg[[s]][[dd]]$thickness,OCP = seg[[s]][[dd]]$ocPriorMat[,ppi],UCP = seg[[s]][[dd]]$ucPriorMat[,ppi])
            #check to see if it's the correct length
            thisSegLength = length(tSim[[dd]]$newThicks)
            if(thisSegLength == newSegLengths[1]){#found it!!!
              core2sim = append(core2sim,dd)
              newSegLengths = append(newSegLengths,thisSegLength)
              break
            }
            if(ddt>1e4){#give up
              #print(paste("gave up trying to get core",as.character(dd),"to be",as.character(newSegLengths[1]),"years long for segment",as.character(s)))
              break
            }
          }
        }
      }
      
      
      
      
      
      #Once same length is matched build new data.frame
      #And assign results into posterior.
      segthickmat=matrix(NA,nrow=max(newSegLengths),ncol=length(core2sim))
      for(c in 1:length(core2sim)){
        #build Thickmat
        segthickmat[1:newSegLengths[c],c] = tSim[[core2sim[c]]]$newThicks
      }
      cm=cor(segthickmat,use="complete.obs")
      cm[which(cm<0)]=0
      meanCorr=mean(cm[lower.tri(cm)])
      
      ###ASSIGN IN POSTERIOR INFORMATION###
      for(c in core2sim){
        #assign posterior info
        #Undercounted
        ucIndex = tSim[[c]]$UCi+min(seg[[s]][[c]]$year)-1
        corePost[[c]]$ucPost[ucIndex,i]=1
        #Overcounted
        ocIndex = tSim[[c]]$OCi+min(seg[[s]][[c]]$year)-1
        corePost[[c]]$ocPost[ocIndex,i]=1
        
        #mean correlation
        corePost[[c]]$meanCorr[min(seg[[s]][[c]]$year) : max(seg[[s]][[c]]$year) ,i] = meanCorr
        
        #store thicknesses for each core
        if(mustMatchLengths[s]){
          if(c==core2sim[1]){#only update if its the first one.
            thisMLyear = lastMLyear+newSegLengths[1]
          }
          corePost[[c]]$mlCum[thisMLyear,i]=ML2sim[s+1]
          #test if sum thicknesses are the same as the original (THEY SHOULD BE!)
          if(abs(sum(tSim[[c]]$newThicks)-sum(seg[[s]][[c]]$thickness))>.001){
            stop("The thicknesses don't sum to the same number!")
          }
          corePost[[c]]$thick2=append(corePost[[c]]$thick2,tSim[[c]]$newThicks)
          corePost[[c]]$thicks[(lastMLyear+1) : thisMLyear ,i] = tSim[[c]]$newThicks
          
          
        }else{
          lastSimYear = max(which(!is.na(corePost[[c]]$thicks[,i])))
          corePost[[c]]$thicks[(lastSimYear+1) : (lastSimYear+newSegLengths[which(c==core2sim)]) ,i] = tSim[[c]]$newThicks
        }
      }
      
      
      #is the total thickness the same as the total thickness should be?
      curThick=sum(corePost[[c]]$thicks[1 : sum(!is.na(corePost[[c]]$thicks[,i])) ,i] )
      lowerMLname= seg[[s]][[c]]$markerLayers[length(seg[[s]][[c]]$markerLayers)]
      
      #assign age into MLpostAge
      if(!is.na(lowerMLname)){
        wML = which(lowerMLname==allMarkerLayers)
        MLpostAge[wML,i] = thisMLyear
      }
      
      origMLi = which(cores[[c]]$markerLayers==lowerMLname)
      if(length(origMLi)<1){
        origMLi = length(cores[[c]]$thickness)
      }
      origThick = sum(cores[[c]]$thickness[1 : origMLi])
      # if(abs(origThick-curThick)> .001){
      #   stop("The current cumulative thickness isn't right")
      # }
      
      
      lastMLyear=thisMLyear#update with new last year
      
      #build Composite
      segCom= rowMeans(scale(segthickmat,center=meanThick[core2sim],scale =sdThick[core2sim]),na.rm = TRUE)
      composite = append(composite,segCom)
    }
    compMat[1:length(composite),i]=composite
    if(i%%10 == 0){
      cat(round(i/nSim * 100),"% \r")
    }
    
  }
  save.image(here(inputFolder,paste0("SkiEkOut.Rdata")))
}
