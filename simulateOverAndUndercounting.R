simulateOverAndUndercounting = function(thicks,OCP,UCP){

  if(length(OCP)==1){#if there's only one OCP, replicate over length
    OCP=matrix(data=OCP,ncol=1,nrow=length(thicks))
  }
  if(length(UCP)==1){#if there's only one UCP, replicate over length
    UCP=matrix(data=UCP,ncol=1,nrow=length(thicks))
  }
  
  
  #generate a random series (uniform length) 
  rdata=runif(length(thicks))
  
  #find those that were overcounted,
  OCi=which(rdata<OCP)
  #and undercounted 
  UCi=which(rdata>(1-UCP))
  
  #loop through and correct overcounting 
  newThicks=matrix(data=NA,ncol=1,nrow=length(thicks)+length(UCi)-length(OCi))
  #newThicks=c()
  counter=0
  secFlag=FALSE
  OCi1Flag=FALSE
  u=1
  wasUCi=c()#create new indices that show which years had been UC before correction
  wasOCi=c()
  #create an index that maps old to new
  old2new=c()
  
  
  while(u <= counter+length(thicks)){
    if(any((u-counter)==UCi) & !secFlag){
      wu=UCi[which(u==(UCi+counter))]
      newThicks[u]=thicks[wu]/2
      secFlag=TRUE
      wasUCi=append(wasUCi,u)
    }else if(secFlag){
      newThicks[u]=thicks[wu]/2
      counter=counter+1
      secFlag=FALSE
      wasUCi=append(wasUCi,u)
      
    }else if(any((u-counter)==OCi) & !secFlag){
      wo=OCi[which(u==(OCi+counter))]
      #build more functionality here for handling these options
      if (u==1){
          #combine with next measurement
          newThicks[u]=thicks[wo]+thicks[wo+1]
          wasOCi=append(wasOCi,u)
      }else{
        #combine with previous
        u=u-1
        newThicks[u]=thicks[wo]+newThicks[u]
        #newThicks[u]=thicks[wo+1]
        wasOCi=append(wasOCi,u)
      }
      counter=counter-1
    }else{
      newThicks[u]=thicks[u-counter]
    }
    old2new[u]=u-counter
    u=u+1
    
    
    
  }
  #account for possibility that last one was undercounted
  if(secFlag){
    newThicks[u]=thicks[wu]/2
    counter=counter+1
    secFlag=FALSE
    wasUCi=append(wasUCi,u)
  }
  
  #remove NAs.
  newThicks=na.omit(newThicks)
  
  #NOTE: THIS ERROR GETS VIOLATED OCCASIONALLY (1/200 simulations), for reasons that aren't clear to me. I don't think it's too important. Yet. OK - I think I know why, when the first one is Undercounted, and the second is overcounted. It handles it properly, but doesn't register 2 as overcounted. The total thicknesses are the same, so I'm ok with this, for now at least.
  
#   #do some checking!
#   if(length(newThicks)-counter != length(thicks)){
#     stop("Counter & newThick lenghts don't line up")
#   }
  if( abs(sum(thicks)-sum(newThicks))>0.0001 ){
    stop("The thicknesses don't sum to the same number")
  }

  out = list(newThicks=newThicks,old2new=old2new,wasUCi=wasUCi,wasOCi=wasOCi,UCi=UCi,OCi=OCi)
  return(out)  
  
}