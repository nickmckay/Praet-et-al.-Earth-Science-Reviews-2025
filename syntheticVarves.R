syntheticVarves = function(n,AR,shape=1.5,mean = 1){
  


#function to create a synthetic varve/sedimentation rate series
#AR=.25
#shape=1.5
#n=10000
gN=n*10
nseq=seq(1,n)  

library(RandomFields)
#generate pool of random data, length = n*10
gPool=rgamma(gN,shape=shape,rate = shape/mean)
gPs=sort(gPool)  

#create autocorrelated fractal brownian motion series
  model_fBm <- RMfbm(alpha=AR, var = 1)
  fBm.sim = RFsimulate(model_fBm, nseq)
  fBm = RFspDataFrame2conventional(fBm.sim)$data
  
  ncdf=pnorm(scale(fBm))
  index=ceiling(ncdf*length(gPool))
  
  varves=gPs[index]
  #plot(varveInit,type="l")
  return(varves)
}



syntheticVarvesHurst = function(n,H,shape=1.5,mean = 1){
  
 
  #function to create a synthetic varve/sedimentation rate series
  #AR=.25
  #shape=1.5
  #n=10000
  gN=n*10
  #nseq=seq_len(n)  
  
  #generate pool of random data, length = n*10
  gPool=rgamma(gN,shape=shape,rate = shape/mean)
  gPs=sort(gPool)  
  
  #create autocorrelated fractal brownian motion series
  fBm = FGN::SimulateFGN(n,H)

  ncdf=pnorm(scale(fBm))
  index=ceiling(ncdf*length(gPool))
  varves=gPs[index]
  #plot(varveInit,type="l")
  return(varves)
}



#second function, that simulates number of years given a depth
simulateVarveCountToThickness = function(thick,H,shape=1.5,mean=1){
  
vs = syntheticVarvesHurst(n = max(2,ceiling(thick/mean*10)),H=H,shape=shape,mean = mean)
nVarves  = min(which(cumsum(vs)>thick))
return(nVarves)
}

simpleVarveModel = function(signal,H,shape=1.5,mean = 1,SNR=.25){
  #simple PSM for varves
  #signal = the climate forcing variable, at annual timescales, (e.g., a vector of temperature)
  #H = Hurst coefficient
  #shape=of the gamma distribution
  #mean of the varve series output
  
  #first, gammify the input
  gamSig = gammify(signal,shape=shape,mean = mean)
  

  #create gammafied autocorrelated fractal brownian motion series
  gamNoise = gammify(FGN::SimulateFGN(length(signal),H),shape = shape, mean = mean)
  
  #combine the signal with the noise, based on the SNR
  varves = (gamSig*SNR + gamNoise*(1/SNR))/(SNR+1/SNR)
  return(varves)
}

