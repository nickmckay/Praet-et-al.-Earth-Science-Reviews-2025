library(MASS)
library(fitdistrplus)
library(tidyverse)
library(here)
library(dplyr)
library(readxl)
library(here)
library(lipdR)
library(geoChronR)
library(Hmisc)


#eklutna priors
EkCorePost <- readRDS("EklutnaPostData.RDS")
ekEns <- EkCorePost[[1]]$thicks
ekDepthEns <- apply(ekEns,2,\(x) cumsum(x) / 10 + 60.19)
skiDataPath <- "Review ages Skilak and Eklutna.xlsx"

ekMls <- read_excel(skiDataPath,sheet = "Eklutna counts and events",skip = 1,range = "A2:F35") |> 
  select(Depth, ML = `CORRELATED EVENT NAME`) |> 
  filter(!is.na(ML), Depth >= 60)

for(m in 1:nrow(ekMls)){
  thisOne <- data.frame(age = unlist(apply(ekDepthEns,2,\(x) which.min(abs(x - ekMls$Depth[m])))) + 21,
                        ML = ekMls$ML[m])
  if(m == 1){
    ekPrior <- thisOne
  }else{
    ekPrior <- bind_rows(ekPrior,thisOne)
  }
}




skiDataPath <- reviewAgesFile <- "Review ages Skilak and Eklutna.xlsx"
getCorePost <- function(file) {
  #' Function for extracting an object from a .RData file created by R's save() command
  #' Inputs: RData file, object name
  load(file=file)
  corePost <<- corePost
  rm(list = ls())
}

#stop("a")
inputFolders <- list.dirs(here("SkilakEklutna","Input"))[str_detect(list.dirs(here("SkilakEklutna","Input")),pattern = "\\d$")]
allMlAges <- data.frame(age = NA,name = NA,counter = NA)


hist.data = data.frame(depth =  c(53.6,180.6,494.7),
                       mid.age = c(38,162,771),
                       min.age = c(36,160,761),
                       max.age = c(40,164,781))

hist.data.EK = data.frame(depth =  c(68),
                       mid.age = c(34),
                       min.age = c(33),
                       max.age = c(35))


passHist <- function(depth,age,hist.data){
  if(all(is.na(depth))){
    return(data.frame(GBT2 = FALSE,GBT6 = FALSE, GT6SD2 = FALSE))
  }
  estage <- approx(depth,age,hist.data$depth)$y |> 
  between(left = hist.data$min.age,right = hist.data$max.age) |> 
 t() |> 
    as.data.frame() |> 
    setNames(c("GBT2","GBT6","GT6SD2"))
  
  return(estage)
}


distHist <- function(depth,age,hist.data){
  if(all(is.na(depth))){
    return(data.frame(GBT2 = NA,GBT6 = NA, GT6SD2 = NA))
  }
  estage <- approx(depth,age,hist.data$depth)$y - hist.data$mid.age |> 
    t() |> 
    as.data.frame() |> 
    setNames(c("GBT2","GBT6","GT6SD2"))
  
  return(estage)
}

distHistEK <- function(depth,age,hist.data){
  if(all(is.na(depth))){
    return(data.frame(E1 = NA))
  }
  estage <- data.frame(E1 = approx(depth,age,hist.data$depth)$y - hist.data$mid.age)
  
  return(estage)
}


for(k in 1:length(inputFolders)){
  thisDir <- file.path("SkilakEklutna/Input",k)
  getCorePost(file = file.path(thisDir,"SkiEkOut.Rdata"))

  depthAt21 <- readxl::read_excel(skiDataPath,sheet = "Skilak counts",skip = 1,.name_repair = "none")[c(1,2) + (k-1)*3] |>
    setNames(c("depth","age")) |> 
    filter(age == 21) |> 
    dplyr::select(depth) |> 
    unlist()
    
  #skilak
  
  cumdepth <- apply(corePost[[2]]$thicks,2,cumsum)/10 + depthAt21
  years <- seq_len(nrow(cumdepth)) + 21
  
  #check historical
  t1 <- array_tree(cumdepth,2) |> 
      map(distHist,age = years, hist.data = hist.data) |> 
    list_rbind()
  
  t1$ensMember <- seq_len(nrow(t1))
  t1$counter <- k
  
  
  #eklutna
  ekCumDepth <-  apply(corePost[[1]]$thicks,2,cumsum)/10 + 60
  t2 <- array_tree(ekCumDepth,2) |> 
    map(distHistEK,age = years, hist.data = hist.data.EK) |> 
    list_rbind()
  
  t1$E1 <- t2$E1
  
  if(k == 1){
    allResults <- t1
  }else{
    allResults <- bind_rows(allResults, t1)
  }
  
}


#get original count ensemble members

  og <- readxl::read_excel(skiDataPath,sheet = "Skilak counts",skip = 1) |> 
    dplyr::select(-starts_with("Events")) |> 
    rename_with(
      .cols = 1:46, 
      .fn = ~ paste0(
        rep(c("Depth_", "Varve_"), times = 23),
        rep(1:23, each = 2)
      )
    ) |> 
    pivot_longer(
      cols = 1:46,
      names_to = c(".value", "counter"),
      names_pattern = "(\\w+)_(\\d+)"
    ) %>%
    
    # Step 3: Rename the final columns to match your desired output
    rename(age = Varve, depth = Depth)

#find the "best" ensemble members

hist(allResults$E1)

strictCut <- allResults |> 
  filter(abs(GBT2) <= 1, 
       abs(GBT6) <= 1, 
       abs(GT6SD2) <= 10)


probCut <-  allResults |> 
  mutate(prob = dnorm(abs(GBT2),sd = 1,log = TRUE) + 
           dnorm(abs(GBT6),sd = 1,log = TRUE) + 
           dnorm(abs(E1),sd = 1,log = TRUE) + 
           dnorm(abs(GT6SD2),sd = 10,log = TRUE)) 


cut1000 <- sort(probCut$prob,decreasing = TRUE)[1000]

probCut <- filter(probCut,prob >= cut1000)
  

allProb <- allResults |> 
  mutate(prob = dnorm(abs(GBT2),sd = 1,log = TRUE) + 
           dnorm(abs(GBT6),sd = 1,log = TRUE) + 
           dnorm(abs(E1),sd = 1,log = TRUE) + 
           dnorm(abs(GT6SD2),sd = 10,log = TRUE)) 




inputFolders <- list.dirs(here("SkilakEklutna","Input"))[str_detect(list.dirs(here("SkilakEklutna","Input")),pattern = "\\d$")]
allMlAges <- data.frame(age = NA,name = NA,counter = NA)
for(k in 1:length(inputFolders)){
  thisDir <- file.path("SkilakEklutna/Input",k)
  load(file.path(thisDir,"SkiEkOut.Rdata"))
  thisMlPrior <- read_csv(file.path(thisDir,"MarkerLayers.csv"))
  thisMlPrior$counter <- basename(thisDir)
  if(k == 1){
    allMlPrior <- thisMlPrior
  }else{
    allMlPrior <- bind_rows(allMlPrior,thisMlPrior)
  }
  mlPriors <- read_csv(file.path(thisDir,"MarkerLayers.csv"))
  
  #only grab "good ensemble members
  theseEns <- filter(probCut,counter == k)
  if(nrow(theseEns) > 0){
    MLpostAge <- as.matrix(MLpostAge[,c(theseEns$ensMember)])
    
    for(mli in 1:length(allMarkerLayers)){
      mlName <- allMarkerLayers[mli]
      theseMlAges <- data.frame(age = as.numeric(MLpostAge[mli,]) + 21, name = mlName, counter = basename(thisDir))
      allMlAges <- bind_rows(allMlAges,theseMlAges)
    }
    
  }
}


#build up grand-weighted ensemble
gc <- unique(probCut$counter)


probCut$thickness <- NA
probCut$depth <- NA
probCut$age <- NA
probCut$thicknessEk <- NA
probCut$depthEk <- NA
probCut$ageEk <- NA
for(k in gc){
  thisDir <- file.path("SkilakEklutna/Input",k)
  getCorePost(file = file.path(thisDir,"SkiEkOut.Rdata"))
  
  depthAt21 <- readxl::read_excel(skiDataPath,sheet = "Skilak counts",skip = 1,.name_repair = "none")[c(1,2) + (k-1)*3] |>
    setNames(c("depth","age")) |> 
    filter(age == 21) |> 
    dplyr::select(depth) |> 
    unlist()
  
  #skilak
  
  cumdepth <- apply(corePost[[2]]$thicks,2,cumsum)/10 + depthAt21
  years <- seq_len(nrow(cumdepth)) + 21
  
  goodCounts <- filter(probCut,counter == k)
  thicks <- depths <- age <- vector(mode = "list",length = nrow(goodCounts))
  for(j in 1:nrow(goodCounts)){
    thicks[[j]] <- corePost[[2]]$thicks[,goodCounts$ensMember[j]]
    depths[[j]] <- cumdepth[,goodCounts$ensMember[j]]
    age[[j]] <- years
  }
  
  tri <- which(probCut$counter == k)
  
  probCut$thickness[tri] <- thicks
  probCut$depth[tri] <- depths
  probCut$age[tri] <- age
  
  #eklutna
  cumdepthEk <- apply(corePost[[1]]$thicks,2,cumsum)/10 + 60
  yearsEk <- seq_len(nrow(cumdepth)) + 21
  
  thicksEk <- depthsEk <- ageEk <- vector(mode = "list",length = nrow(goodCounts))
  for(j in 1:nrow(goodCounts)){
    thicksEk[[j]] <- corePost[[1]]$thicks[,goodCounts$ensMember[j]]
    depthsEk[[j]] <- cumdepthEk[,goodCounts$ensMember[j]]
    ageEk[[j]] <- yearsEk
  }
  
  tri <- which(probCut$counter == k)
  
  probCut$thicknessEk[tri] <- thicksEk
  probCut$depthEk[tri] <- depthsEk
  probCut$ageEk[tri] <- ageEk

}


#get all ensemble members for comparison

allResults$thickness <- NA
allResults$depth <- NA
allResults$age <- NA
allResults$thicknessEk <- NA
allResults$depthEk <- NA
allResults$ageEk <- NA
for(k in unique(allResults$counter)){
  thisDir <- file.path("SkilakEklutna/Input",k)
  getCorePost(file = file.path(thisDir,"SkiEkOut.Rdata"))
  
  depthAt21 <- readxl::read_excel(skiDataPath,sheet = "Skilak counts",skip = 1,.name_repair = "none")[c(1,2) + (k-1)*3] |>
    setNames(c("depth","age")) |> 
    filter(age == 21) |> 
    dplyr::select(depth) |> 
    unlist()
  
  #skilak
  
  cumdepth <- apply(corePost[[2]]$thicks,2,cumsum)/10 + depthAt21
  years <- seq_len(nrow(cumdepth)) + 21
  
  theseCounts <- filter(allResults,counter == k)
  thicks <- depths <- age <- vector(mode = "list",length = nrow(theseCounts))
  for(j in 1:nrow(theseCounts)){
    thicks[[j]] <- corePost[[2]]$thicks[,theseCounts$ensMember[j]]
    depths[[j]] <- cumdepth[,theseCounts$ensMember[j]]
    age[[j]] <- years
  }
  
  tri <- which(allResults$counter == k)
  
  allResults$thickness[tri] <- thicks
  allResults$depth[tri] <- depths
  allResults$age[tri] <- age
  
  
  #eklutna
  
  cumdepthEk <- apply(corePost[[1]]$thicks,2,cumsum)/10 + 60
  yearsEk <- seq_len(nrow(cumdepth)) + 21
  
  thicksEk <- depthsEk <- ageEk <- vector(mode = "list",length = nrow(theseCounts))
  for(j in 1:nrow(theseCounts)){
    thicksEk[[j]] <- corePost[[1]]$thicks[,theseCounts$ensMember[j]]
    depthsEk[[j]] <- cumdepthEk[,theseCounts$ensMember[j]]
    ageEk[[j]] <- yearsEk[1:length(depthsEk[[j]])]
  }
  
  tri <- which(allResults$counter == k)
  
  allResults$thicknessEk[tri] <- thicksEk
  allResults$depthEk[tri] <- depthsEk
  allResults$ageEk[tri] <- ageEk
  
}

allResultsGood <- mutate(allResults,
                         allNa = map_lgl(depth,\(x) all(is.na(x)))) |> 
  filter(!allNa)




#create figure showing counter success
counterSuccess <- probCut |> 
  group_by(counter) |> 
  summarise(`Highest probability given historical events` = n())

counterSuccessAll <- allResultsGood |> 
  group_by(counter) |> 
  summarise(`Modeled with Eklutna` = n())

cs <- left_join(counterSuccessAll,counterSuccess) |> 
  pivot_longer(-counter,names_to = "type",values_to = "counts") |> 
  filter(type == "Highest probability given historical events")

counterBarPlot <- ggplot(cs) + 
  geom_col(aes(x = counter,y = counts,fill = type),position = "dodge") + 
  theme_bw() + 
  scale_y_continuous("Ensemble members",expand = c(0,0),limits = c(0,250)) +
  scale_fill_brewer(NULL,palette = "Set1") +
  theme(legend.position = "inside",legend.position.inside = c(.5,.8))

counterBoxPlot <- ggplot(allProb) + 
  geom_boxplot(aes(x = counter,y = prob,group = counter,fill = factor(counter))) + 
  scale_fill_viridis_d("Counter") + 
  theme_bw() + 
  ylab("Log probability of historical events") + 
  theme(legend.position = "none")

counterPlot <- egg::ggarrange(plots = list(counterBoxPlot,counterBarPlot),ncol = 1)
ggsave("Counter summary.pdf",plot = counterPlot,height = 8,width = 5)

#now use probCut depths and ages to get all the marker layers. 

skiEvents <- readxl::read_excel(reviewAgesFile,sheet = "Skilak events") |> 
  filter(Depth > 50) |> 
  dplyr::select(event = `FINAL EVENT NAME`,depth = Depth)

eventDists <- function(depth,age,hist.data){
  estage <- approx(depth,age,hist.data$depth)$y |> 
    t() |> 
    as.data.frame() |> 
    setNames(hist.data$event)
  
  return(estage)
}



skEventDistsBest <- map2(probCut$depth,probCut$age,eventDists,hist.data = skiEvents) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  mutate(counters = "historically constrained",
         lake = "Skilak")


skEventDistsAll <- map2(allResultsGood$depth,allResultsGood$age,eventDists,hist.data = skiEvents) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  mutate(counters = "all",
         lake = "Skilak")

skiDists <- ggplot() + 
  geom_histogram(data = skEventDistsAll,aes(x = age,y = after_stat(density), fill = "All"),bins = 30) + 
  geom_histogram(data = skEventDistsBest,aes(x = age,y = after_stat(density), fill = "Historically constrained"),alpha = 0.8,bins = 30) + 
  facet_wrap(~ fct_reorder(event,age),ncol = 4,scales = "free") +
  scale_y_continuous(labels = NULL,breaks = NULL) +
  scale_fill_brewer("Ensemble Members",palette = "Set1") +
  theme_bw() +
  xlab("Age (BP)") + ylab(NULL) +
  ggtitle("Skilak Event Age Distributions")


ggsave(plot = skiDists,filename =  "Skilak event distributions.pdf",height = 16,width = 12)


# Eklutna distributions ---------------------------------------------------


ekEvents <- readxl::read_excel(reviewAgesFile,sheet = "Eklutna counts and events",range = "B2:F33") |> 
  filter(Depth > 61) |> 
  dplyr::select(event = `FINAL EVENT NAME`,depth = Depth)

ekEventDistsBest <- map2(probCut$depthEk,probCut$ageEk,eventDists,hist.data = ekEvents) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  mutate(counters = "historically constrained",
         lake = "Eklutna")


ekEventDistsAll <- map2(allResultsGood$depthEk,allResultsGood$ageEk,eventDists,hist.data = ekEvents) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  mutate(counters = "all",
         lake = "Eklutna")

ekDists <- ggplot() + 
  geom_histogram(data = ekEventDistsAll,aes(x = age,y = after_stat(density), fill = "All"),bins = 30) + 
  geom_histogram(data = ekEventDistsBest,aes(x = age,y = after_stat(density), fill = "Historically constrained"),alpha = 0.8,bins = 30) + 
  facet_wrap(~ fct_reorder(event,age),ncol = 4,scales = "free") +
  scale_y_continuous(labels = NULL,breaks = NULL) +
  scale_fill_brewer("Ensemble Members",palette = "Set1") +
  theme_bw() +
  xlab("Age (BP)") + ylab(NULL) +
  ggtitle("Eklutna Event Age Distributions")


ggsave(plot = ekDists,filename =  "Eklutna event distributions.pdf",height = 16,width = 12)


#create a figure for T5
T5Distribution <- ggplot() + 
  geom_histogram(data = filter(ekEventDistsAll,event == "EK-T5"),aes(x = age,y = after_stat(density), fill = "EK - All"),alpha = 0.6,bins = 30) + 
  geom_histogram(data = filter(ekEventDistsBest,event == "EK-T5"),aes(x = age,y = after_stat(density), fill = "EK - Historically constrained"),alpha = 0.6,bins = 30) + 
  geom_histogram(data = filter(skEventDistsAll,event == "SK-T5"),aes(x = age,y = after_stat(density), fill = "SK - All"),alpha = 0.6,bins = 30) + 
  geom_histogram(data = filter(skEventDistsBest,event == "SK-T5"),aes(x = age,y = after_stat(density), fill = "SK - Historically constrained"),alpha = 0.6,bins = 30) +
  scale_y_continuous(labels = NULL,breaks = NULL) +
  scale_fill_brewer("Ensemble Members",palette = "Set1") +
  theme_bw() +
  xlab("Age (BP)") + ylab(NULL) +
  ggtitle("SK & EK T5 Event Age Distributions")
ggsave(plot = T5Distribution,filename =  "T5_distribution.pdf")

# Create a table of event results -----------------------------------------

allEvents <- bind_rows(skEventDistsBest,skEventDistsAll,ekEventDistsBest,ekEventDistsAll)

eventSummaryTable <- allEvents |> 
  group_by(event,counters,lake) |> 
  summarise(medianAge = median(age,na.rm = TRUE),
            as_tibble_row(quantile(age,probs = c(0.025,.25,.75,.975),na.rm = TRUE)))

write_csv(x = eventSummaryTable,file = "Event Summary Table.csv")

# Combine both lakes to come up with best estimates for correlated events
corrEvents <- readxl::read_excel(reviewAgesFile,
                                 sheet = "Preliminar ages",
                                 range = "A2:F43") |> 
  dplyr::select(1,2,6) |>
  setNames(c("Corr","Sk","Ek")) |> 
  filter(!is.na(Corr)) |> 
  pivot_longer(-Corr,values_to = "event")

corrEventData <- allEvents |> 
  filter(event %in% corrEvents$event) |> 
  rowwise() |> 
  mutate(corrEvent = map_chr(event,\(x) corrEvents$Corr[event == corrEvents$event])) |> 
  dplyr::select(corrEvent,age,counters)


corrEventSummaryTable <- corrEventData |> 
  group_by(corrEvent,counters) |> 
  summarise(medianAge = median(age,na.rm = TRUE),
            as_tibble_row(quantile(age,probs = c(0.025,.25,.75,.975),na.rm = TRUE)))

write_csv(x = corrEventSummaryTable,file = "Correllated Event Summary Table.csv")

#load in priors allMlAges (skilak) and ekPrior (eklutna)
EkCorePost <- readRDS("EklutnaPostData.RDS")
ekEns <- EkCorePost[[1]]$thicks
ekDepthEns <- apply(ekEns,2,\(x) cumsum(x) / 10 + 60.19)

ekMls <- readxl::read_excel(skiDataPath,sheet = "Eklutna counts and events",skip = 1,range = "A2:F35") |> 
  dplyr::select(Depth, ML = `CORRELATED EVENT NAME`) |> 
  filter(!is.na(ML), Depth >= 60)

for(m in 1:nrow(ekMls)){
  thisOne <- data.frame(age = unlist(apply(ekDepthEns,2,\(x) which.min(abs(x - ekMls$Depth[m])))) + 21,
                        ML = ekMls$ML[m])
  if(m == 1){
    ekPrior <- thisOne
  }else{
    ekPrior <- bind_rows(ekPrior,thisOne)
  }
}




inputFolders <- list.dirs(here("SkilakEklutna","Input"))[str_detect(list.dirs(here("SkilakEklutna","Input")),pattern = "\\d$")]
allMlAges <- data.frame(age = NA,name = NA,counter = NA)
for(k in 1:length(inputFolders)){
  thisDir <- file.path("SkilakEklutna/Input",k)
  #load(file.path(thisDir,"SkiEkOut.Rdata"))
  thisMlPrior <- read_csv(file.path(thisDir,"MarkerLayers.csv"))
  thisMlPrior$counter <- basename(thisDir)
  if(k == 1){
    allMlPrior <- thisMlPrior
  }else{
    allMlPrior <- bind_rows(allMlPrior,thisMlPrior)
  }
  mlPriors <- read_csv(file.path(thisDir,"MarkerLayers.csv"))
  
  
  theseMlAges <- data.frame(age = as.numeric(mlPriors$SK), name = mlPriors$name, counter = basename(thisDir))
  allMlAges <- bind_rows(allMlAges,theseMlAges)
  
}


#a correlated event plot
corrEventPriors <- allMlAges |> 
  filter(name %in% unique(corrEventData$corrEvent)) |> 
  rename(corrEvent = name) |> 
  mutate(counters = NA,
         lake = "Skilak")

#a correlated event plot
corrEventPriors <- ekPrior |> 
  filter(ML %in% unique(corrEventData$corrEvent)) |> 
  rename(corrEvent = ML) |> 
  mutate(counters = NA,
         lake = "Eklutna") |> 
  bind_rows(corrEventPriors)
# 
corrEventDataPlot <- corrEventData |>
  mutate(lake = NA) |>
  bind_rows(corrEventPriors)
  

corrDists <- ggplot() + 
  geom_histogram(data = corrEventDataPlot,aes(x = age,y = after_stat(density),fill = counters,alpha = counters),bins = 30,position = "identity") +
  geom_density(data = corrEventDataPlot,aes(x = age,y = after_stat(density),fill = NA,color = lake,alpha = counters),position = "identity") +
  facet_wrap(~ fct_reorder(corrEvent,age),ncol = 4,scales = "free") +
  scale_y_continuous(labels = NULL,breaks = NULL) +
  scale_fill_brewer("Posteriors",palette = "Set1", breaks = ~ .x[!is.na(.x)]) +
  scale_color_brewer("Priors",palette = "Dark2",breaks = ~ .x[!is.na(.x)]) +
  scale_alpha_manual(values = c(1,.8,0),guide = "none") + 
  theme_bw() +
  xlab("Age (BP)") + ylab(NULL) +
  ggtitle("Correlated Event Age Distributions")

corrDists
ggsave(plot = corrDists,filename =  "Correlated event distributions.pdf",height = 8,width = 8)

# Recurrence intervals SKilak ----------------------------------------------------

recurrence <- readxl::read_excel(reviewAgesFile,sheet = "Recurrence")

skMega <- dplyr::select(recurrence,`Megathrust Skilak`) |>  filter(`Megathrust Skilak` != "incomplete record") |>  drop_na() 

skRec <- map2(probCut$depth,probCut$age,eventDists,hist.data = skiEvents) |> 
  list_rbind() |> 
  mutate(`SK-EQ2` = -14) |> 
  dplyr::select(!!skMega$`Megathrust Skilak`)

diffs <-  skRec[,-1] - skRec[,-ncol(skRec)]
newNames <- paste(names(skRec)[-ncol(skRec)],"to",names(skRec)[-1])
names(diffs) <- newNames


for(i in 1:nrow(diffs)){
  #fit distribution to each ensemble member in diffs
  data <- as.numeric(na.omit(as.numeric(diffs[i,])))
  
  # Fit distributions using fitdistrplus
  fit_weibull <- fitdist(data, "weibull")
  fit_exponential <- fitdist(data, "exp")  # Exponential distribution
  fit_normal <- fitdist(data, "norm")
  
  # Perform Goodness-of-Fit Tests (Kolmogorov-Smirnov for continuous distributions)
  ks_weibull <- ks.test(data, "pweibull", shape=fit_weibull$estimate["shape"], scale=fit_weibull$estimate["scale"])
  ks_exponential <- ks.test(data, "pexp", rate=fit_exponential$estimate["rate"])
  ks_normal <- ks.test(data, "pnorm", mean=fit_normal$estimate["mean"], sd=fit_normal$estimate["sd"])
  
  fit_params <- data.frame(shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"], rate = fit_exponential$estimate["rate"], mean =fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])
  
  if(i == 1){
    all_fit_params <- fit_params
  }else{
    all_fit_params <- bind_rows(all_fit_params,fit_params)
  }
  
  # Compare distributions using AIC
  if(i == 1){
    aic_values <- data.frame(
      Distribution = c("Weibull", "Exponential", "Normal"),
      AIC = c(
        fit_weibull$aic,
        fit_exponential$aic,
        fit_normal$aic
      )
    )
  }else{
    aic_values <- bind_cols(aic_values,AIC = c(
      fit_weibull$aic,
      fit_exponential$aic,
      fit_normal$aic
    ),.name_repair = "minimal")
  }
  
}

aicQuants <- t(aic_values[,-1]) |> 
  as.data.frame() |> 
  setNames(aic_values[,1]) |> 
  apply(2,quantile,probs = c(.025, .5, .975)) |> 
  as.data.frame()
aicQuants$Lake <- "Skilak"
aicQuants$Event <- "Megathrust"


aicWinners <- t(aic_values[,-1]) |> 
  apply(1,which.min) |> 
  table() |> 
  as.data.frame()

aicWinners$Lake <- "Skilak"
aicWinners$Event <- "Megathrust"



allAicWinners <- aicWinners


allAic <- aicQuants

#make a plot with the distributions and their uncertainties
distQuants <- apply(all_fit_params,2,quantile,probs = c(.025, .5, .975)) |> as.data.frame()

allDatMega <- unlist(diffs) |> na.omit() |> as.numeric()

xr <- range(allDatMega)
xseq <- seq(xr[1],xr[2])

weib <- data.frame(x = xseq,
                   low = dweibull(xseq,shape = distQuants$shape[1],scale = distQuants$scale[1]) ,
                   mid = dweibull(xseq,shape = distQuants$shape[2],scale = distQuants$scale[2]), 
                   hi = dweibull(xseq,shape = distQuants$shape[3],scale = distQuants$scale[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Weibull")

expo <- data.frame(x = xseq,
                   low = dexp(xseq,rate = distQuants$rate[1]) ,
                   mid = dexp(xseq,rate = distQuants$rate[2]), 
                   hi = dexp(xseq,rate = distQuants$rate[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Exponential")


normo <- data.frame(x = xseq,
                    low = dnorm(xseq,mean = distQuants$mean[1],sd = distQuants$sd[1]) ,
                    mid = dnorm(xseq,mean = distQuants$mean[2],sd = distQuants$sd[2]), 
                    hi = dnorm(xseq,mean = distQuants$mean[3],sd = distQuants$sd[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Normal")

distribs <- bind_rows(weib, expo, normo)




skMegathrustModeled <- ggplot() + 
  geom_density(aes(x = allDatMega), fill = "grey60",color = NA) + 
  geom_line(data = distribs,aes(x = x, y = value, linetype = quantile, color = distribution)) + 
  scale_linetype_manual(values = c(2,2,1),guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8)) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Skilak Megathrust Recurrence Intervals") 




recPlot <- diffs |> 
  pivot_longer(everything(),names_to = "Events",values_to = "Duration") |> 
  mutate(Events = fct_relevel(Events,newNames))

skMegathrust <- ggplot(recPlot) + 
  geom_histogram(aes(x = Duration,y = after_stat(density),fill = Events),position = "identity",binwidth = 5,alpha = .7) +
  geom_density(aes(x = Duration,y = after_stat(density) * 5,color = "Total"),fill = NA) +
  scale_color_manual("",values = "black") +
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_viridis_d() +
  theme_bw() +
  xlim(0,600) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Skilak Megathrust Recurrence Intervals") 

#save for table
recurrenceTable <- recPlot |> 
  mutate(lake = "Skilak",
         source = "Megathrust")

# aggregate for distribution
recDistTable <- recPlot |> 
  mutate(window = cut_interval(Duration,length = 5)) |> 
  group_by(window) |> 
  dplyr::summarize(count = n()) |> 
  mutate(lake = "Skilak",
         source = "Megathrust")



# repeat for Skilak intraplate

  
skIntra <- dplyr::select(recurrence,`Intraplate Skilak`) |> drop_na() |> filter(`Intraplate Skilak` != "incomplete record")

skRec <- map2(probCut$depth,probCut$age,eventDists,hist.data = skiEvents) |> 
  list_rbind() |> 
  mutate(`SK-EQ0` = -68,
         `SK-EQ1` = -16,
         `SK-EQ3` = -4) |> 
  dplyr::select(!!skIntra$`Intraplate Skilak`)

diffs <-  skRec[,-1] - skRec[,-ncol(skRec)]
newNames <- paste(names(skRec)[-ncol(skRec)],"to",names(skRec)[-1])
names(diffs) <- newNames


for(i in 1:nrow(diffs)){
  #fit distribution to each ensemble member in diffs
  data <- as.numeric(na.omit(as.numeric(diffs[i,])))
  
  # Fit distributions using fitdistrplus
  fit_weibull <- fitdist(data, "weibull")
  fit_exponential <- fitdist(data, "exp")  # Exponential distribution
  fit_normal <- fitdist(data, "norm")
  
  # Perform Goodness-of-Fit Tests (Kolmogorov-Smirnov for continuous distributions)
  ks_weibull <- ks.test(data, "pweibull", shape=fit_weibull$estimate["shape"], scale=fit_weibull$estimate["scale"])
  ks_exponential <- ks.test(data, "pexp", rate=fit_exponential$estimate["rate"])
  ks_normal <- ks.test(data, "pnorm", mean=fit_normal$estimate["mean"], sd=fit_normal$estimate["sd"])
  
  fit_params <- data.frame(shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"], rate = fit_exponential$estimate["rate"], mean =fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])
  
  if(i == 1){
    all_fit_params <- fit_params
  }else{
    all_fit_params <- bind_rows(all_fit_params,fit_params)
  }
  
  # Compare distributions using AIC
  if(i == 1){
    aic_values <- data.frame(
      Distribution = c("Weibull", "Exponential", "Normal"),
      AIC = c(
        fit_weibull$aic,
        fit_exponential$aic,
        fit_normal$aic
      )
    )
  }else{
    aic_values <- bind_cols(aic_values,AIC = c(
      fit_weibull$aic,
      fit_exponential$aic,
      fit_normal$aic
    ),.name_repair = "minimal")
  }
  
}
aicQuants <- t(aic_values[,-1]) |> 
  as.data.frame() |> 
  setNames(aic_values[,1]) |> 
  apply(2,quantile,probs = c(.025, .5, .975)) |> 
  as.data.frame()
aicQuants$Lake <- "Skilak"
aicQuants$Event <- "Intraplate"


aicWinners <- t(aic_values[,-1]) |> 
  apply(1,which.min) |> 
  table() |> 
  as.data.frame()

aicWinners$Lake <- "Skilak"
aicWinners$Event <- "Intraplate"

allAic <- bind_rows(allAic,aicQuants)
allAicWinners <- bind_rows(allAicWinners,aicWinners)

#make a plot with the distributions and their uncertainties
distQuants <- apply(all_fit_params,2,quantile,probs = c(.025, .5, .975)) |> as.data.frame()

allDat <- unlist(diffs) |> na.omit() |> as.numeric()

xr <- range(allDat)
xseq <- seq(xr[1],xr[2])

weib <- data.frame(x = xseq,
                   low = dweibull(xseq,shape = distQuants$shape[1],scale = distQuants$scale[1]) ,
                   mid = dweibull(xseq,shape = distQuants$shape[2],scale = distQuants$scale[2]), 
                   hi = dweibull(xseq,shape = distQuants$shape[3],scale = distQuants$scale[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Weibull")

expo <- data.frame(x = xseq,
                   low = dexp(xseq,rate = distQuants$rate[1]) ,
                   mid = dexp(xseq,rate = distQuants$rate[2]), 
                   hi = dexp(xseq,rate = distQuants$rate[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Exponential")


normo <- data.frame(x = xseq,
                    low = dnorm(xseq,mean = distQuants$mean[1],sd = distQuants$sd[1]) ,
                    mid = dnorm(xseq,mean = distQuants$mean[2],sd = distQuants$sd[2]), 
                    hi = dnorm(xseq,mean = distQuants$mean[3],sd = distQuants$sd[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Normal")

distribs <- bind_rows(weib, expo, normo)




skIntraplateModeled <- ggplot() + 
  geom_density(aes(x = allDat), fill = "grey60",color = NA) + 
  geom_line(data = distribs,aes(x = x, y = value, linetype = quantile, color = distribution)) + 
  scale_linetype_manual(values = c(2,2,1),guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8)) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Skilak Intraplate Recurrence Intervals") 






recPlot <- diffs |> 
  pivot_longer(everything(),names_to = "Events",values_to = "Duration") |> 
  mutate(Events = fct_relevel(Events,newNames))

skIntraplate <- ggplot(recPlot) + 
  geom_histogram(aes(x = Duration,y = after_stat(density),fill = Events),position = "identity",binwidth = 5,alpha = .7) +
  geom_density(aes(x = Duration,y = after_stat(density) * 5,color = "Total"),fill = NA) +
  scale_color_manual("",values = "black") +
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_viridis_d() +
  theme_bw() +
  xlim(0,600) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Skilak Intraplate Recurrence Intervals") 
  
#save for table
recurrenceTable <- recPlot |> 
  mutate(lake = "Skilak",
         source = "Intraplate") |> 
  bind_rows(recurrenceTable)


# aggregate for distribution
recDistTable <- recPlot |> 
  mutate(window = cut_interval(Duration,length = 5)) |> 
  group_by(window) |> 
  dplyr::summarize(count = n()) |> 
  mutate(lake = "Skilak",
         source = "Intraplate") |> 
  bind_rows(recDistTable)


SkiRecurrModeled <- egg::ggarrange(plots = list(skMegathrustModeled,skIntraplateModeled),ncol = 1)

ggsave(SkiRecurrModeled,filename = "Skilak Recurrence Intervals Modeled.pdf")


SkiRecurr <- egg::ggarrange(plots = list(skMegathrust,skIntraplate),ncol = 1)

ggsave(SkiRecurr,filename = "Skilak Recurrence Intervals.pdf")




# now eklutna

ekMega <- dplyr::select(recurrence,`Megathrust Eklutna`) |>  drop_na()

ekRec <- map2(probCut$depthEk,probCut$ageEk,eventDists,hist.data = ekEvents) |> 
  list_rbind() |> 
  mutate(`EK-EQ1` = -14) |> 
  #select(-`EK-EQ21`) |>  #too deep
  dplyr::select(!!ekMega$`Megathrust Eklutna`)

diffs <-  ekRec[,-1] - ekRec[,-ncol(ekRec)]
newNames <- paste(names(ekRec)[-ncol(ekRec)],"to",names(ekRec)[-1])
names(diffs) <- newNames


for(i in 1:nrow(diffs)){
  #fit distribution to each ensemble member in diffs
  data <- as.numeric(na.omit(as.numeric(diffs[i,])))
  
  # Fit distributions using fitdistrplus
  fit_weibull <- fitdist(data, "weibull")
  fit_exponential <- fitdist(data, "exp")  # Exponential distribution
  fit_normal <- fitdist(data, "norm")
  
  # Perform Goodness-of-Fit Tests (Kolmogorov-Smirnov for continuous distributions)
  ks_weibull <- ks.test(data, "pweibull", shape=fit_weibull$estimate["shape"], scale=fit_weibull$estimate["scale"])
  ks_exponential <- ks.test(data, "pexp", rate=fit_exponential$estimate["rate"])
  ks_normal <- ks.test(data, "pnorm", mean=fit_normal$estimate["mean"], sd=fit_normal$estimate["sd"])
  
  fit_params <- data.frame(shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"], rate = fit_exponential$estimate["rate"], mean =fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])
  
  if(i == 1){
    all_fit_params <- fit_params
  }else{
    all_fit_params <- bind_rows(all_fit_params,fit_params)
  }
  
  # Compare distributions using AIC
  if(i == 1){
    aic_values <- data.frame(
      Distribution = c("Weibull", "Exponential", "Normal"),
      AIC = c(
        fit_weibull$aic,
        fit_exponential$aic,
        fit_normal$aic
      )
    )
  }else{
    aic_values <- bind_cols(aic_values,AIC = c(
      fit_weibull$aic,
      fit_exponential$aic,
      fit_normal$aic
    ),.name_repair = "minimal")
  }
  
}

aicQuants <- t(aic_values[,-1]) |> 
  as.data.frame() |> 
  setNames(aic_values[,1]) |> 
  apply(2,quantile,probs = c(.025, .5, .975)) |> 
  as.data.frame()

aicQuants$Lake <- "Eklutna"
aicQuants$Event <- "Megathrust"


allAic <- bind_rows(allAic,aicQuants)

#make a plot with the distributions and their uncertainties
distQuants <- apply(all_fit_params,2,quantile,probs = c(.025, .5, .975)) |> as.data.frame()

allDatMega <- unlist(diffs) |> na.omit() |> as.numeric()

xr <- range(allDatMega)
xseq <- seq(xr[1],xr[2])

weib <- data.frame(x = xseq,
                   low = dweibull(xseq,shape = distQuants$shape[1],scale = distQuants$scale[1]) ,
                   mid = dweibull(xseq,shape = distQuants$shape[2],scale = distQuants$scale[2]), 
                   hi = dweibull(xseq,shape = distQuants$shape[3],scale = distQuants$scale[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Weibull")

expo <- data.frame(x = xseq,
                   low = dexp(xseq,rate = distQuants$rate[1]) ,
                   mid = dexp(xseq,rate = distQuants$rate[2]), 
                   hi = dexp(xseq,rate = distQuants$rate[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Exponential")


normo <- data.frame(x = xseq,
                    low = dnorm(xseq,mean = distQuants$mean[1],sd = distQuants$sd[1]) ,
                    mid = dnorm(xseq,mean = distQuants$mean[2],sd = distQuants$sd[2]), 
                    hi = dnorm(xseq,mean = distQuants$mean[3],sd = distQuants$sd[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Normal")

distribs <- bind_rows(weib, expo, normo)




ekMegathrustModeled <- ggplot() + 
  geom_density(aes(x = allDatMega), fill = "grey60",color = NA) + 
  geom_line(data = distribs,aes(x = x, y = value, linetype = quantile, color = distribution)) + 
  scale_linetype_manual(values = c(2,2,1),guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8)) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Eklutna Megathrust Recurrence Intervals") 



recPlot <- diffs |> 
  pivot_longer(everything(),names_to = "Events",values_to = "Duration") |> 
  mutate(Events = fct_relevel(Events,newNames))

ekMegathrust <- ggplot(recPlot) + 
  geom_histogram(aes(x = Duration,y = after_stat(density),fill = Events),position = "identity",binwidth = 5,alpha = .7) +
  geom_density(aes(x = Duration,y = after_stat(density) * 5,color = "Total"),fill = NA) +
  scale_color_manual("",values = "black") +
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_viridis_d() +
  theme_bw() +
  xlim(0,600) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Eklutna Megathrust Recurrence Intervals") 

#save for table
recurrenceTable <- recPlot |> 
  mutate(lake = "Eklutna",
         source = "Megathrust") |> 
  bind_rows(recurrenceTable)


# aggregate for distribution
recDistTable <- recPlot |> 
  mutate(window = cut_interval(Duration,length = 5)) |> 
  group_by(window) |> 
  dplyr::summarize(count = n()) |> 
  mutate(lake = "Eklutna",
         source = "Megathrust") |> 
  bind_rows(recDistTable)


# repeat for Eklutna intraplate


ekIntra <- dplyr::select(recurrence,`Intraplate Eklutna`) |> drop_na()

ekRec <- map2(probCut$depthEk,probCut$ageEk,eventDists,hist.data = ekEvents) |> 
  list_rbind() |> 
  mutate(`EK-EQ0` = -68) |> 
  dplyr::select(!!ekIntra$`Intraplate Eklutna`)

diffs <-  ekRec[,-1] - ekRec[,-ncol(ekRec)]
newNames <- paste(names(ekRec)[-ncol(ekRec)],"to",names(ekRec)[-1])
names(diffs) <- newNames

for(i in 1:nrow(diffs)){
#fit distribution to each ensemble member in diffs
data <- as.numeric(na.omit(as.numeric(diffs[i,])))

# Fit distributions using fitdistrplus
fit_weibull <- fitdist(data, "weibull")
fit_exponential <- fitdist(data, "exp")  # Exponential distribution
fit_normal <- fitdist(data, "norm")

# Perform Goodness-of-Fit Tests (Kolmogorov-Smirnov for continuous distributions)
ks_weibull <- ks.test(data, "pweibull", shape=fit_weibull$estimate["shape"], scale=fit_weibull$estimate["scale"])
ks_exponential <- ks.test(data, "pexp", rate=fit_exponential$estimate["rate"])
ks_normal <- ks.test(data, "pnorm", mean=fit_normal$estimate["mean"], sd=fit_normal$estimate["sd"])

fit_params <- data.frame(shape = fit_weibull$estimate["shape"], scale = fit_weibull$estimate["scale"], rate = fit_exponential$estimate["rate"], mean =fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])

if(i == 1){
  all_fit_params <- fit_params
}else{
  all_fit_params <- bind_rows(all_fit_params,fit_params)
}

# Compare distributions using AIC
if(i == 1){
  aic_values <- data.frame(
    Distribution = c("Weibull", "Exponential", "Normal"),
    AIC = c(
      fit_weibull$aic,
      fit_exponential$aic,
      fit_normal$aic
    )
  )
}else{
  aic_values <- bind_cols(aic_values,AIC = c(
    fit_weibull$aic,
    fit_exponential$aic,
    fit_normal$aic
  ),.name_repair = "minimal")
}

}

aicQuants <- t(aic_values[,-1]) |> 
  as.data.frame() |> 
  setNames(aic_values[,1]) |> 
  apply(2,quantile,probs = c(.025, .5, .975)) |> 
  as.data.frame()
aicQuants$Lake <- "Eklutna"
aicQuants$Event <- "Intraplate"


aicWinners <- t(aic_values[,-1]) |> 
  apply(1,which.min) |> 
  table() |> 
  as.data.frame()

aicWinners$Lake <- "Eklutna"
aicWinners$Event <- "Intraplate"

allAic <- bind_rows(allAic,aicQuants)
allAicWinners <- bind_rows(allAicWinners,aicWinners)



allAic$quantiles <- c(.025, .5, .975)

allAic <- dplyr::select(allAic,quantiles, everything())

write_csv(allAic,file = "AIC quantiles.csv")

#make a plot with the distributions and their uncertainties
distQuants <- apply(all_fit_params,2,quantile,probs = c(.025, .5, .975)) |> as.data.frame()

allDat <- unlist(diffs) |> na.omit() |> as.numeric()

xr <- range(allDat)
xseq <- seq(xr[1],xr[2])

weib <- data.frame(x = xseq,
                   low = dweibull(xseq,shape = distQuants$shape[1],scale = distQuants$scale[1]) ,
                   mid = dweibull(xseq,shape = distQuants$shape[2],scale = distQuants$scale[2]), 
                   hi = dweibull(xseq,shape = distQuants$shape[3],scale = distQuants$scale[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Weibull")

expo <- data.frame(x = xseq,
                   low = dexp(xseq,rate = distQuants$rate[1]) ,
                   mid = dexp(xseq,rate = distQuants$rate[2]), 
                   hi = dexp(xseq,rate = distQuants$rate[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Exponential")


normo <- data.frame(x = xseq,
                   low = dnorm(xseq,mean = distQuants$mean[1],sd = distQuants$sd[1]) ,
                   mid = dnorm(xseq,mean = distQuants$mean[2],sd = distQuants$sd[2]), 
                   hi = dnorm(xseq,mean = distQuants$mean[3],sd = distQuants$sd[3])) |> 
  pivot_longer(-x,names_to = "quantile",values_to = "value") |> 
  mutate(distribution = "Normal")

distribs <- bind_rows(weib, expo, normo)




ekIntraplateModeled <- ggplot() + 
  geom_density(aes(x = allDat), fill = "grey60",color = NA) + 
  geom_line(data = distribs,aes(x = x, y = value, linetype = quantile, color = distribution)) + 
  scale_linetype_manual(values = c(2,2,1),guide = "none") +
  scale_color_brewer(palette = "Dark2") +
  theme_bw() + 
  theme(legend.title = element_blank(),legend.position = c(0.8,0.8)) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Eklutna Intraplate Recurrence Intervals") 

  
  
  


recPlot <- diffs |> 
  pivot_longer(everything(),names_to = "Events",values_to = "Duration") |> 
  mutate(Events = fct_relevel(Events,newNames))

ekIntraplate <- ggplot(recPlot) + 
  geom_histogram(aes(x = Duration,y = after_stat(density),fill = Events),position = "identity",binwidth = 5,alpha = .7) +
  geom_density(aes(x = Duration,y = after_stat(density) * 5,color = "Total"),fill = NA) +
  scale_color_manual("",values = "black") +
  #scale_fill_brewer(palette = "Paired") + 
  scale_fill_viridis_d() +
  theme_bw() +
  xlim(0,600) +
  xlab("Duration (yr)") + 
  ylab("Probability density") +
  ggtitle("Eklutna Intraplate Recurrence Intervals") 

#save for table
recurrenceTable <- recPlot |> 
  mutate(lake = "Eklutna",
         source = "Intraplate") |> 
  bind_rows(recurrenceTable)

# aggregate for distribution
recDistTable <- recPlot |> 
  mutate(window = cut_interval(Duration,length = 5)) |> 
  group_by(window) |> 
  dplyr::summarize(count = n()) |> 
  mutate(lake = "Eklutna",
         source = "Intraplate") |> 
  bind_rows(recDistTable)


EkRecurrModeled <- egg::ggarrange(plots = list(ekMegathrustModeled,ekIntraplateModeled),ncol = 1)


EkRecurr <- egg::ggarrange(plots = list(ekMegathrust,ekIntraplate),ncol = 1)

ggsave(EkRecurrModeled,filename = "Eklutna Recurrence Intervals Modeled.pdf")


ggsave(EkRecurr,filename = "Eklutna Recurrence Intervals.pdf")

recurrenceSummaryTable <- recurrenceTable |> 
  group_by(Events,source,lake) |> 
  summarise(medianDuration = median(Duration,na.rm = TRUE),
            as_tibble_row(quantile(Duration,probs = c(0.025,.25,.75,.975),na.rm = TRUE)))

write_csv(x = recurrenceSummaryTable,file = "Recurrence Summary Table.csv")
write_csv(x = recurrenceTable,file = "All recurrence data long.csv")

recurrenceDistributionTable <- recDistTable |> 
  mutate(numericWindow = as.numeric(map_chr(str_extract_all(window,pattern = "\\d{1,}"),pluck,1)) + 12.5) |> 
  pivot_wider(names_from = c("lake","source"),values_from = count,id_cols = numericWindow) |> 
  arrange(numericWindow)

write_csv(x = recurrenceDistributionTable,file = "Recurrence Distribution Table.csv")



# get age distributions to compare with radiocarbon dates -----------------

#skilak 14c Dates 
ski14c <- readxl::read_excel(skiDataPath,sheet = "Skilak 14C") |> 
  dplyr::select(event = `Lab ID`,depth = `Composite depth (cm)`,median14C_age = `Median cal age (yr BP)`) |> 
  mutate(event = as.character(event))


ski14c_comp <- map2(probCut$depth,probCut$age,eventDists,hist.data = ski14c) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  dplyr::summarize(.by = event,medianVarveAge = median(age,na.rm = TRUE)) |> 
  left_join(ski14c,by = "event") |> 
  dplyr::select(event,depth,everything()) |> 
  mutate(ageDifference = medianVarveAge - median14C_age,
         lake = "Skilak")


#now ekutna
ek14c <- readxl::read_excel(skiDataPath,sheet = "Eklutna 14C") |> 
  mutate(`Lab ID` = paste0(`Core section`,"-",`Composite depth at Site 2(cm)`)) |> 
  dplyr::select(event = `Lab ID`,depth = `Composite depth at Site 2(cm)`,median14C_age = `Median cal age (yr BP)`) |> 
  mutate(event = as.character(event))


ek14c_comp <- map2(probCut$depthEk,probCut$ageEk,eventDists,hist.data = ek14c) |> 
  list_rbind() |> 
  pivot_longer(everything(),names_to = "event",values_to = "age") |> 
  dplyr::summarize(.by = event,medianVarveAge = median(age,na.rm = TRUE)) |> 
  left_join(ek14c,by = "event") |> 
  dplyr::select(event,depth,everything()) |> 
  mutate(ageDifference = medianVarveAge - median14C_age,
         lake = "Eklutna")

combined14CEstimates <- bind_rows(ski14c_comp,ek14c_comp) |> rename(identifier = event)
write_csv(combined14CEstimates,file = "Varve-14C comparison.csv")

#depth version
bsd <- 0:3000
bestSkEnsemble <- map2(probCut$depth,probCut$age,\(x,y){
  good <- which(!is.na(x) & !is.na(y))
  approxExtrap(x[good],y[good],xout = bsd)$y}) |>
  list_c() |>
  matrix(ncol = 1000)

fullSkEnsemble <- map2(allResultsGood$depth,allResultsGood$age,\(x,y){
  good <- which(!is.na(x) & !is.na(y))
  approxExtrap(x[good],y[good],xout = bsd)$y}) |>
  list_c() |>
  matrix(ncol = 23000)

bestEkEnsemble <- map2(probCut$depthEk,probCut$ageEk,\(x,y){
  good <- which(!is.na(x) & !is.na(y))
  approx(x[good],y[good],xout = bsd)$y}) |>
  list_c() |>
  matrix(ncol = 1000)

goodRows <- which(rowSums(!is.na(bestEkEnsemble)) > 100)


ekEnsOut <- cbind(bsd[goodRows],bestEkEnsemble[goodRows,]) |> as.data.frame() |> setNames(c("depth (cm)",paste0("age (yr BP) - ens",1:1000)))
write_csv(ekEnsOut,file = "EkBestEns_Core2.csv")

summaryEkOut <- cbind(bsd[goodRows],t(apply(bestEkEnsemble[goodRows,],1,quantile, probs = c(0.025,.5, 0.975),na.rm = TRUE)))  |> as.data.frame() |> setNames(c("depth (cm)","age (yr BP) 2.5%","age (yr BP) 50%","age (yr BP) 97.5%"))
write_csv(summaryEkOut,file = "EkBestSumm_Core2.csv")


# create a lipd model with these probs to use plotModelDistributions()
mod <- list

eml <- unique(ekPrior$ML)
emldepthsinskilak <- c(163.5,272.1,309.1,367.4,463.8,494.7,567.3,627,764.2,1044.4)
distributionTable <- vector(mode = "list",length = length(eml))

for(i in 1:length(eml)){
  thisml <- filter(ekPrior,ML == eml[i])
  r <- range(thisml$age)
  f <- ecdf(thisml$age)
  age <- list(values = seq(r[1],r[2]),
              variableName = "age",
              units= "yr BP")
  probabilityDensity <- list(values = c(diff(f(age$values)),0),
                             variableName = "probabilityDensity")
  
  distributionTable[[i]] <- list(age = age, probabilityDensity = probabilityDensity, depth = emldepthsinskilak[i])
}

L <- list()
L$chronData[[1]]$model[[1]]$distributionTable <- distributionTable


# create a depth distribution table ---------------------------------------



eml <- unique(ekPrior$ML)
emldepthsinskilak <- c(163.5,272.1,309.1,367.4,463.8,494.7,567.3,627,764.2,1044.4)
emlagesinskilak <- c()
distributionTable <- vector(mode = "list",length = length(eml))

for(i in 1:length(eml)){
  thisml <- filter(ekPrior,ML == eml[i])
  r <- range(thisml$age)
  f <- ecdf(thisml$age)
  age <- list(values = seq(r[1],r[2]),
              variableName = "age",
              units= "yr BP")
  
  aged <- f(age$values)
  
  
  probabilityDensity <- list(values = c(diff(f(age$values)),0),
                             variableName = "probabilityDensity")
  
  emlagesinskilak[i] <- mean(r)
  
  distributionTable[[i]] <- list(age = age, probabilityDensity = probabilityDensity, depth = emldepthsinskilak[i])
}

L <- list()
L$chronData[[1]]$model[[1]]$distributionTable <- distributionTable




# create a second distribution table for historical constains -------------

distributionTable <- vector(mode = "list",length = nrow(hist.data))

for(i in 1:nrow(hist.data)){
  ar <- seq(-20,20) + hist.data$mid.age[i]
  age <- list(values = ar,
              variableName = "age",
              units= "yr BP")
  
  aged <- dnorm(ar,mean = hist.data$mid.age[i], sd = (hist.data$mid.age[i] - hist.data$min.age[i])/2)
  
  
  probabilityDensity <- list(values = aged,
                             variableName = "probabilityDensity")
  
  
  
  distributionTable[[i]] <- list(age = age, probabilityDensity = probabilityDensity, depth = hist.data$depth[i])
}

HC <- list()
HC$chronData[[1]]$model[[1]]$distributionTable <- distributionTable




# plotTimeseriesEnsRibbons(add.to.plot = start, Y = 22:3444,X = fullSkEnsemble,probs = c(0.025,0.975),alp = 0.7,color.high = "#E41A1C") |> 
#   plotTimeseriesEnsRibbons( Y = bsy,X = bestSkEnsemble,probs = c(0.025,0.975),alp = 0.7,color.high = "#377EB8") + 
#   coord_flip() + 
#   scale_x_reverse() |> 
#   plotModelDistributions(L)
# ens testing -------------------------------------------------------------
# 
# bestSkEnsRange <- apply(bestSkEnsemble,1,range,na.rm = FALSE)
# fullSkEnsRange <- apply(fullSkEnsemble,1,range,na.rm = FALSE)
bestSkEns95d <- apply(bestSkEnsemble,1,quantile,probs = c(0.025,0.975),na.rm = TRUE) 
fullSkEns95d <- apply(fullSkEnsemble,1,quantile,probs = c(0.025,0.975),na.rm = TRUE) 

plot(bsd, bestSkEns95d[1,],type = "l")
lines(bsd, bestSkEns95d[2,])

plot(bsd, fullSkEns95d[1,],type = "l")
lines(bsd, fullSkEns95d[2,])

bsy <- 22:2500
#now interpolate into age-depth space
bestSkEns95Lo <- approx(bestSkEns95d[1,],bsd,xout = bsy)$y
bestSkEns95Hi<- approx(bestSkEns95d[2,],bsd,xout = bsy)$y

fullSkEns95Lo <- approx(fullSkEns95d[1,],bsd,xout = bsy)$y
fullSkEns95Hi<- approx(fullSkEns95d[2,],bsd,xout = bsy)$y

SC <- lipdR::readLipd("Skilak.Radiocarbon.lpd")
hist.names <- c("GBT2","GBT6","GT6SD2")

start <- ggplot() + 
  geom_ribbon(data = NULL,aes(x = bsy,ymin = fullSkEns95Lo,ymax = fullSkEns95Hi,fill = "Full ensemble"), alpha = 1) + 
  geom_ribbon(data = NULL,aes(x = bsy,ymin = bestSkEns95Lo,ymax = bestSkEns95Hi,fill = "Historically constrained"), alpha = 1) + 
  geom_line(data = og,aes(y = depth,x = age,group = counter, color = "Original varve counters"),linewidth = .1) +
  xlim(c(-25,2600)) +
  geom_ribbon(aes(x = c(-10,-9),ymin = -c(10,9),ymax = -c(8,7), fill = "Age Models (95% HDR)" ),alpha = 0) + #for legend
  geom_ribbon(aes(x = c(-10,-9),ymin = -c(10,9),ymax = -c(8,7), fill = "Age distributions" ),alpha = 0) + #for legend
  geom_ribbon(aes(x = c(-10,-9),ymin = -c(10,9),ymax = -c(8,7), fill = "Radiocarbon ages" ),alpha = 0.4) + #for legend
  geom_ribbon(aes(x = c(-10,-9),ymin = -c(10,9),ymax = -c(8,7), fill = "Historical constraints" ),alpha = 0.8) +#for legend
  geom_ribbon(aes(x = c(-10,-9),ymin = -c(10,9),ymax = -c(8,7), fill = "Eklutna Marker Layers" ),alpha = 0.5) +#for legend
  
  scale_fill_manual(name = NULL,values = c("white","#E41A1C","skyblue","white","coral","darkblue","darkgreen"),
                    breaks = c("Age Models (95% HDR)", "Full ensemble","Historically constrained","Age distributions","Radiocarbon ages","Historical constraints","Eklutna Marker Layers"),
                    label = c(expression(bold("Age Models (95% HDR)")),"Full ensemble","Historically constrained",expression(bold("Age distributions")),"Radiocarbon ages","Historical constraints","Eklutna Marker Layers")) +
  scale_color_manual(name = NULL,values = c("black"))

out1 <- plotModelDistributions(add.to.plot = start, L = SC,alp = 0.4,scale.frac = .025,color = "coral")

out2 <- plotModelDistributions(add.to.plot = out1, L = HC,alp = 0.8,scale.frac = .2,color = "darkblue")

out <- plotModelDistributions(add.to.plot = out2,L = L,alp = 0.5,scale.frac = .10,color = "darkgreen") +
  geom_text(data = NULL,aes(x = emlagesinskilak + 50,y = emldepthsinskilak -100, label = eml), color = "darkgreen") +
  geom_text(data = NULL,aes(x = hist.data$mid.age,y = hist.data$depth + 115, label = hist.names), color = "darkblue",size = 3) +
  
  xlab("Age (BP)") + 
  #coord_cartesian(ylim = c(0,1304),expand = FALSE) +
  ggtitle("Skilak Age Model") +
  scale_y_reverse("Depth (cm)") +
  coord_cartesian(ylim = c(1304,0),expand = FALSE) +
  theme(legend.position = "inside",legend.position.inside = c(0.2,0.3))

out 



ggsave(out,filename = "Skilak Age model figure - updated v2.pdf")  
  





for(i in 1:length(allMarkerLayers)){
  mlName <- allMarkerLayers[i]
  mlAges <- filter(allMlAges,name == mlName)
  
  skilakPrior <- filter(allMlPrior,name == mlName) |> 
    select(age = SK)
  
  eklutnaPrior <- filter(ekPrior,ML == mlName)
  
  mlHist <- ggplot(mlAges) + 
    geom_histogram(aes(age,y = after_stat(density),fill = "Posterior")) + 
    geom_density(data = skilakPrior,aes(age,color = "Skilak"),fill = NA) +
    geom_density(data = eklutnaPrior,aes(age,color = "Eklutna"),fill = NA) +
    
    # geom_vline(data = mlP,aes(xintercept = Year,color = Lake)) +
    xlab("Age (yr BP)") + 
    ggtitle(paste("Posterior age distribution for",mlName)) + 
    scale_color_brewer("Prior",palette = "Set1") +
    scale_fill_manual(NULL,values = "grey30") +
    theme_bw() +
    theme(legend.position = "inside",legend.position.inside = c(.8,.8))
  
  
  ggsave(file.path("SkilakEklutna/figures/",paste0(mlName,".pdf")))
  
}
