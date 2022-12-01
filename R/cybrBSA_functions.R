
#Load packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(cowplot)
library(dplyr)

library(foreach)
library(doParallel)

ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                       "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                       CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                 "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                 "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                 "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))



################################################################################
## Functions


cybrInputGATKTable <- function(rawData, yeast = TRUE){

  require(dplyr)
  require(doParallel)
  require(foreach)

  HNGLCDRXY <- read.table(rawData, header = TRUE)

  #Identify the unique values besides AD/DP/GQ/PL
  gsub(".AD", "",
       gsub(".GQ", "",
            gsub(".DP","",
                 gsub(".PL","",
                      colnames(select(HNGLCDRXY, -CHROM, -POS, -REF, -ALT)))))) %>% unique() -> Samples
  #i <- Samples[1]

  resultscdf <- foreach(i=Samples,.combine=rbind) %dopar% {
    mydf <- HNGLCDRXY %>% select(CHROM, POS, REF, ALT) %>% mutate(Dataset = i)
    AD <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("AD"))
    GQ <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("GQ"))
    DP <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("DP"))
    PL <- select(HNGLCDRXY,matches(c(i), ignore.case = FALSE)) %>% select(., contains("PL"))
    cbind(mydf, AD , GQ , DP, PL) -> mydftotal
    colnames(mydftotal) <- c(colnames(mydf), "AD", "GQ", "DP", "PL")

    mydftotal %>% separate(AD, c('AD.REF','AD.ALT'), extra='drop') %>%
      separate(PL, c('PL.REF','PL.ALT'), extra='drop') -> mycdf

    mycdf
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                         "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                         CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                   "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                   "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                   "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    resultscdf %>% left_join(.,ChromKey) %>% select(-CHROM) %>% mutate(CHROM = chromosomes) %>% select(-chromosomes) -> results
  }else{
    results <- resultscdf
  }
  return(results)

}


################################################################################
### Filter by quality

cybrQualityFilter <- function(gatkdf, GQcutoff = 98, cleandata = TRUE){

  #Filter by quality
  gatkdf %>% filter(GQ > GQcutoff) %>%
    select(-DP, -GQ, -PL.ALT, -PL.REF) %>%
    pivot_longer(c(AD.REF, AD.ALT), names_to = "AltRef_Allele", values_to = "ReadCount") %>%
    na.omit() -> filteredgatkdf

  if(cleandata == TRUE){
    #REMOVE POSITIONS WHERE THERE ISN'T ALL BULKS
filteredgatkdf %>% group_by(CHROM, POS) %>%
  summarize(uniqueDatasets = length(unique(Dataset))) -> testingforglm

filteredgatkdf %>% left_join(testingforglm, by = c("CHROM", "POS")) %>%
  filter(uniqueDatasets == max(uniqueDatasets)) %>% select(-uniqueDatasets) -> cleanedgatkdf
}else{
  cleanedgatkdf <- filteredgatkdf
}

return(cleanedgatkdf)
}

################################################################################
### Convert Parental VCFs to Data Frame

cybrConvertParentalAlleles <- function(ParentFiles = c("Wine_VCF.txt", "Oak_VCF.txt"),
                                       Parents = gsub("_VCF.txt","", ParentFiles), Truncate = TRUE, yeast = TRUE){
  temparent <- list()
  mergeparents <- foreach(i=1:length(ParentFiles), .combine=rbind) %dopar% {
    read.table(ParentFiles[i], header = TRUE) %>% mutate(parent = Parents[i])
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    rbind(mergeparents) %>% arrange(CHROM, POS) %>%
      select(CHROM, POS, REF, ALT, parent) %>%
      merge(ChromKey) %>% select(-CHROM) %>%
      mutate(CHROM = chromosomes) %>% select(-chromosomes) -> ParentalVCF

  }else{
    rbind(mergeparents) %>% arrange(CHROM, POS) %>%
      select(CHROM, POS, REF, ALT, parent) -> ParentalVCF
  }

  ParentalVCF %>% pivot_wider(names_from = parent, values_from = ALT) -> SNPids

  SNPids$Type <- 0
  for(i in Parents){

    #filter rows in which all values of columns of the parent NOT selected are NA
    select(SNPids,-i, -CHROM, -POS, -REF) -> tempdf
    tempdf$Any_NA <- apply(tempdf, 1, function(x) anyNA(x))
    SNPids$Type[which(tempdf$Any_NA)] <- i
    rm(tempdf)
  }


  #Collect it to output
  if(Truncate == TRUE){
    SNPids %>% select(CHROM, POS,  Type) %>% filter(Type != 0) -> SNPids
  }

  return(SNPids)

}


################################################################################
### Combine Parental and Experimental Variants

cybrIDAlleles <- function(BSAdfstart = finaldf, Parentdf = test, yeast = TRUE){

  Parentdf %>% na.omit()
  BSAdf <- left_join(BSAdfstart, Parentdf)
  BSAdf$PAllele <- NA

  Parents <- unique(Parentdf$Type)
  if(length(Parents) == 2){
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[1]] <- Parents[1]
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == Parents[2]] <- Parents[2]

    #Run if only two-parent cross
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[1]] <- Parents[2]
    BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.ALT" & BSAdf$Type == Parents[2]] <- Parents[1]

  }else{
    for(i in Parents){
      BSAdf$PAllele[BSAdf$AltRef_Allele == "AD.REF" & BSAdf$Type == i] <- i
    }
  }

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))
    BSAdf$CHROM <- factor(BSAdf$CHROM, levels = ChromKey$chromosomes)
  }


  #Convert all to factors for glm
  BSAdf$Dataset <- factor(BSAdf$Dataset)
  BSAdf$AltRef_Allele <- factor(BSAdf$AltRef_Allele)
  BSAdf$PAllele <- factor(BSAdf$PAllele)

  return(BSAdf)
}





#### Reformat data so that it has bulk etc included

################################################################################
### Run GLM Function

cybrBSA_GLM <-  function(lrP, chr = "II", windowsize = 100, formula = "PAllele ~ Bulk*Parent",
                         resultscol = c("Intercept", "Bulk", "Parent", "Interaction")){

  lrP <- subset(lrP, CHROM == chr)

  AllResults <- list()
  for(c in unique(lrP$CHROM)){
    lrP <- subset(lrP, CHROM == c)
    #Run the glm on that chromosome
    Results <- foreach (i=unique(lrP$POS), .combine=rbind) %dopar% {
      res <- suppressWarnings(glm(as.formula(formula), weights = ReadCount, family = binomial, data = lrP[lrP$POS == i,]))

      #Output of foreach automatically binds rows of what is printed
      c(c, i, summary(res)$coefficients[((length(summary(res)$coefficients)/2)+1):(length(summary(res)$coefficients) - length(summary(res)$coefficients)/4)])
    }

    resnames <- suppressWarnings(glm(as.formula(formula), weights = ReadCount, family = binomial, data = lrP[lrP$POS == unique(lrP$POS)[1],]))

    #Format the results for the next function
    Results <- as.data.frame(Results)
    colnames(Results) <- c("CHROM", "POS", names(resnames$coefficients))
    #colnames(Results) <- c("CHROM", "POS", resultscol)

    for(i in 2:length(colnames(Results))){
      Results[,i] <- as.numeric(Results[,i])
    }
    Results %>% arrange(POS) -> Results
  }
  return(Results)
}

################################################################################
### New Window Script Function

cybrSmoothBSAWindows <- function(Results, windowsize = 100, chr = unique(Results$CHROM)[1]){

  #extract parameters
  Results %>% select(-CHROM, -POS) %>% names() -> params

  #arrange by position
  Results %>% filter(CHROM == chr) %>% arrange(POS) -> Results

  #run window
  WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {

    smoothedparams <- vector()
    for(p in 1:length(params)){
      smoothedparams[p] <- mean(Results[[params[p]]][(i-windowsize):(i+windowsize)])
    }

    #print CHROM, index, POS, and mean
    unlist(c(chr, i, Results[i,2:(2+length(params))],smoothedparams))
  }

  #rename columns
  WResult <- as.data.frame(WResult)

  colnames(WResult) <- c("CHROM", "Index", "POS",
                         paste(params, "Z", sep = "_"),
                         paste(params, "Zprime", sep = "_"))

  #convert to numeric
  for(i in 2:length(colnames(WResult))){
    WResult[,i] <- as.numeric(WResult[,i])
  }

  return(WResult)
}


cybrSmoothBSAWindows_b <- function(Results, windowsize = 100, chr = unique(Results$CHROM)[1]){

  #extract parameters
  Results %>% select(-CHROM, -POS) %>% names() -> params

  #arrange by position
  Results %>% filter(CHROM == chr) %>% arrange(POS) -> Results

  Results %>% summarise(SNPs = length(unique(POS))) %>%
    distinct() %>%
    mutate(maxW = floor(SNPs/2)) -> tablecounts

  if(windowsize < tablecounts$maxW){
    #run window
    WResult <- foreach(i=windowsize+1:(length(Results$POS) - (2*windowsize)), .combine=rbind) %dopar% {

      smoothedparams <- vector()
      for(p in 1:length(params)){
        smoothedparams[p] <- mean(Results[[params[p]]][(i-windowsize):(i+windowsize)])
      }

      #print CHROM, index, POS, and mean
      unlist(c(chr, i, Results[i,2:(2+length(params))],smoothedparams))
    }

    #rename columns
    WResult <- as.data.frame(WResult)
    colnames(WResult) <- c("CHROM", "Index", "POS",
                           paste(params, "Z", sep = "_"),
                           paste(params, "Zprime", sep = "_"))

    #convert to numeric
    for(i in 2:length(colnames(WResult))){
      WResult[,i] <- as.numeric(WResult[,i])
    }
    return(WResult)
  }else{
    warning(print(paste("Window size of chromosome ", chr, "are larger than number of data points")),
            call. = TRUE, immediate. = FALSE, noBreaks. = FALSE,
            domain = NULL)
  }
}

################################################################################
## Make function for plotting

cybrPlotZPrime <- function(zprimedf,
                           columns = c("Bulk_Zprime", "Parent_Zprime", "Interaction_Zprime"),
                           chromosomes = "All",
                           title = "Smoothed Z Scores",
                           yeast = TRUE,
                           colvalues = c("#345F6F", "#D7335C", "#FFB05C")){

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    zprimedf$CHROM <- factor(zprimedf$CHROM, levels = ChromKey$chromosomes)
  }


  if(chromosomes == "All"){
    zprimedf %>% filter(CHROM != "M") %>%
      pivot_longer(columns, names_to = "Factor", values_to = "Zprime") %>%
      ggplot(aes(x = POS, y = Zprime, color = Factor)) +

      #Add Bulk, Parent, Day, and Interaction Traces
      geom_line(size = 0.75) +
      scale_color_manual(values = colvalues) +

      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zprime")+
      ggtitle(paste("Whole Genome ", title, sep = "")) -> finalplot

  }else{

    zprimedf %>% filter(CHROM %in% chromosomes) %>%
      pivot_longer(columns, names_to = "Factor", values_to = "Zprime") %>%
      ggplot(aes(x = POS, y = Zprime, color = Factor)) +

      #Add Bulk, Parent, Day, and Interaction Traces
      geom_line(size = 0.75) +
      scale_color_manual(values = colvalues) +

      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zprime")+
      ggtitle(title) -> finalplot
  }
  return(finalplot)
}

cybrPlotZScore <- function(zprimedf = CuSO4_wholegenomeBSA,
                           column = "Bulk_Z",
                           chromosomes = "All",
                           title = "Z Scores by Position",
                           yeast = TRUE,
                           dotcolor = "black"){

  if(yeast == TRUE){
    ChromKey <- data.frame(chromosomes = c("I", "II", "III", "IV", "V", "VI", "VII", "VIII",
                                           "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI", "M"),
                           CHROM = c("NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001136.10",
                                     "NC_001137.3", "NC_001138.5", "NC_001139.9", "NC_001140.6",
                                     "NC_001141.2", "NC_001142.9", "NC_001143.9", "NC_001144.5",
                                     "NC_001145.3", "NC_001146.8", "NC_001147.6", "NC_001148.4", "NC_001224.1"))

    zprimedf$CHROM <- factor(zprimedf$CHROM, levels = ChromKey$chromosomes)
  }


  if(chromosomes == "All"){

    zprimedf %>% pivot_longer(column, names_to = "Factor", values_to = "Zscore") %>%
      filter(CHROM != "M") %>% filter(Factor == column) %>%
      ggplot(aes(x = POS, y = Zscore)) +
      geom_point(size = 0.5, alpha = 0.08, color = dotcolor) +
      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zscore") +
      ggtitle(title) -> finalplot

  }else{

    zprimedf %>% pivot_longer(column, names_to = "Factor", values_to = "Zscore") %>%
      filter(CHROM %in% chromosomes) %>% filter(Factor == column) %>%
      ggplot(aes(x = POS, y = Zscore)) +
      geom_point(size = 1, alpha = 0.08, color = dotcolor) +
      #Add the breaks for the facets
      scale_x_continuous(breaks = seq(from = 0, to = max(zprimedf$POS), by = 1e5), name = "Genomic Position") +
      facet_grid(~CHROM, scales = "free_x", space = "free_x") + theme_minimal() +
      theme(legend.position = "bottom", axis.text.x=element_blank(),
            axis.ticks.x=element_blank()) +
      xlab("Genomic Position") + ylab("Zscore") +
      ggtitle(title) -> finalplot
  }
  return(finalplot)
}

################################################################################
# Check number of loci per chromosome for samples to get window sizes

checkWindowLim <- function(Dataset, includechr = TRUE, exceptchr = NULL){
  if(is.null(exceptchr) == FALSE){
    includechr = FALSE
  }
  if(includechr != TRUE){
    if(is.null(exceptchr)){
      Dataset %>% filter(CHROM %in% includechr) -> Dataset
    }else{
      Dataset %>% filter(CHROM %in% exceptchr == FALSE | CHROM %in% includechr) -> Dataset
    }
  }

  Dataset %>% group_by(CHROM) %>%
    summarise(CHROM = CHROM, SNPs = length(unique(POS))) %>%
    distinct() %>%
    mutate(maxW = floor(SNPs/2)) -> tablecounts

  return(tablecounts)
}
