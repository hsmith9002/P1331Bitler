#List of functions:
#getKM.permuted.LRT: returns a list of LRT stats for a "permutation" where you pull random genes only 
#getKM.permuted.LRT.v2: returns a list of LRT stats using a traditional permutation p-value with shuffled expression values and outcome measurements
#getKM.permuted.LRT.v3: returns a list of LRT stats using both a traditional permutation method by shuffling expression values and outcome measurements as well as pulling random expression values (not the genes in your candidate list)
#getPermutedPval: returns a p-value when inputting your true original statistic and your list of permuted statistics


getKM.permuted.LRT = function(expr, arrayListLength, nPerms, outcome, time){
  set.seed(12)
  
  #input needed for function:
    #expr: full array data (gene by subject), including the 97 candidate genes
    #arrayListLength: the 97 candidate genes are gene symbols.  This could correspond to many probes/gene              depending on platform.  Use the same amount of probes used in original cox PH  
    #nPerms: number of permutations you want to do (1,000)
    #outcome: vector of 0 or 1's corresponding to event (survival or disease recurrence or whatever dataset has)
    #time: numeric vector corresponding to time-to-event or time monitored if event didn't occur
    #note: the expr, outcome and time need to be sorted in same subject order (i.e. subject 1's gene expression        corresponds to column 1 in expr, it's outcome corresponds to the first event in outcome and the it's             time-to-event corresponds to the first number in time) 
  
  #STEP 1: get random gene expression data
  #select random genes;
    #make a matrix with ncol = # permutations and nrow = #genes.  The number corresponds to what row number to          select from original expression data
  randomProbes = matrix(nrow=arrayListLength, ncol=nPerms)
  colnames(randomProbes) =  paste("randomPull.", 1:nPerms, sep="")
  
  for(i in 1:nPerms){
   randomProbes[,i] = sample(c(1:nrow(expr)), arrayListLength, replace=FALSE)
  }
  
  #get expr dataset;
    #make a list.  Each object in list is an expression matrix (ncol = subjects and nrow = random genes) dataset       for a single permutation.   
  randomExpr = list(length(nPerms))
  for(i in 1:nPerms){
  randomExpr[[i]] = expr[randomProbes[,i], ]
  }

  #STEP2: Cox PH Regression by Gene;
    #a list of vectors corresponding to the coefficient for each gene in the univariate cox ph
  cox.byGene = lapply(randomExpr, function(aList) apply(aList, 1, function(a) coxph(Surv(time, outcome) ~ a)$coefficients)) 
  
  #STEP3: Get Sample Scores;
   #a list of weighted expression matrices for each permutation  
  expr.weighted = list(nPerms)
  for(i in 1:nPerms){
    expr.weighted[[i]] = apply(randomExpr[[i]], 2, function(a) a*cox.byGene[[i]])
  }
   #calculate sample scores for each subject resulting in a list of vectors 
  sample.scores = lapply(expr.weighted, colSums, na.rm=TRUE)
   #calculate the median sample score for each permutation
  median.scores = lapply(sample.scores, median)
  
   #determine which subjects have a high/low signature for each permutation
  expr.high = matrix(ncol=nPerms, nrow=floor(ncol(expr)/2))
  expr.low = matrix(ncol=nPerms, nrow=ceiling(ncol(expr)/2))
  
  for(i in 1:nPerms){
    expr.high[,i] = names(which(sample.scores[[i]] > median.scores[[i]]))
    expr.low[,i] = names(which(sample.scores[[i]] <= median.scores[[i]]))
  }
  
  ##STEP4: Get Kaplan-Meier Plots;
  bladder.signature = matrix(ncol=nPerms, nrow=ncol(expr))
  rownames(bladder.signature) = colnames(expr)
  colnames(bladder.signature) = paste("perm.", c(1:nPerms), sep="")
  for(c in 1:nPerms){
    for(r in 1:ncol(expr)){
      if(colnames(expr)[r] %in% as.matrix(expr.high[,c])){bladder.signature[r,c]="high"}else{bladder.signature[r,c]="low"}
    }
  }

  expr.cox = list(length(nPerms))
  for(i in 1:nPerms){
   expr.cox[[i]] <- coxph(Surv(time, outcome)~ as.factor(bladder.signature[,i]))
  }
 
  LRTstat = unlist(lapply(expr.cox, function(a) summary(a)$logtest[1]))
  return(LRTstat)

}

#This function will return a vector of LRT statistics (one for each permutation). 


##Method 2: Get the traditional permutation p-value with shuffled expression values and outcome measurements


getKM.permuted.LRT.v2 = function(expr, nPerms, outcome, time){
  set.seed(12)
  
  #get random assortment of outcome (event&time);
  randomOuts = matrix(nrow=ncol(expr), ncol=nPerms)
  colnames(randomOuts) =  paste("randomOut.", 1:nPerms, sep="")
  for(i in 1:nPerms){
   randomOuts[,i] = sample(c(1:ncol(expr)), ncol(expr), replace=FALSE)
  }
  
  randomTimes = matrix(nrow=ncol(expr), ncol=nPerms)
  randomEvents = matrix(nrow=ncol(expr), ncol=nPerms)
  for(i in 1:ncol(randomOuts)){
    randomTimes[,i] = time[randomOuts[,i]]
    randomEvents[,i] = outcome[randomOuts[,i]]
  }
  
    #STEP1: Cox PH Regression by Gene;
  cox.byGene = matrix(nrow=nrow(expr), ncol=nPerms)
  for(i in 1:nPerms){
    cox.byGene[,i] = apply(expr, 1, function(a) coxph(Surv(randomTimes[,i], randomEvents[,i]) ~ a)$coefficients)
  }
  
  #Step2: Get Sample Scores;
  #dim(expr)
  #dim(cox.byGene)
  
  expr.weighted = list(length=nPerms)
  for(i in 1:nPerms){
    expr.weighted[[i]] = apply(expr, 2, function(a) a*cox.byGene[,i])
  }
  sample.scores = lapply(expr.weighted, function(a) colSums(a,na.rm=TRUE ))
  median.scores = lapply(sample.scores, function(a) median(a,na.rm=TRUE))

  expr.high = matrix(ncol=nPerms, nrow=floor(ncol(expr)/2))
  expr.low = matrix(ncol=nPerms, nrow=ceiling(ncol(expr)/2))
  
  for(i in 1:nPerms){
    expr.high[,i] = names(which(sample.scores[[i]] > median.scores[[i]]))
    expr.low[,i] = names(which(sample.scores[[i]] <= median.scores[[i]]))
  }
  
  ##Step3: Get Kaplan-Meier Plots;
  bladder.signature = matrix(ncol=nPerms, nrow=ncol(expr))
  rownames(bladder.signature) = colnames(expr)
  colnames(bladder.signature) = paste("perm.", c(1:nPerms), sep="")
  for(c in 1:nPerms){
    for(r in 1:ncol(expr)){
      if(colnames(expr)[r] %in% as.matrix(expr.high[,c])){bladder.signature[r,c]="high"}else{bladder.signature[r,c]="low"}
    }
  }

  expr.cox = list(length(nPerms))
  for(i in 1:nPerms){
   expr.cox[[i]] <- coxph(Surv(randomTimes[,i], randomEvents[,i])~ as.factor(bladder.signature[,i]))
  }
 
  LRTstat = unlist(lapply(expr.cox, function(a) summary(a)$logtest[1]))
  ScoreStat = unlist(lapply(expr.cox, function(a) summary(a)$sctest[1]))
  WaldStat = unlist(lapply(expr.cox, function(a) summary(a)$waldtest[1]))
  return(list("lltest" =LRTstat, "wtest" = WaldStat))
}

##Method 3: Combine both selecing 97 random genes AND shuffle expression and outcome data;  
nPerms = 1000
getKM.permuted.LRT.v3 = function(expr, arrayListLength, nPerms, outcome, time){
  set.seed(12)
  #input needed for function:
    #expr: full array data (gene by subject), including the 97 candidate genes
    #arrayListLength: the 97 candidate genes are gene symbols.  This could correspond to many probes/gene              depending on platform.  Use the same amount of probes used in original cox PH  
    #nPerms: number of permutations you want to do (1,000)
    #outcome: vector of 0 or 1's corresponding to event (survival or disease recurrence or whatever dataset has)
    #time: numeric vector corresponding to time-to-event or time monitored if event didn't occur
    #note: the expr, outcome and time need to be sorted in same subject order (i.e. subject 1's gene expression        corresponds to column 1 in expr, it's outcome corresponds to the first event in outcome and the it's             time-to-event corresponds to the first number in time) 
  
  #STEP 1: get random gene expression data & shuffle outcomes
  #select random genes;
    #make a matrix with ncol = # permutations and nrow = #genes.  The number corresponds to what row number to          select from original expression data
  randomProbes = matrix(nrow=arrayListLength, ncol=nPerms)
  colnames(randomProbes) =  paste("randomPull.", 1:nPerms, sep="")
  
  for(i in 1:nPerms){
   randomProbes[,i] = sample(c(1:nrow(expr)), arrayListLength, replace=FALSE)
  }

  #get expr dataset;
    #make a list.  Each object in list is an expression matrix (ncol = subjects and nrow = random genes) dataset       for a single permutation.   
  randomExpr = list(length(nPerms))
  for(i in 1:nPerms){
  randomExpr[[i]] = expr[randomProbes[,i], ]
  }

#get random assortment of outcome (event&time);
  randomOuts = matrix(nrow=ncol(expr), ncol=nPerms)
  colnames(randomOuts) =  paste("randomOut.", 1:nPerms, sep="")
  for(i in 1:nPerms){
   randomOuts[,i] = sample(c(1:ncol(expr)), ncol(expr), replace=FALSE)
  }

  randomTimes = matrix(nrow=ncol(expr), ncol=nPerms)
  randomEvents = matrix(nrow=ncol(expr), ncol=nPerms)
  for(i in 1:ncol(randomOuts)){
    randomTimes[,i] = time[randomOuts[,i]]
    randomEvents[,i] = outcome[randomOuts[,i]]
  }
  
#So now we have a list of expression matrices with random genes as well as matrices for shuffled events and times
  
   #STEP2: Cox PH Regression by Gene;
    #a list of vectors corresponding to the coefficient for each gene in the univariate cox ph
  cox.byGene = matrix(nrow=arrayListLength, ncol=nPerms)
  for(i in 1:nPerms){
    cox.byGene[,i] = apply(randomExpr[[i]], 1, function(a) coxph(Surv(randomTimes[,i], randomEvents[,i]) ~ a)$coefficients)
  }

  #STEP3: Get Sample Scores;
   #a list of weighted expression matrices for each permutation  
  expr.weighted = list(nPerms)
  for(i in 1:nPerms){
    expr.weighted[[i]] = apply(randomExpr[[i]], 2, function(a) a*cox.byGene[[i]])
  }
   #calculate sample scores for each subject resulting in a list of vectors 
  sample.scores = lapply(expr.weighted, colSums, na.rm=TRUE)

   #calculate the median sample score for each permutation
  median.scores = lapply(sample.scores, median, na.rm=TRUE)

   #determine which subjects have a high/low signature for each permutation
  expr.high = list(length=nPerms)
  expr.low = list(length=nPerms)

  for(i in 1:nPerms){
    expr.high[[i]] = names(which(sample.scores[[i]] > median.scores[[i]]))
    expr.low[[i]] = names(which(sample.scores[[i]] <= median.scores[[i]]))
  }
  
  ##STEP4: Get Kaplan-Meier Plots;
  bladder.signature = matrix(ncol=nPerms, nrow=ncol(expr))
  rownames(bladder.signature) = colnames(expr)
  colnames(bladder.signature) = paste("perm.", c(1:nPerms), sep="")
  for(c in 1:nPerms){
    for(r in 1:ncol(expr)){
      if(colnames(expr)[r] %in% as.matrix(expr.high[[c]])){bladder.signature[r,c]="high"}else{bladder.signature[r,c]="low"}
    }
  }

  expr.cox = list(length(nPerms))
  for(i in 1:nPerms){
   expr.cox[[i]] <- coxph(Surv(randomTimes[,i], randomEvents[,i])~ as.factor(bladder.signature[,i]))
  }
 
  LRTstat = unlist(lapply(expr.cox, function(a) summary(a)$logtest[1]))
  ScoreStat = unlist(lapply(expr.cox, function(a) summary(a)$sctest[1]))
  WaldStat = unlist(lapply(expr.cox, function(a) summary(a)$waldtest[1]))
  return(list("lltest" =LRTstat, "wtest" = WaldStat))
}
  
getPermutedPval = function(origStat, permStats){
  pval = sum(origStat<permStats)/length(permStats)
  return(pval)
}
  
