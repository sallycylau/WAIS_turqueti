########################
#This is a script published by Bagley et al. (2016) Molecular Ecology (https://doi.org/10.1111/mec.13972), 
#adapted by Sally Lau for maximsiing likelihood distributions for tested trans-Antarctic seaway 4D-SFS models 


###Original script description by Bagley et al. (2016)##

#  Script to analyse the results of fastsimcoal 2.
#  Vitor C. Sousa
#  18/08/2016	
#  used in "History, geography, and host use shape genome-wide patterns of genetic variation in the redheaded pine sawfly (Neodiprion lecontei)"
#  Robin K. Bagley, Vitor C. Sousa, Matthew L. Niemiller, and Catherine R. Linnen

# FILES required:
# - File with the multidimensional observed SFS (*_MSFS.obs) - check folder FilesRequired for the observed joint SFS
# - Files outputted from fastsimcoal, for the best run (out of 100 runs) - check folder FilesRequired for files from run of best model
# - utilFunction3DSFS.r - R script with functions used in this script
# - ParFileInterpreter-VS1.r - R script with functions used to plot a representation of the inferred demographic model

# OPTIONS to change accordingly to the analyses:
# - pop : change this to use as a prefix to all the models (e.g. "ss6")
# - model : change this according to the name of the model (e.g. "_12_0_mig_41")
# - pop.names : vector of strings with the population names
# - obsfilename : file name of multidimensional observed SFS (3D-MSFS) with all-SNPS supposed to be in folder ./ObsSFS

# Load functions
source("utilFunctions4DSFS.r")
source("ParFileInterpreter-VS1.r")

####################
##
## Parameters and function calling (main)
##
##
####################

pop <- "1_WS_AS_RS_EA_" 
model <- "psc_conflow"
pop.names <- c("WS","AS","RS","EA")
obsfilename <- "1_WS_AS_RS_EA_psc_conflow_DSFS.obs" # supposed to be in folder ./ObsSFS
coption <- 10 # -C x option of fastsimcoal2 to pool all SFS entries with less than x SNPs

setwd("./FilesRequired")
  
  # Call script to make the plot of the model and re-scale parameters
  parFileInterpreter(paste(pop, model, "_maxL", sep=""),pop.names)
  
  # Run 100 simulations to get the expected SFS
  system(paste("cp ./ObsSFS/", obsfilename, " .; mv ", obsfilename, " ", pop, model, "_maxL_DSFS.obs;",     
               "./fsc26 -i ", pop, model, "_maxL.par -n 500000 -d -q --multiSFS -c4 -B4 -R100 --logprecision 18;
               rm ", pop, model, "_maxL_DSFS.obs;", sep=""))

  # Read the lhoods and compute the expected SFS
  nsim <- 1
  lhoods <- matrix(, ncol=100, nrow=nsim)
  runFolder <- paste("./", pop, model, "_maxL", sep="")
  fileLhoods <- paste("/",pop, model, "_maxL.lhoods", sep="")  
  lhoods <- scan(paste(runFolder, fileLhoods, sep=""), skip=1)
  # Plot the distribution of the likelihoods for the different replicates
  par(mfrow=c(1,1))
  hist(lhoods, prob=T, main=paste("Lhood"), xlab=paste("EstLhood"))  
  
  # Recompute the likelihoods 
  # read the obs SFS for the full dataset  
  obsSFS <- scan(paste("./ObsSFS/", obsfilename,sep=""),skip=2)
  pop.sizes <- scan(paste("./ObsSFS/", obsfilename,sep=""),skip=1, nlines = 1)[-1]
  numpops <- length(pop.sizes)
  # Get the MaxObsLhood including monomorphic sites
  relobs <- obsSFS/sum(obsSFS)
  monsites <- obsSFS[1]
  polsites <- sum(obsSFS[-1])
  totsites <- monsites + polsites
  maxObsLhoodnomon <- computelhood(obsSFS, relobs, coption)
  maxObsLhoodincmon <- maxObsLhoodnomon+(polsites*log10(polsites))+(monsites*log10(monsites))-(totsites*log10(totsites))
  
  # Read the Expected SFS
  runFolder <- paste("./", pop, model, "_maxL", sep="")
  fileExpSFS <- paste("/",pop, model, "_maxL_DSFS.txt", sep="")
  tmpexpSFS <- matrix(scan(paste(runFolder,fileExpSFS, sep=""), skip=2), nrow=100, byrow=T)
  
  # compute the likelihood for each replicate
  lhoodnomon <- apply(tmpexpSFS, 1, function(row) {computelhood(obsSFS, row, coption)})
  
  # Get the replicate with maxlhood to compare that with the obs SFS
  bestrep <- which(lhoodnomon==max(lhoodnomon))
  
  # Expected lhood of the run with the maximum likelihood
  expSFS <- tmpexpSFS[bestrep,]
  # Mean expected SFS
  mean.expSFS <- colMeans(tmpexpSFS)  
  # Range of values for each entry of the multidimension expected SFS
  expSFS_low95 <- apply(tmpexpSFS, 2, function(x) min(x))
  expSFS_upper95 <- apply(tmpexpSFS, 2, function(x) max(x))
  
  # Plot the 1D Marginal SFS
  pdf(file=paste("1DmarginalSFS", pop, model,".pdf",sep=""), width=16, height=8)
  plotfit1DSFS(obsSFS,expSFS,pop.sizes,numpops,pop.names,ylims = c(1.5,4),expSFS_low95,expSFS_upper95, model)
  plotfit1DSFS(obsSFS,mean.expSFS,pop.sizes,numpops,pop.names,ylims = c(1.5,4),expSFS_low95,expSFS_upper95, model)
  dev.off()
      
  
  # Assess the fit of the multidimensional SFS
  # Read and compute observed SFS, discarding entry for fixed ancestral and fixed derived
  obs.SFS<-obsSFS
  obs.SFS[c(1,length(obs.SFS))]<-0  
  
  # Read and compute expected SFS, by multiplying by number of sites  
  # get the mean exp sfs
  rel.exp.sfs <- mean.expSFS 
  # get the lower and upper range of values in the SFS
  low.exp.sfs <- apply(tmpexpSFS,2,function(x) {quantile(x,0.001)})
  upper.exp.sfs <- apply(tmpexpSFS,2,function(x) {quantile(x,0.999)})
  # get the expected SFS in the same scale as Observed SFS (i.e counts of sites)
  exp.sfs<-rel.exp.sfs*sum(obs.SFS)
  low.exp.sfs <- low.exp.sfs*sum(obs.SFS)
  upper.exp.sfs <- upper.exp.sfs*sum(obs.SFS)
  
  # get the labels with the coordinates of each entry in the SFS
  matrix.sfs <- get_sfscoordinates(numpops, pop.sizes)  
  
  # compute and plot the relative difference between the SFS
  diff.sfs<-((obs.SFS-exp.sfs)/obs.SFS)
  par(mfrow=c(2,1))
  plot(((obs.SFS-exp.sfs)/obs.SFS), xaxt='n', type="l", col="blue", main=model)
  axis(1, at=seq(1,length(diff.sfs), by=100), labels=matrix.sfs[seq(1,length(diff.sfs), by=100)], las=2, cex.axis=0.7)
  
  # compute the likelihood based on the expected sfs, but do this for each entry  
  print(paste("Sum of entries of expected SFS==0 is ", sum(rel.exp.sfs==0)))
  numzeroentries_expsfs <- sum(rel.exp.sfs==0)
  eval <- which(obs.SFS>coption & rel.exp.sfs>0)
  exp.log.lik <- log10(rel.exp.sfs[eval])*obs.SFS[eval]  
  explhoodnomon <- sum(exp.log.lik, na.rm=T)    
  rel.obs.sfs <- obs.SFS/sum(obs.SFS)
  numzeroentries_obssfs <- sum(rel.obs.sfs==0)
  print(paste("Sum of entries of observed SFS==0 is ", numzeroentries_obssfs))
  obs.log.lik <- log10(rel.obs.sfs[eval])*obs.SFS[eval]
  obslhoodnomon <- sum(obs.log.lik, na.rm=T)
  print(sum(exp.log.lik-obs.log.lik))
  
  # compute the difference in likelihood between expected and obs for each entry
  diff.log.lik <- exp.log.lik-obs.log.lik  
  plot(diff.log.lik, main=model)
  
  
  pdf(file=paste("FitJoint4D_SFS_orderentries_", model ,".pdf",sep=""), width=12, height=6)
  # Plot all the SFS that matters to fit (only parts above minentry)
  par(mfrow=c(1,1))
  indexorder <- order(obs.SFS[eval], decreasing = TRUE)
  plot(log10(obs.SFS[eval[indexorder]]),type="l", ylim=c(0,3.5),lwd=2, xaxt='n',xlab="",ylab="log10(#SNPs)")
  lines(log10(exp.sfs[eval[indexorder]]),col=4)
  axis(1, at=seq(1,length(indexorder), by=1), labels=matrix.sfs[eval[indexorder]], las=2, cex.axis=0.7)
  dev.off()
  
  threshold <- as.numeric(quantile(abs(diff.log.lik), (length(diff.log.lik)-30)/length(diff.log.lik)))
  pdf(file=paste("barplots_outlier_entries", model ,"_threshold", threshold,".pdf", sep=""), width=8, height=10)
  par(mfrow=c(4,1), mar=c(10,5,2,1), las=2)
  # Get entries with bad SFS fit 
  # defined as the 30 entries with the largest difference between the likelihood 
  # and the max possible likelihood (i.e. if we had a perfect fit)
  threshold <- as.numeric(quantile(abs(diff.log.lik), (length(diff.log.lik)-30)/length(diff.log.lik)))
  # get the entries with bias
  entries.SFS <- matrix.sfs[eval]
  # plot the diff likelihood for those entries, based on relative SFS and on the observed counts
  rel_expSFS <- exp.sfs[match(entries.SFS, matrix.sfs)]/sum(exp.sfs)
  rel_obsSFS <- obs.SFS[match(entries.SFS, matrix.sfs)]/sum(obs.SFS)
  exp.log.lik <- log10(rel_expSFS)*obs.SFS[match(entries.SFS, matrix.sfs)]
  obs.log.lik <- log10(rel_obsSFS)*obs.SFS[match(entries.SFS, matrix.sfs)]
  diff.log.lik <- exp.log.lik-obs.log.lik  
  barplot(diff.log.lik, names.arg=entries.SFS, ylab="diff. lhood. exp - obs")  
  # plot the barplots comparing counts for the expected and observed SFS (note that the expected SFS was multiplied by the number of SNPs)
  selentries <- match(entries.SFS, matrix.sfs)
  res <- matrix(c(exp.sfs[selentries], obs.SFS[selentries]), ncol=length(entries.SFS), byrow=T)        
  barx <- barplot(res, beside=T, legend.text=c("exp","obs"), names.arg=entries.SFS, ylim=c(0,max(res)*1.5), angle=90, main=model)
  # add the 0.05 and 0.95 quantiles to the expected distribution
  #   par(mfrow=c(1,1))
  #   plot(c(1:30), c(-30:-1), ylim=c(0,6e5))
  #   error.bar(barx[1,], exp.sfs[selentries], upper=upper.exp.sfs[selentries], lower=low.exp.sfs[selentries], length=0.001)
  segments(barx[1,],upper.exp.sfs[selentries], barx[1,], low.exp.sfs[selentries])  
  # relative bias in the entries
  relres <- apply(res, 2, function(col){col[1]/col[2]})  
  barx <- barplot(relres, names.arg=entries.SFS, ylim=c(0,max(relres)*1.15), angle=90, main="Relative fit", ylab="Relative fit=exp/obs")
  abline(h=1)
  segments(barx,upper.exp.sfs[selentries]/res[2,], barx, low.exp.sfs[selentries]/res[2,])
  dev.off()
  
  # save the recunputed likelihood of the 100 reps (maximised likelihood distributions)
  write.table(lhoodnomon, paste(pop, model, "lhoodnomon_maximised.lhoods", sep=""), col.names=F, row.names=F, quote=F)
setwd("../")

