# COMPUTELHOOD
# computes the loglikelihood given an obs.sfs and exp.sfs
# INPUT
#   obs.sfs : numeric vector with the observed SFS
#   exp.sfs : numeric vector with the expected SFS
#   min.entry : numeric value with the minimum entry for the obs SFS. All entries i where obs.sfs[i]<=min.entry are pooled together.
# RETURN
#   log likelihood computed as SUM m_i log10(p_i), 
#   where m_i is the ith entry of obs.sfs and p_i is the ith entry of exp.sfs
# NOTE: for entries where m_i > 0 and p_i=0, we replace p_i by a small value (penalty)
computelhood <- function(obs.sfs, exp.sfs, min.entry) {
  
  lhood <- 0
  
  # remove the first and last entries
  obs.sfs <- obs.sfs[-c(1,length(obs.sfs))]
  exp.sfs <- exp.sfs[-c(1,length(exp.sfs))]
  
  # Get the valid entries, i.e. entries where obs.SFS > 0 
  eval <- which(obs.sfs > 0)
  
  # Calculate expected SFS with the penaltie for entries where obs.SFS > 0 and exp.SFS == 0
  if(sum(exp.sfs[eval]==0) > 0) {
    # Settings (penalty for exp SFS entries with zero)
    penalty <- 1e-10
    minpenalty <- 1e-8
    
    penalty <- min(exp.sfs[exp.sfs>0])/100
    if(penalty > minpenalty) {
      penalty <- minpenalty
    } 
    
    # Get the entries which are zero at the obs SFS to have the penalty
    tmp.exp <- exp.sfs[eval]  # note that the length of tmp.exp is length(eval)
    tmp.exp[tmp.exp==0] <- penalty 
    exp.sfs[eval] <- tmp.exp
  }
  
  # ensure that the sum of exp.sfs is 1
  exp.sfs <- exp.sfs/sum(exp.sfs)
  
  # compute the likelihood
  if(sum(exp.sfs[eval]==0) > 0) {
    stop("ERROR: still entries with exp.sfs==0!!!!")
  } else {

      # Pool together all the entries <= min.entry
      indexLargerMinEntry <- which(obs.sfs[eval] > min.entry)
      indexSmallerMinEntry <- which(obs.sfs[eval] <= min.entry)
      if(length(indexLargerMinEntry)>0) {
        lhood <- sum(obs.sfs[eval[indexLargerMinEntry]]*log10(exp.sfs[eval[indexLargerMinEntry]]))  
      } else {
        stop("ERROR: No entries larger than min.entry (-C x option)")
      }
      if(length(indexSmallerMinEntry)>0) {
        lhood <- lhood+sum(obs.sfs[eval[indexSmallerMinEntry]])*log10(sum(exp.sfs[eval[indexSmallerMinEntry]]))  
      } 
  }
  
  lhood
}




# GET_SFSCOORDINATES
# get the coordinates for the x axis of the plots, where 0 is coded as ancestral, 1 is het derived and 2 is homozygous derived
# INPUT:
#   popsinmodel : number of populations in the model
#   pop.sizes   : sample sizes in number of gene copies for each population
# OUTPUT:
#   matrix.sfs : vector of strings with size prod(pop.sizes+1) with the coordinates for each entry of the multidimensional SFS
# GET_SFSCOORDINATES
# get the coordinates for the x axis of the plots, where 0 is coded as ancestral, 1 is het derived and 2 is homozygous derived
# INPUT:
#   popsinmodel : number of populations in the model
#   pop.sizes   : sample sizes in number of gene copies for each population
# OUTPUT:
#   matrix.sfs : vector of strings with size prod(pop.sizes+1) with the coordinates for each entry of the multidimensional SFS
get_sfscoordinates <- function(popsinmodel, pop.sizes) {
  
  if(popsinmodel == 4) {
    
    matrix.sfs <-c(rep(0, prod(pop.sizes+1)))
    k<-1
    
    for (i in 0:pop.sizes[1]) {
      
      for (j in 0:pop.sizes[2]) {
        
        for (r in 0:pop.sizes[3]) {
          
          for (p in 0:pop.sizes[4]){
            matrix.sfs[k]<-c(paste(i,j,r,p, sep=","))
            k<-k+1
          }
        }
      }  
    }        
  }
  matrix.sfs
}

# PLOTFIT1DSFS
# Plot 1D SFS as a barplot, or as a line plot. Prints several plots.
# INPUT:
#   obsSFS : vector with the multidimensional observed SFS (counts of SNPs for each entry)
#   expSFS : vector with the multidimensional expected SFS (relative probabilities for each entry)
#   pop.sizes : vector with the sample size for each population as 2*number diploid individuals
#   numpops : number of populations
#   pop.names : vector with the populaiton names
#   ylims : vector with two entries corresponding to the minimum and maximum values of the y axis
#   expSFS_low95 : vector with the multidimensional expected SFS (minimum values or lower 0.025% quantile for each entry)
#   expSFS_upper95 : vector with the multidimensional expected SFS (maximum values or lower 0.025% quantile for each entry)
#   models : string with the model tag
plotfit1DSFS <- function(obsSFS,expSFS,pop.sizes,numpops,pop.names,ylims,expSFS_low95,expSFS_upper95,models) {
  
  # discard the monomorphic sites (assess the fit only to polymorphic sites)
  obs.SFS <- obsSFS
  obs.SFS[c(1,length(obs.SFS))] <- 0
  # expected SFS in the same scale as the observed SFS (i.e. number of sites)
  exp.SFS <- expSFS*sum(obs.SFS)
  exp.SFS[c(1,length(obs.SFS))] <- 0
  
  # transform the vectors into multidimensional arrays (3D-SFS)
  sfs.size <- pop.sizes+1  
  dim(obs.SFS)=sfs.size[c(numpops:1)]
  obs.SFS=aperm(obs.SFS, c(numpops:1))
  dim(exp.SFS)=sfs.size[c(numpops:1)]
  exp.SFS=aperm(exp.SFS, c(numpops:1))
  
  # get the 1 D marginals for each pop
  margobs <- list()
  margexp <- list()
  for(i in 1:numpops) {
    margobs[[i]] <- apply(obs.SFS,i,sum)
    str(margobs[[i]])
    
    margexp[[i]] <- apply(exp.SFS,i,sum)
    str(margexp[[i]])
  }
  
  # get the 1D marginal based on the 0.025 quantile
  exp.SFS_low95 <- expSFS_low95*sum(obs.SFS)
  exp.SFS_low95[c(1,length(obs.SFS))] <- 0
  dim(exp.SFS_low95)=sfs.size[c(numpops:1)]
  exp.SFS_low95=aperm(exp.SFS_low95, c(numpops:1))
  margexp_low95 <- list()
  for(i in 1:numpops) {
    margexp_low95[[i]] <- apply(exp.SFS_low95,i,sum)
    str(margexp_low95[[i]])
  }
  
  # get the 1D marginal based on the 0.975 quantile
  exp.SFS_upper95 <- expSFS_upper95*sum(obs.SFS) 
  exp.SFS_upper95[c(1,length(obs.SFS))] <- 0
  dim(exp.SFS_upper95)=sfs.size[c(numpops:1)]
  exp.SFS_upper95=aperm(exp.SFS_upper95, c(numpops:1))
  margexp_upper95 <- list()
  for(i in 1:numpops) {
    margexp_upper95[[i]] <- apply(exp.SFS_upper95,i,sum)
    str(margexp_upper95[[i]])
  }
  
  par(mfrow=c(1,4), oma=c(2,2,5,2))
  for(i in 1:numpops) {
    res <- matrix(c(margobs[[i]], margexp[[i]]), nrow=2, byrow=T)
    barx <- barplot(res, beside=T, legend.text=c("obs","exp"), names.arg=c(0:pop.sizes[i]), main=pop.names[i], ylab="#sites", xlab="minor allele frequency")
    segments(barx[2,],margexp_upper95[[i]], barx[2,], margexp_low95[[i]])
  }
  title(main=models, outer=T)
  
  par(mfrow=c(1,4), oma=c(2,2,5,2))
  for(i in 1:numpops) {
    res <- matrix(c(log10(margobs[[i]]), log10(margexp[[i]])), nrow=2, byrow=T)
    barplot(res, beside=T, legend.text=c("obs","exp"), names.arg=c(0:pop.sizes[i]), main = pop.names[i], ylim=c(-1,4.2), ylab="log10(#sites)", xlab="minor allele frequency")
    segments(barx[2,],log10(margexp_upper95[[i]]), barx[2,], log10(margexp_low95[[i]]))
  }
  title(main=models, outer=T)
  
  par(mfrow=c(1,4), oma=c(2,2,5,2))
  for(i in 1:numpops) {
    plot(0:pop.sizes[i], log10(margobs[[i]]), type="s", ylim=c(1.8,4.0), pch=2, lwd=2, main=pop.names[i], ylab="log10(#sites)", xlab="minor allele frequency")
    lines(0:pop.sizes[i], log10(margexp[[i]]), type="s", pch=3, col=3, lty=1, lwd=2)
    legend("topright", c("obs","exp"), col=c(1,4), lty=c(1,1), lwd=2)    
  }
  title(main=models, outer=T)
  
  par(mfrow=c(1,4), oma=c(2,2,5,2))
  for(i in 1:numpops) {
    plot(0:pop.sizes[i], margobs[[i]], type="s", pch=2, lwd=2, main=pop.names[i], ylab="#sites", xlab="minor allele frequency")
    lines(0:pop.sizes[i], margexp[[i]], type="s", pch=3, col=3, lty=1, lwd=2)
    legend("topright", c("obs","exp"), col=c(1,4), lty=c(1,1), lwd=2)
  }
  title(main=models, outer=T)
  
}









#This function creates a color scale for use with the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
image.scale <- function(z, zlim, col = rainbow(12), breaks, horiz=TRUE, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){ylim<-c(0,1); xlim<-range(breaks)}
  if(!horiz){ylim<-range(breaks); xlim<-c(0,1)}
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)  
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}

