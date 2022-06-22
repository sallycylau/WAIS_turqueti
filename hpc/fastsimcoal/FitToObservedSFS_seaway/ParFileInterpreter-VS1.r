#...........................................................................................
# (c) Laurent Excoffier June 2015
# 
# Attempt at reading and interpreting a par file

# 11.06.15  Added removal of trailing separators within file
#           Handled "keep" keyword
# 12.06.15
#           Fixed radius of max deme size, irrespective of the number of demes
#           Draw pop. size scale at fixed position, on bototm left of the graph
# 14.06.15  
#           Draw segments from top of circles
#           Implement growth rates
#           Draw triangle in growing or contracting populations
# 15.06.15
#           Take growth and resize into account to compute min and max pop sizes
#           Output growth triangles in legned  only if there is growth
# 
# 16.06.15  Added handling of command line:
#              Run: Rscript ParFileInterpreter-v3.r input.par
#              and input.par wil be used as input and interpreted by the program
#           Corrected bug in update of growth rate and migMat number if keep statement
# 
# 05.07.15  Added handling of rescaling factor
#              Output pdf file with rescaled pop sizes, div times, sampling times and migration rates
#              if rescale factor found in file *.rsf.then output _rescaled.pdf file
#              Run: Rscript ParFileInterpreter-v3.r <GenericName>
#                   and r script will look for the GenericName.par and the  GenericName.rsf files
#           Added output of modified par file with rescaled parameters
#              File name : GenericName_rescaled.par
#
# 17.11.15  Vitor Sousa - corrected a bug on the re-scaling of the growth rates. 
#			  Created the vector growthRatesInitial that is not updated when reading the 
#			  historical events.
#
#...........................................................................................


# INPUT
#   args : string with the name of the generic name of the par file

parFileInterpreter <- function(args, pop.names, gentime=6, printPDF=TRUE) {
  
#   args=commandArgs(TRUE)
  print(args)
  
  #Expects generic name on command line  
  if(length(args)) {
    genName=args[1]
    parFileName=paste(genName, ".par",sep="")
    rsfFileName=paste(genName, ".rsf",sep="")
  } else {  
    stop("Arguments for parFileInterpreter not found")
  }
  
  separator=" "
  migrMatCol="coral"
  admixCol="blue"
  popFusionColor="black"
  popCol="lightgrey"
  popBorderCol="black"
  inbreedColor="lightblue"
  oldSampColor="olivedrab4"
  timeCol="tan4"
  growthCol="hotpink3"
  propLastsegment=0.05
  migMatNameProp=0.8
  migMatLineLength=0.3
  timeProp=0.6
  maxRadius=1/40
  minRadius=maxRadius/3
  arrowLength=0.2
  logScaleAxis=""
  timeOffset=0.25
  migrOffset=0.1
  curvedArrowLTY=1
  drawLogPopSize=T
  plotMigrRates=T
  migrRateTextSizeCEX=0.5
  
  #Define plot area size for PDF
  pdf.x.size=7
  pdf.y.size=10
  
  # parFile=readLines(con=parFileName)
  
  ########################### Reading rescaling factor ###########################
  
  rescalingFactor=1.0
  rescaleParams=FALSE
  if (file.exists(rsfFileName)) {
    rescaleParams=TRUE
    rescalingFactor=scan(rsfFileName, skip=1, nlines = 1)
  }
  
  ###########################      Reading par file      #########################
  
  parFile=scan(parFileName, character(0), sep = "\n", strip.white = TRUE) # separate each line
  
  if (rescaleParams) pdfFileName=paste(parFileName, "_rescaled_ZOOM.pdf", sep="") else {
    pdfFileName=paste(parFileName, "_ZOOM.pdf", sep="")
  }
  
  #--- Clean input file by removing consecutive separators and keep ---------------
  
  #--- Function to remove separators within a string
  removeTrailingSep=function(string, sep) {
    temp=strsplit(string, split=sep)
    temp2=temp[[1]][nchar(temp[[1]])>0]
    cleanStr=temp2[1]
    if (length(temp2)>1) {
      for (i in 2:length(temp2)) {
        cleanStr=paste(cleanStr, temp2[i], sep=sep)
      }
    }
    cleanStr
  }
  
  #--- Replace Keep by -9999
  replaceKeep=function(string) {
    if (grepl("keep", string)) {
      return(gsub("keep", "-9999", string))
    }
    return(string)
  }
  
  #Remove both multiple consecutive whitespace and tabs and replace keep keyword
  for (i in 1:length(parFile)) {
    parFile[i]=removeTrailingSep(parFile[i], sep="\t")
    parFile[i]=removeTrailingSep(parFile[i], sep=" ")
    parFile[i]=replaceKeep(parFile[i])
  }
  
  #-------------------------------------------------------------------------------
  
  #--- Get number of samples on line 2 -----
  l.numsamples=parFile[2]
  sp.l=unlist(strsplit(l.numsamples, split=" "))
  tab.l=unlist(strsplit(l.numsamples, split="\t"))
  if (length(sp.l)>length(tab.l)) {
    numSamples=as.numeric(sp.l[1])
    separator=" "
    separator2="\t"
  } else {
    numSamples=as.numeric(tab.l[1])
    separator="\t"
    separator2=" "
  }
  
  #--- Reading numbers on separate lines -----
  getNumbers=function(start, parFile, numSamples) {  
    for (i in 1:numSamples) {
      curnum=as.numeric(unlist(strsplit(parFile[start+i], split=separator))[1])
      if (i==1) {
        num=curnum 
      } else {
        num=c(num, curnum)
      }
    }
    num
  }
  
  #--- Get population sizes -----------------
  start=3
  popSizes=getNumbers(start, parFile, numSamples)
  #Rescaling pop sizes
  popSizes=round(popSizes*rescalingFactor, digits=0)
  
  iniPopSizes=popSizes
  
  #--- Get sample sizes -----------------
  
  readSampleSizesTimesAndInbreedingLevel=function(start, parFile, numSamples) {
    for (i in 1:numSamples) {
      curLine=unlist(strsplit(parFile[start+i], split=separator))
      curSampSize=as.numeric(curLine[1])
      curSampTime=0
      curInbreeding=0
      if (length(curLine)>1) curSampTime=as.numeric(curLine[2])
      if (length(curLine)>2) curInbreeding=as.numeric(curLine[3])
      if (i==1) {
        sampSize=curSampSize
        sampTime=curSampTime
        inbreeding=curInbreeding
      } else {
        sampSize=c(sampSize,curSampSize)
        sampTime=c(sampTime,curSampTime)
        inbreeding=c(inbreeding,curInbreeding)
      }
    }
    list(ss=sampSize, st=sampTime, inb=inbreeding)
  }
  
  start=start+numSamples+1
  # sampSizes=getNumbers(start, parFile, numSamples)
  sampSizesStats=readSampleSizesTimesAndInbreedingLevel(start, parFile, numSamples)
  #Rescaling sample times
  sampSizesStats$st=round(sampSizesStats$st*rescalingFactor, digits=0)
  
  sampSizes=sampSizesStats$ss
  sampTimes=sampSizesStats$st
  inbrCoeff=sampSizesStats$inb
  
  #--- Get growth rates -----------------
  start=start+numSamples+1
  growthRatesInitial=getNumbers(start, parFile, numSamples)
  # save this into growthRates which will be used and updated when printing historical events
  growthRates=growthRatesInitial
  
  #--- Get number of migration matrices -----------------
  start=start+numSamples+1
  numMigMat=as.numeric(unlist(strsplit(parFile[start+1], split=separator))[1])
  
  #--- Read migration matrix
  readMigMat=function(start, parFile, numSamples) {
    for (i in 1:numSamples) {
      tmpseparator=separator
      if(is.na(as.numeric(unlist(strsplit(parFile[start+i], split=separator)))[1])) {
        tmpseparator=separator2
      }
      curmigs=as.numeric(unlist(strsplit(parFile[start+i], split=tmpseparator)))
      if (i==1) {
        migs=curmigs 
      } else {
        migs=rbind(migs, curmigs)
      }
    }
    rownames(migs)=1:numSamples
    migs 
  }
  
  #--- Get migration matrices as a list --------------
  start=start+2
  migMats=list()
  if (numMigMat>0) {
    for (i in 1:numMigMat) {  
      curMigMat=readMigMat(start, parFile, numSamples) 
      migMats[[i]]=curMigMat
      start=start+numSamples+1
    }
  }
  
  #Rescaling migration rates
  if(numMigMat>0) {
    for (i in 1:numMigMat) {
      migMats[[i]]=migMats[[i]]/rescalingFactor
    }      
  }
  
  
  #--- Get number of historical events
  start=start+1
  numHistEvents=as.numeric(unlist(strsplit(parFile[start], split=separator))[1])
  
  ###################### HISTORICAL EVENTS HANDLING ##############################
  
  #..... Read Historical Event .......
  last.he.time=0
  if (numHistEvents>0) {
    for (i in 1:numHistEvents) {
      start=start+1
      curHE=as.numeric(unlist(strsplit(parFile[start], split=separator)))
      #Rescaling time of event
      curHE[1]=round(curHE[1]*rescalingFactor, digits=0)
      if (i==1) {
        histEvents=curHE
        last.he.time=curHE[1]
      } else {
        histEvents=rbind(histEvents, curHE)
        if (histEvents[i,1] > last.he.time) last.he.time=histEvents[i,1] 
      }
    }
    if (numHistEvents>1) rownames(histEvents)=1:numHistEvents
  }
  names(last.he.time)=""
  last.he.time=as.numeric(last.he.time)
  
  yTimeLimit=0
  if (last.he.time!=0) yTimeLimit=last.he.time*(1+propLastsegment) #Add propLastsegment to y axis after last event (to draw stuff)
  
  #Reorder events by their times
  if (numHistEvents>1) histEvents=histEvents[order(histEvents[,1],decreasing=FALSE),] else {
    if (numHistEvents==1) histEvents=matrix(histEvents, nrow=1,  byrow = T)
  }
  
  endReadParFile=start
  
  
  
  ##############################  PLOTTING THE MODELED SCENARIO ######################################
  
  
  #--- Graphical functions ...........................................................................
  
  fullHeadArrow=function(x0, y0, x1, y1, length, angle, color="black",  weight=1) {
    arrows(x0, y0, x1, y1, length, angle, code=2, lty=1, col=color, lwd=weight)
    arrows(x0, y0, x1, y1, length, angle*0.80, code=2, lty=1, col=color, lwd=weight)
    arrows(x0, y0, x1, y1, length, angle*0.60, code=2, lty=1, col=color, lwd=weight)
    arrows(x0, y0, x1, y1, length, angle*0.40, code=2, lty=1, col=color, lwd=weight) 
    arrows(x0, y0, x1, y1, length, angle*0.20, code=2, lty=1, col=color, lwd=weight)
    arrows(x0, y0, x1, y1, length, angle*0.10, code=2, lty=1, col=color, lwd=weight)
  }
  drawTriangle=function(growth, x, y, size, aspRatio, color) {
    if (growth>0) {
      x0=x; y0=y
      x1=x-size/2; y1=y+size/2*aspRatio
      x2=x+size/2; y2=y+size/2*aspRatio
    } else {
      x0=x-size/2; y0=y
      x1=x; y1=y+size/2*aspRatio
      x2=x+size/2; y2=y
    }
    polygon(c(x0, x1, x2, x0), c(y0, y1, y2, y0), col=color)
    return(y+size/2*aspRatio)
  }
  
  
  #--- Computing maximum current pop size ............................................................
  maxPopSize=popSizes[1]
  minPopSize=popSizes[1]
  for (i in 2:length(popSizes)) {
    if (popSizes[i]>maxPopSize) maxPopSize=popSizes[i];
    if (popSizes[i]<minPopSize) minPopSize=popSizes[i];
  }
  
  #--- Find min and max pop sizes over the whole population history ..................................
  ps=popSizes
  isGrowth=FALSE
  if (numHistEvents) {
    for (i in 1:numHistEvents) {
      he=histEvents[i,]
      sink=he[3]+1 
      resize=he[5]
      growth=0
      #Implement possible growth
      if (i>1) {
        curTime=he[1]
        prevHE=histEvents[i-1,]
        prevTime=prevHE[1]
        growth=prevHE[6]
        formerEvent=-1
        if (growth==-9999) { #handle transformed keep keyword
          if (i-2>0) {
            formerEvent=histEvents[i-2,]
            growth=formerEvent[6]
          } else {
            growth=growthRatesInitial[sink]
          }          
        }
        if (growth!=0) isGrowth=TRUE
        ps[sink]=ps[sink]*exp(growth*(curTime-prevTime))
        if (ps[sink]==0) print(paste("Deme ", sink-1, " is zero at time ",  curTime, " due to large negative growth (", growth, ")", sep=""))
      }
      #Implement resize
      ps[sink]=ps[sink]*resize
      if (ps[sink]>maxPopSize) maxPopSize=ps[sink]
      if (ps[sink]<minPopSize) minPopSize=ps[sink]
    }
  }
  
  #-- Function to compute the circle radius for a given pop size ....................................
  interpolRadius=function(curSize, minSize, maxSize, minRadius, maxRadius, logScale) {
    if(logScale) {
      minSize=log10(minSize)
      maxSize=log10(maxSize)
      curSize=log10(curSize)
    }
    curRadius=minRadius+(curSize-minSize)*(maxRadius-minRadius)/(maxSize-minSize)
    curRadius
  }
  
  #==============================                     ========================
  #==============================   BEGIN MAIN PLOT   ========================
  #==============================                     ========================
  
  library("plotrix")
  library("diagram")
  
  if (printPDF) {
    pdf(pdfFileName, width=pdf.x.size, height=pdf.y.size)
  } 
  
  par(xpd=F, mar=c(6,4,3,0.5))
  
  maxRadius=maxRadius*(numSamples+2)
  minRadius=minRadius*(numSamples+2)
  
  title=parFileName
  if (rescaleParams) {  
    title=paste(title, "(rescaled params)")
  }
  
  if(length(histEvents[,1])>1) {
    yTimeLimit=(as.numeric(histEvents[length(histEvents[,1])-1,1])+200)*gentime  
  } else {
    yTimeLimit=(as.numeric(histEvents[1,1])+200)*gentime  
  }
  
  plot(x=1:numSamples, type="n", xlab="", ylab="time (years)", xlim=c(-0.5, numSamples+1.5), cex.main=0.8,
       ylim=c(-5000, yTimeLimit), main=title, xaxt = 'n', log=logScaleAxis, cex.axis=0.9, , cex.lab=0.9)
  axis(side=1, labels=c("Mig Mat",pop.names, " \nTimes"), at=0:(numSamples+1), cex.axis=0.8)
  mtext("demes", side=1,line=2, cex=0.8)
  
  w <- par("pin")[1]/diff(par("usr")[1:2])
  h <- par("pin")[2]/diff(par("usr")[3:4])
  aspRatio <- w/h
  
  # aspRatio=yTimeLimit/(numSamples+2)
  
  #--- Plotting population circles according to their current size and sampling time
  sampTimes_plot <- sampTimes*gentime
  for (i in 1:numSamples) {
    curRadius=interpolRadius(popSizes[i], minPopSize, maxPopSize, minRadius, maxRadius, drawLogPopSize)
    #   print(curRadius)
    if (i==1) topOfCircle=sampTimes_plot[i]+curRadius*aspRatio else topOfCircle=c(topOfCircle, sampTimes_plot[i]+curRadius*aspRatio)
    curColor=popCol
    if (sampSizes[i]==0) curColor="white"
    if (inbrCoeff[i]==0) {
      draw.circle(i, sampTimes_plot[i], radius = curRadius, col=curColor, border=popBorderCol) 
    } else {
      floating.pie(i, sampTimes_plot[i], c(inbrCoeff[i], 1-inbrCoeff[i]), radius = curRadius, 
                   col=c(inbreedColor,curColor), startpos=pi, border=popBorderCol)    
    }
    if (sampTimes_plot[i]>0) {
      text(i, sampTimes_plot[i]-curRadius*aspRatio, labels=sampTimes_plot[i], cex=timeProp, col=oldSampColor, pos=1)
    }
    #Draw a vertical arrow in case of pop growth
    curGrowthRate=growthRates[i]
    
    #--- Draw gtrowing or shrinking triangle on top of pop circle to reflect growth type
    if (curGrowthRate!=0) {
      arLength=0.15
      topOfCircle[i]=drawTriangle(curGrowthRate, x=i, y=topOfCircle[i], arLength, aspRatio, color=growthCol)
    }
  }
  
  #--- Handle first migration matrix .......................
  curMigMatNum=0
  curvature=0.0075*last.he.time*gentime
#   text(0-migrOffset, 0, labels=0, cex=migMatNameProp, col=migrMatCol)
  if (numMigMat) {
    curMigMat=migMats[1][[1]]
    for (sink in 1:numSamples) {
      for (sourc in 1:numSamples) {
        if (sink!=sourc & curMigMat[sourc, sink]>0)  {
          differ=sourc-sink
#           curvedarrow(from=c(sourc, 0), to=c(sink,0), curve=curvature*(abs(differ)*0.55^abs(differ)), arr.adj=1, 
#                       arr.pos=0.5, arr.type="triangle", arr.col=migrMatCol, lwd=1, lty=curvedArrowLTY, 
#                       lcol=migrMatCol, arr.length=arrowLength)
          curvedarrow(from=c(sourc, 0), to=c(sink,0), curve=curvature*(abs(differ)*0.45^abs(differ)), arr.adj=1, 
                      arr.pos=0.5, arr.type="triangle", arr.col=migrMatCol, lwd=1, lty=curvedArrowLTY, 
                      lcol=migrMatCol, arr.length=arrowLength)   
          if (plotMigrRates) {
            curNm=round(curMigMat[sourc, sink]*popSizes[sourc], digits=2)
#             print(paste("source:", sourc))
#             print(paste("sink:",sink))
#             print(paste("2Nm:",curNm))
            if (differ>0) {
              text(sink+abs(differ)/2, aspRatio*0.05*abs(differ), labels=curNm, cex=migrRateTextSizeCEX, col=migrMatCol) 
            } else {
              text(sourc+abs(differ)/2, -aspRatio*0.05*abs(differ), labels=curNm, cex=migrRateTextSizeCEX, col=migrMatCol) 
            } 
          }
        }
      }
    }
  }
  
  #---Draw all events on the population tree
  lastTime=0
  activePops=1:numSamples
  numActivePops=numSamples
  lastSink=-1
  
  if(numMigMat) {
    keepCurNm <- matrix(0,ncol=ncol(curMigMat), nrow=nrow(curMigMat))
  }
  
  
  if (numHistEvents) {
    for (i in 1:numHistEvents) {
#       print(paste("EVENT:", i))
#       if(curMigMatNum>0) {
#         stop("Wrong Mig mat")
#       }
#       print(popSizes)
      #Extract historical event
      he=histEvents[i,]
      he.time=he[1]*gentime
      he.source=he[2]+1 #+1 is due to the use of base 0 for deme number in fsc
      he.sink=he[3]+1
      he.migr=he[4]
      he.resize=he[5]
      he.growth=he[6]
      if(he.growth==-9999) { #handle transformed keep keyword
        #       he.growth=histEvents[i-1,6]
        he.growth=growthRates[he.sink] #Keep current growth rate
      }
      he.migrMat=he[7]
      
      if(he.migrMat==-9999) { #handle transformed keep keyword
        he.migrMat=curMigMatNum        
      }
      
#       print(paste("At event ", i," current mig matrix:", he.migrMat))
      
      #-- Draw event time
      if (i%%2) {
        slide=timeOffset
      } else {
        slide=-timeOffset
      }
      
      #Draw all vertical segments ....................
      for (p in 1:length(activePops)) {
        #Handle population with older sampling times
        if (sampTimes_plot[activePops[p]]>0) {
          minTime=max(lastTime, topOfCircle[activePops[p]]);
        } else {
          minTime=lastTime;
        }
        if (he.time>sampTimes_plot[activePops[p]] & he.time>topOfCircle[activePops[p]]) {
          if (activePops[p]==lastSink | i==1) {
            segments(activePops[p], topOfCircle[activePops[p]], activePops[p], he.time)    
          } else {
            segments(activePops[p], minTime, activePops[p], he.time)        
          }
        }
      }     
      
      #Handle growth rate changes since last event ......................
      #Update pop sizes according to current growth rates
      if (i>1) {
        prev.he=histEvents[i-1,] 
        branchLength=he.time-lastTime
        for (p in 1:length(activePops)) {
          curPop=activePops[p]
          popSizes[curPop]=popSizes[curPop]*exp(growthRates[curPop]*branchLength);
        }
      }
      #Update growth rate
      growthRates[he.sink]=he.growth
      
      
      lastTime=he.time
      #Handle resize of sink pop ........................
      popSizes[he.sink]=popSizes[he.sink]*he.resize
      
      curRadius=interpolRadius(popSizes[he.sink], minPopSize, maxPopSize, minRadius, maxRadius, drawLogPopSize)
      #     popRadius[he.sink]=curRadius
      topOfCircle[he.sink]=he.time+curRadius*aspRatio
      draw.circle(he.sink, he.time, radius = curRadius, col="#0000ff22", border=popBorderCol)
      
      
      #--- Draw gtrowing or shrinking triangle on top of pop circle to reflect growth type
      if (he.growth!=0) {
        arLength=0.15
        topOfCircle[he.sink]=drawTriangle(he.growth, x=he.sink, y=topOfCircle[he.sink], arLength, aspRatio, color=growthCol)
      }
      
      #Handle population fusion .........................
      if (he.migr>=1 & he.sink!=he.source) { #This is a population fusion
        if (numActivePops==numSamples) removedPops=he.source else removedPops=c(removedPops, he.source)
        numActivePops=numActivePops-1
        activePops=(1:numSamples)[-removedPops]
        #Draw connecting arrows from source to sink
        fullHeadArrow(he.source, he.time, he.sink, he.time, length=0.15, angle=20)      
        #Redraw time with the right color
        text(numSamples+1+slide, he.time, labels=he.time, cex=timeProp, col=popFusionColor)
      } else {
        #--- Handle admixture event ........................
        if (he.migr>0 & he.migr<1) { 
          #Draw connecting arrows from source to sink      
          segments(he.source, he.time, he.sink, he.time, col=admixCol, lty=2)
          if (he.sink>he.source) {
            fullHeadArrow(he.sink-0.15, he.time, he.sink, he.time, length=0.15, angle=20, color=admixCol)
            text(he.sink-0.75, he.time+(25*gentime), format(he.migr, digits=3, scientific=F), col=admixCol, cex=0.7)
          } else {
            fullHeadArrow(he.sink+0.15, he.time, he.sink, he.time, length=0.15, angle=20, color=admixCol)
            text(he.sink+0.75, he.time-(25*gentime), format(he.migr, digits=3, scientific=F), col=admixCol, cex=0.7)
          }     
          #Redraw time with the right color
          text(numSamples+1+slide, he.time, labels=he.time, cex=timeProp, col=admixCol)
        }
        else text(numSamples+1+slide, he.time, labels=he.time, cex=timeProp, col=timeCol)
      }    
      
      
      #--- Handle migmat change  ............................
      if (i!=numHistEvents) nextTime=histEvents[i+1,1] else {
        nextTime=yTimeLimit
      }
      time2DrawArrows=(he.time+nextTime)/2
      
      if(numMigMat>0) {
        if (he.migrMat!=curMigMatNum | he.resize < 1.0 | he.resize > 1.0) {
          curMigMat=migMats[he.migrMat+1][[1]]
          for (sink in 1:numSamples) {
            for (sourc in 1:numSamples) {
              if (sink!=sourc & curMigMat[sourc, sink]>0)  {  
                differ=sourc-sink
                #               curvedarrow(from=c(sourc, time2DrawArrows), to=c(sink, time2DrawArrows), curve=curvature*(abs(differ)*0.55^abs(differ)), 
                #                           arr.adj=1, arr.pos=0.5, arr.type="triangle", arr.col=migrMatCol, lwd=1, lty=curvedArrowLTY, 
                #                           lcol=migrMatCol, arr.length=arrowLength)  
                if (plotMigrRates) {
                  #Write Nm values
                  curNm=round(curMigMat[sourc, sink]*popSizes[sourc], digits=2)
                  keepCurNm[sourc, sink] <- curNm
                  if (differ>0) {
                    text(sink+abs(differ)/2, time2DrawArrows+aspRatio*0.005*abs(differ), labels=curNm, cex=migrRateTextSizeCEX, col=migrMatCol) 
                  } else {
                    text(sourc+abs(differ)/2, time2DrawArrows-aspRatio*0.005*abs(differ), labels=curNm, cex=migrRateTextSizeCEX, col=migrMatCol) 
                  } 
                }
              }
            }
            curMigMatNum=he.migrMat
          }
#           print(keepCurNm)
        }
        #--- Draw separation between migration matrices numbers on the left
        segments(-migMatLineLength/2, he.time, migMatLineLength/2, he.time, lty=3, col=migrMatCol)
        #Output current valid migration matrix      
        if (i%%2) {
          slide=migrOffset
        } else {
          slide=-migrOffset
        } 
#         print(paste("Time event", he.time,"Printing the mig matrix event ", i, ", at time", time2DrawArrows,"the he.migmat=", he.migrMat))
#         text(0+slide, time2DrawArrows, labels=he.migrMat, cex=migMatNameProp, col=migrMatCol)    
        text(0+slide, he.time, labels=he.migrMat, cex=migMatNameProp, col=migrMatCol)            
        lastSink=he.sink
      }
      
    }
  }
  
  #-- Draw last branch
  segments(activePops[1],topOfCircle[activePops[1]], activePops[1], yTimeLimit)
  
  #==============================   PLOT LEGENDS IN MARGINS   ========================
  
  #Compute space available in margin
  minY.coo=grconvertY(0, from="nic", to="user")
  
  par(xpd=NA)
  #--- Draw population size scale with circles of different sizes .................
  maxOrder=ceiling(log10(maxPopSize))
  minOrder=floor(log10(minPopSize))
  popSizeRadius=10^(maxOrder:minOrder)
  winWidth=numSamples+2
  ypos=3/4*minY.coo
  text(x=-winWidth/10*1.2, y=ypos, labels="Pop. \nsizes ", cex=.8, pos=2)
  
  for (i in 1:length(popSizeRadius)) {
    curRadius=interpolRadius(popSizeRadius[i], minPopSize, maxPopSize, minRadius, maxRadius, drawLogPopSize)
    #   print(curRadius)
    if (curRadius>0) {
      xpos=-winWidth/10+(i-1)*winWidth/10
      draw.circle(xpos, ypos, radius=curRadius, col=popCol, border=popBorderCol)
    }
    text(xpos, ypos-abs(ypos)*0.1, popSizeRadius[i], cex=0.7, pos=1, col="black")
  }
  
  #--- Legend for growing or shrinking populations ...............................
  if (isGrowth) {
    x=winWidth-1.5*winWidth/10; y=ypos+abs(ypos)*0.2
    text(x, y-abs(ypos)*0.1, labels="Populations", cex=0.8)
    x=winWidth-2*winWidth/10; y=ypos-abs(ypos)*0.1
    drawTriangle(1, x, y, size=0.15, aspRatio, color=growthCol)
    text(x, y-abs(ypos)*0.1, labels="growing", pos=NULL, cex=0.7)
    x=winWidth-winWidth/10
    drawTriangle(-1, x, y, size=0.15, aspRatio, color=growthCol)
    text(x, y-abs(ypos)*0.1, labels="shrinking", pos=NULL, cex=0.7)  
  }
  
  
  if (printPDF) dev.off()
  
  ###########################                            #########################
  ########################### Printing rescaled par file #########################
  ###########################                            #########################
  
  if (rescaleParams) {
    genericName=substr(basename(parFileName), 1, nchar(basename(parFileName)) - 4) 
    rescaledParFileName=paste(genericName,"_rescaled.par",sep="")
    #Redirects output to par file
    sink(file=rescaledParFileName)  
    
    #...... write header .........................................................
    writeLines(paste("//Parameters for the coalescence simulation program",
                     ", edited with R script. Input file:  ",
                     parFileName, ", Rescaling factor: ", rescalingFactor, sep=""))  
    writeLines(paste(numSamples, " samples to simulate"))
    #...... write population sites ...............................................
    writeLines("//Population effective sizes")
    for (i in 1:length(iniPopSizes)) {
      writeLines(paste(iniPopSizes[i], sep=""))
    }
    #...... write sample sizes ...................................................
    writeLines("//Samples sizes")
    for (i in 1:length(sampSizes)) {
      curLine=sampSizes[i]
      if (inbrCoeff[i]>0) {
        curLine=paste(curLine, round(sampTimes[i]), inbrCoeff[i])
      } else {
        if (sampTimes[i]>0) curLine=paste(curLine, round(sampTimes[i]))
      }
      writeLines(paste(curLine, sep=""))
      
    }
    #...... write growth rates ...................................................
    writeLines("//Growth rates : negative growth implies population expansion")
    for (i in 1:length(growthRatesInitial)) {
      writeLines(paste(growthRatesInitial[i], sep=""))
    }
    #...... write migration matrices .............................................
    writeLines("//Number of migration matrices : 0 implies no migration between demes")
    writeLines(paste(numMigMat, sep=""))
    if(numMigMat>0) {
      for (i in 1:length(migMats)) {
        curMigMat=migMats[[i]]
        writeLines(paste("//Matrix ", i-1, sep=""))
        for (r in 1:nrow(curMigMat)) {
          curRow=curMigMat[r,]
          migrs=curRow[1]
          for (c in 2:length(curRow)) {
            migrs=paste(migrs,curRow[c])
          }
          writeLines(paste(migrs, sep=""))
        }
      }
    }
    #...... write historical events...............................................
    writeLines("//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index")
    if(numHistEvents>0) {
      writeLines(paste(nrow(histEvents), " historical events", sep=""))  
      for (i in 1:nrow(histEvents)) {
        curHE=histEvents[i,]
        h.items=curHE[1]
        for (c in 2:5) {
          h.items=paste(h.items, curHE[c])
        }
        for (c in 6:7) {
          if (curHE[c]==-9999) h.items=paste(h.items, "keep") else {
            h.items=paste(h.items, curHE[c])
          }
        }
        writeLines(h.items)
      }        
    } else {
      writeLines("0 historical events")
    }
    
    
    #...... Write rest of the original par file ..................................
    for (i in (endReadParFile+1):length(parFile)) {
      writeLines(parFile[i])
    }
    
    sink(type = "message")
    sink()
  }
}