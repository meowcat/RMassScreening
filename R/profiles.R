
#' @export
fillProfiles <- function(dataDir, files, polarity="+", pattern=NA)
{
  if(is.na(pattern))
  {
    if(polarity == "+")
      pattern <- ".MSlist.pos.RData"
    else
      pattern <- ".MSlist.neg.RData"
  }
  else
    pattern <- paste0(".MSlist.", pattern, ".RData")
  #
  # Set up an enviMass "profiles" container.
  #
  profiles<-list(0)
  profiles[[1]]<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
  colnames(profiles[[1]])<-c("peaks?","agglom?","profiled?","trends?")
  profiles[[2]]<-0  # peaks
  profiles[[3]]<-0  # datetime
  profiles[[4]]<-0  # time
  profiles[[5]]<-0  # place
  profiles[[6]]<-0  # index_agglom
  profiles[[7]]<-0  # index_prof
  profiles[[8]]<-0  # parameters
  profiles[[9]]<-0  # sample type
  names(profiles)<-c("state","peaks","datetime","sampleID","place",
                     "index_agglom","index_prof","parameters","type")
  
  #
  # First load all to determine total matrix size, then fill data into matrix
  #
  #files<-list.files(dataDir, pattern)
  files <- paste0(dataDir, basename(files), pattern)
  leng<-length(files)
  at<-0
  for(i in 1:leng){ 
    load(file=files[i]);
    # A "backwards compatibility" fix:
    if(!exists("MSlist", environment()))
    {
      if(exists("MSlist.pos", environment()))
      {
        MSlist <- MSlist.pos
        rm(MSlist.pos)
      }
      else if(exists("MSlist.neg", environment()))
      {
        MSlist <- MSlist.neg
        rm(MSlist.neg)
      }
    }
    at<-c(at+length(MSlist[[8]][,1]))    
    rm(MSlist)
  }
  peaks<-matrix(nrow=at,ncol=8,0)
  colnames(peaks)<-c("m/z","intensity","RT","peakIDs","componentIDs","sampleIDs","partitionIDs","profileIDs")
  da1<-1;
  for(i in 1:leng){ 
    load(file=files[i]);
    if(!exists("MSlist", environment()))
    {
      if(exists("MSlist.pos", environment()))
      {
        MSlist <- MSlist.pos
        rm(MSlist.pos)
      }
      else if(exists("MSlist.neg", environment()))
      {
        MSlist <- MSlist.neg
        rm(MSlist.neg)
      }
    }
    peaklist<-MSlist[[8]];
    
    peaklist<-cbind(peaklist,peaklist[,c(1,4,5)])
    colnames(peaklist)[12]<-"m/z_corr";
    colnames(peaklist)[13]<-"sum_int_corr";
    colnames(peaklist)[14]<-"RT_corr";
      
    da2<-(da1+length(peaklist[,1])-1)
    that<-c(length(peaklist[,1]))
    peaks[da1:da2,]<-cbind(peaklist[,c(12,13,14,10)],
                           rep(0,that),rep(i,that),
                           rep(0,that),rep(0,that)
    );
    da1<-c(da2+1);
    rm(MSlist,peaklist)
  }
  peaks<-peaks[order(peaks[,1],decreasing=FALSE),]
  profiles[[2]]<-peaks;
  datetime<-as.POSIXct.numeric(seq(1,leng,1),origin = "1960-01-01")
  profiles[[3]]<-datetime;
  rm(peaks)
  return(profiles)
}

#' @export
computeProfiles <- function(profiles, dmass=3, dret=60)
{
  profiles<-agglomer(
    profiles,
    dmass=dmass,
    ppm=TRUE,
    dret=dret
  )
  profiles<-partcluster(
    profiles,
    dmass=dmass,
    ppm=TRUE,
    dret=dret,
    from=FALSE,
    to=FALSE,
    progbar=FALSE,
    plot_it=FALSE
  )
  return(profiles)
}

