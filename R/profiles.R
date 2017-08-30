
#' Fill enviMass profiles from picked files
#' 
#' Reads a set of files containing enviPick results, and fills an enviMass profile container with the data.
#' Does not compute the actual profiles; this can be done with [computeProfiles].
#' 
#' Code is copied and adapted from enviMass.
#' 
#' @param dataDir Folder containing the enviPick files in .MSlist format. May be either MSlist.pos and MSlist.neg files, or files with an arbitrary pattern.
#' 		Typically, filenames should correspond to the names built by [batchPick] or [comprehensiveBatchPick], respectively: 
#' 		* If only MS1 data was peak-picked ([batchPick]), raw file `FILENAME.mzXML` should give filenames `FILENAME.mzXML.MSlist.pos.RData`
#' 			 (or `neg`, respectively) if only MS1 were picked
#'      * If DIA data was picked (\code{\link{comprehensiveBatchPick}}) a raw file `FILENAME.mzXML` should give filenames:
#'          * `FILENAME.mzXML.MSlist.1-pos--` for the MS1 positive
#' 			* e.g. `FILENAME.mzXML.MSlist.2-neg-200-100` for MS2-DIA, negative, center mass 200, isolation width 100 etc.
#' 		However other patterns work too, since `pattern` specifies which file pattern is searched.  
#' @param files The filenames of the raw files, with or without full path. NOT the filenames or path of the enviPick results.
#' @param polarity If `pattern` is not used, either `+` or `-`; this chooses which files to look for.
#' @param pattern The file pattern to search for (i.e. what is appended to the raw file name). 
#' 		If a raw file `FILENAME.mzXML` was picked into `FILENAME.mzXML.MSlist.ABCDE.RData`, the `pattern` must be set to `ABCDE`, see above.
#' 		If this is longer than one, a separate container will be created for each pattern and filled accordingly.
#' @param rtrange Optional RT range of peaks to extract. Use this to cut off e.g. the column dead volume.
#' @return A profile container in enviMass format (for no `pattern` or a single `pattern`, or a list of containers named after the 
#' 	corresponding pattern/scan.
#' 
#' 
#' @author stravsmi, mostly copied from Martin Loos
#' @export
#' @md
fillProfiles <- function(dataDir, files, polarity="+", pattern=NA, rtrange=c(-Inf, Inf)) 
{
  # If this is a multiscan container, fill every one individually
  if(length(pattern) > 1)
  {
	  profiles <- lapply(pattern, function(run) fillProfiles(dataDir, files, pattern = run, rtrange = rtrange))
	  names(profiles) <- pattern
	  return(profiles)
  }
	
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
  # add a slash to dataDir if missing
	regex <- "[/\\]$"
	dataDir <- sub(regex, "", dataDir)
	dataDir <- paste0(dataDir, "/")
	
 
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
	pks <- MSlist[[8]]
	pks <- pks[(pks[,5] >= rtrange[[1]]) & (pks[,5] <= rtrange[[2]]),]
	MSlist[[8]] <- pks
	
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
    
	pks <- MSlist[[8]]
	pks <- pks[(pks[,5] >= rtrange[[1]]) & (pks[,5] <= rtrange[[2]]),]
	peaklist <- pks
	
	
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

#' Compute profiles using enviMass.
#'  
#' Wrapper for [agglomer()] and [partcluster()] from enviMass.
#' 
#' @param profiles A profiles container in enviMass format, as generated by [fillProfiles].
#' 	Alternatively a list of profile containers (for multiple scans), also as generated by
#' [fillProfiles].
#' @param dmass Mass tolerance in ppm
#' @param dret Retention time window in ppm.
#' @return Processed profile container.
#' 
#' @author stravsmi
#' @export
computeProfiles <- function(profiles, dmass=3, dret=60)
{
	# if this is a profile list, compute all
	if(!is.data.frame(profiles[[1]]))
	{
		profiles <- lapply(profiles, computeProfiles,
				dmass=dmass, dret=dret)
		return(profiles)
	}
	
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


#' Cut single profiles
#' 
#' Removes profiles which do not have at least a specified number of peaks. Typically the limit is 1, i.e.
#' requiring a profile to consist of at least two peaks.
#' 
#' @param profiles A profile container or list of profile containers
#' @param minPeaks Minimum peak count required
#' @return Filtered profiles.
#' 
#' @note Check: Should the `profile_ID` `index_prof` also be rewritten?
#' 
#' @author stravsmi
#' @export
cutSingleProfiles <- function(profiles, minPeaks = 1)
{
	# if this is a list of profile containers, process them all independently
	if(!is.data.frame(profiles[[1]]))
	{
		profiles <- lapply(profiles, cutSingleProfiles, minPeaks)
		return(profiles)
	}	
	
	# eliminate all profiles which don't have a minimum number of peaks
	profiles$index_prof <- profiles$index_prof[
			profiles$index_prof[,"number_peaks_total"] > minPeaks,,drop=FALSE
			]
	# rewrite the profile IDs in $peaks to match the new profile row in the index_prof table
	profiles$peaks[,"profileIDs"] <- match(
			profiles$peaks[,"profileIDs"],
			profiles$index_prof[,"profile_ID"],
			0
			)	
	profiles$index_prof[,"profile_ID"] <- seq_len(nrow(profiles$index_prof))
			
			
	# remove peaks with no profile
	profiles$peaks <- profiles$peaks[
			profiles$peaks[,"profileIDs"] != 0,,drop=FALSE]
			
	profiles
}