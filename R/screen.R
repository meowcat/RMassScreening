#' Map a sample list to sample IDs
#' 
#' Merges the `sampleList` with the corresponding enviMass sample ID (which is the position of a name in the `files` array.)
#'  
#' The filename from sampleList is coordinated with a "sampleIDs" ID from profiles via the "files" list,
#' and the rest of the sampleList table is merged into it.
#' 
#' @param files must be the files list in the same order they were used when building profiles.
#' @param sampleList data frame with a `filename` column
#' 
#' @md
#' @export
assignSamples <- function(files, sampleList)
{
  files <- basename(files)
  files.short <- unlist(lapply(files, function(file) strsplit(file, ".", fixed=TRUE)[[1]][[1]]))
  sampleAssigned <- merge(data.frame(sampleIDs = seq_along(files.short), filename=files.short, stringsAsFactors=FALSE), sampleList)
  return(sampleAssigned)
}



#' Generate sample table from filenames
#' 
#' Splits a vector of filenames by a separator to make a corresponding sample table.
#' E.g. 001_file1_15_22.mzXML gets split to 001 - file1 - 15 - 22.
#' Should handle missing entries more or less
#' 
#' @param files Files vector
#' @param sep Separator
#' @param names Names for output columns. Otherwise will be named col1, col2...
#' @return A data frame with results
#' 
#' @author stravsmi
#' @export
generateSampleList <- function(files, sep="_", names=NULL)
{
	files.short <- unlist(lapply(basename(files), function(file) strsplit(file, ".", fixed=TRUE)[[1]][[1]]))
	
	d.classes <- strsplit(basename(files.short), sep, TRUE)
	
	len.max <- max(unlist(lapply(d.classes, length)), na.rm=TRUE)
	
	d.classes <- lapply(d.classes, function(cl)
			{length(cl) <- len.max
				cl}
	)
	d.table <- do.call(rbind, d.classes)
	
	if(!is.null(names))
		colnames(d.table) <- names
	else
		colnames(d.table) <- paste0("col", seq_len(len.max))
	#c("date", "running", "experiment", "time", "sample", "type")
	
	samples <- cbind(data.frame(filename=files.short), d.table)
	return(samples)
}



#' Merge sample list information to a profile container
#' 
#' Adds sample list columns (e.g. timepoint, species, etc as specified in the sample list) to the `peaks` table in the `profiles` container,
#' such that it can later easily be sorted, filtered etc.  
#' 
#' @param profiles EnviMass profile container or list of containers
#' @param sampleList A sample list with attached sampleIDs, as obtained from [assignSamples].
#' @param zerofill If `zerofill` is active, a new zero-intensity result is generated for every possible profile-sample combination
#' 	(otherwise many results are "missing").  
#' @return The updated profile container(s).
#' 
#' @md
#' @author stravsmi
#' @export
assignProfiles <- function(profiles, sampleList, zerofill = FALSE)
{
	# if this is a list of profile containers, process them all independently
	if(!is.data.frame(profiles[[1]]))
	{
		profiles <- lapply(profiles, assignProfiles, sampleList, zerofill)
		return(profiles)
	}	
  
  if(zerofill)
  {
    profiles.complete <- merge(expand.grid(profileIDs = unique(profiles$peaks[,"profileIDs"]),
                                           sampleIDs = unique(sampleList[,"sampleIDs"]))
                               , profiles$peaks, all.x=TRUE, all.y=TRUE)
    profiles.complete[is.na(profiles.complete$intensity), "intensity"] <- 0
    profiles$peaks <- profiles.complete
  }
  
  profiles$peaks <- qmerge(profiles$peaks, sampleList, by="sampleIDs", all.x=TRUE)
  
  return(profiles)
}

#' Calculate ppm distance from a mass.
#' 
#' @note This conflicts with the RMassBank declaration of `ppm`, which does something different!
#' 
#' @param mz mass
#' @param ppmlimit ppm deviation to calculate 
#' 
#' @md
#' @author stravsmi
#' @export
ppm <- function(mz, ppmlimit) (1e-6 * mz * ppmlimit)

	
#' RT filter for suspect screening
#' 
#' @note This conflicts with the RMassBank declaration of `ppm`, which does something different!
#' 
#' @param expected.rt expected retention time
#' @param measured.rt measured retention time 
#' @param rttol retention time tolerance for determining a match
	
#' @md
#' @author schollje
#' @export	
	
RTfilter <- function(expected.rt, measured.rt, rttol){
  
  match <- c()
  if(measured.rt > expected.rt - rttol && measured.rt < expected.rt + rttol)
    match <- TRUE
  else
    match <- FALSE
  
  return(match)
}
				     
#' Screen for suspect masses in a profile container.
#' 
#' Simple screening looking for suspect exact masses in enviMass profiles.
#' 
#' Could be improved.
#' 
#' @param profiles a profiles container
#' @param suspects A data frame that must contain a "mass" column. Additionally a "name" column is useful (though not formally required)
#' @param polarity "+" or "-", determining what mass to look for
#' @param ppmLimit Maximal ppm deviation from target mass
#' 
#' @export
screenProfiles <- function(profiles, suspects, polarity = "+", ppmLimit = getOption("RMassScreening")$screenProfiles$ppmLimit,
			  rttol = 0.25)
{
  if(is.null(ppmLimit))
    stop("No ppm limit specified. Specify either as a parameter or in your settings file.")
  potential <- suspects
  peaklist <- as.data.frame(profiles$index_prof)
  if(polarity == "+")
    potential$mass <- potential$mass + 1.0072
  else 
    potential$mass <- potential$mass - 1.0072
  
  
  potential$X <- NULL
  ppmMax <- ppmLimit
  
  hits <- lapply(1:nrow(potential), function(n)
  {
    mass <- potential$mass[[n]]
    dppm <- ((peaklist$mean_mz / mass) - 1) * 1e6
    hit <- peaklist[abs(dppm) < ppmMax,c("profile_ID", "mean_mz", "mean_RT", "mean_int"),drop=FALSE]
    colnames(hit) <- c("profileID", "mz", "RT", "int")
    #ppmShort <- dppm[abs(dppm) < ppmMax]
    df <- cbind(potential[rep(n, nrow(hit)),,drop=FALSE], hit)
    
  })
  
  hits.total <- do.call(rbind, hits)
  hits.total$ppm <- (hits.total$mz/hits.total$mass -1) * 1e6
  hits.total$RTmatch <- mapply(RTfilter, hits.total$mean_RT, hits.total$RT, rttol)
  return(hits.total)
}
