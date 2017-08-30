#' Make matrix for use with ramclustR.df
#'
#' @param profiles A profiles object from enviMass
#' @param type `ramclust` or `graphviz`, in the latter case the retention time is inserted as a first row (not only in the column names).
#' @param addIndex If `TRUE`, add a sample index column in front of the table, as required by RAMClust.
#' 
#' @return A matrix m x n, where the rows (sic) are the samples and the columns are the profiles (features) without "profile 0".
#' The columns are named profileID_RT. This is a requirement for ramclustR (which usually uses mz_RT, but actually doesn't use
#' mz for calculation and therefore the content is arbitrary. It's easier to use the profile ID for precise lookup.)
#' 
#' @export
#' @md
#'
makeProfileMatrix <- function(profiles, type = c("ramclust", "graphviz"), 
		addIndex = ifelse(type[[1]] == "ramclust",TRUE,FALSE), reindexProfiles = 0)
{
	# if this is a list of profile containers, process them all independently
	if(!is.data.frame(profiles[[1]]))
	{
		# first create the profile matrix without reindexing, just to look up the # columns
		# yes this is a waste of time but it is fast anyway
		profMatrix <- lapply(seq_along(profiles), function(i)
					makeProfileMatrix(profiles[[i]], type, FALSE, 0))
		reindex <- cumsum(c(0, unlist(lapply(profMatrix, ncol))))
		profMatrix <- lapply(seq_along(profiles), function(i)
					 makeProfileMatrix(profiles[[i]], type, FALSE, reindex[[i]]))
		 names(profMatrix) <- names(profiles)
		 if(addIndex)
		 {
			 
			 indexMatrix <-  cbind("SN" = seq_len(nrow(profMatrix[[1]])))
			 profMatrix <- c(list(index = indexMatrix), profMatrix)
			 
		 }
		return(profMatrix)
	}	
	
	# Convert to data frame if necessary
  if(!is.data.frame(profiles$peaks))
	  profiles$peaks <- as.data.frame(profiles$peaks)
  # Drop peaks that are not in a profile
	profiles$peaks <- profiles$peaks[profiles$peaks[,"profileIDs"] != 0,,drop=FALSE]
  profiles$index_prof <- profiles$index_prof[profiles$index_prof[,"profile_ID"] != 0,,drop=FALSE]
  
  # make the sample x intensity table, with acast as this is faster
  allProfiles <- acast(profiles$peaks, sampleIDs ~ profileIDs, value.var="intensity",
                       fill = 0)
  
  featureTitle <- paste0(profiles$index_prof[,"profile_ID"] + reindexProfiles, "_", round(profiles$index_prof[,15],2))
  
  # for graphViz processing, add the RT as first observation
  if(type[[1]]=="graphviz")
  {
    rt <- round(profiles$index_prof[,15],2)
    allProfiles <- rbind(rt, allProfiles)
  }
  dimnames(allProfiles) <- list(
    1:nrow(allProfiles),
    featureTitle
  )
  
  # Add index column at the beginning for sample number, as required by RAMClustR
  if(addIndex)
  	allProfiles <- cbind("SN" = seq_len(nrow(allProfiles)), allProfiles)
  
  allProfiles
}


# Instead of the archaic run indexing "ramclust style", just paste the information table from the profile objects.
#' Extract essential info from profiles
#' 
#' Extracts minimal essential information from profiles table for each profile,
#' and adds profile run index as required
#' 
#' @param profiles A single enviMass profile container or a list of profile containers
#' @param index The run index to append 
#' @return A data frame with profile_ID, RT, mz, intensity and scantype index.
#' 
#' @author stravsmi
#' @export
makeProfileInfo <- function(profiles, index = 1, name = "")
{
	# if this is a list of profile containers, process them all independently
	if(!is.data.frame(profiles[[1]]))
	{
		info <- lapply(seq_along(profiles), function(i)
				{
					makeProfileInfo(profiles[[i]], i, names(profiles)[[i]])
				})
		names(info) <- names(profiles)
		return(info)
	}	
	# extract essential columns from profiles
	info <- as.data.frame(profiles$index_prof[,c("profile_ID", "mean_RT", "mean_mz", "mean_int")])
	info[,"scan"] <- rep(index, nrow(info))
	info[,"scantype"] <- rep(name, nrow(info))
	return(info)
}

#' Process RAMClust output for viewer
#'
#' @param RC Result from RAMClust call
#'
#' @return a dataframe
#' @export
#'
processSpectra <- function(cluster, profileInfo, type=c("ramclust", "graphvis"))
{
	type <- type[[1]]
	if(type == "ramclust")
		stop("Sorry, pure-ramclust result evaluation is currently being refactored. Use largevis-style processing for now.")
	if(type == "graphvis")
	{
		# split the clustering results into spectra 
		splitSpectra <- split(seq_along(cluster$clusters), cluster$clusters)
		# for each neighbor list, retrieve the corresponding profile info and split into MS1 and DIA windows
		spectra <- lapply(names(splitSpectra), function(id)
				{
					neighbors <- splitSpectra[[id]]
					probabilities <- cluster$probabilities[neighbors]
					profiles <- data.frame(neighbors, probabilities, stringsAsFactors = FALSE)
					profiles  <- cbind(profiles, profileInfo[neighbors,])
					
					subspectra <- split(profiles, profiles$scantype)
					
					subspectraRmb <- lapply(subspectra, function(spectrum)
							{
								colnames(spectrum) <- c("infoProfile", "probability", "profile", "RT", "mz",  "i", "scan", "scantype")
								sp <- new("RmbSpectrum2", mz=spectrum$mz, intensity=spectrum$i, rawOK = rep(TRUE, nrow(spectrum)),
										acquisitionNum = unique(spectrum$scan))
								sp@properties <- spectrum[,c("probability", "profile", "infoProfile", "RT", "scan", "scantype")]
								sp
							}
					)
					subspectraRmb <- as(subspectraRmb, "SimpleList")
					names(subspectraRmb) <- names(subspectra)
					
					subspectraRmb
					
				})
	}
	else
		spectra <- type
	w <- newMsmsWorkspace()
	spectra <- lapply(seq_along(spectra), function(i)
			{
				new("RmbSpectraSet",
						found = TRUE,
						complete=TRUE,
						children=spectra[[i]],
						id = as.character(i))
			})
	w@spectra <- as(spectra, "SimpleList")
	w@aggregated <- aggregateSpectra(w)
	w
}


#' Select spectra for profile container
#' 
#' From an msmsWorkspace containing DDA or DIA spectra, select (in index table) those for a specific profile container for use in the viewer.
#' 
#' @param w msmsWorkspace
#' @param scantype Scantype (profile container) to select
#' @return msmsWorkspace
#' 
#' @author stravsmi
#' @export
selectSpectraProfiles <- function(w, scantype)
{
	w@aggregated <- w@aggregated[w@aggregated$scantype == scantype,,drop=FALSE]
	w
}


	
## }

#' Quick RAMClustR application
#' 
#' Applies RAMClustR to a small subset of all profiles that are in the retention time window around a specific profile of interest.
#'
#' @param profileID The profile you want to look up (integer)
#' @param profileMatrix A sample x feature matrix as returned by \link{\code{makeProfileMatrix}}
#' @param profiles A profile list produced by \link{nontarget}
#' @param rttol Retention time window to extract from the profile matrix to cluster on. A value of 60 means +- 30 
#'        (in unit of your retention time, i.e. if you have seconds use seconds)
#' @param raw if TRUE, return the RAMClustR result instead. Useful for development purposes.
#'
#' @return a \code{data.frame} (\code{profileID, mz, RT, int}) with the "pseudospectrum" at that RT 
#'        (using the pseudo-sum-intensity from RAMClustR)
#' @export
#'
quickclust <- function(profileID, profiles, rttol = 60, raw=FALSE, profileMatrix = NULL, sr=0.5 ,st=25)
{
  rt <- profiles$index_prof[profileID , "mean_RT"]
  profiletol <- 0.001
  
  profileIDs <- which(abs(profiles$index_prof[,"mean_RT"] - rt) < rttol/2)
  # If no profile matrix is passed, compute reduced matrix on the fly
  if(is.null(profileMatrix))
  {
	  profileMatrix <- submatrix(rt, profiles, rttol, addIndex = TRUE)
  }
  else
    profileMatrix <- profileMatrix[,profileIDs]
  

  
  # calculate cluster
  clust <- ramclustR(xcmsObj = NULL, ms = profileMatrix, idmsms = NULL,  
                     featdelim = "_", timepos = 2, 
                     st = st, sr = sr, maxt = 20, deepSplit = FALSE, 
                     blocksize = 2000, mult = 5, hmax = 0.3, sampNameCol = 1, 
                     collapse = TRUE, mspout = FALSE, mslev = 1, ExpDes = NULL, 
                     normalize = "TIC", minModuleSize = 2, linkage="average")
  # find where  the feature is in the cluster
  if(raw)
    return(clust)
  featHit <- which(abs(clust$fmz - profileID) < profiletol)
  if(length(featHit) > 1)
    warning("More than one identical profile ID in cluster!")
  feature <- clust$featclus[featHit]
  points <- which(clust$featclus %in% feature)
  m <- data.frame( 
    profileID = as.integer(clust$fmz[points]),
    mz = profiles$index_prof[
      match(as.integer(clust$fmz[points]), profiles$index_prof[,"profile_ID"]), "mean_mz"],
    RT = clust$frt[points],
    int = clust$msint[points]
  )
  m[order(m[,"mz"]),]
  
}


#' Make RAMClust input submatrix for a retention time window
#' 
#' Creates a profile matrix (RAMClust input) from a subset of a profile container (enviMass) in a retention time slice
#' 
#' @param rt Center retention time
#' @param profiles profile container
#' @param rttol Retention time tolerance
#' @param addIndex Add sample index in front of table? Needed for RAMClust 
#' @return Matrix
#' 
#' @author stravsmi
#' @export
submatrix <- function(rt, profiles, rttol = 60, addIndex = TRUE)
  {

    profileIDs <- which(abs(profiles$index_prof[,"mean_RT"] - rt) < rttol/2)
    profiles$index_prof <- profiles$index_prof[profileIDs,,drop=FALSE]
    profiles$peaks <- profiles$peaks[profiles$peaks[,"profileIDs"] %in% profileIDs,,drop=FALSE]
    ti <- system.time(profileMatrix <- makeProfileMatrix(profiles, "ramclust", addIndex))
    cat(paste("Time for building profile matrix:", round(ti[[3]],2), "\r\n"))
    return(profileMatrix)
}