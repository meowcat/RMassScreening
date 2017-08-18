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
#'
#' @examples
makeProfileMatrix <- function(profiles, type = c("ramclust", "graphviz"), addIndex = TRUE)
{
  # Drop peaks that are not in a profile
  profiles$peaks <- profiles$peaks[profiles$peaks[,"profileIDs"] != 0,,drop=FALSE]
  profiles$index_prof <- profiles$index_prof[profiles$index_prof[,"profile_ID"] != 0,,drop=FALSE]
  
  allProfiles <- acast(profiles$peaks, sampleIDs ~ profileIDs, value.var="intensity",
                       fill = 0)
  
  featureTitle <- paste0(profiles$index_prof[,"profile_ID"], "_", round(profiles$index_prof[,15],2))
  
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



#' Process RAMClust output for viewer
#'
#' @param RC Result from RAMClust call
#'
#' @return a dataframe
#' @export
#'
#' @examples
processRc <- function(RC, profile)
{
  rcAssignment <- data.frame(RC$fmz, RC$featclus, RC$frt, RC$msint)
  
  colnames(rcAssignment) <-  c("profileID", "featureID", "RT", "int")
  rcAssignment$mz <- profile$index_prof[match(rcAssignment$profileID, profile$index_prof[,"profile_ID"]), "mean_mz"]
  rcAssignment$profint <- profile$index_prof[match(rcAssignment$profileID, profile$index_prof[,"profile_ID"]), "mean_int"]
  rcAssignment$profRT <- profile$index_prof[match(rcAssignment$profileID, profile$index_prof[,"profile_ID"]), "mean_RT"]
  rcAssignment <- rcAssignment[order(rcAssignment$mz),]
  rcAssignment
}

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
#' @examples
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