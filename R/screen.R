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
  sampleAssigned <- merge(data.frame(sampleIDs = seq_len(files.short), filename=files.short, stringsAsFactors=FALSE), sampleList)
  return(sampleAssigned)
}

#' Merge sample list information to a profile container
#' 
#' Adds sample list columns (e.g. timepoint, species, etc as specified in the sample list) to the `peaks` table in the `profiles` container,
#' such that it can later easily be sorted, filtered etc.  
#' 
#' @param profiles EnviMass profile container
#' @param sampleList A sample list with attached sampleIDs, as obtained from [assignSamples].
#' @param zerofill If `zerofill` is active, a new zero-intensity result is generated for every possible profile-sample combination
#' 	(otherwise many results are "missing").  
#' @return The updated profile container.
#' 
#' @md
#' @author stravsmi
#' @export
assignProfiles <- function(profiles, sampleList, zerofill = FALSE)
{
  
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
screenProfiles <- function(profiles, suspects, polarity = "+", ppmLimit = getOption("RMassScreening")$screenProfiles$ppmLimit)
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
  
  
  hits <- lapply(1:nrow(peaklist), function(n)
  {
    mass <- peaklist[n,"mean_mz"]
    peaklist.hit <- potential[abs(potential$mass - mass) < ppm(mass, ppmMax),,drop=FALSE]
    peaklist.hit$dppm <- ((peaklist.hit$mass / mass)-1)*1e6
    peaklist.hit$profileID <- rep(peaklist[n, "profile_ID"], nrow(peaklist.hit))
    peaklist.hit$mz <- rep(peaklist[n, "mean_mz"], nrow(peaklist.hit))
    peaklist.hit$RT <- rep(peaklist[n, "mean_RT"], nrow(peaklist.hit))
    peaklist.hit$int <- rep(peaklist[n, "mean_int"], nrow(peaklist.hit))
    return(peaklist.hit)
  })
  
  hits.total <- do.call(rbind, hits)
  return(hits.total)
}