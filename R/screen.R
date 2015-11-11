#' 
#' 
#' "files" must be the files list in the same order they were used when building profiles.
#' "sampleList" needs a filename column
#' The filename from sampleList is coordinated with a "sampleIDs" ID from profiles via the "files" list,
#' and the rest of the sampleList table is merged into it.
#' 
#'
#' @export
assignSamples <- function(files, sampleList)
{
  files <- basename(files)
  files.short <- unlist(lapply(files, function(file) strsplit(file, ".", fixed=TRUE)[[1]][[1]]))
  sampleAssigned <- merge(data.frame(sampleIDs = 1:length(files.short), filename=files.short, stringsAsFactors=FALSE), sampleList)
  return(sampleAssigned)
}

#' @export
assignProfiles <- function(profiles, sampleList)
{
  profiles[[2]] <- qmerge(profiles[[2]], sampleList, by="sampleIDs", all.x=TRUE)
  return(profiles)
}

#' @export
ppm <- function(mz, ppmlimit) (1e-6 * mz * ppmlimit)

#'
#' 
#' "profiles" is a profiles container
#' "suspects" is a suspect list and must contain a "mass" column. Additionally a "name" column is useful
#' (though not formally required)
#' @export
screenProfiles <- function(profiles, suspects, polarity = "+", ppmLimit = 2)
{
  potential <- suspects
  peaklist <- as.data.frame(profiles[[7]])
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