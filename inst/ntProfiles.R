library("nontargetData")
profiles <- profiles.pos




getPatterns <- function(profiles)
{
  peaklist <- as.data.frame(profiles[[7]][,c(14,16,15)])
  data(OrbitrapXL_VelosPro_R60000at400_q)
  isos <- make.isos(isotopes, use_isotopes = c("13C","15N","34S","37Cl","81Br"), use_charges=c(1,1,1,1,1))
  data(isotopes)
  isoPatterns <- pattern.search2(peaklist, OrbitrapXL_VelosPro_R60000at400_q,3,TRUE,use_isotopes=c("13C","37Cl","15N","81Br","34S","18O"),
                                 use_charges=c(1),isotopes=isotopes)
  i <- isoPatterns[[1]]
  i <- i[order(i[,"mean_mz"]),]
  i[,"parentMass"] <- !duplicated(i[,"group ID"])
  i <- i[order(i[,"peak ID"]),]
  
  i
}



ntProfiles <- function(profiles)
{
  peaklist <- as.data.frame(profiles[[7]])
  
  hits <- lapply(1:nrow(peaklist), function(n)
  {
    mass <- peaklist[n,"mean_mz"]
    peaklist.hit <- list(mass=mass, name=paste0("TP", round(mass,4)))
    peaklist.hit$dppm <- 0
    peaklist.hit$profileID <- peaklist[n, "profile_ID"]
    peaklist.hit$mz <- peaklist[n, "mean_mz"]
    peaklist.hit$RT <-peaklist[n, "mean_RT"]
    peaklist.hit$int <- peaklist[n, "mean_int"]
    return(peaklist.hit)
  })
  
  hits.total <- do.call(rbind, hits)
  hits.total <- lapply(as.data.frame(hits.total), unlist)
  return(hits.total)
}