
#' @export
getPatterns <- function(profiles)
{
  peaklist <- as.data.frame(profiles$index_prof[,c(14,16,15)])
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
