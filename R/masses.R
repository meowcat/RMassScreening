#' @import RMassBank
NULL

countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

#' Combinatorially apply transformations
#' 
#' `combineReactions` applies `reactions` specified as mass differences to `substances` specified as exact masses.
#' `combineReactions.formula` applies `reactions` specified as formula differences (e.g. `C2H-2O1`) to
#'  `substances` specified as chemical formulas.
#' 
#' For `combineReactions.formula`, all combinations where one of two formulas is not specified are omitted.
#' 
#' @param substances Data frame with column `name` and mandatory column `mass` (`combineReactions`) or `formula`  (`combineReactions.formula`).
#' @param reactions Data frame with column `name` and mandatory column `mass` (`combineReactions`) or `formula`  (`combineReactions.formula`).
#' @param sep Separator between combined reaction names (e.g. `oh` and `deme` with separator `-` get combined to `oh-deme`)
#' @param omit.name A reaction name which should be omitted in the combined reaction, such as "parent"
#' 
#' Expects two data frames with columns mass, name. Combines 
#'
#' @md
#' @export
combineReactions <- function(substances, reactions, sep=" ", omit.name="")
{
  par.peaks <- data.frame(mz = substances$mass, rt=0, into=0)
  rownames(par.peaks) <- substances$name
  
  # Consolidate reactions such that no mass is duplicate:
  # for multiply occurring reactions, the one with least "-" in the name is taken.
  reactions.comb <- consolidateReactions(reactions, sep=sep)
  
  # generate table for Jen's functions
  reactions.comb$massdiff <- NULL
  for(i in 1:nrow(reactions.comb)){
    reactions.comb$massdiff[i] <- reactions.comb$mass[i]
  }
  
  x <- act[,"name"]
  x[x==omit.name] <- ""
  
  
  # calculate!
  options <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      if(x[k] != "")
        options[j,k] <- paste(rownames(par.peaks[j,]), x[k], sep=sep)
      else
        options[j,k] <- rownames(par.peaks[j,])
    }
  }
  
  option.masses <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  rownames(option.masses) <- rownames(par.peaks)
  colnames(option.masses) <- act[,"name"]
  
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      option.masses[j,k] <- par.peaks[j,"mz"] + act[k,"massdiff"]}
  }
  
  potential <- as.vector(as.matrix(option.masses))
  names(potential) <- as.vector(as.matrix(options))
  
  potential.comb <- data.frame(mass=potential, name=names(potential), stringsAsFactors = FALSE)
  
  # Consolidate the reactions again such that multistep reactions are replaced by single-step reactions
  potential.clean <- consolidateReactions(potential.comb)
  return(potential.clean)
}


#' @describeIn combineReactions
#' 
#' @export
combineReactions.formula <- function(substances, reactions, sep=" ", omit.name="", calculateMasses = TRUE, validate=FALSE)
{
  par.peaks <- data.frame(formula = substances$formula, rt=0, into=0, stringsAsFactors = FALSE)
  rownames(par.peaks) <- substances$name
  
  # this is a safeguard for the "parent" reaction where mass is zero, there the formula should be preserved
  # otherwise the "parents" are lost if only formula but not mass is used for combining reactions
  if("mass" %in% colnames(substances))
    par.peaks[substances$mass == 0, "formula"] <- "C0"
  if("mass" %in% colnames(reactions))
    reactions[reactions$mass == 0, "formula"] <- "C0"
  
  # Use only the rows for which there is a formula entry
  par.peaks <- par.peaks[!is.na(par.peaks$formula) & (par.peaks$formula != ""),]
  
  # Consolidate reactions such that no mass is duplicate:
  # for multiply occurring reactions, the one with least "-" in the name is taken.
  
  #reactions.comb <- consolidateReactions.formula(reactions)
  
  # generate table for Jen's functions
  
  act <- reactions[!is.na(reactions$formula)& (reactions$formula != ""),]
  
  x <- act[,"name"]
  x[x==omit.name] <- ""
  
  
  
  
  # calculate!
  options <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      if(x[k] != "")
        options[j,k] <- paste(rownames(par.peaks[j,]), x[k], sep=sep)
      else
        options[j,k] <- rownames(par.peaks[j,])
    }
  }
  
  option.masses <- data.frame(matrix(NA, ncol = nrow(act), nrow = nrow(par.peaks)))
  rownames(option.masses) <- rownames(par.peaks)
  colnames(option.masses) <- act[,"name"]
  
  for(k in 1:nrow(act)){
    for(j in 1:nrow(par.peaks)) {
      option.masses[j,k] <- add.formula(par.peaks[j,"formula"], act[k,"formula"])
    }}
  
  potential <- as.vector(as.matrix(option.masses))
  names(potential) <- as.vector(as.matrix(options))
  
  potential.comb <- data.frame(formula=potential, name=names(potential), stringsAsFactors = FALSE)
  
  # validation: check if any element is below zero
  if(validate)
  {
    invalid <- laply(potential.comb$formula, function(fo)
      any(unlist(formulastring.to.list(fo)) < 0) )
    potential.comb <- potential.comb[!invalid,,drop=FALSE]
  }
  
  if(calculateMasses)
    potential.comb <- convertFormulas(potential.comb)
  
  # Consolidate the reactions again such that multistep reactions are replaced by single-step reactions
  #potential.clean <- consolidateReactions(potential.comb)
  return(potential.comb)
}


#' Merge reactions with identical mass differences
#' 
#' Merges reactions with identical mass difference within `eps` to avoid creating duplicate reaction products.
#' Of multiple reactions, the "simplest" one will be kept, i.e. if one is already combined (contains a separator as specified
#' by `sep` and the other is not, the "shorter" path will be kept. 
#' 
#' @param reactions Data frame with `name` and `mass`.
#' @param eps Merging margin, absolute
#' @param sep Reaction separator token to filter/count for priority
#' 
#' @author stravsmi
#' @md
#' 
#' @export
consolidateReactions <- function(reactions, eps = 0.0002, sep="-")
{
  n <- 1
  while(n < nrow(reactions))
  {
    mz <- reactions$mass[[n]]
    equivalent <- which(abs(reactions$mass- mz) <= eps)
    rNames <- reactions$name[equivalent]
    if(length(rNames) > 1)
    {
      rName <- which.min(countCharOccurrences(sep, rNames))
      reactions <- reactions[-equivalent[-rName],,drop=FALSE]
    }
    n <- n+1
  }
  return(reactions)
}

#' Convert formulas to masses
#' 
#' Calculates the monoisotopic exact masses from formulas. 
#' 
#' @param substances a data.frame with column `formula` and existing column `mass`.
#' @return Data frame with overwritten column `mass`.
#' 
#' @author stravsmi
#' @md
#' @export
convertFormulas <- function(substances)
{
  toConvert <- which(substances$formula != "")
  for(i in toConvert)
  {
    substances[i, "mass"] <- findMz.formula(substances[i, "formula"], "")$mzCenter
  }
  return(substances)
}

#' Merge suspect lists
#' 
#' Merges suspect lists generated from formula calculations and from mass calculations.
#' If formula results are present, they take priority over mass results. Note that formula results must first be
#' converted to masses using [convertFormulas].
#' 
#' @param massSuspects data frame with `name`,`mass` columns
#' @param formulaSuspects data frame with `name`,`mass` columns 
#' 
#' @author stravsmi
#' @export
mergeSuspects <- function(massSuspects, formulaSuspects)
{
  allSuspects <- merge(massSuspects, formulaSuspects, all=TRUE, by="name", suffixes=c("",".f"))
  allSuspects$mass[!is.na(allSuspects$mass.f)]  <- allSuspects$mass.f[!is.na(allSuspects$mass.f)]
  allSuspects$mass.f <- NULL
  allSuspects
}