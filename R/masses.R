countCharOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

#' 
#' Expects two data frames with columns mass, name.
#' Note that the "reactions" dataframe will have dash-separated names,
#' such as "oh-deme"
#'
#' @export
combineReactions <- function(substances, reactions, sep=" ", omit.name="")
{
  par.peaks <- data.frame(mz = substances$mass, rt=0, into=0)
  rownames(par.peaks) <- substances$name
  
  # Consolidate reactions such that no mass is duplicate:
  # for multiply occurring reactions, the one with least "-" in the name is taken.
  reactions.comb <- consolidateReactions(reactions)
  
  # generate table for Jen's functions
  reactions.comb$active <- "x"
  reactions.comb$massdiff <- NULL
  for(i in 1:nrow(reactions.comb)){
    reactions.comb$massdiff[i] <- reactions.comb$mass[i]
  }
  
  act <- subset(reactions.comb,reactions.comb[,"active"]!="")
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


#' Combine reactions using molecular formula, not mass.
#' 
#' Expects two data frames with columns mass, name.
#' Note that the "reactions" dataframe will have dash-separated names,
#' such as "oh-deme"
#'
#' @export
combineReactions.formula <- function(substances, reactions, sep=" ", omit.name="")
{
  par.peaks <- data.frame(formula = substances$formula, rt=0, into=0)
  rownames(par.peaks) <- substances$name
  
  # Use only the rows for which there is a formula entry
  par.peaks <- par.peaks[!is.na(par.peaks$formula) & (par.peaks$formula != ""),]
  
  # Consolidate reactions such that no mass is duplicate:
  # for multiply occurring reactions, the one with least "-" in the name is taken.
  
  #reactions.comb <- consolidateReactions.formula(reactions)
  
  # generate table for Jen's functions
  reactions.comb <- reactions
  reactions.comb$active <- "x"

  act <- subset(reactions.comb,reactions.comb[,"active"]!="")
  act <- reactions.comb[!is.na(reactions.comb$formula)& (reactions.comb$formula != ""),]
  
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
  potential.comb <- convertFormulas(potential.comb)
  
  # Consolidate the reactions again such that multistep reactions are replaced by single-step reactions
  #potential.clean <- consolidateReactions(potential.comb)
  return(potential.comb)
}


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

convertFormulas <- function(substances)
{
  toConvert <- which(substances$formula != "")
  for(i in toConvert)
  {
    substances[i, "mass"] <- findMz.formula(substances[i, "formula"], "")$mzCenter
  }
  return(substances)
}

mergeSuspects <- function(massSuspects, formulaSuspects)
{
  allSuspects <- merge(massSuspects, formulaSuspects, all=TRUE, by="name", suffixes=c("",".f"))
  allSuspects$mass[!is.na(allSuspects$mass.f)]  <- allSuspects$mass.f[!is.na(allSuspects$mass.f)]
  allSuspects$mass.f <- NULL
  allSuspects
}