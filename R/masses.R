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


#' @export
consolidateReactions <- function(reactions)
{
  potential.u <- unique(reactions$mass)
  potential.opt <- unlist(lapply(potential.u, function(mz)
  {
    p.names <- reactions[reactions$mass == mz, "name"]
    p.name <- p.names[[which.min(countCharOccurrences("-", p.names))]]
  }))
  potential.clean <- data.frame(mass=potential.u, name=potential.opt, stringsAsFactors = FALSE)
  return(potential.clean)
}
