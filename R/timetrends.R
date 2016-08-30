
#' 
#' Generalized sample selector which hopefully works automatically like it is supposed to.
#' Takes a sampleList and a sampleAssignment.
#' 
#' The sampleList is the table merged to the filename list via assignSamples, and can contain any arbitrary columns
#' describing "what" the measured sample is. E.g. 
#' filename, sample, timepoint, workup
#' (where "workup" can be different ways of working up a biological sample, e.g. lysed cells vs supernatant)
#' 
#' The sampleAssignment table is the table building relationships between samples, and can contain arbitrary columns
#' describing what group a sample belongs to.
#' E.g.
#' sample, organism, spikedChemicals, stressor
#' 
#' sampleSelector returns the indices of all samples from sampleList which match all specified conditions (be they in 
#' sampleList or sampleAssignemnt). E.g.
#' 
#' sampleSelector(l, a, workup = "supernatant", spikedChemicals = "atrazine")
#' sampleSelector(l, a, species = c("ecoli", "botulinum"), spikedChemicals = "atrazine")
#' 
#' @export
sampleSelector <- function(sampleList, sampleAssignment, ...)
{
  params <- list(...)
  cond.sl <- intersect(colnames(sampleList), names(params))
  cond.sa <- intersect(colnames(sampleAssignment), names(params))
  # warn if there are inconsistencies
  if(length(intersect(cond.sl, cond.sa)) > 0)
    warning("Some conditions are used in both sample list and sample-group assignments!")
  par.unused <- setdiff(names(params), union(cond.sl, cond.sa))
  if(length(par.unused) > 0)
  {
    warning("Some conditions are unused:")
    warning(paste(par.unused, collapse=", "))
  }
  for(cond in cond.sa)
  {
    sampleAssignment <- sampleAssignment[
      sampleAssignment[,cond] %in% params[[cond]],
      ,drop=FALSE]
  }
  if(nrow(sampleList) == 0)
    error("sample list is empty")
  
  sampleList$INDICES <- 1:nrow(sampleList)
  sampleList <- sampleList[sampleList$sample %in% sampleAssignment$sample,,drop=FALSE]
  for(cond in cond.sl)
  {
    sampleList <- sampleList[
      sampleList[,cond] %in% params[[cond]],
      ,drop=FALSE]
  }
  return(sampleList$INDICES)
}

#' @export
summarizeProfiles <- function(profiles, sampleList, sampleIndices, groupName, groupBy = "time", groups = NA)
{
  # select all profiles from samples which are in the selected samples group
  profiles.ss <- profiles[[2]][profiles[[2]][,"sampleIDs"] %in% sampleList[sampleIndices, "sampleIDs"],,drop=FALSE]
  # remove all datapoints that are not associated to a profile
  profiles.ss <- profiles.ss[profiles.ss$profileIDs != 0,,drop=FALSE]
  
  #browser()
  # build a full combination of all sample table entries and all profileIDs so each mean and SD is computed from the same length vectors
  # and zerofill
  
  # Note: the original problem here was thet the zerofill rows didn't have the merged sample information from 
  # sample assignment, and therefore didn't merge with the correct groups!
  # (20160817; problem was noted around 201605 when the zerofill was moved to the profile filling step)
  
  profiles.complete <- merge(expand.grid(profileIDs = unique(profiles.ss$profileIDs),
                                          sampleIDs = unique(sampleList[sampleIndices,"sampleIDs"])),
                            profiles.ss, all.x=TRUE)
  profiles.split <- split(profiles.complete, is.na(profiles.complete$intensity))
  # fill "informative" colums for grouping:
  cols.sub <- intersect(colnames(profiles.complete), colnames(sampleList))
  profiles.split[["FALSE"]][,cols.sub] <- sampleList[match(profiles.split[["FALSE"]]$sampleIDs, sampleList$sampleIDs), cols.sub]
  profiles.complete <- rbind(profiles.split[["TRUE"]], profiles.split[["FALSE"]])
  
  profiles.complete[is.na(profiles.complete$intensity), "intensity"] <- 0
  
  # all sample IDs:
  #sampleIDs <- unique(sampleList[sampleIndices, ])
  
  # group by the summary variable
  groupFormula <- as.formula(paste("intensity ~ profileIDs + ", groupBy))
  
  tt.sd <- as.data.frame(unclass(xtabs(groupFormula, do.call(aggregate, list(groupFormula, data=profiles.ss,
                                                               FUN=sd)), na.action=na.pass)))
  tt.mean <- as.data.frame(unclass(xtabs(groupFormula, do.call(aggregate, list(groupFormula, data=profiles.ss,
                                                                             FUN=mean)))))  
#   
  #tt.sd <- as.data.frame(unclass(xtabs(intensity ~ profileIDs + time, aggregate(intensity ~ profileIDs + time, data=profiles.complete,
   #                                                                             FUN=sd))))
#   
#   tt.mean <- as.data.frame(unclass(xtabs(groupFormula, aggregate(eval(groupFormula), data=profiles.complete,
#                                                                                   FUN= mean))))
  
  # fill up the columns if some are missing and reorder correctly (note: timepoints over unselected list!)$
  if(all(is.na(groups)))
  {
    timepoints <- as.character(unique(sampleList[,groupBy]))
    timepoints <- sort(timepoints)
    timepoints <- timepoints[timepoints != ""]
  }
  else
    timepoints <- groups

  tt.sd[,timepoints[which(!(timepoints %in% colnames(tt.sd)))]] <- NA
  tt.mean[,timepoints[which(!(timepoints %in% colnames(tt.mean)))]] <- NA
  tt.sd <- tt.sd[,timepoints]
  tt.mean <- tt.mean[,timepoints]
  # prepare for merging
  colnames(tt.sd) <- paste0(timepoints, ".sd")
  tt.mean$profileIDs <- rownames(tt.mean)
  tt.sd$profileIDs <- rownames(tt.sd)
  #
  tt <- merge(tt.mean, tt.sd, by="profileIDs", suffixes = c("", ".sd"), all.x=TRUE)
  # here the later values are left out on purpose
  tt
}

#' @export
groupSummaries <- function(
  profiles, sampleList, sampleAssignment, sampleGroups, mock=FALSE, ...
  )
{
  # find columns to use for sample selector
  params <- colnames(sampleGroups)
  cond.sl <- intersect(colnames(sampleList), params)
  cond.sa <- intersect(colnames(sampleAssignment), params)
  cond.total <- union(cond.sl, cond.sa)
  
  tt.extracted <- lapply(1:nrow(sampleGroups), function(ngroup)
  {
    message(sampleGroups[ngroup, "sampleGroup"])
    cond <- sampleGroups[ngroup, cond.total]
    
    cond <- as.vector(as.matrix(cond))
    cond <- strsplit(cond, ",")
    names(cond) <- cond.total
    
    arglist <- list(sampleList=sampleList, sampleAssignment=sampleAssignment)
    arglist <- c(arglist, cond)
    sampleIndices <- do.call(sampleSelector, arglist)
    
    if(!mock)
      summarizeProfiles(profiles, sampleList, sampleIndices, sampleGroups[ngroup, "sampleGroup"], ...)
    else
      sampleList[sampleIndices,,drop=FALSE]
  })
  names(tt.extracted) <- sampleGroups$sampleGroup
  return(tt.extracted)
}

#' @export
mergeGroups <- function(profiles, sampleGroups, summaries)
{
  tt.total <- data.frame( profileIDs = profiles[[7]][,"profile_ID"], stringsAsFactors = FALSE )
  for(ngroup in 1:nrow(sampleGroups))
  {
    colnames(summaries[[ngroup]]) <-   paste(colnames(summaries[[ngroup]]), sampleGroups[ngroup, "sampleGroup"], sep=".")
    gname <- sampleGroups[ngroup, "sampleGroup"]
    tt.total <- qmerge(tt.total, summaries[[ngroup]], all.x=TRUE, by.x = "profileIDs", by.y = paste0("profileIDs.",gname),
                       suffixes = c("", paste0(".", gname)))
  }
  tt.total[is.na(tt.total)] <- 0
  tt.total
}