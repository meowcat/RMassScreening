#' Batch pick raw mzXML files
#' 
#' Peak picking using mostly enviPick. The function handles polarity switching raw files and parallel picking on multiple cores.
#' `batchPickDIA` additionally handles DIA files, including polarity switching DIA files, and picks every DIA scan separately.
#' 
#' Warning: `batchPick` doesn't actually use the `settings` at all! only batchPickDIA does.
#' 
#' @details 
#' `batchPick` picks either positive or negative mode data as specified. For a raw file `RAWFILE.mzXML`, the picked files (`writeData` output) have names
#' `RAWFILE.mzXML.MSlist.pos.RData` for positive mode, or `.neg.RData` for negative mode respectively.  For `writeList`, the files are `.csv` instead.
#' 
#' For `batchPickDIA`, the file names follow the pattern `RAWFILE.mzXML.MSlist.LEVEL-POLARITY-CENTER-WIDTH.RData`. For the MS1 scans, 
#' this is e.g. `1-pos--` because there is no isolation center or isolation width. For the MS2 scans, a positive scan with center 400 and isolation width
#' 200 would have `2-pos-400-200` in the name.
#' 
#' For legacy reasons, `batchPick` produces the objects `MSlist.pos` and `MSlist.neg` respectively, whereas `batchPickDIA` 
#' only makes `MSlist` regardless of polarity.
#' 
#' `batchPick` will soon be deprecated (or rather `batchPickDIA` will take over the name.)
#' 
#' If `multicore` is used, the processing messages ("output") are piped to the `logfile`. For single core, the output is displayed on screen.
#' 
#' @param files Array of paths to the raw files.
#' @param outputDir Folder to write picked file results to.
#' @param polarity For `batchPick`, `+`or `-`. Selects the polarity to pick. On one hand, this is needed to handle polarity switching files correctly; 
#' 	on the other hand, it also determines the file names to use downstream. 
#' @param writeData If `TRUE`, the picked file object (the `MSlist` from enviPick) is saved to disk
#' @param writeList If `TRUE`, the CSV picked peak table is saved to disk
#' @param settings Settings to use for enviPick. Read the settings file template for descriptions, or the enviPick documentation accordingly.
#' 		By default uses the settings as loaded.
#' @param multicore `FALSE` for single core processing, or an integer number for the number of cores to use. 
#' @param log Logfile output for multicore processing (irrelevant for single core).
#'  
#' @note `batchPickDIA` is the new name for `comprehensiveBatchPick`, which should not be used anymore. 
#' 
#' @return For `batchPickDIA`, returns a character vector of all scan types picked (the filename patterns for use in `fillProfiles` etc.).
#' @author stravsmi
#' @md
#' @export
batchPick <- function(files, outputDir, polarity = c("+", "-"), writeData = TRUE, writeList = TRUE,
                      settings=getOption("RMassScreening")$enviPick,
					  multicore = FALSE, log = "logCluster.txt")
{
	if(multicore)
	{
		# Initiate cluster
		cl <- makeCluster(multicore, outfile = log)
		
		clusterEvalQ(cl, library(enviPick))
		clusterEvalQ(cl, library(RMassScreening))
		
		tryCatch(
				parLapply(cl, files, function(file, outputDir, polarity, writeData, writeList, settings)
						{
							batchPick(file, outputDir, polarity, writeData, writeList, settings)
						},outputDir=outputDir, polarity=polarity, writeData=writeData, writeList=writeList, settings=settings)
				,finally=stopCluster(cl))
		return()
	}
	
  for(filepath.mzML in files)
  {
    pos <- ("+" %in% polarity)
    neg <- ("-" %in% polarity)
    message(filepath.mzML)
    message("Step 1 - reading", appendLF = TRUE)
    # (2) Initialize an MSlist object and load this .mzML file into it:
    f <- readMzXmlFile(filepath.mzML)
    if(pos)
    {
      f.pos <- f[which(unlist(lapply(f, function(scan) scan$metaData$polarity == "+")))]
      MSlist.pos <-readMSdata.direct(f.pos, MSlevel=c(1))
    }
    if(neg)
    {
      f.neg <- f[which(unlist(lapply(f, function(scan) scan$metaData$polarity == "-")))]
      MSlist.neg <-readMSdata.direct(f.neg, MSlevel=c(1))
    }
    
    # (3) Partition the measurements now available in MSlist:
    message("Step 2 - agglomeration", appendLF = FALSE)
    if(pos){
      MSlist.pos <- mzagglom(MSlist.pos,dmzgap=10,ppm=TRUE,drtgap=500,minpeak=4,maxint=1E7)
      message(" +", appendLF = FALSE)
    }
    if(neg){
      MSlist.neg <- mzagglom(MSlist.neg,dmzgap=10,ppm=TRUE,drtgap=500,minpeak=4,maxint=1E7)
      message(" -", appendLF = FALSE)
    }
    message("", appendLF = TRUE)
    
    
    message("Step 3 - clustering", appendLF = FALSE)
    if(pos){
      MSlist.pos<-mzclust(MSlist.pos,dmzdens=5,ppm=TRUE,drtdens=120,minpeak=4)      
      message(" +", appendLF = FALSE)
    }
    if(neg){
      MSlist.neg<-mzclust(MSlist.neg,dmzdens=5,ppm=TRUE,drtdens=120,minpeak=4)
      message(" -", appendLF = FALSE)
    }
    message("", appendLF = TRUE)
    
    
    message("Step 4 - picking", appendLF = FALSE)
    if(pos){
      MSlist.pos<-mzpick(MSlist.pos, minpeak = 4, drtsmall = 50, drtfill = 10,  drttotal = 200, recurs = 4, 
                         weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2)
      message(" +", appendLF = FALSE)
    }
    if(neg){
      MSlist.neg<-mzpick(MSlist.neg, minpeak = 4, drtsmall = 50, drtfill = 10,  drttotal = 200, recurs = 4, 
                         weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2)
      message(" -", appendLF = FALSE)
    }
    message("", appendLF = TRUE)
    
    
    message("Step 5 - exporting", appendLF = TRUE)
    
    if(writeData)
    {
      if(pos){
        MSlist <- MSlist.pos
        save(MSlist, file=paste0(outputDir, "/", basename(filepath.mzML), ".MSlist.pos.RData"))
      }
      if(neg){
        MSlist <- MSlist.neg
        save(MSlist, file=paste0(outputDir, "/", basename(filepath.mzML), ".MSlist.neg.RData"))  
      }
    }
    if(writeList)
    {
      if(pos)
      {
        writePeaklist(MSlist.pos,outputDir,
                      filename= paste0(basename(filepath.mzML), ".MSlist.pos.csv", sep = ""))
      }
      if(neg)
      {
        writePeaklist(MSlist.neg,outputDir,
                      filename= paste0(basename(filepath.mzML), ".MSlist.neg.csv", sep = ""))
      }
    }
    message("Done.")
    
  }
  
  
}



#' Call function with specified settings
#' 
#' Wrapper for formals() and do.call() which finds the arguments of \code{fun} and substitutes them
#' with the data in \code{settings} (but \code{settings} can have many more settings stored)
#' Arguments in \code{...} take precedence over those in \code{settings}
#'
#' @param fun (character) Function name to be called
#' @param settings Settings (options) list
#' @param ... Additional settings to substitute into \code{fun} that are not in \code{settings}
#'
#' @export
#'
#' @examples
#' #nope
.call.fun.settings <- function(fun, settings, ...)
{
  funFormals <- formals(fun)
  dotList <- list(...)
  funFormals[intersect(names(funFormals), names(settings))] <- settings[intersect(names(funFormals), names(settings))]
  funFormals[intersect(names(funFormals), names(dotList))] <- dotList[intersect(names(funFormals), names(dotList))]
  do.call(fun, funFormals)
}


#' @describeIn batchPick
#' 
#' @export
batchPickDIA <- function(files, outputDir, writeData = TRUE, writeList = TRUE, settings=getOption("RMassScreening")$enviPick,
                                   multicore = FALSE, log = "logCluster.txt")
{
  if(multicore)
  {
    # Initiate cluster
    cl <- makeCluster(multicore, outfile = log)
    
    clusterEvalQ(cl, library(enviPick))
    clusterEvalQ(cl, library(RMassScreening))
    
    tryCatch(
      runs <- parLapply(cl, files, function(file, outputDir, writeData, writeList, settings)
      {
        batchPickDIA(file, outputDir, writeData, writeList, settings)
      },outputDir=outputDir, writeData=writeData, writeList=writeList, settings=settings),
		error = function(e) {runs <- c()}
		
      ,finally=stopCluster(cl))
    return(unique(unlist(runs)))
  }
  
  for(filepath.mzML in files)
  {
    
    message(filepath.mzML)
    
    message("Step 1 - reading", appendLF = TRUE)
    # (2) Initialize an MSlist object and load this .mzML file into it:
    d <- readMzXmlFile(filepath.mzML)
    
    # read all scantypes from mzxml file headers
    scantype <- unlist(lapply(d, function(scan) 
      paste(scan$metaData$msLevel,
            ifelse(scan$metaData$polarity=="+", "pos", "neg"),
            scan$metaData$precursorMz, 
            scan$metaData$precursor$windowWideness,
            sep="-")))
    # define settings specific for every scantype
    allTypes <- unique(scantype)

    settings.sct <- list() # settings by scantype, named list by scantype name
    settings.general <- settings # global settings
    settings.general$scantypes <- NULL
    for(st in allTypes)
    {
      settings.sct[[st]] <- settings.general
    }
    
    if(!is.null(settings$scantypes))
    {
      types <- names(settings$scantypes)
      types.regex <- glob2rx(types)
      scanSettings <- sapply(types.regex, function(type.regex)
      {
        grepl(type.regex, allTypes)
      })
      rownames(scanSettings) <- allTypes
      for(st in allTypes)
      {
        settingsMatch <- sum(scanSettings[st,])
        # check if there is one and exactly one match for each scan type
        if(settingsMatch == 0)
           warning(paste0("Specific settings exist, but not for scantype ", st, "!"))
        else if(settingsMatch > 1)
          warning(paste0("More than one specific setting type matches for scantype ", st, "!"))
        else
        {
          # overwrite the settings that are scantype-specific, leave the other settings untouched.
          settingScanList <- types[[which(scanSettings[st,])]]
          settingScan <- settings$scantypes[[settingScanList]]
          # for(name in names(settingScan))
          settings.sct[[st]][names(settingScan)] <- settingScan[names(settingScan)]
        }
           
      }
    }


    # split the list of all scans (with data) by scantype
    d.tot <- split(d, scantype)
    # Make sure they are ordered as in the scantype list
    d.tot <- d.tot[allTypes]
    
    MSlist.tot <- lapply(seq_len(length(allTypes)),
                         function(n)
                         {
                           msLevel <- as.numeric(substr(allTypes[[n]],1,1))
                           readMSdata.direct(d.tot[[n]], MSlevel=msLevel)
                         })

    # (3) Partition the measurements now available in MSlist:
    message("Step 2 - agglomeration", appendLF = FALSE)
    MSlist.tot <- lapply(seq_len(length(d.tot)),
                         function(n)
                         {
                           MSlist <- .call.fun.settings("mzagglom", settings.sct[[allTypes[[n]]]], MSlist = MSlist.tot[[n]])
                           message(paste0(" ", names(d.tot)[[n]]), appendLF = FALSE)
                           return(MSlist)
                         })
    message("", appendLF = TRUE)
    
    
    message("Step 3 - clustering", appendLF = FALSE)
    MSlist.tot <- lapply(seq_len(length(d.tot)),
                         function(n)
                         {
                           MSlist <- .call.fun.settings("mzclust", settings.sct[[allTypes[[n]]]], MSlist = MSlist.tot[[n]])
                           message(paste0(" ", names(d.tot)[[n]]), appendLF = FALSE)
                           return(MSlist)
                         })
    message("", appendLF = TRUE)
    
    
    message("Step 4 - picking", appendLF = FALSE)
    MSlist.tot <- lapply(seq_len(length(d.tot)),
                         function(n)
                         {
                           MSlist <- .call.fun.settings("mzpick", settings.sct[[allTypes[[n]]]], MSlist = MSlist.tot[[n]])
                           message(paste0(" ", names(d.tot)[[n]]), appendLF = FALSE)
                           return(MSlist)
                         })
    message("", appendLF = TRUE)
    

    message("Step 5 - exporting", appendLF = TRUE)
    
    if(writeData)
    {
      res <- lapply(seq_len(length(d.tot)),
                    function(n)
                    {
                      MSlist <- MSlist.tot[[n]]
                      save(MSlist, file=paste0(outputDir, "/", basename(filepath.mzML), ".MSlist.", names(d.tot)[[n]],".RData"))
                      message(paste0(" ", names(d.tot)[[n]]), appendLF = FALSE)
                    })
    }
    if(writeList)
    {
        res <- lapply(seq_len(length(d.tot)),
                      function(n)
                      {
                        writePeaklist(MSlist.tot[[n]],outputDir,
                                      filename= paste0(basename(filepath.mzML), ".MSlist.", names(d.tot)[[n]],".csv", sep = ""))
                        message(paste0(" ", names(d.tot)[[n]]), appendLF = FALSE)
                      })
        
    }
    message("Done.")
    
  }
  return(names(d.tot))
  
}

#' @export
#' @describeIn batchPick
#' 
comprehensiveBatchPick <- function(...)
{
	message("This call is deprecated, the function is now called batchPickDIA!")
	batchPickDIA(...)
}