#' @export
mcBatchPick <- function(files, outputDir, polarity = c("+", "-"), writeData = TRUE, writeList = TRUE,
                        log="logCluster.txt", no_cores = 6)
{
  # Initiate cluster
  cl <- makeCluster(no_cores, outfile = log)
  
  clusterEvalQ(cl, library(enviPick))
  
  #   argNames <- c("files", "outputDir", "polarity", "writeData", "writeList")
  #   clusterExport(cl, argNames)
  
  # More ram-saving version: process pos and neg separately
  tryCatch(
  parLapply(cl, files, function(filepath.mzML, outputDir, polarity, writeData, writeList)
  {
    pos <- ("+" %in% polarity)
    neg <- ("-" %in% polarity)
    message(paste0(filepath.mzML,": 0 - reading"))
    # (2) Initialize an MSlist object and load this .mzML file into it:
    f <- readMzXmlFile(filepath.mzML)
    if(pos)
    {
      message(paste0(filepath.mzML,": 1 - initializing +"))
      f.pos <- f[which(unlist(lapply(f, function(scan) scan$metaData$polarity == "+")))]
      MSlist.pos <-readMSdata.direct(f.pos, MSlevel=c(1))
      rm(f.pos)
    }
    if(neg)
    {
      message(paste0(filepath.mzML,": 1 - initializing -"))
      f.neg <- f[which(unlist(lapply(f, function(scan) scan$metaData$polarity == "-")))]
      MSlist.neg <-readMSdata.direct(f.neg, MSlevel=c(1))
      rm(f.neg)
    }
    rm(f)
    if(pos)
    {
      MSlist <- MSlist.pos
      rm(MSlist.pos)
      # (3) Partition the measurements now available in MSlist:
      message(paste0(filepath.mzML,": 2 - agglomeration +"))
      MSlist <- mzagglom(MSlist,dmzgap=10,ppm=TRUE,drtgap=500,minpeak=4,maxint=1E7)
      message(paste0(filepath.mzML,": 3 - clustering +"))
      MSlist<-mzclust(MSlist,dmzdens=5,ppm=TRUE,drtdens=120,minpeak=4)      
      message(paste0(filepath.mzML,": 4 - picking +"))
      MSlist<- mzpick(MSlist, minpeak = 4, drtsmall = 50, drtfill = 10,  drttotal = 200, recurs = 4, 
                           weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2)
      message(paste0(filepath.mzML,": 5 - exporting +"))
      if(writeData)
        save(MSlist, file=paste0(outputDir, "/", basename(filepath.mzML), ".MSlist.pos.RData"))
      if(writeList)
        writePeaklist(MSlist,outputDir,
                      filename= paste0(basename(filepath.mzML), ".MSlist.pos.csv", sep = ""))
      message(paste0(filepath.mzML,": done +"))
    }
    if(neg)
    {
      MSlist <- MSlist.neg
      rm(MSlist.neg)
      # (3) Partition the measurements now available in MSlist:
      message(paste0(filepath.mzML,": 2 - agglomeration -"))
      MSlist <- mzagglom(MSlist,dmzgap=10,ppm=TRUE,drtgap=500,minpeak=4,maxint=1E7)
      message(paste0(filepath.mzML,": 3 - clustering -"))
      MSlist<-mzclust(MSlist,dmzdens=5,ppm=TRUE,drtdens=120,minpeak=4)      
      message(paste0(filepath.mzML,": 4 - picking -"))
      MSlist<- mzpick(MSlist, minpeak = 4, drtsmall = 50, drtfill = 10,  drttotal = 200, recurs = 4, 
                      weight = 2, SB = 3, SN=2, minint = 1E4, maxint = 1e+07, ended = 2)
      message(paste0(filepath.mzML,": 5 - exporting -"))
      if(writeData)
        save(MSlist, file=paste0(outputDir, "/", basename(filepath.mzML), ".MSlist.neg.RData"))
      if(writeList)
        writePeaklist(MSlist,outputDir,
                      filename= paste0(basename(filepath.mzML), ".MSlist.neg.csv", sep = ""))
      message(paste0(filepath.mzML,": done -"))
    }
    
  }, outputDir=outputDir, polarity=polarity, writeData=writeData, writeList=writeList)
  ,finally=stopCluster(cl))
}


#' Multicore batch pick
#'
#' @param files 
#' @param outputDir 
#' @param writeData 
#' @param writeList 
#' @param log 
#' @param no_cores 
#'
#' @return
#' @export
#'
#' @examples
mcComprehensiveBatchPick <- function(files, outputDir, writeData, writeList,
                                             log="logCluster.txt", no_cores = 6)
{
  
  
  # Initiate cluster
  cl <- makeCluster(no_cores, outfile = log)
  
  clusterEvalQ(cl, library(enviPick))
  clusterEvalQ(cl, library(RMassScreening))
  
  tryCatch(
    parLapply(cl, files, function(file, outputDir, writeData, writeList)
    {
      comprehensiveBatchPick(file, outputDir, writeData, writeList)
    },outputDir=outputDir, writeData=writeData, writeList=writeList)
    ,finally=stopCluster(cl))
}
    