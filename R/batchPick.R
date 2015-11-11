#' @export
batchPick <- function(files, outputDir, polarity = c("+", "-"), writeData = TRUE, writeList = TRUE)
{
  
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
    message("", appendLF = TRUE)
    
    
    
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
        save(MSlist, file=paste0(outputDir, basename(filepath.mzML), ".MSlist.pos.RData"))
      }
      if(neg){
        MSlist <- MSlist.neg
        save(MSlist, file=paste0(outputDir, basename(filepath.mzML), ".MSlist.neg.RData"))  
      }
    }
    if(writeList)
    {
      if(pos)
      {
        writePeaklist(MSlist.pos,"results",
                      filename= paste0(outputDir, ".MSlist.pos.csv", sep = ""))
      }
      if(neg)
      {
        writePeaklist(MSlist.neg,"results",
                      filename= paste0(outputDir, ".MSlist.neg.csv", sep = ""))
      }
    }
    message("Done.")
    
  }
  
  
}

