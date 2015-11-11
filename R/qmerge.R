
# profiles[[2]] <- merge(profiles[[2]], sampleAssigned, by="sampleIDs", all.x=TRUE, sort=FALSE)

# Merge function which deals with simple tables (i.e. they have 0-1 y datasets per x dataset and only one column to merge by.)
# 

#' @export
qmerge <- function(x, y, by = NA, by.x = by, by.y = by, all.x = FALSE, suffixes = c(".x",".y"))
{
  # matching
  x.in.y <- match(x[,by.x], y[,by.y])
  # locate column names
  x.cols <- colnames(x)
  y.cols <- colnames(y)
  x.indexcol <- match(by.x, x.cols)
  y.indexcol <- match(by.y, y.cols)
  # suffix column names correctly
  x.cols.edit <- which(x.cols %in% y.cols) 
  y.cols.edit <- which(y.cols %in% x.cols) 
  x.cols[x.cols.edit] <- paste0(x.cols[x.cols.edit], suffixes[[1]])
  y.cols[y.cols.edit] <- paste0(y.cols[y.cols.edit], suffixes[[2]])
  # rename x index column back to original name
  x.cols[x.indexcol] <- by.x
  colnames(x) <- x.cols
  colnames(y) <- y.cols
  # reduce matched set. 
  y.matched <- y[x.in.y,,drop=FALSE]
  # remove y index column
  y.matched <- y.matched[,-y.indexcol]
  # bind and if necessary remove NA hits
  x <- cbind(x, y.matched)
  if(!all.x)
    x <- x[!is.na(x.in.y),,drop=FALSE]
  x
}