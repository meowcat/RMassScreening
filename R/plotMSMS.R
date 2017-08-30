# Various code snippets for plotting MSMS spectra
# E. Schymanski, 11/1/2016, copied from MetFragResultPlotting.R

# plot spectrum with formula annotation
plotSpectrum <- function(mz, int, main, smiles="", labels=NULL, formatFormula=FALSE,
                         ylim_factor=2,max_xlim=0,kekulise=TRUE,color="black", ySpacer = 0.05, srt=45, offset=-0.4,
                         htt = FALSE, htt.rev = FALSE) {
  
  #int <- 1000*int/max(int) 
  rev.factor <- ifelse(htt.rev, -1, 1)
  
  if (max_xlim == 0) {
    max_xlim <- max(mz)*1.2
  }
  abs_ylim <- (max(int)*ylim_factor)
  if(htt)
    ylim <- c(-abs_ylim, abs_ylim)
  else
    ylim <- c(0, abs_ylim)
  
  ySpacer.abs <- ySpacer * abs_ylim
  
  par(yaxs="i")
  if(!htt.rev)
    plot.new()
  plot.window(xlim=c(0, max_xlim), ylim=ylim)
  par(mar=c(3,2,1,1))
  #title(xlab = "m/z", ylab="intensity")
  
  lines(mz, rev.factor*int, type="h", 
       #font=2, #bold the font
       #cex=2, #magnifies text
       cex=1,
       lwd = 2, #line width
       col=color)
  
  if(!htt.rev)
    axis(1)
  axis(2, at=rev.factor * c(0, 250, 500, 750, 1000), labels=c("0%", "25%", "50%", "75%", "100%"))
  #box(bty="l")
  abline(h=0)
  
  #text(mz, int, labels=labels,adj=c(-0.03,0.5),srt=90)#pos=3,offset=0.2, adj=c(0,0))
  if (!is.null(labels) && formatFormula) {
    text(mz, rev.factor*(int + ySpacer.abs), labels=chemistry2expression(labels),srt=srt,pos=4,offset=offset)
  } else if (!is.null(labels)) {
    for(j in seq_along(mz))
    {
      text(mz[[j]], rev.factor * (int[[j]] + ySpacer.abs), labels[[j]], srt=srt,pos=4,offset=offset)
    }
    #text(mz, int + ySpacer.abs, labels=labels,srt=srt,pos=4,offset=offset)
  } 
  #text(mz, int, labels=labels,srt=90,pos=4,offset=0)
  #text(mz, int, labels=labels,pos=3,offset=0.2)
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise=kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol, width = 150, height = 150))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not plotted: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol, width = 150, height = 150))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    # img <- view.image.2d(mol, width = 150, height = 150)
    if (length(img)>2) {
      rasterImage(img, 0,1000, 0.3*max(mz),(max(int)*2.5))
    }
  }
  invisible(NULL)
}



#plot multiple MSMS spectra
plotSpectra <- function(mz, int, main, labels=NULL,labelSpec=1,mz_a=NA,int_a=NA,
                        mz_b=NA,int_b=NA,smiles="",kekulise=TRUE,legend=c("spec1")) {
  #normalise the peak lists before plotting!!!
  n_specs <- 1
  int <- 1000*int/max(int)  
  plot(mz, int, type="h",
       #if(length(mz)>1) {}
       xlim=c(0,max(mz)*1.1), 
       ylim=c(0,max(int)*2.5),
       main=main,
       cex=1, lwd=2, col="black",
       xlab="m/z", ylab="Intensity"
  )
  if (!is.null(labels)&&labelSpec==1) {
    text(mz, int, labels=labels,srt=srt,pos=4,offset=-0.2)
  }
  colour <- c("black")
  line_type <- c("solid")
  if (!is.na(mz_a[1])) {
    int_a <- 1000*int_a/max(int_a)
    lines(mz_a, int_a, type="h",col= "red", lty="dashed",lwd=2)
    n_specs <- n_specs+1
    if (length(legend) < n_specs) legend[n_specs] <- "spec2"
    colour[n_specs] <- "red"
    line_type[n_specs] <- "dashed"
    if (!is.null(labels)&&labelSpec==2) {
      text(mz_a, int_a, labels=labels,srt=srt,pos=4,offset=-0.2)
    }
  }
  if (!is.na(mz_b[1])) {
    int_b <- 1000*int_b/max(int_b)
    lines(mz_b, int_b, type="h",col= "green", lty="dotdash",lwd=1)
    n_specs <- n_specs+1
    if (length(legend) < n_specs) legend[n_specs] <- "spec3"
    colour[n_specs] <- "green"
    line_type[n_specs] <- "dotdash"
    if (!is.null(labels)&&labelSpec==3) {
      text(mz_b, int_b, labels=labels,srt=srt,pos=4,offset=-0.2)
    }
  }
  legend("topright", legend=legend,
         col=colour, lty=line_type, lwd=2, bty="n")
  #plot the molecule onto the spectrum
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise = kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol, width = 150, height = 150))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not plotted: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol, width = 150, height = 150))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    # img <- view.image.2d(mol, width = 150, height = 150)
    if (length(img)>2) {
      rasterImage(img, 0,1000, 0.3*max(mz),(max(int)*2.5))
    }
    #img <- view.image.2d(mol, width = 150, height = 150)
    #rasterImage(img, 0,1000, 75,1500)
  }
  
  #mol <- load.molecules(correct_smiles, aromaticity=TRUE)#, type="smi")
  #aro_mol <- do.aromaticity(mol)
  #mol <- do.typing(mol) #neither of these work
  #get.smiles(mol) #so these smiles look totally different and still have aromaticity dramas
  #minor.tick(nx=10, ny=10, tick.ratio=0.3)
  invisible(NULL)
}

chemistry2expression <- function(formulas, mz="") {
  exprs <- gsub("([\\+\\-])", "^'\\1'", formulas, fixed=FALSE)
  exprs <- gsub("([0-9]+)", "[\\1]*", exprs, fixed=FALSE)
  exprs <- sub("\\*$", "", exprs)
  exprs <- sub("\\*\\^", "\\^", exprs)
  #quote(exprs)
  parse(text=paste0("textstyle(atop(", mz, ",", exprs, "))"))
}

mzformula2expression <- function(mz, formula)
{
  expression(atop(mz,deparse(chemistry2expression(formula))))
}


plotMS <- function(mz, int, main, labels=NULL, formatFormula=TRUE,
                         absInt=FALSE,ylim_factor=2.5) {
  if (!absInt) {
    int <- 1000*int/max(int) 
  }
  plot(mz, int, type="h", 
       xlim=c((min(mz)-5),(max(mz)+5)),
       ylim=c(0,(max(int)*ylim_factor)),
       main=main,
       #font=2, #bold the font
       #cex=2, #magnifies text
       cex=1,
       lwd = 2, #line width
       col="black",
       xlab="m/z", ylab = "Intensity")
  #text(mz, int, labels=labels,adj=c(-0.03,0.5),srt=90)#pos=3,offset=0.2, adj=c(0,0))
  if (!is.null(labels) && formatFormula) {
    text(mz, int, labels=chemistry2expression(labels),srt=srt,pos=4,offset=-0.2)
  } else if (!is.null(labels)) {
    text(mz, int, labels=labels,srt=srt,pos=4,offset=-0.2)
  } 
  #text(mz, int, labels=labels,srt=90,pos=4,offset=0)
  #text(mz, int, labels=labels,pos=3,offset=0.2)
  invisible(NULL)
}

addStructToPlot <- function(smiles,kekulise=FALSE,coords=c(0,1000,100,2500)) {
  if (nchar(smiles)>1) {
    mol <- parse.smiles(smiles,kekulise=kekulise)[[1]]
    img <- tryCatch({
      (view.image.2d(mol, width = 150, height = 150))
    }, error = function(e) {
      img <- ""
      print(paste("Invalid SMILES not plotted: ", smiles, sep=""))
    })
    if (length(img)<=2 && kekulise) {
      mol <- parse.smiles(smiles,kekulise=FALSE)[[1]]
      img <- tryCatch({
        (view.image.2d(mol, width = 150, height = 150))
      }, error = function(e) {
        img <- ""
        print(paste("Invalid SMILES not plotted without kekulise either: ", smiles, sep=""))
      })
    }
    # img <- view.image.2d(mol, width = 150, height = 150)
    if (length(img)>2) {
      rasterImage(img, coords[1],coords[2],coords[3],coords[4])
    }
  }
}
# 
# plot(1,1)
# text(1,1, chemistry2expression("C123H4", 123))
# 
# eval(substitute(
#   expression(textstyle(atop(
#     mz, 
#     formula))
#   ), list(mz=123, formula=chemistry2expression("C6H12O6")))
#   )
# )


# MS plotting function extracted from the data for Biotransformation paper (Stravs et al. 2017)
# Some code is originally from Emma Schymanski
plotTP <- function(sp, SMILES="", main = "", readAnnotation = TRUE, ylim_factor=1.5, max_xlim=0, htt= FALSE, htt.rev = FALSE, ...)
{
	
	merged_peaks <- getData(sp)
	rev.factor <- ifelse(htt.rev, -1, 1)
	
	if ((length(sp@formula)>=1) && (readAnnotation == TRUE)) {
		indexes <- which(merged_peaks$intrel>50) # take all but the smallest
		indexes2 <- which(merged_peaks$mz>100) # take all above mz100, because they dont' get overwritten by struct
		# for lower mz, take those that are below a certain relInt and mz
		indexes3 <- which(merged_peaks$intrel>50) 
		indexes4 <- which(merged_peaks$mz<100)
		indexes5 <- intersect(indexes3, indexes4)
		indexes5 <- intersect(indexes5,indexes) # but don't accept the really small ones
		maxIntpeaks <- order(merged_peaks$intrel,decreasing = TRUE)
		if (length(maxIntpeaks)>10) {
			maxIntpeaks <- maxIntpeaks[1:10]
		} 
		index_final <- intersect(indexes, indexes2) # high intense, high mz
		index_final <- intersect(index_final, maxIntpeaks) # and of those, only the max intense ones
		index_final <- union(index_final,indexes5) # plus some low intense low mz
		index_final <- sort(index_final) # puts these back in index order
		if (length(index_final)<2) {
			index_final <- sort(maxIntpeaks)
		}
		anno_peaks <- merged_peaks[index_final,] # these are the ones to be annotated
		mz1 <- 0
		mz2 <- 0
		relInt1 <- 0
		relInt2 <- 0
		index_rem <- ""
		n_rem <- 0
		if (length(anno_peaks[,1])>2) {
			for (j in 1:(length(anno_peaks[,1])-1)) {
				mz1 <- anno_peaks$mz[j]
				mz2 <- anno_peaks$mz[j+1]
				relInt1 <- anno_peaks$intrel[j]
				relInt2 <- anno_peaks$intrel[j+1]
				if ((mz2-mz1)<10 && (relInt1 > relInt2)) {
					n_rem <- n_rem+1
					index_rem[n_rem] <- j+1
				} else if ((mz2-mz1)<10 && (relInt1 < relInt2)) {
					n_rem <- n_rem+1
					index_rem[n_rem] <- j
				}
			}
		}
		if (is.na(as.numeric(index_rem))) {
			anno_peaks2 <- anno_peaks
		} else {
			anno_peaks2 <- anno_peaks[-as.numeric(index_rem),]
		}
		# plot the spectrum
		
		# construct labels:
		annoLabels <- apply(anno_peaks2, 1, function(pk)
				{
					chemistry2expression(pk[["formula"]],sprintf("%.04f", as.numeric(pk[["mz"]])))
				})
		#browser()
		
		plotSpectrum(mz=anno_peaks2$mz,int=anno_peaks2$intrel,main=main,
				labels=annoLabels,
				#anno_peaks2$formula,
				smiles=SMILES,formatFormula = FALSE, ylim_factor = ylim_factor, 
				max_xlim=max_xlim, htt=htt, htt.rev=htt.rev, ...)
		lines(merged_peaks$mz,rev.factor *merged_peaks$intrel,type="h",lwd=2,col="black")
		# write the peaklists to another directory
		# main <- sub("?","",main,fixed=TRUE)
		#anno_pl_path <- paste(Directory,"/annotated_peaklists/",spec_name,".txt",sep="")
		#write.table(merged_peaks[,c(1,2,6,7)],file=anno_pl_path, quote=FALSE,row.names=FALSE)
	} else {
		
		indexes <- which(merged_peaks$intrel>50)
		
		annoLabels <- apply(merged_peaks[indexes,,drop=FALSE], 1, function(pk)
				{
					sprintf("%.04f", as.numeric(pk[["mz"]]))
				})
		
		plotSpectrum(mz=merged_peaks$mz[indexes],int=merged_peaks$intrel[indexes],main=main,smiles=SMILES,
				labels=annoLabels,formatFormula = FALSE, ylim_factor = ylim_factor, max_xlim=max_xlim
				, htt=htt, htt.rev=htt.rev, ...)
		lines(merged_peaks$mz,rev.factor * merged_peaks$intrel,type="h",lwd=2,col="black")
		# write the peaklists to another directory
		#anno_pl_path <- paste(Directory,"/annotated_peaklists/",spec_name,".txt",sep="")
		#write.table(merged_peaks[,c(1,2,1)],file=anno_pl_path, quote=FALSE,row.names=FALSE)
	}
}

