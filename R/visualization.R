#library(shiny)
#library(shinyAce)
#' @export
erb <- function(x,y,sd,eps,...)
{
  segments(x,y-sd,x,y+sd,...)
  segments(x-eps,y-sd,x+eps,y-sd,...)
  segments(x-eps,y+sd,x+eps,y+sd,...)
}

#' Convert scantype string to readable output
#' 
#' @param x input characters vector 
#' @return Converted vector
#' 
#' @author stravsmi
#' @export
readScantype <- function(x)
{
	xsplit <- strsplit(x, "-")
	xsplit <- lapply(xsplit, function(xx) {
				length(xx) <- 4
				xx
			})
	xtab <- do.call(rbind,xsplit)
	xdf <- data.frame(msLevel = as.integer(xtab[,1]),
				polarity = xtab[,2],
				center = as.numeric(xtab[,3]),
				window = as.numeric(xtab[,4]) )
	xdf$min <- xdf$center - xdf$window/2
	xdf$max <- xdf$center + xdf$window/2
	xdf$string <- sprintf("MS%d %s", xdf$msLevel, xdf$polarity)
	xdf$ms2string <- sprintf("MS%d %s [%.0f - %.0f]", xdf$msLevel, xdf$polarity, xdf$min, xdf$max)
	xdf$string[xdf$msLevel == 2] <- xdf$ms2string[xdf$msLevel==2]
	xdf
}

#' Run the results viewer
#'
#' @param totalTable 
#'
#' @param hits 
#' @param timepoints 
#' @param sampleGroups 
#' @param patterns 
#' @param spectra msmsWorkspace containing extracted spectra for profiles
#' @param addcols Number of additional columns per sample group present in summarized (totalTable) output. This is 1 if regular processing (adding the mean) is used. 
#' @param ... 
#'
#' @export
runViewer <- function(totalTable, hits, timepoints, sampleGroups, patterns = NULL, spectra = NULL, addcols = 1, settings = getOption("RMassScreening"), 
                      profiles = NULL, files = NULL, ...)
{
  fe <- environment()
  values <- reactiveValues()
  values$pinList <- c()
  
  hitsLimit <- settings$viewer$hitsLimit
  if(is.null(hitsLimit))
  {
    warning("Parameter hitsLimit not set. Defaulting to 2000")
    hitsLimit <- 2000
  }
  
  tt.hits.all <- merge(totalTable, hits, by.x="profileIDs", by.y="profileID")
  tt.hits <- tt.hits.all[min(
    seq_len(hitsLimit),
    nrow(tt.hits.all)),]
  
  
#   filterText <- "
# #tt.hits <- tt.hits[tt.hits$t0.mixed < 1e5,]
# tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN)
#                            > (tt.hits$mean.synCTL + tt.hits$mean.chlCTL + tt.hits$mean.mcyCTL) * 10,]
# tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN) > tt.hits$mean.WC * 7,]
# tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN) > 3e5,]
# 
# #tt.hits <- tt.hits[grepl(\"AZY\",tt.hits$name),]
# 
# tt.hits <- tt.hits[order(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN, decreasing = TRUE),]
# 
# "
  
	rv <- reactiveValues(cpd = NA, spectrum = NA)
  
  profileNames <- paste0(tt.hits$name, " (", round(tt.hits$mz,4), "/", round(tt.hits$RT / 60,1),")")
  profileList <- tt.hits$profileIDs
  names(profileList) <- profileNames
  
  updateFilterList <- function(se, session)
  {
    if(length(se$filterList) > 0)
    {
      afList <- 1:length(se$filterList)
      names(afList) <- unlist(lapply(se$filterList, renderFilter))
    } 
    else
      afList <- c()
    updateSelectInput(session, "activeFilters", choices=afList)
  }
  
  server <- function(input, output, session) {
    
    filterList <- list()
    se <- environment()
    
    observe({
      f <- input$profileFilter
      if(f == "")
        updateSelectInput(session, "profile", choices=profileList)
      else
      {
        profileList.m <- which(grepl(f, names(profileList)))
        profileList.f <- profileList[profileList.m]
        updateSelectInput(session, "profile", choices=profileList.f)
      }
    })
    
    
    # On click on the filterApply button: go through the filter list and apply each one to tt.hits. Then update the list of profiles
    observe({
      input$filterApply
      isolate({
        fe$tt.hits <- tt.hits.all
        tt.hits <- fe$tt.hits
        # execute the filter
        for(filter in filterList)
          tt.hits <- applyFilter(tt.hits, filter)
        
        # safety measure to prevent super slowdowns
        if(nrow(tt.hits) > hitsLimit)
          fe$tt.hits <- tt.hits[seq_len(hitsLimit),,drop=FALSE]
        else
          fe$tt.hits <- tt.hits
        
        fe$profileNames <- paste0(fe$tt.hits$name, " (", round(fe$tt.hits$mz,4), "/", round(fe$tt.hits$RT / 60,1),")")
        fe$profileList <- fe$tt.hits$profileIDs
        
        names(fe$profileList) <- fe$profileNames
        updateSelectInput(session, "profile", choices=fe$profileList)
        updateTextInput(session, "profileFilter", value="")
      })
    })
    
    ### Pinned profiles list: 
    # Pin profile button adds / removes profile ID from pin list
    # pinnedProfiles data table lists mz/rt/profile ID of pinned profiles
    # exportPinned exports that table as csv (WIP)
    #### Pin profile button
    observe({
      input$pin
      if(input$pin == 0)
        return(c())
      values$pinList <- isolate({
        
        prof <- input$profile
        matchPos <- match(prof, values$pinList)
        if(!is.na(matchPos))
          li <- values$pinList[-matchPos]
        else
          li <- c(values$pinList, prof)
        li
      })
    })
    
    observe({
      hit <- which(tt.hits.all$profileIDs %in% values$pinList)
      values$pinnedProfiles <- tt.hits.all[hit,c("profileIDs", "mass", "name", "dppm", "mz", "RT", "int"),drop=FALSE]
    })
    
    # pinned profiles download
    output$exportPinned <- reactive({values$pinnedProfiles})
    # pinned profiles table
    output$pinnedProfiles <- renderDataTable(values$pinnedProfiles)
    
  #   values$pinnedProfiles <- reactive({
  #     # find profile index in hits
  #     values$pinList
  #     message(paste(values$pinList, collapse=" "))
  #     hits <- match(values$pinList, tt.hits.all$profileIDs)
  #     message(paste(hits, collapse=" "))
  #     tt.hits.all[hits,c("mass", "name", "dppm", "mz", "RT", "int"),drop=FALSE]
  #   })
  #   
  #   output$pinnedProfiles <- renderDataTable(reactive({
  #     values$pinList
  #     message(paste(values$pinList, collapse=" "))
  #     hits <- match(values$pinList, tt.hits.all$profileIDs)
  #     message(paste(hits, collapse=" "))
  #     tt.hits.all[hits,c("mass", "name", "dppm", "mz", "RT", "int"),drop=FALSE]
  # }))
    #output$pinList <- renderText(tt.hits.all[values$pinList,])


	#######################################
	# MS2 integration
	########################################

	# get the RMassBank compound associated with the selected profile,
	# if there is one
    cpd <- reactive({
				p <- input$profile
				if(is.null(spectra)) return(NULL)
				specrow <- which(spectra@aggregated$profile == p)
				if(length(specrow) > 0)
				{
					cpdID <- spectra@aggregated$cpdID[specrow]
					return(spectra@spectra[[as.numeric(cpdID)]])
				}
				else
					return(NULL)
			})
	
	# update the spectrum selection dropdown with the spectra
	# present for this compound
	observe({
				cp <- cpd()
				if(is.null(cp))
					updateSelectInput(session, "spectraSelect", choices = c("."))
				else
				{
					scans <- names(cp@children)
					names(scans) <- readScantype(scans)$string
					updateSelectInput(session, "spectraSelect", choices = scans)
				}
				# set the compound in the reactive container, which is used for RMB processing
				rv$cpd <- cp
			})
	
	# read out the selected spectrum
	spectrum <- reactive({
				cp <- rv$cpd
				if(is.null(cp))
					return(NULL)
				scan <- input$spectraSelect
				sp <- cp@children[[scan]]
				sp
			})
	
	# simply plot MS
	output$spectrumPlot <- renderPlot({
				sp <- getData(spectrum())
				plot(intensity ~ mz, data = sp, type="h", col="red")
			})
	
	# formatted data table output
	output$spectrumTable <- renderDataTable({
				sp <- getData(spectrum())
				spFormatted <- data.frame(
						mz = round(sp$mz, 5),
						intensity = signif(sp$intensity, 3),
						prob = signif(sp$probability, 3) * 100,
						profile = sp$profile, 
						RT = round(sp$RT  / 60, 2)
						)
						sp
			},
			options=list(rownames = FALSE))
	
	observeEvent(input$targetAnalyze,
			{
				withProgress(message = "Processing...",
						value = 0,
						{
				formula <- input$targetFormula
				# get the unprocessed compound
				cpd <- cpd()
				cpd@formula <- formula
				cpd@mz <- 0
				cpdPos <- cpd
				cpdNeg <- cpd
				scans <- readScantype(names(cpd@children))
				cpdPos@children <- cpdPos@children[
						(scans$msLevel == 2) &
						(scans$polarity == "pos")
						]
				cpdNeg@children <- cpdNeg@children[
						(scans$msLevel == 2) &
								(scans$polarity == "neg")
				]
				filterSettings <- list(
						fineCut = 0,
						fineCutRatio = 0,
						specOkLimit = 0,
						dbeMinLimit = 0,
						satelliteMzLimit = 0.5,
						satelliteIntLimit = 0.05,
						ppmFine = 5)
				
				if(length(cpdPos@children) > 1)
				{
					cpdPos <- .depropSpectra(cpdPos)
					cpdPos$new <- analyzeMsMs.formula(cpdPos$new, "pH", TRUE, "",
							filterSettings)
					cpdPos <- .repropSpectra(cpdPos)
				}
					
				if(length(cpdNeg@children) > 1)
				{
					cpdNeg <- .depropSpectra(cpdNeg)
					cpdNeg$new <- analyzeMsMs.formula(cpdNeg$new, "mH", TRUE, "",
							filterSettings)
					cpdNeg <- .repropSpectra(cpdNeg)
				}
				cpd@children[names(cpdPos@children)] <- cpdPos@children
				cpd@children[names(cpdNeg@children)] <- cpdNeg@children
				rv$cpd <- cpd
			})
			})
	
	# update the formula selection if there are formulas for the suspect hits
	observe({
				p <- input$profile
				profileHits <- tt.hits[tt.hits$profile == p,]
				if(is.null(profileHits$formula))
					formulas <- c()
				else
				{
					profileHits <- profileHits[!is.na(profileHits$formula) & (profileHits$formula != ""),,drop=FALSE]
					formulas <- profileHits$formula
					names(formulas) <- c(sprintf("%s [%s]", profileHits$name, profileHits$formula)) 
				}
				formulas <- c(formulas, "Custom formula" = "")
				updateSelectInput(session, "targetSelect", choices=formulas)
			})
	# update the formula field if a target is chosen from the dropdown 
	observe({
				updateTextInput(session, "targetFormula", value=input$targetSelect)
			})
	
	
#    #####################################
#    # RAMClust data
#    clusterTable <- reactive({
#      p <- input$profile
#      if(is.null(clusters)) return(data.frame())
#      isolate({
#        row <- match(p, rcAssignment$profileID)
#        fid <- rcAssignment[row,"featureID"]
#        if(fid != 0)
#        {
#          sp <- (spectra[[as.character(fid)]])
#          sp$int <- signif(sp$int, 3)
#          sp$dmz <- sp$mz - rcAssignment[row, "mz"]
#          sp$mz <- round(sp$mz, 4)
#          sp$dmz <- round(sp$dmz, 4)
#          sp$profint <- signif(sp$profint, 3)
#          sp$RT <- round(sp$RT/60,2)
#          sp$profRT <- round(sp$profRT,1)
#          sp <- sp[,c("featureID", "profileID", "name", "mz", "dmz", "int", "profint", "RT", "profRT")]
#          sp
#        }
#          
#        else
#          return(data.frame())
#      })
#    })
    
#    output$clusterTable <- renderDataTable(clusterTable())
#    
#    output$clusterPlot <- renderPlot(
#      {
#        spec <- clusterTable()
#        if(nrow(spec) == 0)
#        {
#          plot.new()
#          return()
#        }
#        plot.new()
#        plot.window(xlim=range(spec$mz)+c(-10,10), ylim=range(0, spec$int))
#        axis(1)
#        axis(2)
#        isolate({
#          p <- input$profile
#          peak <- spec[spec$profileID == p,]
#          abline(v=peak$mz, col="red")
#        })
#        
#        lines(spec$mz, spec$int, type="h", col="blue", lwd=2)
#        
#      }
#    )
    
    
    
    #####################################
    # Filter buttons and filter list
    
    # Add, delete
    observe({
      input$filterAdd
      if(input$filterAdd == 0) return()
      se$filterList <- c(se$filterList, list(list(type="empty")))
      #print(se$filterList)
      updateFilterList(se, session)
      })
    
    observe({
      input$filterDel
      if(input$filterDel == 0) return()
      isolate({
        #print(input$activeFilters)
        se$filterList <- se$filterList[-as.numeric(input$activeFilters)]
        #print(se$filterList)
        updateFilterList(se, session)
        
      })
    })
    
    # download and upload filter files
    output$filterSave <- downloadHandler(
      filename = function() {
        paste('filter-', Sys.Date(), '.RData', sep='')
      },
      content = function(con) {
        filterList <- se$filterList
        save(filterList, file=con)
      }
    )
    
    observe({
        input$filterLoad
        
        inFile <- input$filterLoad
        
        if (is.null(inFile))
          return(NULL)
        
        load(inFile$datapath)
        se$filterList <- filterList
        updateFilterList(se, session)
        
      })
    
    # clicking on the filter list updates the options set in the right panel
    observe({
      input$activeFilters
      
      if(length(input$activeFilters) == 0) return()
      
      isolate({
        i <- as.numeric(input$activeFilters)
        #print(se$filterList[[i]])
        editFilter(session, as.list(se$filterList[[i]]))
      })
    })
    
    # clicking on "insert" in the literal filter tab inserts the column name to the end of the shinyAce editor
    observe({
      input$insertLiteral
      if(input$insertLiteral == 0) return()
      isolate({
      updateAceEditor(session, "literalFilter", value=paste0(input$literalFilter, input$litFilterInsert))
      })
    })
    
    
    # clicking on any "save filter" options overwrites the currently selected filter with the new info    
    observe({
      input$updateOrder
      input$updateRatioFilter
      input$updateNameFilter
      input$updateLiteralFilter
      if(input$updateOrder+
           input$updateRatioFilter+
           input$updateNameFilter+
          input$updateLiteralFilter
         == 0) return()
      
      isolate({
        if(is.null(input$activeFilters))
          return()
        i <- as.numeric(input$activeFilters)
        if(i == 0)
          return()
        #print(getFilter(input))
        se$filterList[[i]] <- as.list(getFilter(input))
        updateFilterList(se, session)
      })
    })
    
    
    # Plot
    output$profilePlot <- renderPlot({
      
      # find profile index in hits
      hit <- match(input$profile,tt.hits$profileIDs)
      #samples.s <- input$samples
      plotProfile(tt.hits, tt.hits.all, hit, sampleGroups, input$samples, timepoints, addcols)
    })
    
    # Plot
    output$hitTable <- renderDataTable({
      
      # find profile index in hits
      hit <- which(tt.hits.all$profileIDs == input$profile)
      #samples.s <- input$samples
      tt.hits.all[hit,c("mass", "name", "dppm", "mz", "RT", "int"),drop=FALSE]
    })
    output$patternTable <- renderDataTable({
      p <- input$profile
      if(is.null(patterns))
        return(data.frame())
      gid <- patterns[patterns[,"peak ID"] == p,"group ID"]
      if(as.numeric(as.character(gid)) == 0)
        return(patterns[c(),,drop=FALSE])
      gp <- patterns[patterns[,"group ID"] == gid,,drop=FALSE]
      gp
    })
    
    # show files including the profile
    output$peaksTable <- renderDataTable({
      if(is.null(profiles))
         return(data.frame())
      if(is.null(files))
        return(data.frame())
      
      hitPeaks <- with(profiles, peaks[
        peaks$profileIDs == input$profile,
        ])
      hitPeaks <- hitPeaks[order(hitPeaks$intensity, decreasing = TRUE),]
      hitPeaks <- hitPeaks[,c("filename", "intensity", "m/z",  "RT")]
      #files[grepl(hitPeaks$filename)]
      
      hitPeaks$file <- files[
        unlist(lapply(hitPeaks$filename, 
                      function(f) grep(f, files)))]
      hitPeaks$file <- shQuote(hitPeaks$file)
      hitPeaks
    }, options=list(rownames=FALSE))
    
  }
  
  runApp(shinyApp(ui = visualization.ui(colnames(tt.hits.all), profileList, sampleGroups), server = server), ...)
}


# in profiles, col 1 is the profile ID, 
# cols 2:n are the datapoints: timepoints, sd(timepoints), additionalcalculatedcolumns like mean
# finally the cols from the merging(mass, name, dppm, mz, RT, int)
# addcols is the number of additional columns per sample group (default: mean, = 1 additional column)
#' @export
plotProfile <- function(profiles, profiles.all, hit, sampleGroups, selectedGroups, timepoints, addcols = 1)
{
  par(lwd=1.5)
  
  if(is.na(hit)) return()
  if(is.null(hit)) return()
  
  
  col.samples <- sampleGroups$color
  names(col.samples) <- sampleGroups$sampleGroup
  
  if("pch" %in% colnames(sampleGroups))
  pch.samples <- sampleGroups$pch
  else
    pch.samples <- rep(1, nrow(sampleGroups))
  names(pch.samples) <- sampleGroups$sampleGroup
  
  if("lwd" %in% colnames(sampleGroups))
    lwd.samples <- sampleGroups$lwd
  else
    lwd.samples <- rep(1, nrow(sampleGroups))
  names(lwd.samples) <- sampleGroups$sampleGroup
  
  if("lty" %in% colnames(sampleGroups))
    lty.samples <- sampleGroups$lty
  else
    lty.samples <- rep(1, nrow(sampleGroups))
  names(lty.samples) <- sampleGroups$sampleGroup
  
  mat.rownames <- c(
    timepoints$name,
    paste0(timepoints$name, ".sd"),
    "mean")
  
  # extract the entire row from the profiles table containing data for all samples,
  # then reformat as a DF (from original format a.t1, a.t2, a.tN, a.t1.sd, a.t2.sd,
  # a.tN.sd, a.mean, b.t1 etc)
  row <- as.numeric(profiles[hit,2:(nrow(sampleGroups)*(length(mat.rownames))+addcols)])
  mat <- as.data.frame(matrix(row, nrow=length(mat.rownames)))
  rownames(mat) <- mat.rownames
  colnames(mat) <- sampleGroups$sampleGroup
  # split into DF with mean values (top) and SDs (bottom)
  mat.mean <- mat[1:nrow(timepoints),]
  mat.sd <- mat[nrow(timepoints) + (1:nrow(timepoints)),]
  
  plot.new()
  plot.window(xlim=range(0,timepoints$t), ylim=range(0 - mat.sd,mat.mean + mat.sd))
  axis(1)
  axis(2)
  #yshift <- c(sample=0, mix1=100, mix2=200, ctl=300, wc=400)
  
  # The groups are now rows in the mat.mean and mat.sd matrix; plot lines and SD bars by row.
  for(spl in selectedGroups)
  {
    lines(timepoints$t, mat.mean[,spl], col=col.samples[spl], pch=pch.samples[spl],
          lwd = lwd.samples[spl], lty = lty.samples[spl], type="b")
    erb(timepoints$t, mat.mean[,spl], mat.sd[,spl], 0.1, col=col.samples[spl])
    # mark inexistent SDs!
    fsd <- which(mat.sd[,spl] == -1)
    points(timepoints[fsd,"t"], rep(0,length(fsd)), col=col.samples[spl], pch=pch.samples[spl])
  }
  
  # collect all corresponding suspect hits, since more than one suspect can be a hit for a profile
  titles <- profiles.all[which(profiles.all$profileIDs == profiles[hit,1]), "name"]
  titles <- paste(titles, collapse=", ")
  title(main=titles, sub=paste0(round(profiles[hit, "mz"],4), " / ", round(profiles[hit, "RT"] / 60,1)))
  if(length(selectedGroups) > 0)
    legend("topleft", legend=selectedGroups, fill=col.samples[selectedGroups],
           bty="n")
}


# UI is generated at runtime when running the app
# The sidebar has two conditional panels that are shown depending on the tab 
# selected in the main panel,
# the main panel is separated into tabs:
# * profile display (with subtabs: plot, data, RAMClust results)
# * filter building (with sub-UIs for the filter types)
visualization.ui <- function(filterCols, profileList, sampleGroups)
  
  fluidPage(
  sidebarLayout(
    sidebarPanel(
      
      conditionalPanel("input.tab =='visualize'",
                       textInput("profileFilter", "Filter profiles:", ""),
                       selectInput("profile", "Profile:", profileList, size=20, selectize=FALSE),
                       actionButton("pin", "pin profile"),
                       checkboxGroupInput("samples", "Sample groups:", as.list(sampleGroups$sampleGroup), as.list(sampleGroups$sampleGroup))
                       
      ),
      
      conditionalPanel("input.tab == 'buildFilter'",
                       selectInput("activeFilters", "Active filter list", c(), size=20, selectize=FALSE),
                       fluidRow(
                         column(3, actionButton("filterDel", "delete")),
                         column(2, actionButton("filterAdd", "add")),
                         column(3, actionButton("filterApply", "apply")),
                         column(3, downloadButton("filterSave", "save"))
                         
                       ),
                       fluidRow(
                         column(12, fileInput("filterLoad", "load"))
                       )
      )
    ),
    mainPanel(
      tabsetPanel(id="tab",
                  tabPanel("Visualize", value="visualize",
                           tabsetPanel(
                             tabPanel("time series", 
                                      plotOutput("profilePlot", width="600px", height="600px",
                                                 click="profilePlotClick")
                             ),
                             tabPanel("data", 
                                      verticalLayout(
                                        dataTableOutput("hitTable"),
                                        dataTableOutput("patternTable"),
                                        dataTableOutput("peaksTable")                                      
										)
                                      )
                             ,
                             tabPanel("spectra",
                                      verticalLayout(
										fluidRow(column(width=4,
												selectInput("spectraSelect", label="Spectrum:", choices=NULL)),
										column(width=6,
												selectInput("targetSelect", label="Target:", choices=NULL)
												)),
										fluidRow(
												column(width=4, actionButton("targetAnalyze", "Analyze spectra")),
												column(width=6, textInput("targetFormula", label="Formula", value=""))
										),
								
                                        plotOutput("spectrumPlot"),
                                        dataTableOutput("spectrumTable")
                             )
                            )
                             )
                            )
                           
                           ,
                  tabPanel("Build filter", value="buildFilter",
                           tabsetPanel(id="filterType",
                                       tabPanel("Ratio filter", value="ratioFilter", ratioFilterTab(filterCols)),
                                       tabPanel("Name filter", value="nameFilter", nameFilterTab),
                                       tabPanel("Order", value="order", orderTab(filterCols)),
                                       tabPanel("Literal", value="literalFilter", literalFilterTab(filterCols))
                           )),
                tabPanel("Pinned", value="pinned",
                         dataTableOutput("pinnedProfiles"),
                         #textOutput("pinList"),
                         downloadButton("exportPinned", "Export")#,
                         #actionButton("importPinned", "Import")
                         )
                           ))
  )
)
