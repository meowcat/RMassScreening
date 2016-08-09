#library(shiny)
#library(shinyAce)
#' @export
erb <- function(x,y,sd,eps,...)
{
  segments(x,y-sd,x,y+sd,...)
  segments(x-eps,y-sd,x+eps,y-sd,...)
  segments(x-eps,y+sd,x+eps,y+sd,...)
}

#' @export
runViewer <- function(totalTable, hits, timepoints, sampleGroups, patterns = NULL, addcols = 1, ...)
{
  fe <- environment()
  
  tt.hits.all <- merge(totalTable, hits, by.x="profileIDs", by.y="profileID")
  tt.hits <- tt.hits.all[1:2000,]
  
  
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
        if(nrow(tt.hits) > 2000)
          fe$tt.hits <- tt.hits[1:2000,,drop=FALSE]
        else
          fe$tt.hits <- tt.hits
        
        fe$profileNames <- paste0(fe$tt.hits$name, " (", round(fe$tt.hits$mz,4), "/", round(fe$tt.hits$RT / 60,1),")")
        fe$profileList <- fe$tt.hits$profileIDs
        
        names(fe$profileList) <- fe$profileNames
        updateSelectInput(session, "profile", choices=fe$profileList)
        updateTextInput(session, "profileFilter", value="")
      })
    })
    
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
  
  pch.samples <- sampleGroups$pch
  names(pch.samples) <- sampleGroups$sampleGroup
  
  mat.rownames <- c(
    timepoints$name,
    paste0(timepoints$name, ".sd"),
    "mean")
  
  row <- as.numeric(profiles[hit,2:(nrow(sampleGroups)*(length(mat.rownames))+addcols)])
  mat <- as.data.frame(matrix(row, nrow=length(mat.rownames)))
  rownames(mat) <- mat.rownames
  colnames(mat) <- sampleGroups$sampleGroup
  mat.mean <- mat[1:nrow(timepoints),]
  mat.sd <- mat[nrow(timepoints) + (1:nrow(timepoints)),]
  
  plot.new()
  plot.window(xlim=range(0,timepoints$t), ylim=range(0 - mat.sd,mat.mean + mat.sd))
  axis(1)
  axis(2)
  #yshift <- c(sample=0, mix1=100, mix2=200, ctl=300, wc=400)
  for(spl in selectedGroups)
  {
    lines(timepoints$t, mat.mean[,spl], col=col.samples[spl], pch=pch.samples[spl], type="b")
    erb(timepoints$t, mat.mean[,spl], mat.sd[,spl], 0.1, col=col.samples[spl])
    # mark inexistent SDs!
    fsd <- which(mat.sd[,spl] == -1)
    points(timepoints[fsd,"t"], rep(0,length(fsd)), col=col.samples[spl], pch=pch.samples[spl])
  }
  
  titles <- profiles.all[which(profiles.all$profileIDs == profiles[hit,1]), "name"]
  titles <- paste(titles, collapse=", ")
  title(main=titles, sub=paste0(round(profiles[hit, "mz"],4), " / ", round(profiles[hit, "RT"] / 60,1)))
  if(length(selectedGroups) > 0)
    legend("topleft", legend=selectedGroups, fill=col.samples[selectedGroups],
           bty="n")
}



visualization.ui <- function(filterCols, profileList, sampleGroups)
  
  fluidPage(
  sidebarLayout(
    sidebarPanel(
      conditionalPanel("input.tab =='visualize'",
                       textInput("profileFilter", "Filter profiles:", ""),
                       selectInput("profile", "Profile:", profileList, size=20, selectize=FALSE),
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
                                      plotOutput("profilePlot", width="600px", height="600px" )),
                             tabPanel("data", 
                                      verticalLayout(
                                        dataTableOutput("hitTable"),
                                        dataTableOutput("patternTable")
                                      )
                                      )
#                              ,
#                              tabPanel("cluster",
#                                       verticalLayout(
#                                         actionButton("findClusters", "Find clusters"),
#                                         numericInput("rttolCluster", "RT tolerance", 30),
#                                         dataTableOutput("clusterTable")
#                                       )
                               
                             )
                            )
                           
                           ,
                  tabPanel("Build filter", value="buildFilter",
                           tabsetPanel(id="filterType",
                                       tabPanel("Ratio filter", value="ratioFilter", ratioFilterTab(filterCols)),
                                       tabPanel("Name filter", value="nameFilter", nameFilterTab),
                                       tabPanel("Order", value="order", orderTab(filterCols)),
                                       tabPanel("Literal", value="literalFilter", literalFilterTab(filterCols))
                                       
                                       
                           ))
      ))
  )
)
