library(shinyAce)

profileList <- c(1,2,3)
sampleGroups <- data.frame(sampleGroup = c("a","b", "c"))
filterText <- ""
filterCols <- cols[2:28,2]
filterCols <- as.character(filterCols)

server <- function(input, output, session)
{
  filterList <- list()
  se <- environment()
  
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
  
  observe({
      input$filterAdd
      if(input$filterAdd == 0) return()
      se$filterList <- c(se$filterList, list(list(type="empty")))
      #print(se$filterList)
      afList <- 1:length(se$filterList)
      names(afList) <- unlist(lapply(se$filterList, renderFilter))
      updateSelectInput(session, "activeFilters", choices=afList)    })
  
  observe({
    input$filterDel
    if(input$filterDel == 0) return()
    isolate({
      #print(input$activeFilters)
      se$filterList <- se$filterList[-as.numeric(input$activeFilters)]
      #print(se$filterList)
      afList <- 1:length(se$filterList)
      names(afList) <- unlist(lapply(se$filterList, renderFilter))
      updateSelectInput(session, "activeFilters", choices=afList)
    })
  })
  
  observe({
    input$activeFilters
    isolate({
      i <- as.numeric(input$activeFilters)
      #print(se$filterList[[i]])
      editFilter(session, as.list(se$filterList[[i]]))
    })
  })
  
  observe({
    input$updateOrder
    input$updateRatioFilter
    input$updateNameFilter
    if(input$updateOrder+
       input$updateRatioFilter+
       input$updateNameFilter
       == 0) return()
    
    isolate({
      i <- as.numeric(input$activeFilters)
      #print(getFilter(input))
      se$filterList[[i]] <- as.list(getFilter(input))
      afList <- 1:length(se$filterList)
      names(afList) <- unlist(lapply(se$filterList, renderFilter))
      updateSelectInput(session, "activeFilters", choices=afList)
    })
  })
}



