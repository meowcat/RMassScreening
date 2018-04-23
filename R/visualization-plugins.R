

pinPlugin <- list(
  server = function(env)
  {
    with(env, {
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
    output$exportPinned <- reactive({
      values$pinnedProfiles
    })
    
    output$exportPinned <- downloadHandler(
      filename = function() {
        paste("pinned-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(values$pinnedProfiles, file)
      })
    # pinned profiles table
    output$pinnedProfiles <- renderDataTable(values$pinnedProfiles)
    })
    
  },
  ui = function()
  {
    message("Pinned UI was run")
    return(
    tabPanel("Pinned", value="pinned",
             dataTableOutput("pinnedProfiles"),
             #textOutput("pinList"),
             downloadButton("exportPinned", "Export")#,
             #actionButton("importPinned", "Import")
    )
    )
  }

)