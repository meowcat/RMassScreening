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
runViewer <- function(totalTable, hits, timepoints, sampleGroups)
{
  fe <- environment()
  
  tt.hits.all <- merge(totalTable, hits, by.x="profileIDs", by.y="profileID")
  tt.hits <- tt.hits.all[1:2000,]
  
  
  filterText <- "
#tt.hits <- tt.hits[tt.hits$t0.mixed < 1e5,]
tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN)
                           > (tt.hits$mean.synCTL + tt.hits$mean.chlCTL + tt.hits$mean.mcyCTL) * 10,]
tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN) > tt.hits$mean.WC * 7,]
tt.hits <- tt.hits[(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN) > 3e5,]

#tt.hits <- tt.hits[grepl(\"AZY\",tt.hits$name),]

tt.hits <- tt.hits[order(tt.hits$mean.synSN + tt.hits$mean.chlSN + tt.hits$mean.mcySN, decreasing = TRUE),]

"
  
  profileNames <- paste0(tt.hits$name, " (", round(tt.hits$mz,4), "/", round(tt.hits$RT / 60,1),")")
  profileList <- tt.hits$profileIDs
  names(profileList) <- profileNames
  
  
  
  server <- function(input, output, session) {
    
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
    
    
    observe({
      i <- input$filterApply
      if(i>0)
      {
        fe$tt.hits <- tt.hits.all
        tt.hits <- fe$tt.hits
        # execute the filter
        isolate(eval(parse(text=input$filter)))
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
      }
    })
    
    output$profilePlot <- renderPlot({
      
      # find profile index in hits
      hit <- match(input$profile,tt.hits$profileIDs)
      #samples.s <- input$samples
      
      
      plotProfile(tt.hits, tt.hits.all, hit, sampleGroups, input$samples, timepoints)
      
      
      
    })
  }
  
  ui <- fluidPage(
      sidebarLayout(
        sidebarPanel(
          textInput("profileFilter", "Filter profiles:", ""),
          selectInput("profile", "Profile:", profileList, size=20, selectize=FALSE),
          checkboxGroupInput("samples", "Sample groups:", as.list(sampleGroups$sampleGroup))
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Visualize", plotOutput("profilePlot", width="600px", height="600px" )),
            tabPanel("Filter", aceEditor("filter", filterText), actionButton("filterApply", "apply"))
          ))
      )
    )
  shinyApp(ui = ui, server = server)
}


# in profiles, col 1 is the profile ID, 
# cols 2:n are the datapoints: timepoints, sd(timepoints), additionalcalculatedcolumns like mean
# finally the cols from the merging(mass, name, dppm, mz, RT, int)
#' @export
plotProfile <- function(profiles, profiles.all, hit, sampleGroups, selectedGroups, timepoints)
{
  par(lwd=1.5)
  
  
  col.samples <- sampleGroups$color
  names(col.samples) <- sampleGroups$sampleGroup
  
  pch.samples <- sampleGroups$pch
  names(pch.samples) <- sampleGroups$sampleGroup
  
  mat.rownames <- c(
    timepoints$name,
    paste0(timepoints$name, ".sd"),
    "mean")
  
  row <- as.numeric(profiles[hit,2:(nrow(sampleGroups)*(length(mat.rownames))+1)])
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
  
  titles <- profiles.all[which(profiles.all$profileIDs == row[[1]]), "name"]
  titles <- paste(titles, collapse=", ")
  title(main=titles, sub=paste0(round(profiles[hit, "mz"],4), " / ", round(profiles[hit, "RT"] / 60,1)))
  legend("topleft", legend=selectedGroups, fill=col.samples[selectedGroups],
         bty="n")
}

