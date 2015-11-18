
getRatioFilter <- function(input)
{
  fLHS <- input$filterLHS
  fRHS <- input$filterRHS
  fRatio <- input$filterRatio
  fOperator <- input$filterOperator
  list(type="ratio", LHS = fLHS, RHS = fRHS, ratio = fRatio, operator = fOperator)
}

editRatioFilter <- function(session, filter)
{
  updateTabsetPanel(session, "filterType", selected="ratioFilter")
  updateSelectInput(session, "filterLHS", selected=filter$LHS)
  updateSelectInput(session, "filterRHS", selected=filter$RHS)
  updateNumericInput(session, "filterRatio", value=filter$ratio)
  updateSelectInput(session, "filterOperator", selected=filter$operator)
}

applyRatioFilter <- function(df, filter)
{
  if(filter$operator == ">")
    factor <- 1
  else
    factor <- -1
  df[
    factor * ((rowSums(df[,filter$LHS,drop=FALSE])) - 
      (filter$ratio * rowSums(df[,filter$RHS,drop=FALSE])) )
    > 0
    ,,drop=FALSE]
}

getNameFilter <- function(input)
{
  list(type="name", name=input$filterName)
}

editNameFilter <- function(session, filter)
{
  updateTabsetPanel(session, "filterType", selected="nameFilter")
  updateTextInput(session, "filterName", value=filter$name)
}

applyNameFilter <- function(df, filter)
{
  df[grepl(filter$name, df$name),,drop=FALSE]
}

getOrder <- function(input)
{
  list(type="order", cols=input$orderCol, dir=as.numeric(input$orderDir))
}

applyOrder <- function(df, filter)
{
  #print(filter)
  #print(rowSums(df[,filter$cols, drop=FALSE])[1:10])
  df[order(rowSums(df[,filter$cols, drop=FALSE]) * filter$dir)
  ,,drop=FALSE]
}

editOrder <- function(session, filter)
{
  updateTabsetPanel(session, "filterType", selected="order")
  updateSelectInput(session, "orderCol", selected=filter$cols)
  updateSelectInput(session, "orderDir", selected=filter$dir)
}

renderOrder <- function(filter)
{
  paste("Order",
        ifelse(filter$dir == 1, "/\\", "\\/"),
        "Sum",
        paste(filter$cols, collapse="+")
  )
}

renderRatioFilter <- function(filter)
{
  paste("Ratio",
        paste(filter$LHS, collapse="+"),
        filter$operator,
        filter$ratio,
        "x",
        paste(filter$RHS, collapse="+"))
}


renderNameFilter <- function(filter)
{
  paste("Name contains",
        filter$name)
}


renderFilter <- function(filter)
{
  if(filter$type == "name") return(renderNameFilter(filter))
  if(filter$type == "order") return(renderOrder(filter))
  if(filter$type == "ratio") return(renderRatioFilter(filter))
  if(filter$type == "empty") return("[empty]")
}

applyFilter <- function(df, filter)
{
  print(filter)
  if(filter$type == "name") return(applyNameFilter(df, filter))
  if(filter$type == "order") return(applyOrder(df, filter))
  if(filter$type == "ratio") return(applyRatioFilter(df, filter))
  if(filter$type == "empty") return(df)
  stop("filter type is undefined")
}

getFilter <- function(input)
{
  if(input$filterType == "ratioFilter") return(getRatioFilter(input))
  if(input$filterType == "order") return(getOrder(input))
  if(input$filterType == "nameFilter") return(getNameFilter(input))
  return(list(type="empty"))
}

editFilter <- function(session, filter)
{
  print(filter$type)
  if(filter$type == "name") editNameFilter(session, filter)
  if(filter$type == "order") editOrder(session, filter)
  if(filter$type == "ratio") editRatioFilter(session, filter)
}




ratioFilterTab <- function(filterCols)
  fluidRow(
  column(5, selectInput("filterLHS", "Column sum", filterCols, multiple=TRUE, selectize=FALSE)),
  column(2, selectInput("filterOperator", "Operator", c(">", "<"), selectize=FALSE),
         numericInput("filterRatio", "Filter ratio", 1,0.001, 1000),
         actionButton("updateRatioFilter", "Save filter")
  ),
  column(5, selectInput("filterRHS", "Column sum", filterCols, multiple=TRUE, selectize=FALSE))
)

nameFilterTab <- verticalLayout(
  textInput("filterName", "Contains:"),
  actionButton("updateNameFilter", "Save filter"))

orderTab <- function(filterCols)
  verticalLayout(
  selectInput("orderCol", "Column sum", filterCols, multiple=TRUE, selectize=FALSE),
  selectInput("orderDir", "Direction", c("ascending" = 1, "descending" = -1)),
  actionButton("updateOrder", "Save filter")
)

