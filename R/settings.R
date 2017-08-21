#' Load RMassScreening settings
#'
#' @param file Settings to load
#'
#' @export
#'
loadRmsSettings <- function(file)
{
  set <- yaml.load_file(file)
  # check for scientific notation and convert to numbers
  options("RMassScreening" = set)
}

#' Create RMassScreening settings template file
#'
#' @param target Target filename
#'
#' @export
#'
RmsSettingsTemplate <- function(target = "settings.ini")
{
  file.copy(
    system.file("settings.ini", "RMassScreening"),
    target
  )
}

#' Load default RMassScreening settings 
#'
#' @export
#'
RmsDefaultSettings <- function()
{
 set <- yaml.load_file(system.file("settings.ini", "RMassScreening"))
 options("RMassScreening" = set)
}