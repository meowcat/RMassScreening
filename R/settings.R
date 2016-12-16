#' Load RMassScreening settings
#'
#' @param file Settings to load
#'
#' @return
#' @export
#'
#' @examples
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
#' @return
#' @export
#'
#' @examples
RmsSettingsTemplate <- function(target = "settings.ini")
{
  file.copy(
    system.file("settings.ini", "RMassScreening"),
    target
  )
}

RmsDefaultSettings <- function()
{
 set <- yaml.load_file(system.file("settings.ini", "RMassScreening"))
 options("RMassScreening" = set)
}