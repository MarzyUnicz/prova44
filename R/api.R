#' @title Get MMRF-COMPASS Project Summary from GDC
#' @export
#' @examples
#' \dontrun{
#' getProjectSummary()
#' }
MMRF_getProjectSummary <- function(){

  project<-"MMRF-COMMPASS"
  baseURL <- "https://api.gdc.cancer.gov/projects/"
  url <- paste0(baseURL, project,"?expand=summary,summary.data_categories&pretty=true")
  summary<-fromJSON(url,simplifyDataFrame = TRUE)$data$summary
  return(summary)
}


#' @title Get Number of cases of MMRF-COMMPASS project from GDC
#' @param data.category A  GDC project data category
#' @export
#' @examples
#' \dontrun{
#' getNbCases("Clinical")
#' }
getNbCases <- function(data.category){
  summary <- MMRF_getProjectSummary()
  if(data.category %in% summary$data_categories$data_category){
    summary <- getProjectSummary()$data_categories
    nb <- summary[summary$data_category == data.category,"case_count"]
  } else {
    nb <- summary$case_count
  }
  return(nb)
}






#' @title Get Number of files in MMRF-COMMPASS project
#' @param data.category A  GDC project data category
#' @param legacy Select between Harmonized or Legacy database
#' @export
#' @examples
#' \dontrun{
#' getNbFiles("Clinical")
#' }
getNbFiles <- function(data.category){
  summary <- getProjectSummary()
  if(data.category %in% summary$data_categories$data_category){
    summary <- getProjectSummary()$data_categories
    nb <- summary[summary$data_category == data.category,"file_count"]
  } else {
    nb <- summary$file_count
  }
  return(nb)
}
