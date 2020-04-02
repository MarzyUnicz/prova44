


#' @title Retrieve barcode filtered by tissue or sample type.
#' @description
#'   MMRFfilterBarcodesByDrugs returns the barcode of samples filtered on theses types.
#' @param clin.bio is data frame with the clinical information (A valid type is Biospecimen)
#' @param sampletype a character vector indicating sample type
#' @param tissuetype a character vector indicating tissue type
#' Example:
#' \tabular{ll}{
#'TB \tab   Primary Blood Derived Cancer-Peripheral Blood \cr
#'TRBM \tab Recurrent Blood Derived Cancer-Bone Marrow \cr
#'TBM \tab  Primary Blood Derived Cancer-Bone Marrow \cr
#'NB \tab   Blood Derived Normal \cr
#'TRB \tab Recurrent Blood Derived Cancer - Peripheral Blood \cr
#'}
#' @export
#' @examples
#' clin.bio<-MMRFquery_clinic(type = "biospecimen")
#'  selection of Normal samples
#'  MMRFfilterBarcodesByTypes(clin.bio,tissuetype="Normal")
#'  selection of Tumor samples
#'  MMRFfilterBarcodesByTypes(clin.bio)
#'  selection of barcodes with Tumor tissue and sample type "Recurrent Blood Derived Cancer - Bone Marrow" (TRBM)
#'  MMRFfilterBarcodesByTypes(clin.bio,"TRBM","Tumor")
#' @return a list of barcode filtered by  sample or  tissue type selected



MMRFfilterBarcodesByDrugs<-function(clin, sampletype){

  s_type <- c("Recurrent Blood Derived Cancer - Bone Marrow",
              "Primary Blood Derived Cancer - Bone Marrow",
              "Blood Derived Normal",
              "Primary Blood Derived Cancer - Peripheral Blood",
              "Recurrent Blood Derived Cancer - Peripheral Blood")
  names(s_type) <- c("TRBM","TBM","NB","TB","TRB")




  if(missing(sampletype)) stop("A saple type need to be provided")



  if (!is.null(sampletype) && !is.element(sampletype, names(s_type)) ){
    return("Error message: one or more sample types do not exist")
  }

  aux<-data.frame(clin.bio$submitter_id, clin.bio$sample_type)
  colnames(aux)<-c('bcr_patient_barcode','sample_type')




  switch (sampletype,
          "TRBM"= {type.sub <- aux[aux$sample_type %in% s_type[[1]],]},
          "TBM"=  {type.sub <- aux[aux$sample_type %in% s_type[[2]],]},
          "NB"=   {type.sub <- aux[aux$sample_type %in% s_type[[3]],]},
          "TB"=   {type.sub <- aux[aux$sample_type %in% s_type[[4]],]},
          "TRB"=  {type.sub <- aux[aux$sample_type %in% s_type[[5]],]}

  )




  if(nrow(type.sub)==0) stop("there exist no matched samples")

  barcode<-unique(substr(type.sub$bcr_patient_barcode,1,9))



  return(barcode)

}

















#' @title Retrieve barcode filtered by tissue or sample type.
#' @description
#'   MMRFfilterBarcodesByTypes returns the barcode of samples filtered on theses types.
#' @param clin.bio is data frame with the clinical information (A valid type is Biospecimen)
#' @param sampletype a character vector indicating sample type
#' @param tissuetype a character vector indicating tissue type
#' Example:
#' \tabular{ll}{
#'TB \tab   Primary Blood Derived Cancer-Peripheral Blood \cr
#'TRBM \tab Recurrent Blood Derived Cancer-Bone Marrow \cr
#'TBM \tab  Primary Blood Derived Cancer-Bone Marrow \cr
#'NB \tab   Blood Derived Normal \cr
#'TRB \tab Recurrent Blood Derived Cancer - Peripheral Blood \cr
#'}
#' @export
#' @examples
#' clin.bio<-MMRFquery_clinic(type = "biospecimen")
#'  selection of Normal samples
#'  MMRFfilterBarcodesByTypes(clin.bio,tissuetype="Normal")
#'  selection of Tumor samples
#'  MMRFfilterBarcodesByTypes(clin.bio)
#'  selection of barcodes with Tumor tissue and sample type "Recurrent Blood Derived Cancer - Bone Marrow" (TRBM)
#'  MMRFfilterBarcodesByTypes(clin.bio,"TRBM","Tumor")
#' @return a list of barcode filtered by  sample or  tissue type selected



MMRFfilterBarcodesByTypes<-function(clin.bio, sampletype){

  s_type <- c("Recurrent Blood Derived Cancer - Bone Marrow",
              "Primary Blood Derived Cancer - Bone Marrow",
              "Blood Derived Normal",
              "Primary Blood Derived Cancer - Peripheral Blood",
              "Recurrent Blood Derived Cancer - Peripheral Blood")
  names(s_type) <- c("TRBM","TBM","NB","TB","TRB")




  if(missing(sampletype)) stop("A saple type need to be provided")



   if (!is.null(sampletype) && !is.element(sampletype, names(s_type)) ){
    return("Error message: one or more sample types do not exist")
  }

  aux<-data.frame(clin.bio$submitter_id, clin.bio$sample_type)
  colnames(aux)<-c('bcr_patient_barcode','sample_type')




  switch (sampletype,
          "TRBM"= {type.sub <- aux[aux$sample_type %in% s_type[[1]],]},
          "TBM"=  {type.sub <- aux[aux$sample_type %in% s_type[[2]],]},
          "NB"=   {type.sub <- aux[aux$sample_type %in% s_type[[3]],]},
          "TB"=   {type.sub <- aux[aux$sample_type %in% s_type[[4]],]},
          "TRB"=  {type.sub <- aux[aux$sample_type %in% s_type[[5]],]}

         )




  if(nrow(type.sub)==0) stop("there exist no matched samples")

   barcode<-unique(substr(type.sub$bcr_patient_barcode,1,9))



  return(barcode)

}

#' @title Get GDC clinical data
#' @description
#' MMRFquery_clinical will download and prepare all clinical information from the API
#' as the one with using the button from each project
#' @param "MMRF-COMMPASS" project
#' \itemize{
#' \item{ MMRF-COMMPASS }
#' }
#' @param type A valid type. Options "clinical", "Biospecimen"  (see list with getGDCprojects()$project_id)]
#' @param save.csv Write clinical information into a csv document
#' @export
#' @importFrom data.table rbindlist as.data.table
#' @importFrom jsonlite fromJSON
#' @examples
#' clin.mm<-MMRFquery_clinic(type = "clinical", save.csv = TRUE)
#' clin.mm<-MMRFquery_clinic(type = "biospecimen", save.csv = TRUE)

#' \dontrun{
#' clin <- MMRFquery_clinical(type = "clinical")
#' }
#' @return A data frame with the clinical information




MMRFquery_clinic<- function(type = "clinical", save.csv = FALSE){

  clin<-TCGAbiolinks::GDCquery_clinic(project="MMRF-COMMPASS", type, save.csv)

  #names(clin)[names(clin) == 'submitter_id'] <- 'bcr_patient_barcode'

  #if clinical

  if (type!="clinical" & type!="Biospecimen" ){
    return("Error message: type does not exist")
  }

  if (type=="clinical"){
    names(clin)[1] <- "bcr_patient_barcode"

  }
  else if  (type=="Biospecimen"){
    bcr_patient_barcode.sub<-clin[colnames(clin)=='submitter_id']
    colnames(bcr_patient_barcode.sub)[1]<-"bcr_patient_barcode"
    bcr_patient_barcode.sub$bcr_patient_barcode<-substr(bcr_patient_barcode.sub$bcr_patient_barcode,1,9)
    clin<-cbind(bcr_patient_barcode.sub,clin)

  }


  return(clin)

}








