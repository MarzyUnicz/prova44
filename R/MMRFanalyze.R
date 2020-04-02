

#' @title Array Array Intensity correlation (AAIC) and correlation boxplot to define outlier
#' @description TCGAanalyze_Preprocessing perform Array Array Intensity correlation (AAIC).
#' It defines a square symmetric matrix of spearman correlation among samples.
#' According this matrix and boxplot of correlation samples by samples it is possible
#' to find samples with low correlation that can be identified as possible outliers.
#' @param object of gene expression of class RangedSummarizedExperiment from TCGAprepare
#' @param cor.cut is a threshold to filter samples according their spearman correlation in
#' samples by samples. default cor.cut is 0
#' @param filename Filename of the image file
#' @param width Image width
#' @param height Image height
#' @param datatype is a string from RangedSummarizedExperiment assay
#' @importFrom grDevices dev.list
#' @importFrom SummarizedExperiment assays
#' @export
#' @return Plot with array array intensity correlation and boxplot of correlation samples by samples

MMRFanalyze_Preprocessing <- function(object,
                                      cor.cut = 0,
                                      filename = NULL,
                                      width = 1000,
                                      height = 1000,
                                      datatype = names(assays(object))[1]){


  TCGAbiolinks::TCGAanalyze_Preprocessing(object,cor.cut, filename, width,height,datatype)

}








#' @title survival analysis (SA) univariate with Kaplan-Meier (KM) method.
#' @description MMRFanalyze_SurvivalKM perform an univariate Kaplan-Meier (KM) survival analysis (SA).
#' It performed Kaplan-Meier survival univariate using complete follow up with all days
#' taking one gene a time from Genelist of gene symbols.
#' For each gene according its level of mean expression in cancer samples,
#' defining two thresholds for quantile
#' expression of that gene in all samples (default ThreshTop=0.67,ThreshDown=0.33) it is possible
#' to define a threshold of intensity of gene expression to divide the samples in 3 groups
#' (High, intermediate, low).
#' MMRFanalyze_SurvivalKM performs SA between High and low groups using following functions
#' from survival package
#' \enumerate{
#' \item survival::Surv
#' \item survival::survdiff
#' \item survival::survfit
#' }
#' @param clinical_patient is a data.frame using function 'clinic' with information
#' related to barcode / samples such as bcr_patient_barcode, days_to_death ,
#' days_to_last_follow_up , vital_status, etc
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from TCGAprepare
#' @param Genelist is a list of gene symbols where perform survival KM.
#' @param Survresult is a parameter (default = FALSE) if is TRUE will show KM plot and results.
#' @param ThreshTop is a quantile threshold to identify samples with high expression of a gene
#' @param ThreshDown is a quantile threshold to identify samples with low expression of a gene
#' @param p.cut p.values threshold. Default: 0.05
#' @param group1 a string containing the barcode list of the samples in in control group
#' @param group2 a string containing the barcode list of the samples in in disease group
#' @importFrom survival Surv survdiff survfit
#' @export
#' @return table with survival genes pvalues from KM.
#' @examples
#'  # Selecting only 20 genes for example
#'  dataMMRFcomplete <- log2(DataMMRF[1:20,] + 1)
#'
#'  # clinical_patient_Cancer <- MMRFquery_clinic("clinical")
#'  clinical_patient_Cancer <- data.frame(
#'       bcr_patient_barcode = substr(colnames(dataMMRFcomplete),1,12),
#'       vital_status = c(rep("alive",3),"dead",rep("alive",2),rep(c("dead","alive"),2)),
#'       days_to_death = c(NA,NA,NA,172,NA,NA,3472,NA,786,NA),
#'       days_to_last_follow_up = c(3011,965,718,NA,1914,423,NA,5,656,1417)
#'  )
#'
#'  group1 <- TCGAquery_SampleTypes(colnames(dataBRCAcomplete), typesample = c("NT"))
#'  group2 <- TCGAquery_SampleTypes(colnames(dataBRCAcomplete), typesample = c("TP"))
#'
#'  tabSurvKM <- TCGAanalyze_SurvivalKM(clinical_patient_Cancer,
#'                                      dataMMRFcomplete,
#'                                      Genelist = rownames(dataMMRFcomplete),
#'                                      Survresult = FALSE,
#'                                      p.cut = 0.4,
#'                                      ThreshTop = 0.67,
#'                                      ThreshDown = 0.33,
#'                                      group1 = group1, # Control group
#'                                      group2 = group2) # Disease group
#'
#'  # If the groups are not specified group1 == group2 and all samples are used
#'  \dontrun{
#'  tabSurvKM <- MMRFanalyze_SurvivalKM(clinical_patient_Cancer,
#'                                      dataMMRFcomplete,
#'                                      Genelist = rownames(dataMMRFcomplete),
#'                                      Survresult = TRUE,
#'                                      p.cut = 0.2,
#'                                      ThreshTop = 0.67,
#'                                      ThreshDown = 0.33)
#' }
MMRFanalyze_SurvivalKM <- function(clinical_patient,
                                   dataGE,
                                   Genelist,
                                   Survresult = FALSE,
                                   ThreshTop = 0.67,
                                   ThreshDown = 0.33,
                                   p.cut = 0.05,
                                   group1,
                                   group2){

  colnames(dataGE) <- substr(colnames(dataGE),1,9)
  tabSurvKM<-TCGAbiolinks::TCGAanalyze_SurvivalKM(clinical_patient,dataGE,Genelist, Survresult, ThreshTop, ThreshDown, p.cut, group1, group2)

  return(tabSurvKM)


}







#' #' @title MMRF_get_barcodeByGenes.
#' @description
#'   MMRF_get_barcodeByGenes returns a subset of the Gene expression matrix (genes in rows, samples in cols) from MMRFprepare.
#' @param dataGE is a matrix of Gene expression (genes in rows, samples in cols) from MMRFprepare
#' @param geneList a character vector of the genes symbols.
#' @return a subset of the Gene expression matrix
#' @export

MMRF_FPKM_filterByGenes<- function(dataGE,genes){

  ENSEMBL<-row.names(dataGE)
  ENSEMBL<-as.data.frame(geneListGE)


 # colnames(geneListGE)[colnames(geneListGE) == "geneListGE"] <- "ENSEMBL"

  dataGE.geneSymbols<-cbind(geneListGE,dataGE)


  eg = as.data.frame(bitr(geneListGE$ENSEMBL,
                          fromType="ENSEMBL",
                          toType="SYMBOL",
                          OrgDb="org.Hs.eg.db"))


  dataGE.geneSymbols<-merge(x = eg, y = dataGE.geneSymbols, by = "ENSEMBL", all = TRUE)

  dataGE.geneSymbols.filtered<-NULL

  for(genes.i in 1:length(genes)){
    aux<-genes[genes.i]
    if (is.element(aux,dataGE.geneSymbols$SYMBOL)==TRUE){
      gene.sub<-subset(dataGE.geneSymbols,SYMBOL==aux)
      dataGE.geneSymbols.filtered<-rbind(dataGE.geneSymbols.filtered,gene.sub)
    }

    else {
      return("Error message: one or more  genes do not exist")
    }

  }





return(dataGE.geneSymbols.filtered)




}





