## code to prepare `DataMM` dataset goes here

## code to prepare `dataMM` dataset goes here

file.path <- system.file(
  "extdata",
  "dataMM.rda",
  package = "MMRFBiolinksX"
)

dataMM<-get(load(file = file.path))


usethis::use_data(dataMM, overwrite = TRUE)





