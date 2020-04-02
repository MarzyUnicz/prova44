## code to prepare `mm.exp` dataset goes here

file.path <- system.file(
  "extdata",
  "mm.exp.rda",
  package = "MMRFBiolinksX"
)

mm.exp<-get(load(file = file.path))


usethis::use_data(mm.exp, overwrite = TRUE)


