# BSure

To install this R package

library(devtools)

install_github("magStra/BSure",build_vignettes = T,force=T)

If install_github fails on your system, download zip from github repository, extract, and set your R working directory to the unzipped folder named folderName. You may need to remove all files in the subfolder src for sucessful install (but do not delete the folder). 
Then use install.packages("../folderName",repos=NULL,type="source")
 
To get started after installing

library(BSure)

vignette("BSure")
