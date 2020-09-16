# BSure

To install this R package

if(!require(BiocManager))

{install.packages("BiocManager")}

if(!require(BiocStyle))

{BiocManager::install("BiocStyle")}

library(devtools)

install_github("magStra/BSure",build_vignettes = T,force=T,dependencies=T)

On Windows: source installing for the stringi package currently does not work. If this leads to an error installing BSure install stringi from binaries before: install.packages("stringi", type = "win.binary") and then proceed with the BSure installation.

If install_github fails on your system, download zip from github repository, extract, and set your R working directory to the unzipped folder named folderName. You may need to remove all files in the subfolder src for sucessful install (but do not delete the folder). 
Then use install.packages("../folderName",repos=NULL,type="source")
 
To get started after installing

library(BSure)

vignette("BSure")
