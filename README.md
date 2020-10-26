# BSure

To install this R package

if(!require(BiocManager))

{install.packages("BiocManager")}

if(!require(BiocStyle))

{BiocManager::install("BiocStyle")}

library(devtools)

install_github("magStra/BSure",build_vignettes = T,force=T,dependencies=T)

For problems with rstan dependence, install the versions package and use install.versions("rstan",versions="2.19.3") before proceeding to install BSure. 

On Windows: source installing for the stringi package currently does not work. If this leads to an error installing BSure install stringi from binaries before: install.packages("stringi", type = "win.binary")-or choose installation from binary for the stringi package during BSure installation- and then proceed with the BSure installation. Proceed in the same way, in case this problem occurs with other packages than stringi. 

For similar problems on Mac, replace "win.binary" by "mac.binary" or "mac.binary.el-capitan".

If install_github fails on your system, download zip from github repository, extract, and set your R working directory to the unzipped folder named folderName. You may need to remove all files in the subfolder src for sucessful install (but do not delete the folder). 
Then use install.packages("../folderName",repos=NULL,type="source",build_vignettes = T)

If this fails, you may need to, after deleting the contents of the src folder, rebuild and recompile the package as follows (the working directory needs to be the folder containing the package - folderName defined above):  

pkgbuild::compile_dll()

roxygen2::roxygenize()
 
To get started after installing

library(BSure)

vignette("BSure")

Further troubleshooting concerning rstan - advice from a college (Windows), see also https://cran.r-project.org/web/packages/tlrmvnmvt/readme/README.html and https://github.com/stan-dev/stan/issues/1613:

dotR <- file.path(Sys.getenv("HOME"), ".R")

if (!file.exists(dotR)) 

dir.create(dotR)

M <- file.path(dotR, "Makevars.win")

if (!file.exists(M)) 

file.create(M)

cat("\nCXX14FLAGS=-O3 -Wno-unused-variable -Wno-unused-function",
"CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y",
"CXX11FLAGS=-O3 -Wno-unused-variable -Wno-unused-function",
file = M, sep = "\n", append = TRUE)
