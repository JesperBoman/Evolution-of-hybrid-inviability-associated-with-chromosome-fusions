suppressPackageStartupMessages(is_installed <- require(devtools))
if (!is_installed) {
  install.packages("devtools")  
  suppressPackageStartupMessages(require(devtools))  
}

suppressPackageStartupMessages(is_installed <- require(visPedigree))
if (!is_installed) {
  install_github("luansheng/visPedigree")  
  suppressPackageStartupMessages(require(visPedigree))
}

require(visPedigree)

f0_f2_pedigree <- read.csv(file = file.choose() , header = T, sep = ";", dec =",", stringsAsFactors = F, na.strings=c(""," ", "NA"))
f0_f2_pedigree_tidy <- tidyped(f0_f2_pedigree[,1:4])
visped(f0_f2_pedigree_tidy, compact=F, cex=1)
