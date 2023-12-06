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
