pgks <- c("optparse","gridExtra")
pgks_bioc <- c("dada2","DECIPHER")

for(i in seq_along(pgks)){
  install.packages(pgks[i], repos='http://cran.us.r-project.org')
}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

for(i in seq_along(pgks_bioc)){
  BiocManager::install(pgks_bioc[i])
}



# End 
