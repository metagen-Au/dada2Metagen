pgks <- c("dada2","optparse","DECIPHER","ggplot2","gridExtra")

for(i in seq_along(pgks)){
  install.packages(pgks, repos='http://cran.us.r-project.org')
}

# End 
