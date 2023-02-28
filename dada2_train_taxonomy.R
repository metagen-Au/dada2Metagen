


library(optparse)
option_list = list(
  optparse::make_option(c("-r", "--run"), type="character", default=NULL,
                        help="dataset run number", metavar="character"),
  optparse::make_option(c("-g", "--group"), type="character", default=NULL,
                        help="taxonomy database to use", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# TMUX
# Load packages
library(dada2);#library(DECIPHER) # load packages for sequence variants, taxomic assignments


run = opt$run
group = opt$group

stopifnot(group %in% c("16S","18S","NEM","ITS"))

print(1)

seq_path = file.path(getwd(),"/",paste0(run),"/",group, "/")
output_path =file.path(paste0(seq_path,"output/"))
sv_path =file.path(paste0(output_path,"SV_tables/"))

print(2)
print(sv_path)
# Move run data


fns<- list.files(sv_path)
print(fns)
print(3)
nc<- fns[grepl("NC", fns)]
print(nc)
print(fns)
nochim.seqtab<- readRDS(paste0(sv_path,nc))

if(group=="18S"){
tax<- assignTaxonomy(getSequences(colnames(nochim.seqtab)),
                     "18S_database.fasta.gz",multithread=TRUE,
                     minBoot = 60 ,
                     tryRC=TRUE,
                     taxLevels =c("Domain",  "Kingdom", "Phylum" , "Class"  , "Order" ,  "Family" , "Genus" , "Species") )

}else if(group=="NEM"){
  tax<- assignTaxonomy(getSequences(colnames(nochim.seqtab)),
                       "18S_database.fasta.gz",multithread=TRUE,
                       minBoot = 60 ,
                       tryRC=TRUE,
                       taxLevels =c("Domain",  "Kingdom", "Phylum" , "Class"  , "Order" ,  "Family" , "Genus" , "Species") )

}else if(group=="16S"){

  tax<- assignTaxonomy(getSequences(colnames(nochim.seqtab)),
                       "16S_database.fa.gz",
                       minBoot = 80 ,tryRC=TRUE , multithread = TRUE)

}else if(group=="ITS"){

  tax<- assignTaxonomy(getSequences(colnames(nochim.seqtab)),
                       "its_database.fasta",
                       minBoot = 60 ,tryRC=TRUE , multithread = TRUE)

}




saveRDS(tax,paste0(sv_path,"taxonomy_",group,"_",run,'.RDS'))
