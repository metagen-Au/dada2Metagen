

library("optparse")

option_list = list(
  make_option(c("-r", "--run"), type="character", default=NULL, 
              help="dataset run number", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
# TMUX
# Load packages
library(dada2);#library(DECIPHER) # load packages for sequence variants, taxomic assignments
library(ggplot2)
library(gridExtra) # packages for visualiations
getN <- function(x) sum(getUniques(x))

run = opt$run 

seq_path = file.path(getwd(),"/",paste0(run),"/16S/")
dir.create(paste0(seq_path,"output"))
output_path =file.path(paste0(seq_path,"output/"))
dir.create(paste0(output_path,"SV_tables"))
sv_path =file.path(paste0(output_path,"SV_tables/"))

# Move run data
dir.create(paste0(output_path,"run_data"))
run_path = file.path(output_path,"run_data/")
dir.create(paste0(seq_path,"cutadapt"))
path.cut <- file.path(paste0(seq_path, "cutadapt/"))

print(getwd())
print(output_path)
print(run_path)
print(path.cut)
#Make a datestamped ID for files of this run
#Create file path, get files

ID<- paste0(Sys.Date(),"_",run,"_16S_" )
#seq_path<-  "P:/soil_health_reports/bioinformatics/Report_Fastqs/16122020/metagen11/16S/"
#ps_path<- "P:/soil_health_reports/phyloseq_objects/"
fns <- sort(list.files(seq_path, full.names = TRUE)) # or fns <- sort(list.files( full.names = TRUE))


# Sort files by fwd and rvs

fnFs <- fns[grepl("r1.fq.gz",fns)]
# Check fastq naming convention.
if(length(fnFs)==0){

  fns2<- gsub("fastq","fq",basename(fns))
  file.rename(list.files(seq_path,full.names = TRUE),paste0(seq_path,fns2))

  fnFs <- fns[grepl("r1.fq.gz",fns)]
  fnRs <- fns[grepl("r2.fq.gz",fns)]

}else{

  fnRs <- fns[grepl("r2.fq.gz",fns)]

}
print(fnFs)

# Fix the names up, make sure there is matching damples
namesF<- sapply(strsplit(basename(fnFs), "_16S"), `[`, 1)
namesR<- sapply(strsplit(basename(fnRs), "_16S"), `[`, 1)
match(namesF,namesR)

names(fnFs)<-names(fnRs)<- namesF
if(any(is.na(match(namesF,namesR)))){
  warning("File Names Don't Match")
}

message(paste0("You have: ", length(fnFs)," samples"))
toplot<- sample(length(fnFs),ifelse(length(fnFs)<10,5,10))

qual_plots_fwd<-qual_plots_rvs<-  vector("list",length=length(toplot))
for(k in seq_along(toplot)){
  print(k)
  qual_plots_fwd[[k]]<- plotQualityProfile(fnFs[toplot[k]])
  qual_plots_rvs[[k]]<-  plotQualityProfile(fnRs[toplot[k]])
}



ggsave(arrangeGrob(grobs= qual_plots_fwd[1:9],ncol=3),
       file=paste0(run_path,"FWD_qc",ID,".png"),
       dpi=250,units="cm",height=30,width=30)


ggsave(arrangeGrob(grobs= qual_plots_rvs[1:9],ncol=3),
       file=paste0(run_path,"RVS_qc",ID,".png"),
       dpi=250,units="cm",height=30,width=30)

message("Plots made. ")

FWD <- "GACTACNVGGGTATCTAATCC"
REV<- "CCTACGGGNBGCASCAG"


## Make a list of primers at all orientations
library(ShortRead )
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
message("loadedshortreads. ")


FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

fnFs.filtN <- file.path(seq_path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(seq_path, "filtN", basename(fnRs))


filterAndTrim(fnFs, fnFs.filtN,
              fnRs, fnRs.filtN, maxN = 0, multithread =TRUE)

message("MERGED RAW reads. ")


primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

testsample<- sample(1:length(fnFs),1)
pre_trim<- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[testsample]]),
                 FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[testsample]]),
                 REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[testsample]]),
                 REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[testsample]]))


#cutadapt <- "~/.local/lib/python2.7/site-packages/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2("cutadapt", args = "--version") # Run shell commands from R



fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files

}


message("cut adapted. ")

post_trim<- rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[testsample]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[testsample]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[testsample]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[testsample]])
)

write(post_trim,paste0(run_path,ID,"_post_trim_check.txt"))
write(pre_trim,paste0(run_path,ID,"_pre_trim_check.txt"))

if(sum(grepl("fq",list.files(path.cut)))==0){
fq<- "fastq"
}else{
  fq<- "fq"
}

cutFs <- sort(list.files(path.cut, pattern = paste0("r1.",fq ,".gz"), full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern =paste0("r2.",fq ,".gz"), full.names = TRUE))



# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

namesF<- sapply(strsplit(basename(cutFs), "r1"), `[`, 1)
namesR<- sapply(strsplit(basename(cutRs), "r2"), `[`, 1)




out<- filterAndTrim(cutFs,filtFs ,
                    cutRs,filtRs ,
                    maxEE=c(2,3),
                    multithread = TRUE,
                    truncLen=c(270,240),maxN = 0,
                    compress=TRUE,verbose=TRUE,matchIDs = TRUE,
                    rm.phix = TRUE,  minLen = 50)

saveRDS(out,paste0(run_path,ID,"ReadsInOut.RDS",sep=""))

# Need to remove samples with low counts
# Usually they end up with 0 reads. This causes bugs downstream

remove <- out[,"reads.out"] < 100 # Or other cutoff
if(sum(remove==TRUE)==0){
  Fs<- filtFs
  Rs<- filtRs

}else{
  Fs <- file.path(filtFs[-(which(remove==TRUE))])
  Rs <- file.path(filtRs[-(which(remove==TRUE))])
}




samNames<- sapply(strsplit(basename(Fs),"_r1"),"[",1)
samNamesR<- sapply(strsplit(basename(Rs),"_r2"),"[",1)
match(samNames,samNamesR)

names(Fs)<-samNames
names(Rs)<-samNames
# Learn base pair transition error rates
errF <- learnErrors(Fs,nbases=1e8,randomize=TRUE,multithread = TRUE)
errR <- learnErrors(Rs,nbases=1e8,randomize=TRUE,multithread = TRUE)
pef<- plotErrors(errF, nominalQ=TRUE)
per<- plotErrors(errR, nominalQ=TRUE)
ggsave(pef,file=paste0(run_path,ID,"errorsF.png"),dpi=300)
ggsave(per,file=paste0(run_path,ID,"errorsR.png"),dpi=300)

errors<- list(fwd=errF,rev=errR)
saveRDS(errors,file=paste0(run_path,ID,"errors.RDS"))
## Made inference, dereplicate seqs and call SVs

## Made inference, dereplicate seqs and call SVs
reads_passed<-mergers <- vector("list", length(Fs))
names(mergers) <- samNames
for(sam in samNames) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(Fs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(Rs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)


  no_reads<- data.frame(denoisedF= getN(ddF),
                        denoisedR= getN(ddR),
                        merged= getN(merger))

  mergers[[sam]] <- merger
  reads_passed[[sam]]<- no_reads
}

seqtab<- makeSequenceTable(mergers)
saveRDS(seqtab,paste0(sv_path,ID,"seqtab.RDS"))
richness<- function(seq){apply(seq,1,function(x)sum(x>0))}
nochim.seqtab<-  dada2::removeBimeraDenovo(seqtab)

saveRDS(nochim.seqtab,paste0(sv_path,ID,"NCseqtab.RDS"))
chimera_sum<- 1-rowSums(nochim.seqtab)/rowSums(seqtab)
chimera_svs<- 1- richness(nochim.seqtab)/richness(seqtab)

chimeric_reads<- data.frame(chimera_sum= chimera_sum, chimeric_sv=chimera_svs,ID=rownames(nochim.seqtab) )
write.csv(chimeric_reads,paste0(run_path,ID,"chimeric_read_log.csv"))



#getN <- function(x) sum(getUniques(x))
track <- data.frame(out, plyr::ldply(reads_passed,data.frame), rowSums(nochim.seqtab),chimera_sum= chimera_sum, chimeric_sv=chimera_svs)

write.csv(track,paste0(run_path,ID,"dada2_read_log.csv"))
