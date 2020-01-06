### modified from: https://benjjneb.github.io/dada2/tutorial.html
### Start script by running Rscript dada2_16s_workflow.R [Input folder with fastq files] [ID for run] [output folder]
### The input folder must contain forward and reverse read for the samples 
### if no RunID is given, it is generated from the standard folder structure
### if no output folder is given, standard folder defined in "outbase_std" is used

args <- commandArgs(trailingOnly = TRUE)

### standard folder for output
outbase_std="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/Runs_v.1.10"

path=gsub("/$","",args[1])
if(is.na(args[2])){runid=strsplit(rev(strsplit(path, split="/")[[1]])[2],split="_")[[1]][1]}else{runid=args[2]}
if(is.na(args[3])){outbase=outbase_std}else{outbase=gsub("/$","",args[3])}

numbases=1e+09
threads=8

outdir=paste0(outbase,"/",runid)

cat(paste0("All output will be saved in: ",outdir,"\n"))

dir.create(outdir,recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/plots"),recursive=T,showWarnings=F)
dir.create(paste0(outdir,"/errors"),recursive=T,showWarnings=F)

### inclued R library folder
.libPaths("/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/R_Librarires/")

library(dada2)
version=packageVersion("dada2")
set.seed(666)
cat(paste0("Running DADA2 version ",version,"\n"))

fns <- list.files(path)
#fns

####################
### load your data
####################

fastqs <- fns[grepl(".fastq.gz$", fns)]
fastqs <- sort(fastqs) # Sort ensures forward/reverse reads are in same order
### make sure that R1 is for forward read and R2 for reverse

fnFs <- fastqs[grepl("R1_001.fastq.gz", fastqs)] ## Just the forward read files
fnRs <- fastqs[grepl("R2_001.fastq.gz", fastqs)] ## Just the reverse read files

## Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, "(-L1)*_((S[0-9]+)|([ACTG-]+))_L001_R1_001.fastq.gz"), `[`, 1)
cat(paste0("Starting processing for ",length(sample.names)," samples\n"))

## Fully specify the path for the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)


###########################################
## Examine quality profiles of F & R reads
###########################################
pdf(paste0(outdir,"/plots/plotQualityProfile.pdf"), onefile=T)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
dev.off()

##################################
## Perform filtering and trimming
##################################
filt_path <- file.path(outdir, "filtered") # Place filtered files in filtered/subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

## Filter  the forward and reverse reads:
## Important to remove primers and low quality regions
## adjusted for V1-V2 (27f/338r)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,180), 
                     trimLeft=c(5, 5),
                     maxN=0, maxEE=c(2,2), truncQ=5, rm.phix=TRUE,
                     compress=TRUE, multithread=threads) #


## Examine quality profiles of filtered reads
pdf(paste0(outdir,"/plots/plotQualityProfile.filt.pdf"), onefile=T)
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])
dev.off()

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

#########################
## Learn the Error Rates
#########################
## Learn forward error rates
errF <- learnErrors(filtFs, nbases=numbases, multithread=threads,randomize=T)
## Learn reverse error rates
errR <- learnErrors(filtRs, nbases=numbases, multithread=threads,randomize=T)

### save error profiles
saveRDS(errF, paste0(outdir,"/errors/errF.Rds"))
saveRDS(errR, paste0(outdir,"/errors/errR.Rds"))


## Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

## Plot estimated error as sanity check
pdf(paste0(outdir,"/plots/plotErrors_F.pdf"), onefile=T)
plotErrors(errF, nominalQ=TRUE)
dev.off()

pdf(paste0(outdir,"/plots/plotErrors_R.pdf"), onefile=T)
plotErrors(errR, nominalQ=TRUE)
dev.off()


#########################
## Perform dereplication
#########################
## Dereplicate the filtered fastq files
derepRs <- derepFastq(filtRs, verbose=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)

sample.names=sample.names[exists]

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


####################
## Sample Inference
####################
## Apply the core sequence-variant inference algorithm to the dereplicated data
## Infer the sequence variants in each sample
dadaFs <- dada(derepFs, err=errF, multithread=threads)
dadaRs <- dada(derepRs, err=errR, multithread=threads)


## Inspect the dada-class object returned by dada
#dadaFs[[1]]

## Merge the denoised forward and reverse reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=10)
## Inspect the merger data.frame from the first sample
#head(mergers[[1]])


############################
## Construct sequence table
############################
seqtab <- makeSequenceTable(mergers)
## Get dimensions
#dim(seqtab)

## Inspect distribution of sequence lengths
#table(nchar(getSequences(seqtab)))


saveRDS(seqtab, paste0(outdir,"/seqtab.Rds"))

###################
## Remove chimeras
###################
## Remove chimeric sequences:
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)
#dim(seqtab.nochim)
#sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim, paste0(outdir,"/seqtab_nochim.Rds"))

####################################
## Track reads through the pipeline
####################################
getN <- function(x) sum(getUniques(x))
track <- cbind(out[exists,], sapply(dadaFs, getN), sapply(m2, getN), rowSums(seqtab), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
#head(track)

write.table(track, paste0(outdir,"/track_reads.txt"),sep="\t",quote=FALSE)

###################
## Assign taxonomy
###################
database_file="/ifs/data/nfs_share/sukmb276/Microbiome/clean_data_from_dada2/reference_data/rdp_train_set_16.fa.gz"
taxHS <- assignTaxonomy(seqtab.nochim, database_file, multithread=threads) ## CHANGE to directory and pertinent database
colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
#unname(head(taxHS))
#unname(tail(taxHS))


write.table(taxHS, file = paste0(outdir,"/taxa_SV.tsv"), quote=FALSE)
