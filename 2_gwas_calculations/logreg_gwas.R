### include library path
.libPaths("/ifs/data/nfs_share/sukmb276/GWAS/180718_German_mGWAS/190208_MicroGWAS/R_libraries/")

args <- commandArgs(trailingOnly = TRUE)
### define what feature to test
tax=args[1]
### level="otu97" or "otu99" or "asv" or taxonomic level
level=args[2]
### cohort="focus" or "popgen" or "kora" or "ship" or "shipt"
cohort=args[3]

### number of cores to use
ncores=4

library(snpStats)
library(MASS)
library(parallel)
library(phyloseq)

### folder holding all data
datafolder="/ifs/data/nfs_share/sukmb276/GWAS/180718_German_mGWAS/data/"

## output folder
outfolder=paste0("/home/sukmb276/Isilon/GWAS/180718_German_mGWAS/",cohort,"_logreg/LR_",tax)
dir.create(,recursive=T)

### load covariate data: age, gender, BMI and genetic data ID (genid)
covariate_rds=paste0(datafolder,cohort,"_covars.rds")
gwasdata<-readRDS(covariate_rds)

### load asv table phyloseq object
phylo_rds=paste0(datafolder,cohort,"_phyloseq.rds")
phylo<-readRDS(phylo_rds)

### load taxonomic assigments; also includes asv- and otu-level labels
tax_rds=paste0(datafolder,cohort,"_tax.rds")
tax_full<-readRDS(tax_rds)

### load genetic PCs
pca_file=paste0(datafolder,cohort,"_pca.rds")
pca<-read.table(pca_file,head=T)
if(cohort in c("popgen", "focus")){rownames(pca)<-pca$IID}
if(cohort in c("ship", "shipt")){pca<-pca[pca$FID!="missing",];rownames(pca)<-pca$FID}
if(cohort=="kora"){rownames(pca)<-paste0("KG",pca$IID)}


### extract asv table from phyloseq
asv=as.matrix(otu_table(phylo))

alltax=cbind(tax_full[colnames(asv),],tax_table(phylo)[match(colnames(asv), rownames(tax_table(phylo))),])

### get columns belonging to the specified feature 
taxcount=rowSums(otu[,alltax[,level]==tax & !is.na(alltax[,level]),drop=F])

### add to covariate data
gwasdata[,tax]<-taxcount[rownames(gwasdata)]

### data frame for GWAS
testdata=gwasdata[,c("age","gender","bmi","genid",tax)]
colnames(testdata)<-c("age","gender","bmi","genid","tax")

### make abundance binary
testdata$tax<-ifelse(testdata$tax>0,1,0)
testdata<-testdata[rowSums(is.na(testdata))==0,]

### add principal component data
if(cohort=="kora"){testdata=testdata[testdata$genid %in% rownames(pca),]}
pca_sub<-pca[testdata$genid,]
testdata<-data.frame(testdata, pca_sub)

### iterate over choromosomes
for(chrom in 1:22){
print(chrom)
files<-c(paste(paste0(datafolder,cohort,"_genotype_data/chr",chrom,".maf.clean"),c("bed","bim"),sep="."),paste0(datafolder,cohort,"_genotype_data/",cohort,".fam"))
plink<-read.plink(files[1],files[2],files[3])
if(cohort in c("popgen", "focus")){rownames(plink$genotypes)<-plink$fam$member}
if(cohort in c("ship", "shipt")){rownames(plink$genotypes)<-plink$fam$pedigree}
if(cohort=="kora"){rownames(plink$genotypes)<-paste0("KG",plink$fam$member)}

### convert genotypes to numerics (0, 1, 2)
gts<-as(plink$genotypes[testdata$genid,],"numeric")

#filter out all genotypes with MAF <= 2.5%; later in processing, only MAF > 5% is included
sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
gts<-gts[,sub]

### create output dataframe
out<-data.frame(plink$map[sub,],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,Beta=NA,StdErr=NA,Z=NA,P=NA)
colnames(out)<-c("chr","snp.name","cM","position","A","B","tax","n","AA","AB","BB","Beta","StdErr","Z","P")
out$tax<-tax

### do calculations

system.time(f<-mclapply(as.list(sub), function(snp){
df<-data.frame(testdata,gt=gts[,snp])
df<-df[is.na(df$gt)==F & is.na(df$tax)==F,]
if(nrow(df)>100 & length(unique(df$gt))>1){
  res<-c(nrow(df),table(factor(df$gt,levels=c(0,1,2))),tryCatch(summary(glm(tax ~ age + gender + bmi + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + gt , data=df,family = binomial(link = "logit")))$coefficients[15,], error=function(x) return(rep(NA,4))))
  names(res)<-c("n","AA","AB","BB","Beta","StdErr","Z","P")
}else{
  res<-rep(NA,8)
  names(res)<-c("n","AA","AB","BB","Beta","StdErr","Z","P")
}
return(res)},mc.cores=getOption("mc.cores", ncores)))

### collect results
out[,c("n","AA","AB","BB","Beta","StdErr","Z","P")]<-data.frame(do.call(rbind, f))
print(head(out))

### save results
saveRDS(out,paste0(outfolder,"/chr",chrom,".pca.rds"))

}
