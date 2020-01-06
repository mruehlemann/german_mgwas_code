### This code was run on a machine with GPUs available for calculation
### 
### Starting R as "CUDA_VISIBLE_DEVICES=[device_number] R" selects which device to use, if multiple are available

### arguments are used to split the dataset into subsets for parallel calculation

args = commandArgs(trailingOnly=TRUE)
## chromosome
chrom = args[1]
## cohort
cohort = args[2]
## "set" defines from which CPUs should be selected to minimize overlap in calculation using multiple GPUs 
set = args[3]
## which chunk
chunk = as.numeric(args[4])
## total number of chunks for the chromosome
nchunks = as.numeric(args[5])
## dist; should be "bray" or "unifrac"
dist=args[6]

usecores=1L

setwd(paste0("/home/mruehlemann/DBF_GPU/",cohort))

library(snpStats) ### über Bioconductor
library(combinat)
library(snpStats)
library(MASS)
library(parallel)

#source("/home/mruehlemann/DBF_GPU_test/DBF_functions.Rscript")

"trace.matrix"<- function(x)
{
        return(sum(diag(x)))
}


dm<-function(genotype, distmat){
require(gmatrix)
#genotype[genotype==0]<-NA
allgroups<-unique(na.omit(genotype))
no.grps<-length(allgroups)
ind.grps <- lapply(1:no.grps, function(i) which(genotype==allgroups[i]))
genotype<-g(as.numeric(genotype))
distmat<-g(distmat)
n<-length(genotype)
n.h<-as.numeric(n)
d2 <- -(1/2)*distmat^2
ones <- gmatrix(data=1/n.h,ncol=n.h,nrow=n.h, dup=FALSE)
centmat <- (diag(n.h)-ones)
G <- centmat%*%d2%*%centmat
L <- gmatrix(0,nrow=n.h,ncol=no.grps)
for (i in 1:no.grps){
                L[as.numeric(ind.grps[[i]]),i] <- 1
        }
H<-g(scale(L%*%solve(t(L)%*%L)%*%t(L),center=TRUE,scale=FALSE))
dbf0<-h(1/((trace.matrix(G)/trace.matrix(H%*%G))-1))
Am2 <- H%*%H
Am3 <- H%*%Am2
dA <- diag(H)
Ap2 <- H^2
Ap3 <- Ap2*H
A1 <- h(sum(dA))
A2 <- h(trace.matrix(Am2))
A3 <- h(trace.matrix(Ap2))
A4 <- h(trace.matrix(Am3))
A5 <- h(trace.matrix(Ap3))
A6 <- h(sum(Ap3))
A7 <- h(as.gvector(dA%*%diag(Am2)))
A8 <- h(as.gvector(dA%*%H%*%dA))
Bm2 <- G%*%G
Bm3 <- G%*%Bm2
dB <- diag(G)
Bp2 <- G^2
Bp3 <- Bp2*G
B1 <- h(sum(dB))
B2 <- h(trace.matrix(Bm2))
B3 <- h(trace.matrix(Bp2))
B4 <- h(trace.matrix(Bm3))
B5 <- h(trace.matrix(Bp3))
B6 <- h(sum(G^3))
B7 <- h(as.gvector(dB%*%diag(Bm2)))
B8 <- h(as.gvector(dB%*%G%*%dB))

n<-as.numeric(n)
n2 <- n.h^2
n3 <- n.h^3
n4 <- n.h^4


mom1<-h(trace.matrix(H)*trace.matrix(G)/(n-1))
mom2<-((n-1)*A3*B3+(B1^2-B3)*(A1^2-A3)+2*(A2-A3)*(B2-B3)+4*A3*B3)/(n*(n-1))+(4*(n-3)*(2*A3-A2)*(2*B3-B2)+2*(2*A3-A1^2)*(2*B3-B1^2)*(n-3)+(2*A2+A1^2-6*A3)*(2*B2+B1^2-6*B3))/(n*(n-1)*(n-2)*(n-3))
mom3<-(n2*(n+1)*(n2+15*n-4)*A5*B5 + 4*(n4-8*n3+19*n2-4*n-16)*A6*B6 + 24*(n2-n-4)*(A6*B8+B6*A8)+ 6*(n4-8*n3+21*n2-6*n-24)*A8*B8 + 12*(n4-n3-8*n2+36*n-48)*A7*B7 + 12*(n3-2*n2+9*n-12)*(A1*A3*B7 + A7*B1*B3) + 3*(n4 - 4*n3 - 2*n2+9*n-12)*A1*B1*A3*B3 + 24*((n3 - 3*n2 - 2*n+8)*(A7*B6 + A6*B7) + (n3 - 2*n2 - 3*n+12)*(A7*B8 + A8*B7)) + 12*(n2 - n + 4)*(A1*A3*B6 + B1*B3*A6) + 6*(2*n3 - 7*n2 - 3*n + 12)*(A1*A3*B8 + A8*B1*B3) - 2*n*(n-1)*(n2-n+4)*((2*A6+3*A8)*B5+(2*B6+3*B8)*A5) - 3*n*(n-1)*(n-1)*(n+4)*((A1*A3+4*A7)*B5+(B1*B3+4*B7)*A5) + 2*n*(n-1)*(n-2)*((A1^3 + 6*A1*A2 + 8*A4)*B5+(B1^3 + 6*B1*B2 + 8*B4)*A5) + (A1^3)*((n3-9*n2+23*n-14)*(B1^3)+6*(n-4)*B1*B2+8*B4) + 6*A1*A2*((n-4)*(B1^3)+(n3-9*n2+24*n-14)*B1*B2 + 4*(n-3)*B4) + 8*A4*((B1^3)+3*(n-3)*B1*B2+(n3-9*n2+26*n-22)*B4) - 16*((A1^3)*B6+A6*(B1^3)) - 6*(A1*A2*B6 + A6*B1*B2)*(2*n2-10*n+16) - 8*(A4*B6+A6*B4)*(3*n2-15*n+16)-((A1^3)*B8+A8*(B1^3))*(6*n2-30*n+24)-6*(A1*A2*B8+A8*B1*B2)*(4*n2-20*n+24) - 8*(A4*B8+A8*B4)*(3*n2-15*n+24) - (n-2)*(24*((A1^3)*B7+A7*(B1^3))+6*(A1*A2*B7+A7*B1*B2)*(2*n2-10*n+24)+8*(A4*B7+A7*B4)*(3*n2-15*n+24)+(3*n2-15*n+6)*((A1^3)*B1*B3+A1*A3*(B1^3))+6*(A1*A2*B1*B3+A1*A3*B1*B2)*(n2-5*n+6) + 48*(A4*B1*B3+A1*A3*B4)))/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
mean.T <- mom1
variance.T <- (mom2 - mom1^2)
skewness.T <- (mom3 - 3 * variance.T * mom1 - mom1^3)/(variance.T^(3/2))
tr.G <- h(trace.matrix(G))
aa<-list(mean.T=mean.T,variance.T=variance.T,skewness.T=skewness.T,tr.G=tr.G)
p.dbf0 <- 1- cdf.F(aa,dbf0)
return(c(dbf.statistic=dbf0,dbf.p.value=p.dbf0))
}



### load covariate data: age, gender, BMI and genetic data ID (genid)
covariate_rds=paste0(datafolder,cohort,"_covars.rds")
gwasdata<-readRDS(covariate_rds)


### load covariate data: age, gender, BMI and genetic data ID (genid)
residuals_rds=paste0(datafolder,cohort,"_",dist,"_residuals.rds")
cap.resid<-readRDS(residuals_rds)

#cap.resid<-readRDS("/home/mruehlemann/DBF_GPU/focus/FOCUS_weightedUniFrac.residuals.Rds")
files<-c(paste(paste0(datafolder,cohort,"genotype_data/chr",chrom,".maf.clean"),c("bed","bim","fam"),sep="."))
plink<-read.plink(files[1],files[2],files[3])
if(cohort in c("popgen", "focus")){rownames(plink$genotypes)<-plink$fam$member}
if(cohort in c("ship", "shipt")){rownames(plink$genotypes)<-plink$fam$pedigree}
if(cohort=="kora"){rownames(plink$genotypes)<-paste0("KG",plink$fam$member)}

### convert genotypes to numerics (0, 1, 2)
gts<-as(plink$genotypes[gwasdata$genid,],"numeric")

map=plink$map

sub<-colnames(gts)[colMeans(gts,na.rm=T)/2>0.025 & colMeans(gts,na.rm=T)/2<0.975 & is.na(colMeans(gts,na.rm=T))==F]
print(head(sub))

gts=gts[,sub]

sub<-sub[1:length(sub)%%nchunks==(chunk-1)]

print(length(sub))

gts<-gts[,sub]

out<-data.frame(map[sub,],tax=NA,n=NA,AA=NA,AB=NA,BB=NA,stat=NA,P=NA)
colnames(out)<-c("chr","snp.name","cM","position","A","B","tax","n","AA","AB","BB","stat","P")
out$tax=dist

mcaffinity((1:32)[1:32%%4==set])

system.time(f<-mclapply(as.list(sub), function(snp){
thisgt<-as(gts[,snp],"numeric")
thisgt.sub<-thisgt[!is.na(thisgt)]
if(length(unique(thisgt.sub))>1 & length(thisgt.sub)>100){
thisdmat<-cap.resid[names(thisgt.sub),names(thisgt.sub)]
thisn<-nrow(thisdmat)
thisgtdist<-as.numeric(table(factor(thisgt,levels=c(0,1,2))))

#thisdbf<-DBF.test(thisdmat,thisgt.sub,thisn)
thisdbf<-dm(thisgt.sub,thisdmat)

res<-c(thisn,as.numeric(thisgtdist),as.numeric(thisdbf))
}else{res<-rep(NA,6)}
names(res)<-c("n","AA","AB","BB","stat","P");return(res)},mc.cores=getOption("mc.cores", usecores)))


out[,c("n","AA","AB","BB","stat","P")]<-data.frame(do.call(rbind, f))
print(head(out))

dir.create("results/BetaUnifrac",recursive=T)
saveRDS(out,paste0("results/",dist,"/",cohort,"/chr",chrom,".chunk.",chunk,".of.",nchunks,".rds"))
