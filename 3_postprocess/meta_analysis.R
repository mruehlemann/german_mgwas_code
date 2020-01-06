args <- commandArgs(trailingOnly = TRUE)
tax2=args[1]
### either NB (abundance) or LR (logistic regression)
model=args[2]

model2=ifelse(model=="LR","logreg","negbin")

calcmaf<-function(x){
n=x[1]
het=x[3]
if(x[2]>x[4]){hom_maj=x[2];hom_min=x[4]}else{hom_maj=x[4];hom_min=x[2]}
maf=(het+2*hom_min)/(2*n)
return(maf)
}

allbsp<-data.frame(read.table(paste("popgen_",model2,"/popgen_all_traits")),cohort="POPGEN",num=1)
allfocus<-data.frame(read.table(paste("focus_",model2,"/focus_all_traits")),cohort="FOCUS",num=1)
allkora<-data.frame(read.table(paste("kora_",model2,"/kora_all_traits")),cohort="KORA",num=1)
allship<-data.frame(read.table(paste("ship_",model2,"/ship_all_traits")),cohort="SHIP",num=1)
allshipt<-data.frame(read.table(paste("shipt_",model2,"/shipt_all_traits")),cohort="SHIPT",num=1)


library(reshape2)

cc<-rbind(allbsp,allfocus,allkora,allship,allshipt)
cctab<-data.frame(dcast(cohort ~ V1, data=cc,fill=0),row.names=1)
print(cctab)

#for(tax2 in colnames(cctab)){
print(tax2)
tax=paste0(model,"_",tax2)

outfolder=paste0("/ifs/data/nfs_share/sukmb276/GWAS/180718_German_mGWAS/META_",model2,"/",tax)

print(cctab[,tax2,drop=F])

for(chrom in seq(1,22)){
print(chrom)

alllist=list()
hits=c()
cohorts=c()

if(cctab["POPGEN",tax2,drop=F]==1){
bspspcfile=paste0("popgen_",model2,"/",tax,"/chr",chrom,".pca.rds")
if(file.exists(bspspcfile)==F){print("POPGEN missing");next}}

if(cctab["FOCUS",tax2,drop=F]==1){
focusfile=paste0("focus_",model2,"/",tax,"/chr",chrom,".pca.rds")
if(file.exists(focusfile)==F){print("FOCUS missing");next}}

if(cctab["KORA",tax2,drop=F]==1){
korafile=paste0("kora_",model2,"/",tax,"/chr",chrom,".pca.rds")
if(file.exists(korafile)==F){print("KORA missing");next}}

if(cctab["SHIP",tax2,drop=F]==1){
shipfile=paste0("ship_",model2,"/",tax,"/chr",chrom,".pca.rds")
if(file.exists(shipfile)==F){print("SHIP missing");next}}

if(cctab["SHIPT",tax2,drop=F]==1){
shiptfile=paste0("shipt_",model2,"/",tax,"/chr",chrom,".pca.rds")
if(file.exists(shiptfile)==F){print("SHIPT missing");next}}


if(cctab["POPGEN",tax2,drop=F]==1){
if(file.exists(bspspcfile)){
bspspc<-readRDS(bspspcfile)
cohorts<-c(cohorts,"bspspc")
bspspc=bspspc[is.na(bspspc$n)==F,]
bspspc$MAF<-apply(bspspc[,c("n","AA","AB","BB")],1,calcmaf)
bspspc<-bspspc[bspspc$MAF>=0.05,]
bspspc$BETA2<-ifelse(bspspc$A > bspspc$B, -bspspc$Beta, bspspc$Beta)
bspspc$w<-1/(bspspc$StdErr^2)
bspspc$BETA.w<-bspspc$w*bspspc$BETA2
bspspc.lambda<-read.table(paste0("bspspc_full/",tax,"/lambda.pca.txt"))
if(bspspc.lambda$lambda>1){bspspc$Z.adj<-(bspspc$Z/bspspc.lambda$lambda)*(sign(bspspc$Beta)*sign(bspspc$BETA2))}else{bspspc$Z.adj<-bspspc$Z*(sign(bspspc$Beta)*sign(bspspc$BETA2))}
bspspc$P.adj<-2*pnorm(-abs(bspspc$Z.adj))
bspspc$wN<-sqrt(bspspc$n)
bspspc$wZ<-bspspc$wN*bspspc$Z.adj
hits<-unique(c(hits,bspspc$snp.name))
colnames(bspspc)<-paste0("bspspc.",colnames(bspspc))
alllist[[length(alllist)+1]]<-bspspc
}else{print("POPGEN missing");next}}

if(cctab["FOCUS",tax2,drop=F]==1){
if(file.exists(focusfile)){
focus<-readRDS(focusfile)
cohorts<-c(cohorts,"focus")
focus=focus[is.na(focus$n)==F,]
focus$MAF<-apply(focus[,c("n","AA","AB","BB")],1,calcmaf)
focus<-focus[focus$MAF>=0.05,]
focus$BETA2<-ifelse(focus$A > focus$B, -focus$Beta, focus$Beta)
focus$w<-1/(focus$StdErr^2)
focus$BETA.w<-focus$w*focus$BETA2
focus.lambda<-read.table(paste0("focus_full/",tax,"/lambda.pca.txt"))
if(focus.lambda$lambda>1){focus$Z.adj<-(focus$Z/focus.lambda$lambda)*(sign(focus$Beta)*sign(focus$BETA2))}else{focus$Z.adj<-focus$Z*(sign(focus$Beta)*sign(focus$BETA2))}
focus$P.adj<-2*pnorm(-abs(focus$Z.adj))
focus$wN<-sqrt(focus$n)
focus$wZ<-focus$wN*focus$Z.adj
hits<-unique(c(hits,focus$snp.name))
colnames(focus)<-paste0("focus.",colnames(focus))
alllist[[length(alllist)+1]]<-focus
}else{print("FOCUS missing");next}}

if(cctab["KORA",tax2,drop=F]==1){
if(file.exists(korafile)){
kora<-readRDS(korafile)
cohorts<-c(cohorts,"kora")
kora=kora[is.na(kora$n)==F,]
kora$MAF<-apply(kora[,c("n","AA","AB","BB")],1,calcmaf)
kora<-kora[kora$MAF>=0.05,]
kora$BETA2<-ifelse(kora$A > kora$B, -kora$Beta, kora$Beta)
kora$w<-1/(kora$StdErr^2)
kora$BETA.w<-kora$w*kora$BETA2
kora.lambda<-read.table(paste0("kora_v12_full/",tax,"/lambda.pca.txt"))
if(kora.lambda$lambda>1){kora$Z.adj<-(kora$Z/kora.lambda$lambda)*(sign(kora$Beta)*sign(kora$BETA2))}else{kora$Z.adj<-kora$Z*(sign(kora$Beta)*sign(kora$BETA2))}
kora$P.adj<-2*pnorm(-abs(kora$Z.adj))
kora$wN<-sqrt(kora$n)
kora$wZ<-kora$wN*kora$Z.adj
hits<-unique(c(hits,kora$snp.name))
colnames(kora)<-paste0("kora.",colnames(kora))
alllist[[length(alllist)+1]]<-kora
}else{print("KORA missing");next}}

if(cctab["SHIP",tax2,drop=F]==1){
if(file.exists(shipfile)){
ship<-readRDS(shipfile)
cohorts<-c(cohorts,"ship")
ship=ship[is.na(ship$n)==F,]
ship$MAF<-apply(ship[,c("n","AA","AB","BB")],1,calcmaf)
ship<-ship[ship$MAF>=0.05,]
ship$BETA2<-ifelse(ship$A > ship$B, -ship$Beta, ship$Beta)
ship$w<-1/(ship$StdErr^2)
ship$BETA.w<-ship$w*ship$BETA2
ship$snp.name<-paste0(ship$chr,":",ship$position)
ship.lambda<-read.table(paste0("ship_full/",tax,"/lambda.pca.txt"))
if(ship.lambda$lambda>1){ship$Z.adj<-(ship$Z/ship.lambda$lambda)*(sign(ship$Beta)*sign(ship$BETA2))}else{ship$Z.adj<-ship$Z*(sign(ship$Beta)*sign(ship$BETA2))}
ship$P.adj<-2*pnorm(-abs(ship$Z.adj))
ship$wN<-sqrt(ship$n)
ship$wZ<-ship$wN*ship$Z.adj
hits<-unique(c(hits,ship$snp.name))
colnames(ship)<-paste0("ship.",colnames(ship))
alllist[[length(alllist)+1]]<-ship
}else{print("SHIP missing");next}}

if(cctab["SHIPT",tax2,drop=F]==1){
if(file.exists(shiptfile)){
shipt<-readRDS(shiptfile)
cohorts<-c(cohorts,"shipt")
shipt=shipt[is.na(shipt$n)==F,]
shipt$MAF<-apply(shipt[,c("n","AA","AB","BB")],1,calcmaf)
shipt<-shipt[shipt$MAF>=0.05,]
shipt$BETA2<-ifelse(shipt$A > shipt$B, -shipt$Beta, shipt$Beta)
shipt$w<-1/(shipt$StdErr^2)
shipt$BETA.w<-shipt$w*shipt$BETA2
shipt$snp.name<-paste0(shipt$chr,":",shipt$position)
shipt.lambda<-read.table(paste0("shipt_full/",tax,"/lambda.pca.txt"))
if(shipt.lambda$lambda>1){shipt$Z.adj<-(shipt$Z/shipt.lambda$lambda)*(sign(shipt$Beta)*sign(shipt$BETA2))}else{shipt$Z.adj<-shipt$Z*(sign(shipt$Beta)*sign(shipt$BETA2))}
shipt$P.adj<-2*pnorm(-abs(shipt$Z.adj))
shipt$wN<-sqrt(shipt$n)
shipt$wZ<-shipt$wN*shipt$Z.adj
hits<-unique(c(hits,shipt$snp.name))
colnames(shipt)<-paste0("shipt.",colnames(shipt))
alllist[[length(alllist)+1]]<-shipt
}else{print("SHIPT missing");next}}

for(i in 1:length(alllist)){
if(i==1){
all<-alllist[[1]][match(hits, alllist[[1]][,paste0(cohorts[1],".snp.name")]),]
}else{
all<-cbind(all,alllist[[i]][match(hits, alllist[[i]][,paste0(cohorts[i],".snp.name")]),])
}}

rownames(all)<-hits
### inverse-variance weighted meta analysis
all$META.se<-sqrt(1/rowSums(all[,paste0(cohorts,".w")],na.rm=T))
all$META.BETA<-rowSums(all[,paste0(cohorts,".BETA.w")],na.rm=T)/rowSums(all[,paste0(cohorts,".w")],na.rm=T)
all$META.Z<-all$META.BETA/all$META.se
all$META.P<-2*pnorm(-abs(all$META.Z))
all$META.N<-rowSums(all[,paste0(cohorts,".n")],na.rm=T)
all$META.numcohorts<-rowSums(is.na(all[,paste0(cohorts,".n")])==F,na.rm=T)
all$META.numsigcohorts<-rowSums(all[,paste0(cohorts,".P")] < 0.05,na.rm=T)
all$META.sigcohorts<-apply(all[,paste0(cohorts,".P")], 1, function(x) paste(cohorts[which(x<0.05)],collapse=","))

### sample size weighted meta analysis (not used, as all Lambda_GC were < 1.05)
all$META2.wZ<-rowSums(all[,paste0(cohorts,".wZ")],na.rm=T)
all$META2.wN<-sqrt(rowSums(all[,paste0(cohorts,".wN")]^2,na.rm=T))
all$META2.Z<-all$META2.wZ/all$META2.wN
all$META2.P<-2*pnorm(-abs(all$META2.Z))
all$META2.numsigcohorts<-rowSums(all[,paste0(cohorts,".P.adj")] < 0.05,na.rm=T)
all$META2.sigcohorts<-apply(all[,paste0(cohorts,".P.adj")], 1, function(x) paste(cohorts[which(x<0.05)],collapse=","))

all<-all[all$META.numcohorts>0,]

dir.create(outfolder,recursive=T)
saveRDS(all,paste0(outfolder,"/chr",chrom,".pca.rds"))

}

