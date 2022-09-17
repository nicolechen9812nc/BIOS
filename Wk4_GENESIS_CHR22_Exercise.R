## Install GENESIS

## With R version 3.5 or greater, need BiocManager to intall Bioconductor packages
## The current release of Bioconductor is version 3.15; it works with R version 4.2. 
## Install R version 4.2.1 (Windows: https://cran.r-project.org/bin/windows/base/; Mac: https://cran.r-project.org/bin/macosx/base/) 
## Then install BiocManager using the code below
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")
BiocManager::install("GENESIS")
## You might also need to install these
BiocManager::install("SeqArray")
BiocManager::install("SeqVarTools")
BiocManager::install("Biobase")
install.packages("SPAtest")
install.packages("qqman")

## With R version < 3.5, install GENESIS using the code below
## Try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GENESIS")

## Load GENESIS 
## Manual can be downloaded here: https://bioconductor.org/packages/release/bioc/manuals/GENESIS/man/GENESIS.pdf
library(GENESIS)


################################START##############################
library(GENESIS)
library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(SPAtest)
library(qqman)

## Set working directory, clean work space
setwd("C:/Users/nicole/Desktop/NINI")
rm(list=ls())

## Open the data, preview the contents
data<-seqOpen("CHR22.gds")
data

## Get sample IDs, know the sample size N=195, preview sample IDs
sample.id<-seqGetData(data, "sample.id")
length(sample.id)
head(sample.id)

## Get SNP IDs, know the total number of SNPs is 5199, a unique integer ID is assigned to each variant
variant.id<-seqGetData(data, "variant.id")
length(variant.id)
head(variant.id)

## Get SNP IDs - rs#, aligned/ordered the same
rs<-seqGetData(data,"annotation/id")
length(rs)
head(rs)

## Get SNP positions, aligned/ordered the same
pos<-seqGetData(data,"position")
length(pos)
head(pos)

## Get alleles, aligned/ordered the same
allele<-seqGetData(data,"allele")
length(allele)
head(allele)

## Get alternative/effect alleles, aligned/ordered the same
alt<-seqGetData(data,"$alt")
length(alt)
head(alt)

## Construct a SNP annotation data frame
annot<-data.frame(variant.id=variant.id,rs=rs,position=pos,allele=allele,alt=alt,stringsAsFactors=F)

## Get genotype, know the dimension, view subset of the data
## Genotype data is 3-dimensional
## The first dimension is always 2 for diploid genotypes
## The second and third dimensions are samples and variants respectively
geno<-seqGetData(data,"genotype")
dim(geno)
geno[,1:5,20:25]

## Get missing rate of genotype across 5199 SNPs, draw a histogram
## get the SNPs for exclusion if missing >10%
byvariant<-missingGenotypeRate(data,margin="by.variant")
length(byvariant)
hist(byvariant,breaks=50)
missingbyvariant<-byvariant>0.1 #set missing rate>0.1
table(missingbyvariant) #675 TRUE
byvariant_exclude<-variant.id[missingbyvariant]

## Get minor allele frequency of each variant, draw a histogram
## Do you observe rare variants? 
## Let's get the rare variants for exclusion as QC later
maf<-seqAlleleFreq(data,minor=T)
hist(maf,breaks=10)
rare<-maf<0.1 #set MAF<0.1
table(rare) #1777 rare
maf_exclude<-variant.id[rare]

## Run a Hardy-Weinberg Equilibrium test on each variant
## variant.id is unique identifier for the variant
## nAA is number of reference homozygotes
## nAa is number of heterozygotes
## naa is number of alternate homozygotes
## afreq is reference allele frequency
## p is p values for the exact test
## f is inbreeding coefficient, 1 - observed heterozygosity / expected heterozygosity
hwe.res<-hwe(data)
dim(hwe.res)
head(hwe.res)

## Find the SNPs with low HWE p value
## Let's get the variants with p<1e-4 for exclusion as a QC  ##choose?
lowp<-!is.na(hwe.res$p)&hwe.res$p<1e-6 #set HWE p-value<1e-6
table(lowp) #1
hwe_exclude<-variant.id[lowp]

## Integrate all SNPs together for exclusion
excludeSNP<-c(byvariant_exclude,maf_exclude,hwe_exclude)
excludeSNP<-unique(excludeSNP)
length(excludeSNP) #2219

## Get the SNPs to keep after quality control
keepSNP<-setdiff(variant.id,excludeSNP)
length(keepSNP) #2980=5199-2219

## Define a filter, use the SNPs passed QC to conduct sample QC
## If you would like to reset filter, use seqResetFilter(data)
seqSetFilter(data,variant.id=keepSNP)

## Get missing rate of genotype across 195 samples, draw a histogram
## get the samples for exclusion if missing >10%
bysample<-missingGenotypeRate(data,margin="by.sample")
length(bysample)
hist(bysample,breaks=50)
missingbysample<-bysample>0.1 #set missing rate>0.1
table(missingbysample) #1
bysample_exclude<-sample.id[missingbysample]

## Get the samples to keep after quality control
keepSample<-setdiff(sample.id,bysample_exclude)
length(keepSample) #194=195-1

## Read the phenotype, know the variables available
## Our outcome1 is triglyceride (continuous), covariates are age and sex
## Our outcome2 is Obese (yes or no), covariates are age and sex
## We would like to restrict the analysis to Caucasians
pheno<-read.csv("CHR22_pheno.csv",stringsAsFactors=F)
dim(pheno)
pheno[1:5,]
names(pheno)
hist(pheno$Trig)
hist(pheno$Age)
table(pheno$sex)
hist(pheno$BMI)
pheno$Trig_high<-ifelse(pheno$Trig>5.65,1,0)
table(pheno$Trig_high)

######### Single-variant test 1: BMI ~ SNP + Age + Sex #########
######### Single-variant test 2: Trig_high ~ SNP + Age + Sex #########
## First fit null model Trig ~ Age + Sex, Obese ~ Age + Sex 
## Then pass the null model results to variants testing
pheno<-as(pheno,"AnnotatedDataFrame") #Create an annotated dataframe
varMetadata(pheno)
head(pData(pheno))
summary(pData(pheno))
table(pData(pheno)$inferred_population)

##Get the IDs for Caucasians
CaucID<-pData(pheno)[pData(pheno)$inferred_population=="Cauc" & !is.na(pData(pheno)$inferred_population),]$sample.id
keepSample<-intersect(keepSample,CaucID) #178 samples to keep

## Single-variant test 1: BMI ~ SNP + Age + Sex
seqSetFilter(data,variant.id=keepSNP,sample.id=keepSample) #Define a filter, analyze data passed QC. To reset filter, use seqResetFilter(data)
seqData<-SeqVarData(data,sampleData=pheno) #Create seqData - a combination of geno & pheno
iterator<-SeqVarBlockIterator(seqData) #Create SeqVarIterator object for input of testing

nullmod1<-fitNullModel(iterator,outcome="BMI",covars=c("Age","sex"),family="gaussian") #Fit null model
assoc1<-assocTestSingle(iterator,nullmod1,test="Score") #Fit final model
head(assoc1)
summary(assoc1$Score.pval)

assoc1<-merge(annot,assoc1,by="variant.id",all.y=T)
assoc1<-assoc1[order(assoc1$Score.pval),]

## Single-variant test 2: Obese ~ SNP + Age + Sex
seqSetFilter(data,variant.id=keepSNP,sample.id=keepSample) #Define a filter, analyze data passed QC. To reset filter, use seqResetFilter(data)
seqData<-SeqVarData(data,sampleData=pheno) #Create seqData - a combination of geno & pheno
iterator<-SeqVarBlockIterator(seqData) #Create SeqVarIterator object for input of testing

nullmod2<-fitNullModel(iterator,outcome="Trig_high",covars=c("Age","sex"),family="binomial") #Fit null model
assoc2<-assocTestSingle(iterator,nullmod2,test="Score.SPA",recalc.pval.thresh=1) #Fit final model
head(assoc2)
summary(assoc2$SPA.pval)

assoc2<-merge(annot,assoc2,by="variant.id",all.y=T)
assoc2<-assoc2[order(assoc2$SPA.pval),]

## Draw QQ plot to visualize the results
qq(assoc1$Score.pval)
qq(assoc2$SPA.pval)

## Draw Manhattan plot to visualize the results
assoc1$chr<-as.numeric(assoc1$chr)
manhattan(assoc1,chr="chr",bp="pos",p="Score.pval",snp="rs",suggestiveline=-log10(1e-2),genomewideline=-log10(5e-08),xlim=c(min(assoc1$pos)-1e6,max(assoc1$pos)+1e6))
assoc2$chr<-as.numeric(assoc2$chr)
manhattan(assoc2,chr="chr",bp="pos",p="SPA.pval",snp="rs",suggestiveline=-log10(1e-2),genomewideline=-log10(5e-08),xlim=c(min(assoc1$pos)-1e6,max(assoc1$pos)+1e6))

## Calculate the inflation factor lambda
lambda1<-dchisq(qchisq(median(assoc1$Score.pval),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
lambda1

lambda2<-dchisq(qchisq(median(assoc2$SPA.pval),df=1),df=1)/dchisq(qchisq(0.5,df=1),df=1)
lambda2

## On a side note: if you would like to save the data passed QC as a new GDS file, use the code below
# seqSetFilter(data,variant.id=keepSNP,sample.id=keepSample)
# seqExport(data,out.fn="NewData.gds")

## Close the GDS file
seqClose(data)

################################END##############################
