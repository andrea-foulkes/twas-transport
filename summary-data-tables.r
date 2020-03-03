#############################################
# The code in this file was used to generate:
# Table 1 and Supplement Table 1
# Foulkes et al, Scientific Reports (2020)
#############################################

require(Hmisc)
require(dplyr)

gene <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/GENE_LPS_pheno_data.csv")

cric <- read.csv("/Users/afoulkes/GTeX_Functions/Data/CRiC/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]

nhanes.bmx <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/DEMO_I.XPT")
nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")

pheno <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Isoform-Obesity/data/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"

# CRIC
bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
print(quantile(bmi.cric,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.cric.white.women = c(round(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.cric,breaks=c(0,25,30,100)))),3),length(bmi.cric))
length(bmi.cric)

bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
print(quantile(bmi.cric,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.cric.white.men = c(round(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))),3),length(bmi.cric))
length(bmi.cric)

bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$RACE_CAT_1==2 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
print(quantile(bmi.cric,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.cric.aa.women = c(round(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))),3),length(bmi.cric))
length(bmi.cric)

bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$RACE_CAT_1==2 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
print(quantile(bmi.cric,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.cric.aa.men = c(round(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.cric,breaks=c(0,18.5,25,30,100)))),3),length(bmi.cric))


#GENE
bmi.gene = gene[gene$race==1 & gene$sex == 2 & gene$time==0,]$bmi
print(quantile(bmi.gene,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gene.white.women = c(round(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gene))
length(bmi.gene)

bmi.gene = gene[gene$race==1 & gene$sex == 1 & gene$time==0,]$bmi
print(quantile(bmi.gene,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gene.white.men = c(round(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gene))
length(bmi.gene)

bmi.gene = gene[gene$race==2 & gene$sex == 2 & gene$time==0,]$bmi
print(quantile(bmi.gene,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gene.aa.women = c(round(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gene))
length(bmi.gene)

bmi.gene = gene[gene$race==2 & gene$sex == 1 & gene$time==0,]$bmi
print(quantile(bmi.gene,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gene.aa.men = c(round(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gene,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gene))


# GTEX
bmi.gtex.sub = pheno$BMI[pheno$SEX==2 & pheno$RACE==3 & pheno$ETHNCTY==0]
print(quantile(bmi.gtex.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gtex.white.women = c(round(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gtex.sub))
length(bmi.gtex.sub)

bmi.gtex.sub = pheno$BMI[pheno$SEX==1 & pheno$RACE==3 & pheno$ETHNCTY==0]
print(quantile(bmi.gtex.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gtex.white.men = c(round(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gtex.sub))
length(bmi.gtex.sub)

bmi.gtex.sub = pheno$BMI[pheno$SEX==2 & pheno$RACE==2 & pheno$ETHNCTY==0]
print(quantile(bmi.gtex.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gtex.aa.women = c(round(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gtex.sub))
length(bmi.gtex.sub)

bmi.gtex.sub = pheno$BMI[pheno$SEX==1 & pheno$RACE==2 & pheno$ETHNCTY==0]
print(quantile(bmi.gtex.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.gtex.aa.men = c(round(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.gtex.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.gtex.sub))


# NHANES
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
print(quantile(bmi.nhanes.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.nhanes.white.women = c(round(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.nhanes.sub))

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
print(quantile(bmi.nhanes.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.nhanes.white.men = c(round(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.nhanes.sub))

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
print(quantile(bmi.nhanes.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.nhanes.aa.women = c(round(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.nhanes.sub))

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
print(quantile(bmi.nhanes.sub,prob=c(0,0.2,0.4,0.6,0.8,1),na.rm=TRUE))
tab1.nhanes.aa.men = c(round(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))/sum(table(cut(bmi.nhanes.sub,breaks=c(0,18.5,25,30,100)))),3),length(bmi.nhanes.sub))

tab1.results = rbind(tab1.cric.white.women,tab1.cric.white.men,tab1.cric.aa.women,tab1.cric.aa.men,
    tab1.gene.white.women,tab1.gene.white.men,tab1.gene.aa.women,tab1.gene.aa.men,
    tab1.gtex.white.women,tab1.gtex.white.men,tab1.gtex.aa.women,tab1.gtex.aa.men,
    tab1.nhanes.white.women,tab1.nhanes.white.men,tab1.nhanes.aa.women,tab1.nhanes.aa.men)

write.csv(tab1.results,"/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/results/table1.results.csv")

# PHS
round(c(2789,5545,4893,2006,1099)/sum(c(2789,5545,4893,2006,1099)),3)

# WHS
round(c(6792,7421,6849,4160,7478)/sum(c(6792,7421,6849,4160,7478)),3)

# COMPARING NHANES AND GTEX

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
bmi.gtex.sub = pheno$BMI[pheno$SEX==2 & pheno$RACE==3 & pheno$ETHNCTY==0]
ks.gtex.white.women = c(ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")
w.gtex.white.women = c(wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
bmi.gtex.sub = pheno$BMI[pheno$SEX==1 & pheno$RACE==3 & pheno$ETHNCTY==0]
ks.gtex.white.men =c(ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")
w.gtex.white.men =c(wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==4 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
bmi.gtex.sub = pheno$BMI[pheno$SEX==2 & pheno$RACE==2 & pheno$ETHNCTY==0]
ks.gtex.aa.women = c(ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")
w.gtex.aa.women = c(wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")

bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
bmi.gtex.sub = pheno$BMI[pheno$SEX==1 & pheno$RACE==2 & pheno$ETHNCTY==0]
ks.gtex.aa.men = c(ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", ks.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")
w.gtex.aa.men = c(wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gtex.sub)$"p.value")

# COMPARING CRIC AND NHANES

bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
ks.cric.white.women =  c(ks.test(bmi.nhanes.sub,bmi.cric)$"statistic", ks.test(bmi.nhanes.sub,bmi.cric)$"p.value")
w.cric.white.women =  c(wilcox.test(bmi.nhanes.sub,bmi.cric)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.cric)$"p.value")

bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==3 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
ks.cric.white.men =   c(ks.test(bmi.nhanes.sub,bmi.cric)$"statistic", ks.test(bmi.nhanes.sub,bmi.cric)$"p.value")
w.cric.white.men =   c(wilcox.test(bmi.nhanes.sub,bmi.cric)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.cric)$"p.value")

bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$RACE_CAT_1==2 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
ks.cric.aa.women =   c(ks.test(bmi.nhanes.sub,bmi.cric)$"statistic", ks.test(bmi.nhanes.sub,bmi.cric)$"p.value")
w.cric.aa.women =   c(wilcox.test(bmi.nhanes.sub,bmi.cric)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.cric)$"p.value")

bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$RACE_CAT_1==2  & cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==4 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
ks.cric.aa.men =  c(ks.test(bmi.nhanes.sub,bmi.cric)$"statistic", ks.test(bmi.nhanes.sub,bmi.cric)$"p.value")
w.cric.aa.men =  c(wilcox.test(bmi.nhanes.sub,bmi.cric)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.cric)$"p.value")

# COMPARING GENE AND NHANES

bmi.gene = gene[gene$race==1 & gene$sex == 2 & gene$time==0,]$bmi
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
ks.test(bmi.gene,bmi.nhanes.sub)
ks.gene.white.women =  c(ks.test(bmi.nhanes.sub,bmi.gene)$"statistic", ks.test(bmi.nhanes.sub,bmi.gene)$"p.value")
w.gene.white.women =  c(wilcox.test(bmi.nhanes.sub,bmi.gene)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gene)$"p.value")

bmi.gene = gene[gene$race==1 & gene$sex == 1 & gene$time==0,]$bmi
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==3 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
ks.gene.white.men =  c(ks.test(bmi.nhanes.sub,bmi.gene)$"statistic", ks.test(bmi.nhanes.sub,bmi.gene)$"p.value")
w.gene.white.men =  c(wilcox.test(bmi.nhanes.sub,bmi.gene)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gene)$"p.value")

bmi.gene = gene[gene$race==2 & gene$sex == 2 & gene$time==0,]$bmi
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <=70]
ks.gene.aa.women =  c(ks.test(bmi.nhanes.sub,bmi.gene)$"statistic", ks.test(bmi.nhanes.sub,bmi.gene)$"p.value")
w.gene.aa.women =  c(wilcox.test(bmi.nhanes.sub,bmi.gene)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gene)$"p.value")

bmi.gene = gene[gene$race==2 & gene$sex == 1 & gene$time==0,]$bmi
bmi.nhanes.sub = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==4 & nhanes$ridageyr >=21& nhanes$ridageyr <=70]
ks.gene.aa.men =  c(ks.test(bmi.nhanes.sub,bmi.gene)$"statistic", ks.test(bmi.nhanes.sub,bmi.gene)$"p.value")
w.gene.aa.men =  c(wilcox.test(bmi.nhanes.sub,bmi.gene)$"statistic", wilcox.test(bmi.nhanes.sub,bmi.gene)$"p.value")


tab1.ks.results = rbind(ks.gtex.white.women,ks.gtex.white.men,ks.gtex.aa.women,ks.gtex.aa.men,
ks.cric.white.women,ks.cric.white.men,ks.cric.aa.women,ks.cric.aa.men,
ks.gene.white.women,ks.gene.white.men,ks.gene.aa.women,ks.gene.aa.men,
w.gtex.white.women,w.gtex.white.men,w.gtex.aa.women,w.gtex.aa.men,
w.cric.white.women,w.cric.white.men,w.cric.aa.women,w.cric.aa.men,
w.gene.white.women,w.gene.white.men,w.gene.aa.women,w.gene.aa.men)

write.csv(tab1.ks.results,"/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/results/table1.ks.results.csv")

