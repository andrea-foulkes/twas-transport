#############################################
# The code in this file was used to generate:
# Figure 1 and Supplement Figure 1
# Foulkes et al, Scientific Reports (2020)
#############################################

### GENERATING FIGURE 1 (Whites)

# reading in gtex and nhanes data - MEN
pheno <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Isoform-Obesity/data/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"
bmi.gtex = pheno[pheno$RACE==3 & pheno$SEX==1 & pheno$ETHNCTY==0,]$BMI
gtex.rna <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/Inf5kb/Tx/transcriptome.csv")
il1b.exp = log2(gtex.rna$IL1B[gtex.rna$IL1B!=0])

nhanes.bmx <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/DEMO_I.XPT")
nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")
bmi.nhanes = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
bmi.nhanes = bmi.nhanes[!is.na(bmi.nhanes)]

gene <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/GENE_LPS_pheno_data.csv")
bmi.gene = gene[gene$race==1 & gene$sex == 1 & gene$time==0,]$bmi

cric <- read.csv("/Users/afoulkes/GTeX_Functions/Data/CRiC/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]
bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]

# generating figure of bmi

dta.bmi = data.frame(bmi=c(bmi.nhanes, bmi.gtex,bmi.cric,bmi.gene),
cohort=c(rep("NHANES",length(bmi.nhanes)),rep("GTEx",length(bmi.gtex)),rep("CRIC",length(bmi.cric)),rep("GENE",length(bmi.gene))))

dta.bmi$cohort <- factor(dta.bmi$cohort , levels = c("NHANES", "GTEx", "CRIC", "GENE"))

ggplot(data=dta.bmi,aes(x=bmi,colour=cohort,fill=cohort)) + geom_density(alpha=.1,size=.6) +
scale_fill_manual(values=c("black","steelblue2","#009E73","slateblue4"))+
scale_colour_manual(values=c("black","steelblue2","#009E73","slateblue4"))+
xlim(15,70)+ylim(0,.155)+ theme_bw()+
theme(legend.position="bottom",legend.text=element_text(size=14),legend.title = element_blank())

require(viridis)
pdf("bmi_men_white.pdf")
ggplot(dta.bmi, aes(x = `bmi`, y = `cohort`, fill = cohort)) +
geom_density_ridges(scale=5,alpha=.9,quantile_lines = TRUE, quantiles = 2) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
xlim(15,70)+theme_minimal() + theme(axis.title.y = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=14))+theme(legend.position="none")
dev.off()

# reading in gtex and nhanes data - WOMEN
pheno <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Isoform-Obesity/data/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"
bmi.gtex = pheno[pheno$RACE==3 & pheno$SEX==2 & pheno$ETHNCTY==0,]$BMI
gtex.rna <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/Inf5kb/Tx/transcriptome.csv")
il1b.exp = log2(gtex.rna$IL1B[gtex.rna$IL1B!=0])

nhanes.bmx <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/DEMO_I.XPT")
nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")
bmi.nhanes = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
bmi.nhanes = bmi.nhanes[!is.na(bmi.nhanes)]

gene <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/GENE_LPS_pheno_data.csv")
bmi.gene = gene[gene$race==1 & gene$sex == 2 & gene$time==0,]$bmi

cric <- read.csv("/Users/afoulkes/GTeX_Functions/Data/CRiC/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]
bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]

# generating figure of bmi

dta.bmi = data.frame(bmi=c(bmi.nhanes, bmi.gtex,bmi.cric,bmi.gene),
cohort=c(rep("NHANES",length(bmi.nhanes)),rep("GTEx",length(bmi.gtex)),rep("CRIC",length(bmi.cric)),rep("GENE",length(bmi.gene))))

dta.bmi$cohort <- factor(dta.bmi$cohort , levels = c("NHANES", "GTEx", "CRIC", "GENE"))

ggplot(data=dta.bmi,aes(x=bmi,colour=cohort,fill=cohort)) + geom_density(alpha=0.4) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
scale_colour_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+xlim(15,70)+ylim(0,.155)+
theme(legend.position="bottom",legend.text=element_text(size=20),legend.title = element_blank())

pdf("bmi_women_white.pdf")
ggplot(dta.bmi, aes(x = `bmi`, y = `cohort`, fill = cohort)) +
geom_density_ridges(scale=5,alpha=.9,quantile_lines = TRUE, quantiles = 2) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
xlim(15,70)+theme_minimal() + theme(axis.title.y = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=14))+theme(legend.position="none")
dev.off()



### GENERATING Supplement FIGURE 1 (AAs)

# reading in gtex and nhanes data - MEN
pheno <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Isoform-Obesity/data/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"
bmi.gtex = pheno[pheno$RACE==2 & pheno$SEX==1 & pheno$ETHNCTY==0,]$BMI
gtex.rna <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/Inf5kb/Tx/transcriptome.csv")
il1b.exp = log2(gtex.rna$IL1B[gtex.rna$IL1B!=0])

nhanes.bmx <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/DEMO_I.XPT")
nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")
bmi.nhanes = nhanes$bmxbmi[nhanes$riagendr==1 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
bmi.nhanes = bmi.nhanes[!is.na(bmi.nhanes)]

gene <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/GENE_LPS_pheno_data.csv")
bmi.gene = gene[gene$race==2 & gene$sex == 1 & gene$time==0,]$bmi

cric <- read.csv("/Users/afoulkes/GTeX_Functions/Data/CRiC/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]
bmi.cric = cric.first$BMI[cric.first$SEX==1 & cric.first$RACE_CAT_1==2 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]

# generating figure of bmi

dta.bmi = data.frame(bmi=c(bmi.nhanes, bmi.gtex,bmi.cric,bmi.gene),
cohort=c(rep("NHANES",length(bmi.nhanes)),rep("GTEx",length(bmi.gtex)),rep("CRIC",length(bmi.cric)),rep("GENE",length(bmi.gene))))

dta.bmi$cohort <- factor(dta.bmi$cohort , levels = c("NHANES", "GTEx", "CRIC", "GENE"))

pdf("bmi_men_aa.pdf")
ggplot(dta.bmi, aes(x = `bmi`, y = `cohort`, fill = cohort)) +
geom_density_ridges(scale=5,alpha=.9,quantile_lines = TRUE, quantiles = 2) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
xlim(15,70)+theme_minimal() + theme(axis.title.y = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=14))+theme(legend.position="none")
dev.off()

# reading in gtex and nhanes data - WOMEN
pheno <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Isoform-Obesity/data/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"
bmi.gtex = pheno[pheno$RACE==2 & pheno$SEX==2 & pheno$ETHNCTY==0,]$BMI
gtex.rna <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/Inf5kb/Tx/transcriptome.csv")
il1b.exp = log2(gtex.rna$IL1B[gtex.rna$IL1B!=0])

nhanes.bmx <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Transport/NHANES data/DEMO_I.XPT")
nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")
bmi.nhanes = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==4 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
bmi.nhanes = bmi.nhanes[!is.na(bmi.nhanes)]

gene <- read.csv("/Volumes/GoogleDrive/My Drive/Projects/TWAS-Longitudinal/data/GENE_LPS_pheno_data.csv")
bmi.gene = gene[gene$race==2 & gene$sex == 2 & gene$time==0,]$bmi

cric <- read.csv("/Users/afoulkes/GTeX_Functions/Data/CRiC/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]
bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$RACE_CAT_1==2 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]

# generating figure of bmi

dta.bmi = data.frame(bmi=c(bmi.nhanes, bmi.gtex,bmi.cric,bmi.gene),
cohort=c(rep("NHANES",length(bmi.nhanes)),rep("GTEx",length(bmi.gtex)),rep("CRIC",length(bmi.cric)),rep("GENE",length(bmi.gene))))

dta.bmi$cohort <- factor(dta.bmi$cohort , levels = c("NHANES", "GTEx", "CRIC", "GENE"))

pdf("bmi_women_aa.pdf")
ggplot(dta.bmi, aes(x = `bmi`, y = `cohort`, fill = cohort)) +
geom_density_ridges(scale=5,alpha=.9,quantile_lines = TRUE, quantiles = 2) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
xlim(15,70)+theme_minimal() + theme(axis.title.y = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=14))+theme(legend.position="none")
dev.off()

