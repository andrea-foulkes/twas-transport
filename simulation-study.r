#############################################
# The code in this file was used to generate:
# Figures 3a, 3b and 4; Table 3
# Foulkes et al, Scientific Reports (2020)
#############################################

require(ggplot2)
require(dplyr)
require(Hmisc)
require(jtools)
require(ggridges)
require(gtools)

# reading in gtex and nhanes data - WOMEN
pheno <- read.csv("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/GTEx_subject_phenotype.csv")
colnames(pheno)[2] <- "ID"
bmi.gtex = pheno[pheno$RACE==3 & pheno$SEX==2 & pheno$ETHNCTY==0,]$BMI
gtex.rna <- read.csv("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/transcriptome.csv")
il1b.exp = log2(gtex.rna$IL1B[gtex.rna$IL1B!=0])

nhanes.bmx <- sasxport.get("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/BMX_I.XPT")
nhanes.demog <- sasxport.get("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/DEMO_I.XPT")
#nhanes.diab <- sasxport.get("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/DIQ_I.XPT")

nhanes <- left_join(nhanes.bmx,nhanes.demog, by="seqn")
bmi.nhanes = nhanes$bmxbmi[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70 ]
bmi.nhanes = bmi.nhanes[!is.na(bmi.nhanes)]


gene <- read.csv("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/GENE_LPS_pheno_data.csv")
bmi.gene = gene[gene$race==1 & gene$sex == 2 & gene$time==0,]$bmi

cric <- read.csv("/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/dta/clinical.csv")
cric.first = cric[!duplicated(cric$PID),]
bmi.cric = cric.first$BMI[cric.first$SEX==2 & cric.first$WHITE==1 &  cric.first$HISPANIC==0 & cric.first$AGE_INTEGER <=70]
bmi.cric = bmi.cric[!is.na(bmi.cric)]

# fit model to estimate weights for gtext data
q.nhanes.age = quantile(nhanes[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70,]$ridageyr,probs=c(0,0.2,0.4,0.6,0.8,1))
mod.weights = lm(log(bmxbmi)~cut(ridageyr,q.nhanes.age)+bmxwt,data=nhanes[nhanes$riagendr==2 & nhanes$ridreth1==3 & nhanes$ridageyr >=21 & nhanes$ridageyr <= 70,])

mod.weights.predict = predict(mod.weights,data.frame(ridageyr=pheno[pheno$RACE==3 & pheno$SEX==2 & pheno$ETHNCTY==0 ,]$AGE, bmxwt = pheno[pheno$RACE==3 & pheno$SEX==2 & pheno$ETHNCTY==0 ,]$WGHT*0.453592),se.fit=TRUE)

n.gtex = length(pheno[pheno$RACE==3 & pheno$SEX==2 & pheno$ETHNCTY==0,]$AGE)
gtex.bmi.predict = exp(mod.weights.predict$fit+rnorm(n.gtex,0,sqrt(var(residuals(mod.weights)))))

# Simulating population data
n = 1000000
allele1 <- sample(c("A","a"),size=n, prob=c(0.8,0.2),replace=TRUE)
allele2 <- sample(c("A","a"),size=n, prob=c(0.8,0.2),replace=TRUE)
z <- as.numeric(allele1=="a") + as.numeric(allele2=="a")
    
allele1 <- sample(c("A","a"),size=n, prob=c(0.8,0.2),replace=TRUE)
allele2 <- sample(c("A","a"),size=n, prob=c(0.8,0.2),replace=TRUE)
zb <- as.numeric(allele1=="a") + as.numeric(allele2=="a")

alpha = 0.2
alphab = 0.2
sigma2 = sqrt(var(log(bmi.nhanes)))
            
alpha0 = mean(il1b.exp) - 0.4*(alpha+alphab)
sigma1 = sqrt(var(il1b.exp))
x =  alpha0 + z*alpha + zb*alphab +  rnorm(n,0,sigma1)
alpha.pop = summary(lm(x~z+zb))$coefficients[2,1]

gamma = 0.15
gamma0 = mean(log(bmi.nhanes)) - mean(x*gamma)
bmi.sim <- gamma0 + x*gamma + rnorm(n,0,sqrt(sigma2^2-gamma^2*var(x)))
gamma.pop = summary(lm(bmi.sim~x))$coefficients[2,1]

# determing sample weight for gtex
q.bmi = exp(quantile(bmi.sim,probs = seq(0, 1, 0.1)))
w.length = length(q.bmi)-1

weights <- rep(0,w.length)
bias.prob = rep(0,n)
for (j in 1:w.length){
    weights[j] = mean((q.bmi[j] <= bmi.gtex) & (bmi.gtex < q.bmi[j+1]),na.rm=TRUE)
    bias.prob[(q.bmi[j] <= exp(bmi.sim)) & (exp(bmi.sim) <= q.bmi[j+1])] = weights[j]
}

weights.pred <- rep(0,w.length)
bias.prob.pred = rep(0,n)
for (j in 1:w.length){
    weights.pred[j] = mean((q.bmi[j] <= gtex.bmi.predict) & (gtex.bmi.predict < q.bmi[j+1]),na.rm=TRUE)
    bias.prob.pred[(q.bmi[j] <= exp(bmi.sim)) & (exp(bmi.sim) <= q.bmi[j+1])] = weights.pred[j]
}

weights.cric <- rep(0,w.length)
bias.prob.cric = rep(0,n)
for (j in 1:w.length){
    weights.cric[j] = mean((q.bmi[j] <= bmi.cric) & (bmi.cric < q.bmi[j+1]),na.rm=TRUE)
    bias.prob.cric[(q.bmi[j] <= exp(bmi.sim)) & (exp(bmi.sim) <= q.bmi[j+1])] = weights.cric[j]
}

weights.gene <- rep(0,w.length)
bias.prob.gene = rep(0,n)
for (j in 1:w.length){
    weights.gene[j] = mean((q.bmi[j] <= bmi.gene) & (bmi.gene < q.bmi[j+1]),na.rm=TRUE)
    bias.prob.gene[(q.bmi[j] <= exp(bmi.sim)) & (exp(bmi.sim) <= q.bmi[j+1])] = weights.gene[j]
}


# Simulation with sampling from population data
n1=750
n2=1500

nsims = 2000
result = matrix(nrow=nsims,ncol=2)
predict.bias = matrix(nrow=n2,ncol=nsims)
predict.sim = matrix(nrow=n2,ncol=nsims)
predict.ipw = matrix(nrow=n2,ncol=nsims)
predict.ipw.est = matrix(nrow=n2,ncol=nsims)
predict.cric = matrix(nrow=n2,ncol=nsims)
predict.gene = matrix(nrow=n2,ncol=nsims)

result.test = matrix(nrow=nsims,ncol=5)
est.new1 = matrix(nrow=nsims,ncol=3)
est.new2 = matrix(nrow=nsims,ncol=3)
est.new2.ipw = matrix(nrow=nsims,ncol=3)
est.new2.ipw.est = matrix(nrow=nsims,ncol=3)
est.new.cric = matrix(nrow=nsims,ncol=3)
est.new.gene = matrix(nrow=nsims,ncol=3)

x.sim = matrix(nrow=n2,ncol=nsims)
y.sim = matrix(nrow=n2,ncol=nsims)

est.unbiased = matrix(nrow=nsims,ncol=3)
est.biased = matrix(nrow=nsims,ncol=3)
est.biased.ipw = matrix(nrow=nsims,ncol=3)
est.biased.ipw.est = matrix(nrow=nsims,ncol=3)
est.biased.cric = matrix(nrow=nsims,ncol=3)
est.biased.gene = matrix(nrow=nsims,ncol=3)

for (i in 1:nsims){
    print(i)
    gwas.samp = base::sample(seq(1:n),size=n2,replace=TRUE)
    
    # Unbiased sampling
    samp1 = base::sample(seq(1:n),size=n1,replace=TRUE)
    sum.lm = summary(lm(x[samp1]~z[samp1]+zb[samp1]))$coefficients
    coef.est = summary(lm(x[samp1]~z[samp1]+zb[samp1]))$coefficients[,1]
    est.unbiased[i,] = coef.est
    
    pred.exp.sim = cbind(1,z[gwas.samp],zb[gwas.samp])%*%coef.est
    sum.lm = summary(lm(bmi.sim[gwas.samp]~pred.exp.sim))$coefficients
    result.test[i,1] = sum.lm[2,3]
    est.new1[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])

    # Biased sampling
    samp = base::sample(seq(1:n),size=n1,prob=bias.prob,replace=TRUE)
    sum.lm = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients
    coef.est = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients[,1]
    est.biased[i,] =  coef.est
    
    pred.exp.bias = cbind(1,z[gwas.samp],zb[gwas.samp])%*%coef.est
    sum.lm = summary(lm(bmi.sim[gwas.samp]~pred.exp.bias))$coefficients
    result.test[i,2] = sum.lm[2,3]
    est.new2[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])

    # IPW
    w = 1/bias.prob[samp]
    fit = lm(x[samp]~z[samp]+zb[samp],weights=w)
    sum.lm = jtools::summ(fit,robust=TRUE)$coeftable # robut standard error estimates
    coef.est = summary(lm(x[samp]~z[samp]+zb[samp],weights=w))$coefficients[,1]
    est.biased.ipw[i,] = coef.est

    pred.exp.ipw = cbind(1,z[gwas.samp],zb[gwas.samp])%*%coef.est
    sum.lm = summary(lm(bmi.sim[gwas.samp]~pred.exp.ipw))$coefficients
    result.test[i,3] = sum.lm[2,3]
    est.new2.ipw[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])

# IPW est
w = 1/bias.prob.pred[samp]
fit = lm(x[samp]~z[samp]+zb[samp],weights=w)
sum.lm = jtools::summ(fit,robust=TRUE)$coeftable # robut standard error estimates
coef.est = summary(lm(x[samp]~z[samp]+zb[samp],weights=w))$coefficients[,1]
est.biased.ipw.est[i,] = coef.est

pred.exp.ipw.est = cbind(1,z[gwas.samp],zb[gwas.samp])%*%coef.est
sum.lm = summary(lm(bmi.sim[gwas.samp]~pred.exp.ipw.est))$coefficients
result.test[i,3] = sum.lm[2,3]
est.new2.ipw.est[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])


    # Doubly biased sampling with CRIC
    gwas.samp.cric = base::sample(seq(1:n),size=n2,prob=bias.prob.cric,replace=TRUE)
    samp = base::sample(seq(1:n),size=n1,prob=bias.prob,replace=TRUE)
    sum.lm = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients
    coef.est = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients[,1]
    est.biased.cric[i,] = coef.est
    
    pred.exp.bias.cric = cbind(1,z[gwas.samp.cric],zb[gwas.samp.cric])%*%coef.est
    sum.lm = summary(lm(bmi.sim[gwas.samp.cric]~pred.exp.bias.cric))$coefficients
    result.test[i,4] = sum.lm[2,3]
    est.new.cric[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])


# Doubly biased sampling with GENE
gwas.samp.gene = base::sample(seq(1:n),size=n2,prob=bias.prob.gene,replace=TRUE)
samp = base::sample(seq(1:n),size=n1,prob=bias.prob,replace=TRUE)
sum.lm = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients
coef.est = summary(lm(x[samp]~z[samp]+zb[samp]))$coefficients[,1]
est.biased.gene[i,] = coef.est

pred.exp.bias.gene = cbind(1,z[gwas.samp.gene],zb[gwas.samp.gene])%*%coef.est
sum.lm = summary(lm(bmi.sim[gwas.samp.gene]~pred.exp.bias.gene))$coefficients
result.test[i,5] = sum.lm[2,3]
est.new.gene[i,] = c(sum.lm[2,1],sum.lm[2,1]-1.96*sum.lm[2,2], sum.lm[2,1]+1.96*sum.lm[2,2])

    predict.bias[,i] = pred.exp.bias
    predict.sim[,i] = pred.exp.sim
    predict.ipw[,i] = pred.exp.ipw
    predict.cric[,i] = pred.exp.bias.cric
    predict.gene[,i] = pred.exp.bias.gene

    x.sim[,i] = x[gwas.samp]
    y.sim[,i] = bmi.sim[gwas.samp]
}

pred.df = data.frame(predict= (2^c(apply(predict.bias,2,mean),apply(predict.cric,2,mean),apply(predict.gene,2,mean),apply(predict.sim,2,mean))),type=c(rep("Scenario 2 (GTEx-RS)",nsims),rep("Scenario 3 (GTEx-CRIC)",nsims),rep("Scenario 4 (GTEx-GENE)",nsims),rep("Scenario 1 (RS)",nsims)))
pred.df$type <- factor(pred.df$type , levels = c("Scenario 1 (RS)","Scenario 2 (GTEx-RS)","Scenario 3 (GTEx-CRIC)","Scenario 4 (GTEx-GENE)"))

saveRDS(pred.df,"predict.rds")
saveRDS(predict.ipw,"predict-ipw.rds")

# Fig 3b
pdf("prediction.pdf")
ggplot(pred.df, aes(x = `predict`, y = `type`, fill = type)) +
geom_density_ridges(scale=5,alpha=.9,quantile_lines = TRUE, quantiles = 2) +
scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
xlab("Average predicted normalized expression (RPKM)") + theme_bw() + theme(axis.title.y = element_blank(),axis.text=element_text(size=14),axis.title=element_text(size=14))+geom_vline(xintercept=2^mean(apply(x.sim,2,mean)),linetype="dotted")+ annotate(geom="text", x=1475, y=5.75, label="True mean expression",color="black",size=3)+theme(legend.position="none")
dev.off()

tapply(pred.df$predict,pred.df$type,mean)
sqrt(tapply(pred.df$predict,pred.df$type,var))
2^mean(apply(x.sim,2,mean))

dta = data.frame(estimate = c(est.new2[,1],est.new1[,1],est.new.cric[,1],est.new.gene[,1],est.new2.ipw.est[,1],est.new2.ipw[,1])-gamma,type=c(rep("2: GTEx-RS",nsims),rep("1: RS",nsims),rep("3: GTEx-CRIC",nsims),rep("4: GTEx-GENE",nsims),rep("IPW (est)",nsims),rep("IPW (known)",nsims)),scenario=c(rep("Unadjusted",nsims),rep("Unadjusted",nsims),rep("Unadjusted",nsims),rep("Unadjusted",nsims),rep("IPW (Scenario 2)",nsims),rep("IPW (Scenario 2)",nsims)))
dta$type <- factor(dta$type , levels = c("1: RS", "2: GTEx-RS", "3: GTEx-CRIC", "4: GTEx-GENE","IPW (known)","IPW (est)"))
dta$scenario <- factor(dta$scenario , levels = c("Unadjusted", "IPW (Scenario 2)"))

saveRDS(dta,"dta.rds")

# Fig 4
pdf("Expr_Trait_association_redo.pdf")
p <- ggplot(data=dta,aes(x=type,y=estimate/gamma, color=type,fill=type)) + geom_boxplot(alpha=.9)+ scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4","deeppink4","gold"))+
scale_colour_manual(values=c("#999999","steelblue2","#009E73","slateblue4","deeppink4","gold")) + ylab("Relative bias in expression-trait association estimate") + theme_bw()+theme(axis.text=element_text(size=12),
axis.title=element_text(size=12)) + #ylim(-450,450) +
xlab("Sampling scenario") +
theme(legend.position="none",axis.text.x = element_text(size=8)) + scale_y_continuous(breaks=c(-3,-2,-1,0, 1,2, 3),limits=c(-3.5,3.5),labels = scales::percent) +
geom_hline(yintercept=0, linetype="dotted")+stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black",
fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))})+facet_wrap(~scenario,scales = "free_x")
gp <- ggplotGrob(p)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
x.var <- sapply(ggplot_build(p)$layout$panel_scales_x,
function(l) length(l$range$range))
gp$widths[facet.columns] <- gp$widths[facet.columns] * x.var
grid::grid.draw(gp)

dev.off()

dta1 = data.frame(estimate = c(c(est.unbiased[,1]-alpha0,est.biased[,1]-alpha0,est.biased.cric[,1]-alpha0,est.biased.gene[,1]-alpha0)/alpha0,
    c(est.unbiased[,2]-alpha,est.biased[,2]-alpha,est.biased.cric[,2]-alpha,est.biased.gene[,2]-alpha)/alpha,
c(est.unbiased[,3]-alphab,est.biased[,3]-alphab,est.biased.cric[,3]-alphab,est.biased.gene[,3]-alphab)/alphab),type=rep(c(rep("1:RS",nsims),rep("2:GTEx-RS",nsims),rep("3:GTEx-CRIC",nsims),rep("4:GTEx-GENE",nsims)),3),parameter=c(rep("intercept",4*nsims),rep("SNP1",4*nsims),rep("SNP2",4*nsims)))

dta1$type <- factor(dta1$type , levels = c("1:RS", "2:GTEx-RS", "3:GTEx-CRIC", "4:GTEx-GENE"))

saveRDS(dta1,"dta1.rds")

# Fig 3a
pdf("Gene-Expression_association.pdf")
ggplot(data=dta1,aes(x= type,y=estimate, color=type,parameter,fill=type)) + geom_boxplot(alpha=.9)+ scale_fill_manual(values=c("#999999","steelblue2","#009E73","slateblue4"))+
scale_colour_manual(values=c("#999999","steelblue2","#009E73","slateblue4")) +
ylab("Relative bias in coefficient estimates for genotype-expression model") +
theme(legend.position="bottom",axis.title.x=element_blank(),
axis.text.x=element_blank()) +
scale_y_continuous(labels = scales::percent) +
geom_hline(yintercept=0, linetype="dotted")+stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black",
fun.data = function(x){c(y=median(x), ymin=median(x), ymax=median(x))}) +facet_wrap(~parameter, scales = "free")
dev.off()

sum.results = rbind(c(round(mean(est.new1[,1]),3), 100*round(mean(est.new1[,1]-gamma)/gamma,3),100*round(quantile(est.new1[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new1[,2] >gamma  | est.new1[,3] < gamma)), round(mean(est.new1[,3] - est.new1[,2]),3)),
c(round(mean(est.new2[,1]),3), 100*round(mean(est.new2[,1]-gamma)/gamma,3),100*round(quantile(est.new2[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new2[,2] >gamma  | est.new2[,3] < gamma)),round(mean(est.new2[,3] - est.new2[,2]),3)),
c(round(mean(est.new.cric[,1]),3), 100*round(mean(est.new.cric[,1]-gamma)/gamma,3),100*round(quantile(est.new.cric[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new.cric[,2] >gamma  | est.new.cric[,3] < gamma)), round(mean(est.new.cric[,3] - est.new.cric[,2]),3)),
c(round(mean(est.new.gene[,1]),3), 100*round(mean(est.new.gene[,1]-gamma)/gamma,3),100*round(quantile(est.new.gene[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new.gene[,2] >gamma  | est.new.gene[,3] < gamma)), round(mean(est.new.gene[,3] - est.new.gene[,2]),3)),
c(round(mean(est.new2.ipw[,1]),3), 100*round(mean(est.new2.ipw[,1]-gamma)/gamma,3),100*round(quantile(est.new2.ipw[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new2.ipw[,2] >gamma  | est.new2.ipw[,3] < gamma)), round(mean(est.new2.ipw[,3] - est.new2.ipw[,2]),3)),
c(round(mean(est.new2.ipw.est[,1]),3), 100*round(mean(est.new2.ipw.est[,1]-gamma)/gamma,3),100*round(quantile(est.new2.ipw.est[,1]-gamma,c(0.5, 0.25,0.75))/gamma,3), 1-mean( (est.new2.ipw.est[,2] >gamma  | est.new2.ipw.est[,3] < gamma)), round(mean(est.new2.ipw.est[,3] - est.new2.ipw.est[,2]),3)))


colnames(sum.results) = c("Estimate","Relative Bias","Q50","Q25","Q75","Coverage","CI length")

row.names(sum.results) = c("RS","2: GTEx-RS", "3: GTEx-CRIC", "4: GTEx-GENE", "GTEx-RS + IPW (known)", "GTEx-RS + IPW (est)")

# Table 3
write.csv(sum.results,"/Users/afoulkes/OneDrive - Partners HealthCare/Projects/TWAS-Transport/results/sum.results.csv")

