

#ÒýÓÃ°ü
library(rms)
library(rmda)

data=read.table("merge.txt", header=T, sep="\t", check.names=F, row.names=1)

row.names(data)=gsub("-", "_", row.names(data))

#Read the gene list file and extract the expression amount of disease characteristic genes
geneRT=read.table("interGenes.txt", header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#Obtain sample grouping information
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=cbind(as.data.frame(data), Type=group)
paste(colnames(data), collapse="+")


ddist=datadist(rt)
options(datadist="ddist")

#Build a model and draw a nomogram
lrmModel=lrm(Type~ Sdhb+Pdha1+Ndufb2+Ndufb6, data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
              fun.at=c(0.00001,0.1,0.3,0.5,0.7,0.9,0.99),
              lp=F, funlabel="Risk of Disease")
#Output nomogram
pdf("Nomo.pdf", width=12, height=6)
plot(nomo)
dev.off()


#Draw calibration curve
cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)
dev.off()



#Plot decision curve
rt$Type=ifelse(rt$Type=="Control", 0, 1)
dc=decision_curve(Type ~ Sdhb+Pdha1+Ndufb2+Ndufb6, data=rt, 
                  family = binomial(link ='logit'),
                  thresholds= seq(0,1,by = 0.01),
                  confidence.intervals = 0.95)
#Output DCA graph
pdf(file="DCA.pdf", width=5.5, height=5.5)
plot_decision_curve(dc,
                    curve.names="Model",
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    col="#FC4CC9",
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

#drawing clinical impact curve 
pdf(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(dc,
                     confidence.intervals=T,
                     col = c("#FC4CC9", "#FE9601"))
dev.off()

