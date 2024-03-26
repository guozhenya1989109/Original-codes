



library(glmnet)     #引用包
inputFile="diffGeneExp(CRG).txt"       #输入文件


#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#构建模型
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)", "\\2", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 100)
#绘制Lasso回归的图形
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()
#绘制交叉验证的图形
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)



