
library(pROC)                  #引用包
expFile="GSE111828.txt"        #表达数据文件
geneFile="interGenes.txt"      #疾病特征基因的文件
setwd("N:\\生信文章\\APAP肝损伤\\1.bioinformatics\\27.ROC")    #设置工作目录

#读取表达数据文件，并对文件进行整理

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)


#如果数据没有取log2,会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#对交集基因进行循环，绘制ROC曲线
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)    #定义图形的颜色
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  #绘制ROC曲线
  roc1=roc(y, as.numeric(rt[x,]))     #得到ROC曲线的参数
  if(k==1){
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

#构建逻辑模型
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #得到模型的打分

#绘制模型的ROC曲线
roc1=roc(y, as.numeric(pred))      #得到模型ROC曲线的参数
ci1=ci.auc(roc1, method="bootstrap")     #得到ROC曲线下面积的波动范围
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="#FC4CC9", legacy.axes=T, main="Model")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="#FC4CC9")
dev.off()


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

