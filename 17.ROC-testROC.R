######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("pROC")


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

#对疾病特征基因进行循环，绘制ROC曲线
for(x in as.vector(geneRT[,1])){
  #绘制ROC曲线
  roc1=roc(y, as.numeric(rt[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  #输出ROC曲线
  pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

