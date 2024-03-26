######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("pROC")


library(pROC)                  #���ð�
expFile="GSE111828.txt"        #���������ļ�
geneFile="interGenes.txt"      #��������������ļ�
setwd("N:\\��������\\APAP������\\1.bioinformatics\\27.ROC")    #���ù���Ŀ¼

#��ȡ���������ļ��������ļ���������

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)


#�������û��ȡlog2,��������Զ�ȡlog2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#��ȡ�����б��ļ�
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#�Լ��������������ѭ��������ROC����
for(x in as.vector(geneRT[,1])){
  #����ROC����
  roc1=roc(y, as.numeric(rt[x,]))
  ci1=ci.auc(roc1, method="bootstrap")
  ciVec=as.numeric(ci1)
  #���ROC����
  pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056
