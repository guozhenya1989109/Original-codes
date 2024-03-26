
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

#�Խ����������ѭ��������ROC����
bioCol=rainbow(nrow(geneRT), s=0.9, v=0.9)    #����ͼ�ε���ɫ
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
  k=k+1
  #����ROC����
  roc1=roc(y, as.numeric(rt[x,]))     #�õ�ROC���ߵĲ���
  if(k==1){
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}
#����ͼ�����õ�ROC�����µ����
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:(ncol(rt)-1)])
dev.off()

#�����߼�ģ��
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #�õ�ģ�͵Ĵ��

#����ģ�͵�ROC����
roc1=roc(y, as.numeric(pred))      #�õ�ģ��ROC���ߵĲ���
ci1=ci.auc(roc1, method="bootstrap")     #�õ�ROC����������Ĳ�����Χ
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc1, print.auc=TRUE, col="#FC4CC9", legacy.axes=T, main="Model")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="#FC4CC9")
dev.off()


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056
