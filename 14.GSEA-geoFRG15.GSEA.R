

#���ð�
library(limma)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)





rm(list = ls())

expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="gene_go.backgroud.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}




#========================
expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="m5.go.bp.v2023.2.Mm.symbols.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_BP.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO_BP.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}








rm(list = ls())
##==================================================
######Video source: https://ke.biowolf.cn
expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="m5.go.cc.v2023.2.Mm.symbols.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_CC.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO_CC.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





rm(list = ls())
##==================================================
######Video source: https://ke.biowolf.cn
expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="m5.go.mf.v2023.2.Mm.symbols.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_mf.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO_mf.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





rm(list = ls())
######Video source: https://ke.biowolf.cn
expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="mh.hallmark .v2023.2.Mm.symbols.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_hallmark.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_hallmark.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





#=====================================
rm(list = ls())
######Video source: https://ke.biowolf.cn
expFile="merge.txt"     #���������ļ�
gene="Ndufb6"                 #��������
gmtFile="gene_kegg.backgroud.gmt"     #�����ļ�

#��ȡ�ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#ȥ�����������Ʒ
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#����Ŀ��������������Ʒ���з��飬�õ��ߵͱ������logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #�ͱ����������
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #�߱����������
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#����logFC�Ի����������
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#��ȡ�����ļ�
gmt=read.gmt(gmtFile)

#GSEA��������
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_KEGG.txt",sep="\t",quote=F,row.names = F)

#����GSEA������ͼ��
termNum=6     #չʾͨ·����Ŀ
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_KEGG.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}