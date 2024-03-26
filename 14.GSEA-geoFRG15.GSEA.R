

#引用包
library(limma)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)





rm(list = ls())

expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="gene_go.backgroud.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}




#========================
expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="m5.go.bp.v2023.2.Mm.symbols.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_BP.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
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
expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="m5.go.cc.v2023.2.Mm.symbols.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_CC.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
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
expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="m5.go.mf.v2023.2.Mm.symbols.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_GO_mf.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_GO_mf.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}





rm(list = ls())
######Video source: https://ke.biowolf.cn
expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="mh.hallmark .v2023.2.Mm.symbols.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_hallmark.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
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
expFile="merge.txt"     #表达数据文件
gene="Ndufb6"                 #基因名称
gmtFile="gene_kegg.backgroud.gmt"     #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#根据目标基因表达量对样品进行分组，得到高低表达组的logFC
dataL=data[,data[gene,]<median(data[gene,]),drop=F]     #低表达组的数据
dataH=data[,data[gene,]>=median(data[gene,]),drop=F]    #高表达组的数据
meanL=rowMeans(dataL)
meanH=rowMeans(dataH)
meanL[meanL<0.00001]=0.00001
meanH[meanH<0.00001]=0.00001
logFC=log2(meanH)-log2(meanL)
#根据logFC对基因进行排序
logFC=sort(logFC, decreasing=T)
genes=names(logFC)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$pvalue<0.05,]
#kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_Ndufb6_KEGG.txt",sep="\t",quote=F,row.names = F)

#绘制GSEA富集的图形
termNum=6     #展示通路的数目
if(nrow(kkTab)>=termNum){
  showTerm=row.names(kkTab)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title=gene)
  pdf(file="GSEA_Ndufb6_KEGG.pdf", width=7.5, height=5.5)
  print(gseaplot)
  dev.off()
}