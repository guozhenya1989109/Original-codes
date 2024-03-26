


#引用包
library(limma)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="diff.txt"      #输入文件
gmtFile="gene_kegg.backgroud.gmt"      #基因集文件

#读取文件,并对输入文件进行整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)


rt1=read.table("interGene.txt", header=F, sep="\t", check.names=F)
colnames(rt1)[1]<-"id"

rt2<-merge(rt1,rt,by="id")



logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)

#读取基因集文件
gmt=read.gmt(gmtFile)

#GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result_KEGG.txt",sep="\t",quote=F,row.names = F)

#输出实验组富集的图形
termNum=6      #设置展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")
	pdf(file="GSEA.treat_KEGG.pdf", width=6.5, height=5)
	print(gseaplot)
	dev.off()
}

#输出对照组富集的图形
termNum=5      #设置展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA.con_KEGG.pdf", width=6.5, height=5)
	print(gseaplot)
	dev.off()
}




