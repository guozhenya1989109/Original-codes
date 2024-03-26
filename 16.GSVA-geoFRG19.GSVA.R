


#引用包
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

gene="Ndufb6"      #基因名称
expFile="merge.txt"              #表达数据文件
gmtFile="gene_kegg.backgroud.gmt"     #基因集文件
 

#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]

#根据目标基因的表达量对样品进行分组
lowName=colnames(data)[data[gene,]<median(data[gene,])]       #低表达组的样品
highName=colnames(data)[data[gene,]>=median(data[gene,])]     #高表达组的样品
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#基因差异分析
outTab=data.frame()
for(i in row.names(data)){
	test=t.test(data[i,] ~ Type)
	pvalue=test$p.value
	t=test$statistic
	Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
	outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#绘制柱状图
pdf(file="barplot_Ndufb6_KEGG.pdf", width=20, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
		palette=c("#fc8d62","#8da0cb","#66c2a5"), sort.val = "asc", sort.by.groups = T,
		rotate=TRUE, legend="right", title=gene,
		xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()


#=======================================


gene="Ndufb6"      #基因名称
expFile="merge.txt"              #表达数据文件
gmtFile="gene_go.backgroud.gmt"     #基因集文件
#读取表达输入文件,并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#去除对照组的样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#GSVA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
ssgseaScore=ssgseaScore[order(apply(ssgseaScore,1,sd),decreasing=T),]
ssgseaScore=ssgseaScore[1:50,]

#根据目标基因的表达量对样品进行分组
lowName=colnames(data)[data[gene,]<median(data[gene,])]       #低表达组的样品
highName=colnames(data)[data[gene,]>=median(data[gene,])]     #高表达组的样品
lowScore=ssgseaScore[,lowName]
highScore=ssgseaScore[,highName]
data=cbind(lowScore, highScore)
conNum=ncol(lowScore)
treatNum=ncol(highScore)
Type=c(rep("Control",conNum), rep("Treat",treatNum))

#基因差异分析
outTab=data.frame()
for(i in row.names(data)){
  test=t.test(data[i,] ~ Type)
  pvalue=test$p.value
  t=test$statistic
  Sig=ifelse(pvalue>0.05, "Not", ifelse(t>0,"Up","Down"))
  outTab=rbind(outTab, cbind(Pathway=i, t, pvalue, Sig))
}

#绘制柱状图
pdf(file="barplot_Ndufb6_GO.pdf", width=20, height=9)
outTab$t=as.numeric(outTab$t)
outTab$Sig=factor(outTab$Sig, levels=c("Down", "Not", "Up"))
gg1=ggbarplot(outTab, x="Pathway", y="t", fill = "Sig", color = "white",
              palette=c("#8da0cb","#66c2a5","#fc8d62"), sort.val = "asc", sort.by.groups = T,
              rotate=TRUE, legend="right", title=gene,
              xlab="Term", ylab="t value of GSVA score",  legend.title="Group", x.text.angle=60)
print(gg1)
dev.off()






