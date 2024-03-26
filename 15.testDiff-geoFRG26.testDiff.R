

#引用包
library(limma)
library(ggpubr)
library(reshape)

expFile="GSE111828.txt"      #表达数据文件
geneFile="interGenes.txt"     #基因列表文件



#读取表达数据文件
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



#读取基因列表文件,提取疾病特征基因的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)


data<- as.data.frame(data)

data=data[as.vector(geneRT[,1]),,drop=F]


exp=data

#提取样品的分组信息
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
conNum=length(Type[Type=="Control"])
treatNum=length(Type[Type=="Treat"])


#把表达数据转换成ggplot2输入文件

exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Type", 
            xlab="",
            ylab="Gene expression",
            legend.title="Type",
            palette = c("#FE9601", "#FC4CC9"),
            add="jitter",add.params=list(color = "Type",size=0.6),
            width=0.8)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="t.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出箱线图
pdf(file="boxplot_GSE15480qqqq.pdf", width=8, height=5)
print(p1)
dev.off()

