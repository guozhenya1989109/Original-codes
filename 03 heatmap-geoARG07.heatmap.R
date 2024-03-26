


#引用包
library(limma)
library(pheatmap)


#读取表达数据文件,并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因列表文件, 提取交集的表达量
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

write.table(data, file="diffGeneExp(CRG).txt", sep="\t", quote=F, col.names=T)

#提取样品的分组信息(对照组和实验组)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
names(Type)=colnames(data)
Type=as.data.frame(Type)
Project=gsub("(.+)\\_(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))     #提取GEO研究的id
Type=cbind(Project, Type)
#绘制交集基因的热图
pdf(file="heatmap.pdf", width=8, height=5)
pheatmap(data, 
         annotation=Type, 
         color = colorRampPalette(c(rep("#FE9601",3), "white", rep("#FC4CC9",3)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 7,
         fontsize_row=5,
         fontsize_col=7)
dev.off()


