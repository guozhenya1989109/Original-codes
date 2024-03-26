#install.packages("VennDiagram")


library(VennDiagram)      #引用包


#读取差异分析的结果文件
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取差异基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对差异基因取unique
geneList[["DEG"]]=uniqGene

#读取自噬基因的列表文件
rt=read.table("gene_37个基因转换成小写了.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取自噬基因的名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #对自噬基因取unique
geneList[["Cuproptosis"]]=uniqGene

#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#FC4CC9", "#FE9601"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集基因的列表
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGene.txt", sep="\t", quote=F, col.names=F, row.names=F)


