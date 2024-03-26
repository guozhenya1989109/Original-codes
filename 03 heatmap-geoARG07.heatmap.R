


#���ð�
library(limma)
library(pheatmap)


#��ȡ���������ļ�,���������ļ���������
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#��ȡ�����б��ļ�, ��ȡ�����ı�����
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

write.table(data, file="diffGeneExp(CRG).txt", sep="\t", quote=F, col.names=T)

#��ȡ��Ʒ�ķ�����Ϣ(�������ʵ����)
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
names(Type)=colnames(data)
Type=as.data.frame(Type)
Project=gsub("(.+)\\_(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))     #��ȡGEO�о���id
Type=cbind(Project, Type)
#���ƽ����������ͼ
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

