#install.packages("VennDiagram")


library(VennDiagram)      #���ð�


#��ȡ��������Ľ���ļ�
rt=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #��ȡ������������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #�Բ������ȡunique
geneList[["DEG"]]=uniqGene

#��ȡ���ɻ�����б��ļ�
rt=read.table("gene_37������ת����Сд��.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #��ȡ���ɻ��������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #�����ɻ���ȡunique
geneList[["Cuproptosis"]]=uniqGene

#����vennͼ
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#FC4CC9", "#FE9601"),scaled=FALSE,cat.pos=c(-1,1),cat.col = c("cornflowerblue", "darkorchid1"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#�������������б�
interGenes=Reduce(intersect, geneList)
write.table(interGenes, file="interGene.txt", sep="\t", quote=F, col.names=F, row.names=F)

