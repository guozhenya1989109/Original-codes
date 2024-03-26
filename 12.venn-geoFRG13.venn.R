######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056

#install.packages("VennDiagram")


library(VennDiagram)      #���ð�
lassoFile="LASSO.gene.txt"      #lasso�ع�Ļ����б��ļ�
RFFile="rfGenes.txt"      #����ѧϰ�Ļ����б��ļ�
svmFile="SVM-RFE.gene.txt"      #����ѧϰ�Ļ����б��ļ�


setwd("N:\\��������\\APAP������\\1.bioinformatics\\13.venn")    #���ù���Ŀ¼
geneList=list()

#��ȡlasso�ع�Ļ����б��ļ�
rt=read.table(lassoFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #��ȡ��������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #����ȡunique
geneList[["LASSO"]]=uniqGene             #��lasso�ع��ҵ�����������ŵ�geneList����

#��ȡ����ѧϰ�Ļ����б��ļ�
rt=read.table(RFFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #��ȡ��������
geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
uniqGene=unique(geneNames)               #����ȡunique
geneList[["RF"]]=uniqGene               #��SVM�����ҵ�����������ŵ�geneList����


# 
# # #��ȡ����ѧϰ�Ļ����б��ļ�
# rt=read.table(svmFile, header=F, sep="\t", check.names=F)
# geneNames=as.vector(rt[,1])              #��ȡ��������
# geneNames=gsub("^ | $","",geneNames)     #ȥ��������β�Ŀո�
# uniqGene=unique(geneNames)               #����ȡunique
# geneList[["SVM"]]=uniqGene               #��SVM�����ҵ�����������ŵ�geneList����



#����vennͼ
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#FC4CC9", "#FE9601" ),
                       scaled=FALSE,cat.col = c("#FC4CC9", "#FE9601"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#�������������б�
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)


######Video source: https://ke.biowolf.cn
######������ѧ��: https://www.biowolf.cn/
######΢�Ź��ںţ�biowolf_cn
######�������䣺biowolf@foxmail.com
######����΢��: 18520221056
