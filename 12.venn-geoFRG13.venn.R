######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

#install.packages("VennDiagram")


library(VennDiagram)      #引用包
lassoFile="LASSO.gene.txt"      #lasso回归的基因列表文件
RFFile="rfGenes.txt"      #机器学习的基因列表文件
svmFile="SVM-RFE.gene.txt"      #机器学习的基因列表文件


setwd("N:\\生信文章\\APAP肝损伤\\1.bioinformatics\\13.venn")    #设置工作目录
geneList=list()

#读取lasso回归的基因列表文件
rt=read.table(lassoFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["LASSO"]]=uniqGene             #把lasso回归找到的特征基因放到geneList里面

#读取机器学习的基因列表文件
rt=read.table(RFFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])              #提取基因名称
geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
uniqGene=unique(geneNames)               #基因取unique
geneList[["RF"]]=uniqGene               #把SVM方法找到的特征基因放到geneList里面


# 
# # #读取机器学习的基因列表文件
# rt=read.table(svmFile, header=F, sep="\t", check.names=F)
# geneNames=as.vector(rt[,1])              #提取基因名称
# geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
# uniqGene=unique(geneNames)               #基因取unique
# geneList[["SVM"]]=uniqGene               #把SVM方法找到的特征基因放到geneList里面



#绘制venn图
venn.plot=venn.diagram(geneList,filename=NULL,fill=c("#FC4CC9", "#FE9601" ),
                       scaled=FALSE,cat.col = c("#FC4CC9", "#FE9601"),cat.cex = 1.2)
pdf(file="venn.pdf", width=5, height=5)
grid.draw(venn.plot)
dev.off()

#输出交集基因的列表
interGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)


######Video source: https://ke.biowolf.cn
######生信自学网: https://www.biowolf.cn/
######微信公众号：biowolf_cn
######合作邮箱：biowolf@foxmail.com
######答疑微信: 18520221056

