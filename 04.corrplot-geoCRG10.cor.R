

#引用包
library(corrplot)
library(circlize)

inputFile="diffGeneExp(CRG).txt"    #输入文件


#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

#去除对照组样品
group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=data[,group=="Treat",drop=F]
rt=t(data)

#计算基因间相关系数
cor1=cor(rt)



#设置图形颜色
col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,1,0,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,1,0,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

#绘制圈图
pdf(file="circos.pdf", width=7, height=7)
par(mar=c(2,2,2,4))
circos.par(gap.degree=c(3,rep(2, nrow(cor1)-1)), start.degree = 180)
chordDiagram(cor1, grid.col=rainbow(ncol(rt)),  transparency = 0.5, symmetric = T)
par(xpd=T)
#绘制图例
colorlegend(col, vertical = T,labels=c(1,0,-1),xlim=c(1.1,1.3),ylim=c(-0.4,0.4))
dev.off()
circos.clear()

#绘制相关性图形
pdf(file="corrplot.pdf", width=8, height=8)
corrplot(cor1,
         method = "circle",
         order = "hclust",
         tl.col="black", addCoef.col = "black",
         type = "upper",insig = "pch",tl.cex=0.8,number.cex = 0.5,
         col=colorRampPalette(c("#FE9601", "white", "#FC4CC9"))(50)
         )
dev.off()

