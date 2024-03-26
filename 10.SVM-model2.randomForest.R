

#引用�?
library(randomForest)
set.seed(123)




#读取输入文件
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(..-*)\\_(.*.)", "\\2", row.names(data))

#随机森林�?
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差�?小的�?
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要�??
importance=importance(x=rf2)

#绘制基因的重要�?�图
pdf(file="geneImportance.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

#挑�?�疾病特征基�?
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>0.5])      #挑�?�重要�?�评分大�?2的基�?
#rfGenes=names(rfGenes[1:5])           #挑�?�重要�?�评分最高的5个基�?
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

#输出重要基因的表达量
sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)



