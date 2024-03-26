
inputFile="merge.txt"     #表达数据文件
conFile="s1.txt"               #对照组的样品信息文件
treatFile="s2.txt"             #实验组的样品信息文件


#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric,(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

#如果数据没有取log2, 会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx [6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  (rt[rt<0]=0)
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

#读取样品信息的文件(对照组和实验组)
sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#设置比较组，进行差异分析
Type=c(rep("con",conNum), rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#输出所有基因的差异情况
allDiff=topTable(fit2, adjust='fdr', number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut, file="all.txt", sep="\t", quote=F, col.names=F)

#输出所有基因矫正后的表达量
Type=c(rep("Control",conNum),rep("Treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

#对差异结果进行过滤,输出显著的差异基因
diffSig=allDiff[with(allDiff, (P.Value< 0.05)), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="diff.txt",sep="\t",quote=F,col.names=F)

#输出显著的差异基因表达量
result1 <- merge(diffSig, rt, by = "row.names")

result1=rt[row.names(diffSig),]

write.table(result1, file="diffGeneExp.txt", sep="\t", quote=F, col.names=T)



#绘制差异基因热图
geneNum=10000    #设置基因的数目
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=data[hmGene,]
#定义注释文件
Type=c(rep("Control",conNum),rep("Treat",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
#绘制热图
pdf(file="heatmap.pdf", width=9, height=6.5)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("#FE9601", "white", "#FC4CC9"))(50),
         cluster_cols =F,
         show_colnames = F,show_rownames=F,
         scale="row",
         fontsize = 7,
         fontsize_row=5,
         fontsize_col=7)
dev.off()


#定义显著性
allDiff$logFC[allDiff$logFC>20]=20
allDiff$logFC[allDiff$logFC< -20]=-20
Significant=ifelse((allDiff$adj.P.Val<adj.P.Val.Filter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Significant))+
  scale_color_manual(values=c("#FE9601", "black", "#FC4CC9"))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
#输出火山图
pdf(file="vol.pdf", width=5.5, height=5)
print(p)
dev.off()


