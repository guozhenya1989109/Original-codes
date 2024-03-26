

#���ð�
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

inputFile="diffGeneExp(CRG).txt"      #���������ļ�


#��ȡ�����ļ�
data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

#���ɭ����ģ��
control=trainControl(method="repeatedcv", number = 5, savePredictions=TRUE)
mod_rf = train(Type ~ .,data = data, method='rf', trControl = control)

#����ѧϰģ��
mod_svm=train(Type ~., data = data, method = "svmRadial", prob.model = TRUE, trControl=control)

#����Ԥ�⺯��
p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(data$Type=="Control", 0, 1)

#���ɭ����ģ��Ԥ����
explainer_rf=explain(mod_rf, label = "RF",
                     data = data, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)
#����ѧϰģ��Ԥ����
explainer_svm=explain(mod_svm, label = "SVM",
                      data = data, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)

#���Ʋв�ķ����ۼƷֲ�ͼ
pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm)
print(p1)
dev.off()


#���Ʋв������ͼ
pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, geom = "boxplot")
print(p2)
dev.off()


#����ROC����
pred1=predict(mod_rf, newx=data, type="prob")
pred2=predict(mod_svm, newx=data, type="prob")
roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
ci1=ci.auc(roc1, method="bootstrap")
ci2=ci.auc(roc2, method="bootstrap")
pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="#FC4CC9")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="#FE9601", add=T)
legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc))),
       col=c("#FC4CC9","#FE9601"), lwd=2, bty = 'n')
dev.off()

