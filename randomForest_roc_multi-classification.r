library(tidyverse)
library(randomForest)
library(caret)
library(pROC)
library(data.table)

input_abundance_file<-"after_remove_batch_effect_table_Genus.txt"
input_target_genus_file<-"target_genus.txt"
order<-c("HC","HUA","Gout")
title_name<-paste(order,collapse = " vs ")

#Prepare the relative abundance table of genus levels
table_Genus<-data.frame(fread(input_abundance_file,
sep = '\t',check.names=F,header=TRUE),check.names=F)
dim(table_Genus) 
table_Genus<-aggregate(table_Genus[,2:ncol(table_Genus)],
by=list(genus=table_Genus[,1]),FUN=sum)
rownames(table_Genus)<-table_Genus[,1]
table_Genus<-table_Genus[,-1]
dim(table_Genus)
table_Genus<-table_Genus[,names(which(colSums(table_Genus)>3000))]
a<-table_Genus
a[a>0]<-1
table_Genus_filter<-table_Genus[which(apply(a,1,sum)>ncol(a)*0.1),]
abundance=0.0001
table_Genus_relative<-data.frame(apply(table_Genus_filter,2,function(x) {x/sum(x)}),check.names=F)
table_Genus_relative<-table_Genus_relative[which(rowMeans(table_Genus_relative)>=abundance),]
table_Genus_relative<-data.frame(apply(table_Genus_relative,2,function(x) {x/sum(x)}))

metadata<-read.table("sample-info.txt",sep = '\t',
header=TRUE,check.names=F)
metadata<-subset(metadata,Group%in%order)
common_name<-intersect(colnames(table_Genus_relative),metadata$sample_id)
metadata<-subset(metadata,sample_id%in%common_name)
metadata<-metadata[order(factor(metadata$Group,levels=order)),]
table_Genus_relative<-table_Genus_relative[,metadata$sample_id]

#7:3 Divide the training set and the test set
set.seed(123)
data_split <- createDataPartition(metadata$Group, p = 0.7, list = FALSE)
target_genus<-data.frame(fread(input_target_genus_file,sep = '\t',check.names=F,header=F),check.names=F)
Genus_order<-table_Genus_relative[target_genus[,1],]
Genus_order1<-t(Genus_order)
Genus_order1<-data.frame(Genus_order1[metadata$sample_id,])
colnames(Genus_order1)<-rownames(Genus_order)

train_data <- Genus_order1[data_split, ]
train_metadata<-metadata[data_split,]
train_metadata$Group<-factor(train_metadata$Group,levels=order)
test_data <- Genus_order1[-data_split, ]
test_metadata<-metadata[-data_split,]
test_metadata$Group<-factor(test_metadata$Group,levels=order)

#Build a random forest model
set.seed(123)
rf <- randomForest(x = train_data,
                       y = factor(train_metadata$Group),
                       ntree = 1000,
                       do.trace=TRUE,
					   importance = TRUE,
					   proximity = TRUE)
saveRDS(rf, "trained_model.rds")
importance<-as.matrix(sort(importance(rf)[,"MeanDecreaseGini"], decreasing = TRUE))
colnames(importance)<-"MeanDecreaseGini"


#Draw the MeanDecreaseGini scatter plot
library(ggplot2)
Genus_classified<-c()
for(i in unique(metadata$Group)){
group<-rowMeans(Genus_order[,which(metadata$Group==i)])
Genus_classified<<-cbind(Genus_classified,group)
}
colnames(Genus_classified)<-as.character(unique(metadata$Group))
Genus_importance<-Genus_classified[rownames(importance),order]
Genus_importance<-data.frame(Genus_importance)

Genus_importance1<-Genus_importance
for(i in 1:nrow(Genus_importance)){
for(j in colnames(Genus_importance)){
if(max(Genus_importance[i,])==Genus_importance[i,j]){
Genus_importance1[i,3]<-j}
}}
colnames(Genus_importance1)[ncol(Genus_importance1)]<-"Group"
importance_plot<-data.frame(Genus=rownames(importance),importance,
Group=Genus_importance1$Group)
importance_plot<-importance_plot[order(importance_plot$MeanDecreaseGini),]
title_name<-paste(order,collapse = " vs ")
write.csv(importance_plot,paste0(title_name,"_importance.csv"),row.names=F)

importance_plot$Genus<-factor(importance_plot$Genus,levels=importance_plot$Genus)
importance_plot$Group<-factor(importance_plot$Group,levels=order)
mycolors<-c("#BC3C29FF","#0072B5FF","#E18727FF")
importance_dotplot<-
ggplot(importance_plot,aes(x=MeanDecreaseGini,y=Genus,colour=Group,
size=MeanDecreaseGini))+
theme_bw()+
theme(text=element_text(size=14,family="sans",color="black"))+
theme(axis.title.x = element_text(size = 14,family="sans",color="black"),
axis.text.y =element_text(size=14,family="sans",color="black"),
axis.text.x =element_text(size=14,family="sans",color="black"),
legend.text=element_text(size=14))+
geom_point()+
ylab("")+
scale_color_manual(values=mycolors)

number=nrow(importance_plot)
if(number<=30){
pdf_height=5}else{
pdf_height=5+2/10*(number-30)
}

len<-max(nchar(as.character(importance_plot$Genus)))
if(len>42){
pdf_width=13
}else{pdf_width=8}

pdf(paste0(title_name,"_importance_dotplot.pdf"),width=pdf_width,height=pdf_height,family="sans")
print(importance_dotplot)
dev.off()


#5-fold cross-validation
inner_predict_result<-c()
set.seed(123)
folds<-createFolds(1:nrow(train_data),k=5)
 for(i in 1:5){
 set.seed(123)
 fold_test<-train_data[folds[[i]],] 
  fold_train<-train_data[-folds[[i]],]
  set.seed(123)
  inner_rf<-randomForest(x=fold_train,
                             y =factor(train_metadata[-folds[[i]],"Group"]),
                               ntree = 1000,
                               do.trace=TRUE)
			temp_predict<-predict(inner_rf,fold_test,type = "prob")
			inner_predict_result<-rbind(inner_predict_result,temp_predict)
}
inner_predict_result<-inner_predict_result[train_metadata$sample_id,]
inner_roc_result <- multiclass.roc(train_metadata$Group,inner_predict_result,levels=order)
inner_roc_result_auc<-inner_roc_result$auc
inner_roc_result_auc<-gsub("Multi-class area under the curve: ","",inner_roc_result_auc)
inner_roc_result_auc<-round(as.numeric(inner_roc_result_auc),digits=3)
write.table(inner_roc_result_auc,paste0(title_name,"_inner_roc_result_auc.txt"),quote=F,sep="\t",row.names=F,col.names=F)

#The model constructed using the training set makes predictions on the test set
predict_result<-data.frame(predict(rf,test_data,type = "prob"))
predict_result1<-cbind(sample_id=rownames(predict_result),predict_result)
write.table(predict_result1,paste0(title_name,"_predict_result.txt"),quote=F,sep="\t",row.names=F)
library(pROC)
roc_objects <- lapply(1:ncol(predict_result),function(cls) {
  class_pred <- predict_result[,cls]
  class_labels <- as.numeric(test_metadata$Group == colnames(predict_result)[cls])
  roc_obj <- roc(class_labels, class_pred)
  return(roc_obj)
})
for(i in 1:length(roc_objects)){
tmp_auc<-gsub("Area under the curve: ","",roc_objects[[i]]$auc)
tmp_auc<-round(as.numeric(tmp_auc),digits=3)
tmp_legend_name<-paste0(colnames(predict_result)[i]," (AUC=",tmp_auc,")")
names(roc_objects)[i]<-tmp_legend_name
}
mycolors<-c("#BC3C29FF","#0072B5FF","#E18727FF")
roc_plot<-ggroc(roc_objects,size = 1.2,alpha=1,
,legacy.axes=T)+
labs(x="False Positive Rate (1-Specificity)",y="True Positive Rate (Sensitivity)")+
theme(panel.grid =element_blank(),panel.background=
element_rect(fill="transparent",colour="black"))+
theme(text=element_text(size=16,family="sans"))+
theme(axis.title.x = element_text(size = 16,family="sans"),
axis.title.y =element_text(size=16,family="sans"),
axis.text.x = element_text(color="black"),
axis.text.y = element_text(color="black"))+
scale_color_manual(values = mycolors)+
labs(colour = "test set")+
theme(legend.position=c(0.55,0.20),
legend.text=element_text(size=16,family="sans"))+
theme(legend.background=element_rect(fill="white", colour="black"))
pdf(paste0(title_name,"_roc_plot.pdf"),width=5,height=5,family="sans")
print(roc_plot)
dev.off()

#confusionMatrix
predict_result_class<-predict(rf,test_data)
predict_result_class1<-data.frame(predict_result_class)
predict_result_class1<-cbind(sample_id=rownames(predict_result_class1),predict_result_class1)
write.table(predict_result_class,paste0(title_name,"_predict_result_class.txt"),quote=F,sep="\t",row.names=F)
cm <- confusionMatrix(predict_result_class,test_metadata$Group)
write.table(cbind(rownames(cm$table),cm$table),paste0(title_name,"_confusionMatrix.txt"),quote=F,sep="\t",row.names=F)
write.table(cbind(rownames(cm$byClass),cm$byClass),paste0(title_name,"_confusionMatrix_Accuracy.txt"),quote=F,sep="\t",row.names=F)
cm_d <- as.data.frame(cm$table)
if(length(order)>2){
cm_st <-data.frame(t(cm$byClass))
colnames(cm_st)<-order
cm_st<-cm_st[c(1,2,3,4,8,11),]
}else{
cm_st <-data.frame(cm$byClass)
row_names<-rownames(cm_st)[c(1,2,3,4,8,11)]
cm_st<-data.frame(cm_st[c(1,2,3,4,8,11),])
rownames(cm_st)<-row_names
colnames(cm_st)<-"value"}

library(ggplot2)  
library(gridExtra)  
library(grid)  		 

# plotting the matrix
cm_d_p <- ggplot(data = cm_d, aes(x = Prediction, y = Reference,fill= Freq)) +
  geom_tile() +
  scale_fill_gradient(low="#D6EAF8",high = "#2E86C1")+
  geom_text(aes(label = paste("", Freq)), color = 'red', size = 8) +
  theme_light() +
  guides(fill=FALSE) +
  theme(text=element_text(size=14, family="sans"),
        axis.text.x = element_text(size=14,color="black"),
        axis.text.y = element_text(size=14,color="black"),
		axis.title.x = element_text(size = 14, color = "black"),
         axis.title.y = element_text(size = 14, color = "black"))
 
# plotting the stats
 cm_st_p <-tableGrob(cm_st)
pdf(paste0(title_name,"_confusionMatrix_plot.pdf"),width=14,height=5,family="sans")
grid.arrange(cm_d_p, cm_st_p,nrow = 1, ncol = 2, 
             top=textGrob("Confusion Matrix and Statistics",gp=gpar(fontsize=18,font=1,fontfamily = "sans")))
dev.off()

