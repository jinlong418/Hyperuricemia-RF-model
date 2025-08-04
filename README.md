# Hyperuricemia-RF-model

The files and the code in this repository correspond to the following manuscript:

Multi-cohort analysis unveils novel microbial targets for the treatment of hyperuricemia and gout

To determine if microbial-related biomarkers could distinguish HC from the HUA and gout groups, we used the R package randomForest v4.6.14 to construct a classification model for disease diagnosis.

Data preparing:after_remove_batch_effect_table_Genus.txt/sample-info.txt/target_genus.txt

code:randomForest_roc_multi-classification.r

RF model:trained_model.rds

result:HC vs HUA vs Gout_confusionMatrix_Accuracy.txt/HC vs HUA vs Gout_importance.csv
