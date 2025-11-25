
[1] The R version is used for the method development is :
R version 4.5.1 (2025-06-13 ucrt) -- "Great Square Root"

[2] Install proper R library before loading the function

[3] Execute the "Flexitau Paper Source - utility function.R" file to preload the R utility functions.

[4] We used different types of files : "PTM", "cleavage", and "FLEXITau" data, respectively, to build the model and conduct the analysis. 

[5] The source code "Flexitau Paper Source - data preprocessing.R" contain modularized code to read each training-testing paired samples, based on the 3 types of data discussed in step [4].

[6] After loading the proper training data set and test data set in step [5], run the code "Flexitau Paper Source - downstream analysis.R" in this step for model construction . This script does the following things

[6.1] : For each diseases (in total, we have 7 different tauopathies: "AD"   "CBD"  "CTE"  "CTRL" "PiD"  "PSP" "DLB")

  		For each machine learning model ( in total, we chose 3 different models: "LASSO" "RandomForest" "Boruta" )

			save a disease-model specific machine learning classifier : write down the model into a object list : recorder the model's feature importance and write them into a table

	Return the model important attributes as a table for each disease (i.e. AD_Sampled_DataMatrix.txt)
	Return the optimal machine learning model as a R object (i.e. AD_Random_Forest_Model.RData), which could be read and loaded for future prediction given any test dataset

[6.2] : the model will also generate a detailed list of figures for performance evaluation for each disease type, along with details for each machine learning model. For example, for each disease (i.e. AD disease) model, we will have a figure shows the AUC values with respect to lambda cutoff in the LASSO model, the feature importance figure, the OOB error vs mTry in the random forest plot, and the Random Forest importance dot-plot, another figure will be generated on the feature importance based on the MeanDecreaseGini. In total , 12 figures will be generated for each disease. For 7 tauopathies, it will be 84 figures.  


[6.3] : the model will also generate AUC and random forest cross-validation accuracy figure in a file ending with "ROC_7_tauopathies.pdf". A total of 3 panels will be generated, representing the three machine learning methods. Each panel shows the 7 AUC curve, one curve for each disease type.

