# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''


from numpy import genfromtxt
from sklearn import cross_validation
from sklearn import svm
from sklearn import linear_model
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import BaggingClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
import time
import os
from sklearn.preprocessing import Imputer
from sklearn.feature_selection import VarianceThreshold 
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import VotingClassifier
    
def load_files(NumDataset): 
    #carrega os ficheiros 
    Input = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset_mixed.csv", delimiter=",")   #X
    Output = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Out_Attributes_Dataset_mixed.csv", delimiter=",")   #Y
    return Input, Output
    
  
if __name__ == '__main__':
    NumDataset=input("Dataset number?")

    directory="../Data/BestModels/Using_Dataset"+str(NumDataset)+"mixed"
 
    if not os.path.exists(directory):
        os.makedirs(directory)
    start_time1 = time.time()
    Input,Output=load_files(NumDataset)
    Input = Imputer().fit_transform(Input)
   
    
    #Variabilidade
    sel    =    VarianceThreshold()    
    filt    =    sel.fit_transform (Input)    
    print(filt.shape)
    print("VarianceThreshold filter used")
    print("StandardScaler used for KNN, LogisticRegression and SVM Models")

    directory2="../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis"
    if not os.path.exists(directory2):
        os.makedirs(directory2)
    
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "w")
    file.close()
    

    NB_model= BaggingClassifier(GaussianNB())
    tree_model = BaggingClassifier(ExtraTreesClassifier())
    knn = BaggingClassifier(KNeighborsClassifier())
    logistic = BaggingClassifier(linear_model.LogisticRegression())
    GB_model= GradientBoostingClassifier()
    RF_model= BaggingClassifier(RandomForestClassifier())
    svm_model = svm.SVC()
    
    pipeknn=make_pipeline(preprocessing.StandardScaler(),knn)
    pipelogistic=make_pipeline(preprocessing.StandardScaler(),logistic)
    pipesvm=make_pipeline(preprocessing.StandardScaler(),svm_model)
    
    for model,label in zip([NB_model,tree_model,knn,logistic,GB_model,RF_model,svm_model],["Naive Bayes","ExtraTreeClassifier","K nearest neighbours","Logistic Regression","GradientBoosting","RandomForest","SVM"]):
        file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "a")
        scores = cross_validation.cross_val_score(model, filt,Output, cv=5)
        file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),label))
        print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),label)) 
        scoresf1 = cross_validation.cross_val_score(model, filt,Output, cv=5,scoring="f1")
        file.write("F1: %0.2f (+/- %0.2f) [%s]\n" % (scoresf1.mean(),scoresf1.std(),label))
        print("F1: %0.2f (+/- %0.2f) [%s]" % (scoresf1.mean(),scoresf1.std(),label))
        scoresAUC= cross_validation.cross_val_score(model, filt, Output, cv=5, scoring="roc_auc")
        file.write("ROC-AUC: %0.2f (+/- %0.2f) [%s]\n" % (scoresAUC.mean(),scoresAUC.std(),label))
        print("ROC-AUC: %0.2f (+/- %0.2f) [%s]" % (scoresAUC.mean(),scoresAUC.std(),label))
        file.close
    
    HardVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',pipesvm)],voting='hard',weights=[1,2,2,3,3,2,3])
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "a")
    file.write("Hard Vote\n")
    print("Hard Vote")
    scores = cross_validation.cross_val_score(HardVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"HardVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"HardVote")) 
    scoresf1 = cross_validation.cross_val_score(HardVote, filt,Output, cv=5,scoring="f1")
    file.write("F1: %0.2f (+/- %0.2f) [%s]\n" % (scoresf1.mean(),scoresf1.std(),"HardVote"))
    print("F1: %0.2f (+/- %0.2f) [%s]" % (scoresf1.mean(),scoresf1.std(),"HardVote"))
    file.close
  
    WeightedHardVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',pipesvm)],voting='hard',weights=[1,6,3,10,10,7,4])
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "a")
    file.write("Weighted Hard Vote weights [1,6,3,10,10,7,4]\n")
    print("WeightedHard Vote")
    scores = cross_validation.cross_val_score(WeightedHardVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"WeightedHardVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"WeightedHardVote")) 
    scoresf1 = cross_validation.cross_val_score(WeightedHardVote, filt,Output, cv=5,scoring="f1")
    file.write("F1: %0.2f (+/- %0.2f) [%s]\n" % (scoresf1.mean(),scoresf1.std(),"WeightedHardVote"))
    print("F1: %0.2f (+/- %0.2f) [%s]" % (scoresf1.mean(),scoresf1.std(),"WeightedHardVote"))
    file.close
    
    bestHardVote=VotingClassifier(estimators=[('logistic',pipelogistic),('GB',GB_model),('SVM',pipesvm)],voting='hard')
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "a")
    file.write("Hard Vote with Logistic Regression, GB and SVM models\n")
    print("Hard Vote")
    scores = cross_validation.cross_val_score(WeightedHardVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"WeightedHardVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"WeightedHardVote")) 
    scoresf1 = cross_validation.cross_val_score(WeightedHardVote, filt,Output, cv=5,scoring="f1")
    file.write("F1: %0.2f (+/- %0.2f) [%s]\n" % (scoresf1.mean(),scoresf1.std(),"WeightedHardVote"))
    print("F1: %0.2f (+/- %0.2f) [%s]" % (scoresf1.mean(),scoresf1.std(),"WeightedHardVote"))
    file.close
    
    bestWeightedHardVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',pipesvm)],voting='hard',weights=[1,3,3,10,10,3,10])
    
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/BestModels.txt", "a")   
    file.write("Soft Vote\n")  
    print("Soft Vote")
    SoftVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',svm_model)],voting='soft',weights=[1,2,1,3,3,2,3])
    scores = cross_validation.cross_val_score(SoftVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"SoftVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"SoftVote"))
    scoresf1 = cross_validation.cross_val_score(SoftVote, filt,Output, cv=5,scoring="f1")
    file.write("F1: %0.2f (+/- %0.2f) [%s]\n" % (scoresf1.mean(),scoresf1.std(),"SoftVote"))
    print("F1: %0.2f (+/- %0.2f) [%s]" % (scoresf1.mean(),scoresf1.std(),"SoftVote"))   
    file.close
 
    #Saving the NB model 
    NBM=NB_model.fit(filt,Output)
    learned_NB=directory+"/NB_model.pkl"
    joblib.dump(NBM,learned_NB)
        
         
    #Saving the Tree model
    treeM=tree_model.fit(filt,Output)
    learned_tree=directory+"/Tree_model.pkl"
    joblib.dump(treeM,learned_tree)
             
    #Saving the knn model 
    knnM = knn.fit(filt,Output)
    learned_knn=directory+"/Knn_model.pkl"
    joblib.dump(knnM,learned_knn)
          
    #Saving the logistic regression model
    logisticM=logistic.fit(filt,Output)
    learned_logistic=directory+"/Logistic_model.pkl"
    joblib.dump(logisticM,learned_logistic)
        
          
    #Saving the GB model 
    GBM=GB_model.fit(filt,Output)
    learned_GB=directory+"/GB_model.pkl"
    joblib.dump(GBM,learned_GB)
           
    #Saving the RF model 
    RFM=RF_model.fit(filt,Output)
    learned_RF=directory+"/RF_model.pkl"
    joblib.dump(RFM,learned_RF)
             
    #Saving the SVM model
    svmM=svm_model.fit(filt,Output)
    learned_svm=directory+"/SVM_model.pkl"
    joblib.dump(svmM,learned_svm)
    
    #Saving the best HardVote model
    BHVM=bestHardVote.fit(filt,Output)
    learned_HVM=directory+"/bestHardVote_model.pkl"
    joblib.dump(BHVM,learned_HVM)
    
    #Saving the best weigthedHardVote model
    WHVM=bestWeightedHardVote.fit(filt,Output)
    learned_WHVM=directory+"/bestWeightedHardVote_model.pkl"
    joblib.dump(WHVM,learned_WHVM)
        
    #Saving the HardVote model
    HVM=HardVote.fit(filt,Output)
    learned_HVM=directory+"/HardVote_model.pkl"
    joblib.dump(HVM,learned_HVM)
     

           
    