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
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import VotingClassifier
    
def load_files(NumDataset, numIDs): 
    #carrega os ficheiros 
    Input = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Features/Features_Dataset"+numIDs+"_mixed.csv", delimiter=",")   #X
    Output = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Out_attributes/Out_Attributes_Dataset_"+numIDs+"_mixed.csv",delimiter=",")   #Y
    return Input, Output
    
  
if __name__ == '__main__':
    NumDataset=input("Dataset number?")
    lvlTC=input("Level of TC system to reach?")
    if lvlTC=="1":
        numIDs="TCDBIDs1"
    elif lvlTC=="2":
        numIDs="TCDBIDs2"
    else:
        print("Error")

    directory="../Data/BestModels/Using_TransporterDataset"
 
    if not os.path.exists(directory):
        os.makedirs(directory)
    start_time1 = time.time()
    Input,Output=load_files(NumDataset, numIDs)
    Input = Imputer().fit_transform(Input)
   
    
    #Variabilidade
    sel    =    VarianceThreshold()    
    filt    =    sel.fit_transform (Input)    
    print(filt.shape)
    print("VarianceThreshold filter used")
    print("StandardScaler used for KNN, LogisticRegression and SVM Models")

    directory2="../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis"
    if not os.path.exists(directory2):
        os.makedirs(directory2)
    
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis/BestModels.txt", "w")
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
        file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis/BestModels.txt", "a")
        scores = cross_validation.cross_val_score(model, filt,Output, cv=5)
        file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),label))
        print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),label)) 
        file.close
      
    bestWeightedHardVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',pipesvm)],voting='hard',weights=[1,2,2,3,3,2,3])
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis/BestModels.txt", "a")
    file.write("Weighted Hard Vote [1,8,3,10,9,6,3]\n")
    print("WeightedHard Vote")
    scores = cross_validation.cross_val_score(bestWeightedHardVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"bestWeightedHardVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"bestWeightedHardVote")) 
    file.close
    
    bestHardVote=VotingClassifier(estimators=[('NB',NB_model),('tree',tree_model),('knn',pipeknn),('logistic',pipelogistic),('GB',GB_model),('RF',RF_model),('SVM',pipesvm)],voting='hard')
    
    HardVote=VotingClassifier(estimators=[('logistic',pipelogistic),('GB',GB_model),('tree',tree_model)],voting='hard')
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis/BestModels.txt", "a")
    file.write("Hard Vote with Logistic, GB and ExtraTree models\n")
    print("Hard Vote")
    scores = cross_validation.cross_val_score(HardVote, filt,Output, cv=5)
    file.write("PECC: %0.2f (+/- %0.2f) [%s]\n" % (scores.mean(),scores.std(),"HardVote"))
    print("PECC: %0.2f (+/- %0.2f) [%s]" % (scores.mean(),scores.std(),"HardVote")) 
    file.close
    
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Statistical_analysis/BestModels.txt", "a")   
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
     
    #Saving the weigthedHardVote model
    WHVM=bestWeightedHardVote.fit(filt,Output)
    learned_WHVM=directory+"/bestWeightedHardVote_model.pkl"
    joblib.dump(WHVM,learned_WHVM)
      
    #Saving the HardVote model
    HVM=bestHardVote.fit(filt,Output)
    learned_HVM=directory+"/bestHardVote_model.pkl"
    joblib.dump(HVM,learned_HVM)
  
    #===========================================================================
    # Saving the SoftVote model
    # SVM=SoftVote.fit(filt,Output)
    # learned_SVM=directory+"/SoftVote_model.pkl"
    # joblib.dump(SVM,learned_SVM)
    #===========================================================================
         
