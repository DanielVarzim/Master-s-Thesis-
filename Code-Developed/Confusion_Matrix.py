# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import numpy as np
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
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn import preprocessing
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import VotingClassifier
from sklearn.metrics.classification import confusion_matrix
from pandas_ml import ConfusionMatrix

def load_files(NumDataset): 
    #carrega os ficheiros 
    Input = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset_mixed.csv", delimiter=",")   #X
    Output = genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Out_Attributes_Dataset_mixed.csv", delimiter=",")   #Y
    return Input, Output

def load_files_transporter():
    Input = genfromtxt("../Data/Datasets/Dataset1/TransporterDataset/Features/Features_DatasetTCDBIDs1_mixed.csv", delimiter=",")   #X
    Output = genfromtxt("../Data/Datasets/Dataset1/TransporterDataset/Out_attributes/Out_Attributes_Dataset_TCDBIDs1_mixed.csv", delimiter=",")   #Y
    return Input, Output

def cross_val(Input, Output, prop = 0.3): 
    #efectua cross validation dos dados para determinar os datasets de treino e de teste
    X_train,X_test,y_train,y_test = cross_validation.train_test_split(Input, Output, test_size=prop)
    return X_train,X_test,y_train,y_test

if __name__ == '__main__':
    start_time1 = time.time()
    Trans=input("""Perform Confusion matrix on Dataset:\n 1-Dataset1
2-Dataset2
3-Dataset12345
4-Transporter Dataset
    """)
    if Trans =="4":
        directory="../Data/Datasets/Dataset1/TransporterDataset/Statistical_analysis/"
        
        if not os.path.exists(directory):
            os.makedirs(directory)
        start_time1 = time.time()
        Input,Output=load_files_transporter()
        Input = Imputer().fit_transform(Input) 
        
    if Trans=="1" or Trans=="2" or Trans=="3":
        NumDataset=input("Dataset number?")

        directory="../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        start_time1 = time.time()
        Input,Output=load_files(NumDataset)
        Input = Imputer().fit_transform(Input) 

   
    
    #Variabilidade
    sel    =    VarianceThreshold()    
    filt    =    sel.fit_transform (Input)    
    print(filt.shape)
    
    ################################################################
    if Trans=="1" or Trans=="2":
        filt=np.delete(filt,27575 ,axis=0)
        print(filt.shape)
        
        Output=np.delete(Output, 27575, axis=0)
        print(Output.shape)
    
    if Trans =="3":
        filt=np.delete(filt, [23429,23430,23431,23432,23434] ,axis=0)
        print(filt.shape)
         
        Output=np.delete(Output, [23429,23430,23431,23432,23434], axis=0)
        print(Output.shape)
    
    if Trans=="4":
        filt=np.delete(filt, [11715,11716,11717] ,axis=0)
        print(filt.shape)
         
        Output=np.delete(Output, [11715,11716,11717], axis=0)
        print(Output.shape)
    
    SplitedX=np.split(filt, 5, axis=0)
    SplitedY=np.split(Output,5, axis=0)
    
    X1=SplitedX[0]
    print(X1.shape)
    X2=SplitedX[1]
    print(X2.shape)
    X3=SplitedX[2]
    print(X3.shape)    
    X4=SplitedX[3]
    print(X4.shape)
    X5=SplitedX[4]
    print(X5.shape)
        
    Y1=SplitedY[0]
    print(Y1.shape)
    Y2=SplitedY[1]
    print(Y2.shape)
    Y3=SplitedY[2]
    print(Y3.shape)
    Y4=SplitedY[3]
    print(Y4.shape)
    Y5=SplitedY[4]
    print(Y5.shape)
    ################################################################
    
    XTrain1=np.concatenate((X1,X2,X3,X4),axis=0)
    YTrain1=np.concatenate((Y1,Y2,Y3,Y4),axis=0)
    
    XTrain2=np.concatenate((X1,X2,X3,X5),axis=0)
    YTrain2=np.concatenate((Y1,Y2,Y3,Y5),axis=0)
    
    XTrain3=np.concatenate((X1,X2,X4,X5),axis=0)
    YTrain3=np.concatenate((Y1,Y2,Y4,Y5),axis=0)
    
    XTrain4=np.concatenate((X1,X3,X4,X5),axis=0)
    YTrain4=np.concatenate((Y1,Y3,Y4,Y5),axis=0)
    
    XTrain5=np.concatenate((X2,X3,X4,X5),axis=0)
    YTrain5=np.concatenate((Y2,Y3,Y4,Y5),axis=0)
    
    elapsed_time1 = time.time() - start_time1
    print("Time spent Filtering and Dividing the Dataset: %s" % elapsed_time1)
    
    start_time2 = time.time()
    
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
    
    elapsed_time2 = time.time() - start_time2
    print("Time spent creating the models: %s" % elapsed_time2)
    
    for X_train,y_train, X_test,y_test, label in zip([XTrain1,XTrain2,XTrain3,XTrain4,XTrain5],[YTrain1,YTrain2,YTrain3,YTrain4,YTrain5],
                                                     [X5,X4,X3,X2,X1],[Y5,Y4,Y3,Y2,Y1],["Test5","Test4","Test3","Test2","Test1"]):
        start_time3 = time.time()
        if Trans=="1"or Trans=="2" or Trans=="3":
            file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/Confusion_Matrices_test_"+label+".txt", "w")
            file.close()
        elif Trans=="4":
            file=open("../Data/Datasets/Dataset1/TransporterDataset/Statistical_analysis/Confusion_Matrices_"+label+".txt", "w")
            file.close()
        
        #Saving the NB model 
        NBM=NB_model.fit(X_train,y_train)
        #Saving the Tree model
        treeM=tree_model.fit(X_train,y_train)       
        #Saving the knn model 
        knnM = pipeknn.fit(X_train,y_train)
        #Saving the logistic regression model
        logisticM=pipelogistic.fit(X_train,y_train) 
        #Saving the GB model 
        GBM=GB_model.fit(X_train,y_train)    
        #Saving the RF model 
        RFM=RF_model.fit(X_train,y_train)     
        #Saving the SVM model
        svmM=pipesvm.fit(X_train,y_train)
          
        
        elapsed_time3 = time.time() - start_time3
        print("Time spent training the models: %s" % elapsed_time3)
        
        ####################################################################################
        start_time4 = time.time()
        
        predicted_NBM=NBM.predict(X_test)
        predicted_treeM=treeM.predict(X_test)
        predicted_knnM=knnM.predict(X_test)
        predicted_logisticM=logisticM.predict(X_test)
        predicted_GBM=GBM.predict(X_test)
        predicted_RFM=RFM.predict(X_test)
        predicted_svmM=svmM.predict(X_test)
        

        
        #=======================================================================
        # matrixNB=confusion_matrix(y_test,predicted_NBM)
        # matrixTree=confusion_matrix(y_test, predicted_treeM)
        # matrixKnn= confusion_matrix(y_test, predicted_knnM)
        # matrixLogistic= confusion_matrix(y_test, predicted_logisticM)
        # matrixGB= confusion_matrix(y_test, predicted_GBM)
        # matrixRF= confusion_matrix(y_test, predicted_RFM)
        # matrixSVM=confusion_matrix(y_test, predicted_svmM)
        #=======================================================================
        
        #=======================================================================
        # #Com labels=[0,1]
        # matrixNB=confusion_matrix(y_test,predicted_NBM, labels)
        # matrixTree=confusion_matrix(y_test, predicted_treeM,labels)
        # matrixKnn= confusion_matrix(y_test, predicted_knnM,labels)
        # matrixLogistic= confusion_matrix(y_test, predicted_logisticM,labels)
        # matrixGB= confusion_matrix(y_test, predicted_GBM,labels)
        # matrixRF= confusion_matrix(y_test, predicted_RFM,labels)
        # matrixSVM=confusion_matrix(y_test, predicted_svmM,labels)
        #=======================================================================
        
        #Com pandas_ml
        matrixNB=ConfusionMatrix(y_test,predicted_NBM)
        matrixTree=ConfusionMatrix(y_test, predicted_treeM)
        matrixKnn= ConfusionMatrix(y_test, predicted_knnM)
        matrixLogistic= ConfusionMatrix(y_test, predicted_logisticM)
        matrixGB= ConfusionMatrix(y_test, predicted_GBM)
        matrixRF= ConfusionMatrix(y_test, predicted_RFM)
        matrixSVM=ConfusionMatrix(y_test, predicted_svmM)
        
        print(matrixNB)
        print(matrixTree)
        print(matrixKnn)
        print(matrixLogistic)
        print(matrixGB)
        print(matrixRF)
        print(matrixSVM)
        
        elapsed_time4 = time.time() - start_time4
        print("Time spent making prediction on the test models %s:" % elapsed_time4)
        
        
        if Trans=="1"or Trans=="2" or Trans=="3":
            file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Statistical_analysis/Confusion_Matrices_test_"+label+".txt", "a")
        elif Trans=="4":
            file=open("../Data/Datasets/Dataset1/TransporterDataset/Statistical_analysis/Confusion_Matrices_"+label+".txt", "a")
        file.write("Confusion Matrix: NB model\n")
        file.write(str(matrixNB))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: ExtraTree model\n")
        file.write(str(matrixTree))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: KNN model\n")
        file.write(str(matrixKnn))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: Logistic Regression model\n")
        file.write(str(matrixLogistic))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: GradientBoosting model\n")
        file.write(str(matrixGB))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: RandomForest model\n")
        file.write(str(matrixRF))
        file.write("\n-----------------------------------------------------\n")
        file.write("Confusion Matrix: SVM model\n")
        file.write(str(matrixSVM))
        file.write("\n-----------------------------------------------------\n")
        file.close()
        
    