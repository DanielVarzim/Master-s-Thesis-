# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

from numpy import *
import numpy as np
import pandas
import os


def insert_features_OUTattribute(features,OUTattributes):
    """
    Inserts all the features created and the Out attribute into the final dataset
    Returns the final dataset
    """
    
    finalDataset=np.array((features))
    finalDataset=np.hstack((finalDataset,OUTattributes))
    return finalDataset

def mix_dataset(NumDataset,numIDs):
    dataset=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numIDs+".csv", sep=",", header=None)
    #print(dataset.shape)
    MD=dataset.iloc[np.random.permutation(len(dataset))]
    mixedDataset=np.array(MD)
    print("mixed dataset shape")
    print(mixedDataset.shape)
    #print(mixedDataset)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numIDs+"mixed.csv",mixedDataset,delimiter=",")

def split_features_OUTattributes(NumDataset,features,numIDs):
    mixed_dataset=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numIDs+"mixed.csv",delimiter=",")
    mixed_features_dataset= np.delete(mixed_dataset, s_[len(features[0]):], axis=1) 
    mixed_OUTattributes_dataset=np.delete(mixed_dataset,s_[0:len(features[0])],axis=1) 
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Features/Features_Dataset"+numIDs+"_mixed.csv",mixed_features_dataset, delimiter=",")
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Out_attributes/Out_Attributes_Dataset_"+numIDs+"_mixed.csv",mixed_OUTattributes_dataset,delimiter=",")
    print("mixed feature shape")
    print(mixed_features_dataset.shape)
    print("mixed out shape")
    print(mixed_OUTattributes_dataset.shape)


if __name__ == '__main__':
    NumDataset=input("Qual o número do dataset?")
    directory="../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset"
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    directory2="../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Features"
    if not os.path.exists(directory2):
        os.makedirs(directory2)
    
    directory3="../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Out_attributes"
    if not os.path.exists(directory3):
        os.makedirs(directory3)
    
    features=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset.csv", delimiter=",")
    index=list(range(11717,27577))

    
    numID1="TCDBIDs1"
    Ot1=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB1ID.csv",delimiter=",")
    OUT1=Ot1[:,np.newaxis]
    dataset1=insert_features_OUTattribute(features, OUT1)
    cutedDataset1=np.delete(dataset1, index, axis=0)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numID1+".csv",cutedDataset1,delimiter=",")
    mix_dataset(NumDataset, numID1)
    split_features_OUTattributes(NumDataset, features, numID1)
    
    numID12="TCDBIDs12"
    Ot12=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB12ID.csv",delimiter=",")
    OUT12=Ot12[:,np.newaxis]
    dataset12=insert_features_OUTattribute(features, OUT12)
    cutedDataset12=np.delete(dataset12, index, axis=0)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numID12+".csv",cutedDataset12,delimiter=",")
    mix_dataset(NumDataset, numID12)
    split_features_OUTattributes(NumDataset, features, numID12)
    
    numID123="TCDBIDs123"
    Ot123=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB123ID.csv",delimiter=",")
    OUT123=Ot123[:,np.newaxis]
    dataset123=insert_features_OUTattribute(features, OUT123)
    cutedDataset123=np.delete(dataset123, index, axis=0)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numID123+".csv",cutedDataset123,delimiter=",")
    mix_dataset(NumDataset, numID123)
    split_features_OUTattributes(NumDataset, features, numID123)
    
    numID1234="TCDBIDs1234"
    Ot1234=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB1234ID.csv",delimiter=",")
    OUT1234=Ot1234[:,np.newaxis]
    dataset1234=insert_features_OUTattribute(features, OUT1234)
    cutedDataset1234=np.delete(dataset1234, index, axis=0)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/TransporterDataset/Dataset"+numID1234+".csv",cutedDataset1234,delimiter=",")
    mix_dataset(NumDataset, numID1234)
    split_features_OUTattributes(NumDataset, features, numID1234)
    
    #===========================================================================
    # Ot=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Out_Attributes_Dataset.csv",delimiter=",")
    # 
    # OUTattributes=Ot[:,np.newaxis]
    # 
    # #Usar esta parte se o dataset ainda não estiver criado
    # dataset=insert_features_OUTattribute(features, OUTattributes)
    # #
    # print("dataset shape")
    # print (dataset.shape)
    # np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Dataset"+str(NumDataset)+".csv",dataset,delimiter=",")
    # 
    # mix_dataset(NumDataset)
    # split_features_OUTattributes(NumDataset,features)
    #===========================================================================