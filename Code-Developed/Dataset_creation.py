# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

from numpy import *
import numpy as np
import pandas



def insert_features_OUTattribute(features,OUTattributes):
    """
    Inserts all the features created and the Out attribute into the final dataset
    Returns the final dataset
    """
    
    finalDataset=np.array((features))
    finalDataset=np.hstack((finalDataset,OUTattributes))
    return finalDataset

def mix_dataset(NumDataset):
    dataset=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Dataset"+str(NumDataset)+".csv", sep=",", header=None)
    #print(dataset.shape)
    MD=dataset.iloc[np.random.permutation(len(dataset))]
    mixedDataset=np.array(MD)
    print("mixed dataset shape")
    print(mixedDataset.shape)
    #print(mixedDataset)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Dataset"+str(NumDataset)+"mixed.csv",mixedDataset,delimiter=",")

def split_features_OUTattributes(NumDataset,features):
    mixed_dataset=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Dataset"+str(NumDataset)+"mixed.csv",delimiter=",")
    mixed_features_dataset= np.delete(mixed_dataset, s_[len(features[0]):], axis=1) 
    mixed_features_dataset.astype(float64)
    mixed_OUTattributes_dataset=np.delete(mixed_dataset,s_[0:len(features[0])],axis=1) 
    mixed_OUTattributes_dataset.astype(float64)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset_mixed.csv",mixed_features_dataset, delimiter=",")
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Out_Attributes_Dataset_mixed.csv",mixed_OUTattributes_dataset,delimiter=",")
    print("mixed feature shape")
    print(mixed_features_dataset.shape)
    print("mixed out shape")
    print(mixed_OUTattributes_dataset.shape)
    
if __name__ == '__main__':
    NumDataset=input("Qual o número do dataset?")
    
    features=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset.csv", delimiter=",")
    #print("features shape")
    #print(features.shape)
    
    #usar quando só quisermos saber se é transportadora ou não
    Ot=genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Is_Transporter.csv", delimiter=",")
    OUTattributes=Ot[:,np.newaxis]
    dataset=insert_features_OUTattribute(features, OUTattributes)

    #print("dataset shape")
    #print (dataset.shape)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Dataset"+str(NumDataset)+".csv",dataset,delimiter=",")
    
    mix_dataset(NumDataset)
    split_features_OUTattributes(NumDataset,features)
