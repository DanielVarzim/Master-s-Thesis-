# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import os
import numpy as np
import pandas


def create_IsTransporter_dataset(data):
    """
    Creates a numpy array with 1 or 0´s depending if the protein is a transporter(1) or 
    a negative case
    """
    i=True
    index = range(len(data))
    for row in index:
        if i:
            Out_Atribute=np.array((data.values[row][4]))
            i=False
        else:
            Out_Atribute=np.vstack((Out_Atribute,data.values[row][4]))
    return Out_Atribute

def __letter_to_index(letter):
    alphabet=[';','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
    for i in range(len(alphabet)):
        if letter == alphabet[i]:
            id2=i
    return id2

def __divide_TCDB_ID(TCDB_ID):
    
    division= TCDB_ID.split(".")
    id1=int(division[0])
    id2letter=division[1]
    id2=__letter_to_index(id2letter)
    id3=int(division[2])
    id4=int(division[3])
    return id1,id2,id3,id4
    

def create_TCDB_ID_dataset(data):
    i=True
    index = range(len(data))
    for row in index:
        if i:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                Out_attributes=np.array((id1,id2,id3,id4))
            else:
                Out_attributes=np.array((0,0,0,0))
            i=False
        else:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                OA=np.array((id1,id2,id3,id4))
                Out_attributes=np.vstack((Out_attributes,OA))
            else:
                OA=np.array((0,0,0,0))
                Out_attributes=np.vstack((Out_attributes,OA))
    return Out_attributes

def create_TCDB_1ID(data):
    i=True
    index = range(len(data))
    for row in index:
        print(row)
        if i:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                Out_attributes=np.array((id1))
            else:
                Out_attributes=np.array((0))
            i=False
        else:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                OA=np.array((id1))
                Out_attributes=np.vstack((Out_attributes,OA))
            else:
                OA=np.array((0))
                Out_attributes=np.vstack((Out_attributes,OA))
    return Out_attributes

def create_TCDB_12ID(data):
    i=True
    index = range(len(data))
    for row in index:
        if i:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2))
                Out_attributes=np.array((joinedIDs))
            else:
                Out_attributes=np.array((0))
            i=False
        else:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2))
                OA=np.array((joinedIDs))
                Out_attributes=np.vstack((Out_attributes,OA))
            else:
                OA=np.array((0))
                Out_attributes=np.vstack((Out_attributes,OA))
    return Out_attributes

def create_TCDB_123ID(data):
    i=True
    index = range(len(data))
    for row in index:
        if i:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2)+"00"+str(id3))
                Out_attributes=np.array((joinedIDs))
            else:
                Out_attributes=np.array((0))
            i=False
        else:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2)+"00"+str(id3))
                OA=np.array((joinedIDs))
                Out_attributes=np.vstack((Out_attributes,OA))
            else:
                OA=np.array((0))
                Out_attributes=np.vstack((Out_attributes,OA))
    return Out_attributes

def create_TCDB_1234ID(data):
    i=True
    index = range(len(data))
    for row in index:
        if i:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2)+"00"+str(id3)+"00"+str(id4))
                Out_attributes=np.array((joinedIDs))
            else:
                Out_attributes=np.array((0))
            i=False
        else:
            TCDB_ID=(data.values[row][5])
            if TCDB_ID != "0":
                id1,id2,id3,id4=__divide_TCDB_ID(TCDB_ID)
                joinedIDs=int(str(id1)+"00"+str(id2)+"00"+str(id3)+"00"+str(id4))
                OA=np.array((joinedIDs))
                Out_attributes=np.vstack((Out_attributes,OA))
            else:
                OA=np.array((0))
                Out_attributes=np.vstack((Out_attributes,OA))
    return Out_attributes


def create_OutAttributes_dataset(OutAttributes):
    """
    Inserts all the out attributes created into a dataset 
    Returns the dataset
    """
    i=True
    for OA in OutAttributes:
        if i:
            finalDataset=np.array((OA[:,np.newaxis]))
            i=False
        else:
            finalDataset=np.hstack((finalDataset,OA))
    return finalDataset

if __name__ == '__main__':
    NumDataset=input("Qual o número do dataset?")
    directory="../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes"
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    #sempre 2 proteinas cujo TCDB ID corresponde ao Acession number.(P0C851-->8.A.64.1.1)e (P81694-->8.B.19.2.2 ) alterar manualmente senão ocorre erro    
    data=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/DataDoms.csv",sep=",")
    
    IsTransporter=create_IsTransporter_dataset(data)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Is_Transporter.csv",IsTransporter,delimiter=",")
      
    TCDB_IDs=create_TCDB_ID_dataset(data)
    ID1=create_TCDB_1ID(data)
    ID12=create_TCDB_12ID(data)
    ID123=create_TCDB_123ID(data)
    ID1234=create_TCDB_1234ID(data)
    
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB_IDs.csv",TCDB_IDs,delimiter=",")
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB1ID.csv",ID1,delimiter=",") 
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB12ID.csv",ID12,delimiter=",")
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB123ID.csv",ID123,delimiter=",")
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB1234ID.csv",ID1234,delimiter=",")
      
    #===========================================================================
    # Is_Transporter=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Is_Transporter.csv",delimiter=",")
    # TCDB_id=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/TCDB_IDs.csv",delimiter=",")
    # OutAttributes=[Is_Transporter,TCDB_id]
    # dataset=create_OutAttributes_dataset(OutAttributes)
    # np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Out_attributes/Out_Attributes_Dataset.csv",dataset,delimiter=",")
    #===========================================================================
     
     
