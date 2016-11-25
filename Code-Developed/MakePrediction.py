# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import xlrd
import numpy as np
from Bio import Entrez
import joblib
import time
import PhobiusAPI as P_API
import BOMP_API as B_API
import LocTreeAPI as L_API
import PfamAPI as Pfam_API


def count_aminoacid(aminoacid,sequence):
    """
    Returns the number of a given aminoacid in a sequence
    """
    sequence=sequence
    count=sequence.count(aminoacid)
    return count

def create_aminoacid_occurrence(seq): #utilização desta feature em: "Functional discrimination of membrane proteins using machine learning techniques"
    """
    Returns an array with the aminoacid occurence of each example in the dataset
    """
    sequence=seq
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19

    feature=np.array((count_aminoacid(aminoacids[0],sequence),count_aminoacid(aminoacids[1],sequence),
                      count_aminoacid(aminoacids[2],sequence),count_aminoacid(aminoacids[3],sequence),
                      count_aminoacid(aminoacids[4],sequence),count_aminoacid(aminoacids[5],sequence),
                      count_aminoacid(aminoacids[6],sequence),count_aminoacid(aminoacids[7],sequence),
                      count_aminoacid(aminoacids[8],sequence),count_aminoacid(aminoacids[9],sequence),
                      count_aminoacid(aminoacids[10],sequence),count_aminoacid(aminoacids[11],sequence),
                      count_aminoacid(aminoacids[12],sequence),count_aminoacid(aminoacids[13],sequence),
                      count_aminoacid(aminoacids[14],sequence),count_aminoacid(aminoacids[15],sequence),
                      count_aminoacid(aminoacids[16],sequence),count_aminoacid(aminoacids[17],sequence),
                      count_aminoacid(aminoacids[18],sequence),count_aminoacid(aminoacids[19],sequence)))

    return feature

    
def create_aminoacid_composition(seq): #utilização desta feature em:"Functional discrimination of membrane proteins using machine learning techniques"
    """
    Returns an array with the aminoacid composition of each example in the dataset
    """
    
    sequence=seq
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19

    feature=np.array((count_aminoacid(aminoacids[0],sequence)/len(sequence),
                      count_aminoacid(aminoacids[1],sequence)/len(sequence),
                      count_aminoacid(aminoacids[2],sequence)/len(sequence),
                      count_aminoacid(aminoacids[3],sequence)/len(sequence),
                      count_aminoacid(aminoacids[4],sequence)/len(sequence),
                      count_aminoacid(aminoacids[5],sequence)/len(sequence),
                      count_aminoacid(aminoacids[6],sequence)/len(sequence),
                      count_aminoacid(aminoacids[7],sequence)/len(sequence),
                      count_aminoacid(aminoacids[8],sequence)/len(sequence),
                      count_aminoacid(aminoacids[9],sequence)/len(sequence),
                      count_aminoacid(aminoacids[10],sequence)/len(sequence),
                      count_aminoacid(aminoacids[11],sequence)/len(sequence),
                      count_aminoacid(aminoacids[12],sequence)/len(sequence),
                      count_aminoacid(aminoacids[13],sequence)/len(sequence),
                      count_aminoacid(aminoacids[14],sequence)/len(sequence),
                      count_aminoacid(aminoacids[15],sequence)/len(sequence),
                      count_aminoacid(aminoacids[16],sequence)/len(sequence),
                      count_aminoacid(aminoacids[17],sequence)/len(sequence),
                      count_aminoacid(aminoacids[18],sequence)/len(sequence),
                      count_aminoacid(aminoacids[19],sequence)/len(sequence)))
    
    return feature



def create_aminoacid_physico_chemical_occurrence(seq):
    """
    Returns an array with the aminoacid composition based on the physico-chemical properties of each example in the dataset
    """
    sequence=seq
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
   
    feature=np.array((count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[1], sequence),
                      count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[19], sequence),
                      count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[17], sequence)+count_aminoacid(aminoacids[18], sequence),
                      count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[5], sequence)+count_aminoacid(aminoacids[2], sequence),
                      count_aminoacid(aminoacids[0], sequence)+count_aminoacid(aminoacids[7], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[14], sequence)+count_aminoacid(aminoacids[15], sequence)+count_aminoacid(aminoacids[16], sequence)+count_aminoacid(aminoacids[18], sequence),
                      count_aminoacid(aminoacids[4], sequence)+count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[12], sequence)+count_aminoacid(aminoacids[19], sequence)+count_aminoacid(aminoacids[17], sequence),
                      count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[8], sequence),
                      count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence),
                      count_aminoacid(aminoacids[0], sequence)+count_aminoacid(aminoacids[4], sequence)+count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[7], sequence)+count_aminoacid(aminoacids[15], sequence)+count_aminoacid(aminoacids[16], sequence),
                      count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[12], sequence)+count_aminoacid(aminoacids[2], sequence)+count_aminoacid(aminoacids[14], sequence)+count_aminoacid(aminoacids[5], sequence)+count_aminoacid(aminoacids[19], sequence),
                      count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[17], sequence)+count_aminoacid(aminoacids[18], sequence)))

    return feature

def create_aminoacid_physico_chemical_composition(seq): 
    """
    Returns an array with the aminoacid composition based on the physico-chemical properties of each example in the dataset
    """
    
    sequence=seq
    
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
    
    feature=np.array(((count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[1], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[19], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[17], sequence)+count_aminoacid(aminoacids[18], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[5], sequence)+count_aminoacid(aminoacids[2], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[0], sequence)+count_aminoacid(aminoacids[7], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[14], sequence)+count_aminoacid(aminoacids[15], sequence)+count_aminoacid(aminoacids[16], sequence)+count_aminoacid(aminoacids[18], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[4], sequence)+count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[12], sequence)+count_aminoacid(aminoacids[19], sequence)+count_aminoacid(aminoacids[17], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[8], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[6], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[0], sequence)+count_aminoacid(aminoacids[4], sequence)+count_aminoacid(aminoacids[3], sequence)+count_aminoacid(aminoacids[7], sequence)+count_aminoacid(aminoacids[15], sequence)+count_aminoacid(aminoacids[16], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[6], sequence)+count_aminoacid(aminoacids[8], sequence)+count_aminoacid(aminoacids[9], sequence)+count_aminoacid(aminoacids[10], sequence)+count_aminoacid(aminoacids[11], sequence)+count_aminoacid(aminoacids[12], sequence)+count_aminoacid(aminoacids[2], sequence)+count_aminoacid(aminoacids[14], sequence)+count_aminoacid(aminoacids[5], sequence)+count_aminoacid(aminoacids[19], sequence))/len(sequence),
                      (count_aminoacid(aminoacids[13], sequence)+count_aminoacid(aminoacids[1], sequence)+count_aminoacid(aminoacids[17], sequence)+count_aminoacid(aminoacids[18], sequence))/len(sequence)))
    
    return feature


def count_dipeptide(dipeptide,sequence):
    """
    Returns the number of a given aminoacid in a sequence
    """
    sequence=sequence
    count=sequence.count(dipeptide)
    return count



def dipeptide_composition(sequence):
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    Dipeptides=[]
    for item in aminoacids:
        pep1=item
        for item in aminoacids:
            pep2=str(pep1)+str(item)
            Dipeptides.append(pep2)
    #print(Dipeptides)
    dipcompo=[]
    for dip in Dipeptides:
        count=count_dipeptide(dip, sequence)
        dipcomposition=count/(len(sequence)-1)
        dipcompo.append(dipcomposition)
    
    dipcompotuple=tuple(dipcompo)
    return dipcompotuple

def create_dipeptide_composition(seq): 
    """
    Returns an array with the dipeptide composition of each example in the dataset
    """
    
    sequence=seq
    feature=np.array(dipeptide_composition(sequence))
    
    return feature

def create_PhobiusJobIds(seq):
    seq=seq
    try:
        PhobiusId= (P_API.run(email="danielvarzim@hotmail.com",
                   format="short",
                   stype="protein",
                   sequence=seq))
        PhobiusId=str(PhobiusId)
        PhobiusId=PhobiusId.strip("b'")
               
    except:
        print("Error on Phobius API")   
        
    return PhobiusId
    
def create_num_alphahelices_signalpeptide(PhobiusJobId): 
    PhobiusJobId=PhobiusJobId
    time.sleep(25)
    try:
        res=P_API.result(job_id=PhobiusJobId,result_type="out")
         
        AlphaHelices=str(res["TM"])
        numberAlphaHelices=int(AlphaHelices.strip("b'"))
         
        SignalPeptide=str(res["SP"])
        if SignalPeptide.strip("b'") == "Y" or "found":
            HaveSignalPeptide=1
        else:
            HaveSignalPeptide=0
        

        feature=np.array((numberAlphaHelices,HaveSignalPeptide))     
    except:
        feature=None
        print("Error on Phobius feature")
          
    return feature

def create_BOMPJobIds(seq):
    seq=seq
    try:
        seq=">SeqName\n"+seq
        BOMPJobId= (B_API.run(seqs=seq))
        #print(BOMPJobId)        
    except:
        print("Error on BOMP API")    
        
    return BOMPJobId

def create_betabarrels(BOMPJobId): 
    BOMPJobId=BOMPJobId
    time.sleep(25)
    try:
        res=B_API.result(job_id=BOMPJobId)              
        feature=np.array((res))
    except:
        feature=None
        print("Error on BOMP feature")
        
    feature=feature.flatten()             
    return feature


def create_LocTree3JobIds(seq,domain):
    Entrez.email = "danielvarzim@hotmail.com"
    seq=seq
    dom=domain
    try:
        if dom =="1":
            dom ="euka"
        elif dom =="2":
            dom="bact"
        elif dom =="3":
            dom="arch"
        else:
            dom="arch"
        seq=">SeqName\n"+seq
        reqid,LocTree3JobId=(L_API.run(domain=dom, email="danielvarzim@hotmail.com",sequence=seq))
     
    except:
        print("Error on LocTree3 API")
    
    return reqid, LocTree3JobId    
            

def create_location_prediction(Loctree3JobId, reqid):
    try:
        if L_API.LoadingFinished(reqid):
        
            res=L_API.result(result_url=Loctree3JobId)
            res=res.lower()
            #print(res)   
            if res=="chloroplast":
                resfinal=1
            elif res =="chloroplast membrane":
                resfinal=2
            elif res =="cytosol":
                resfinal=3
            elif res =="endoplasmic reticulum":
                resfinal=4
            elif res =="endoplasmic reticulum membrane":
                resfinal=5
            elif res =="extra-cellular":
                resfinal=6
            elif res =="ﬁmbrium":
                resfinal=7
            elif res =="Golgi apparatus":
                resfinal=8
            elif res =="Golgi apparatus membrane":
                resfinal=9
            elif res =="mitochondrion":
                resfinal=10
            elif res =="mitochondrion membrane":
                resfinal=11
            elif res =="nucleus":
                resfinal=12
            elif res =="nucleus membrane":
                resfinal=13
            elif res =="outer membrane":
                resfinal=14
            elif res =="periplasmic space":
                resfinal=15
            elif res =="peroxisome":
                resfinal=16
            elif res =="peroxisome membrane":
                resfinal=17
            elif res =="plasma membrane":
                resfinal=18
            elif res =="plastid":
                resfinal=19
            elif res =="vacuole":
                resfinal=20  
            elif res =="vacuole membrane":
                resfinal=21
            elif res =="secreted":
                resfinal=22
            elif res =="cytoplasm":
                resfinal=23
            
            feature=np.array((resfinal))
    except:
        feature=None
        print("Error on LocTree3 feature")
    
    feature=feature.flatten()
    return feature

def create_transporter_related_pfam_domains(seq, Transport_Pfam_domains):
    seq=seq
    numPfamDomains=0
    try:
        job_id=Pfam_API.run(seq)   
        time.sleep(25)     
        AllPfams=Pfam_API.result(job_id)
        for PfamId in AllPfams:
            if PfamId in Transport_Pfam_domains:
                numPfamDomains+=1                  
        feature=np.array((numPfamDomains))
        feature=feature.flatten()
    except:
        feature=None
        print("Error on Pfam feature")
            
    return feature

def create_single_pfam_features(Transport_Pfam_domains,seq):
    seq=seq
    job_id=Pfam_API.run(seq)   
    time.sleep(25)     
    AllPfams=Pfam_API.result(job_id)
    
    t=True
    for Pfam in Transport_Pfam_domains:
        if t:
            FEATURES=np.array((str(Pfam)))
            t=False
        else:
            FEATURES=np.hstack((FEATURES,str(Pfam)))       
    
    indexx= range(len(FEATURES))
    tam=len(FEATURES)
    FEATURES=np.vstack((FEATURES,np.zeros(tam)))
    for column in indexx:
        feature=FEATURES[0][column]
        for item in AllPfams:
            if item == feature:
                FEATURES[1][column]=1
    
    finalFeature=np.delete(FEATURES, 0, axis=0)
    finalFeature1D=finalFeature.flatten()
    
    return (finalFeature1D)


def Pfam_domains():

    workbook= xlrd.open_workbook("../Data/Pfam_domains.xlsx")
    worksheet= workbook.sheet_by_name("Pfam_domains")
    Pfam_domains=[]
    for row in range(worksheet.nrows):
        EntryText=str(worksheet.row(row))
        Entry=EntryText[7:(len(EntryText)-2)]
        Pfam_domains.append(Entry)
            
    return(Pfam_domains)

def create_features_dataset(seq,domain, Transport_Pfam_domains):
    aminoacid_composition=create_aminoacid_composition(seq)
    aminoacid_occurrence=create_aminoacid_occurrence(seq)
    aminoacid_physico_chemical_composition=create_aminoacid_physico_chemical_composition(seq)
    aminoacid_physico_chemical_occurrence=create_aminoacid_physico_chemical_occurrence(seq)
    Dipeptide_composition=create_dipeptide_composition(seq)
    Transporter_Pfam_domains=create_transporter_related_pfam_domains(seq, Transport_Pfam_domains)
    single_Transporter_Pfam_domains=create_single_pfam_features(Transport_Pfam_domains, seq)
    PhobiusJobId=create_PhobiusJobIds(seq)
    Ahelices_Signalpeptide=create_num_alphahelices_signalpeptide( PhobiusJobId)
    BOMPJobId=create_BOMPJobIds(seq)
    Bbarrels=create_betabarrels(BOMPJobId)
    reqid,Loctree3JobId=create_LocTree3JobIds(seq, domain)
    LocPred=create_location_prediction(Loctree3JobId, reqid)
     
    features=np.hstack((aminoacid_composition,aminoacid_occurrence))
    features=np.hstack((features,aminoacid_physico_chemical_composition))
    features=np.hstack((features,aminoacid_physico_chemical_occurrence))
    features=np.hstack((features,Dipeptide_composition))
    features=np.hstack((features,Transporter_Pfam_domains))
    features=np.hstack((features,single_Transporter_Pfam_domains))
    features=np.hstack((features,Ahelices_Signalpeptide))
    features=np.hstack((features,Bbarrels))
    features=np.hstack((features,LocPred))
    
    return features

def load_models():
    directory1="../Data/BestModels/Using_Dataset1mixed/bestHardVote_model.pkl"
    learned_HVM=joblib.load(directory1)
    return learned_HVM


if __name__ == '__main__':
    ProteinSeq=input("Insert the protein sequence to be predicted:")
    domain= input("Insert the protein Domain to be predicted: 1-euka, 2-bact, 3-arch, 4-other")
    
    Transport_Pfam_domains=Pfam_domains()
    
    featureDataset=create_features_dataset(ProteinSeq, domain,Transport_Pfam_domains)
    
    print(featureDataset.shape)
    
    best_model=load_models()
    
    prediction=best_model.predict(featureDataset)
    print(prediction)
    if prediction=="1":
        learned_HVM=joblib.load("../Data/BestModels/Using_TransporterDataset/bestHardVote_model.pkl")
        prediction2=learned_HVM.predict(featureDataset)
        print(prediction2)
        
    