# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import pickle
import xlrd
import numpy as np
import pandas
import os
from Bio import SeqIO
from Bio import Entrez
import re
import requests
from bs4 import BeautifulSoup
import time
import PhobiusAPI as P_API
import BOMP_API as B_API
import LocTreeAPI as L_API

    

def count_aminoacid(aminoacid,sequence):
    """
    Returns the number of a given aminoacid in a sequence
    """
    sequence=sequence
    count=sequence.count(aminoacid)
    return count

def create_aminoacid_occurrence(data): #utilização desta feature em: "Functional discrimination of membrane proteins using machine learning techniques"
    """
    Returns an array with the aminoacid occurence of each example in the dataset
    """
    data=data
    FEATURES=np.array(("Alanine-A","Arginine-R","Asparagine-N","Aspartiv Acid-D","Cysteine-C","Glumanine-Q",
                       "Glutamic acid-E","Glycine-G","Histidine-H","Isoleucine-I","Leucine-L","Lysine-K",
                       "Methionine-M","Phenylalanine-F","Proline-P","Serine-S","Threonine-T","Tryptophan-W",
                       "Tyrosine-Y","Valine-V"))
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
    index = range(len(data))
    for row in index:
        sequence=data.values[row][3]
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
        FEATURES=np.vstack((FEATURES,feature))
    return FEATURES

    
def create_aminoacid_composition(data): #utilização desta feature em:"Functional discrimination of membrane proteins using machine learning techniques"
    """
    Returns an array with the aminoacid composition of each example in the dataset
    """
    
    data=data
    FEATURES=np.array(("Alanine-A","Arginine-R","Asparagine-N","Aspartiv Acid-D","Cysteine-C","Glumanine-Q",
                       "Glutamic acid-E","Glycine-G","Histidine-H","Isoleucine-I","Leucine-L","Lysine-K",
                       "Methionine-M","Phenylalanine-F","Proline-P","Serine-S","Threonine-T","Tryptophan-W",
                       "Tyrosine-Y","Valine-V"))
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
    index = range(len(data))
    for row in index:
        sequence=data.values[row][3]
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
        FEATURES=np.vstack((FEATURES,feature))
    return FEATURES



def create_aminoacid_physico_chemical_occurrence(data):
    """
    Returns an array with the aminoacid composition based on the physico-chemical properties of each example in the dataset
    """
    
    data=data
    FEATURES=np.array(("Charged(DEKHR)","Aliphatic(ILV)","Aromatic(FHWY)","Polar(DERKQN)","Neutral(AGHPSTY)","Hydrophobic(CFILMVW)","+charged(KRH)","-charged(DE)","Tiny(ACDGST)","Small(EHILKMNPQV)","Large(FRWY)"))
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
    index = range(len(data))
    for row in index:
        sequence=data.values[row][3]
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
        FEATURES=np.vstack((FEATURES,feature))
    return FEATURES

def create_aminoacid_physico_chemical_composition(data): 
    """
    Returns an array with the aminoacid composition based on the physico-chemical properties of each example in the dataset
    """
    
    data=data
    FEATURES=np.array(("Charged(DEKHR)","Aliphatic(ILV)","Aromatic(FHWY)","Polar(DERKQN)","Neutral(AGHPSTY)","Hydrophobic(CFILMVW)","+charged(KRH)","-charged(DE)","Tiny(ACDGST)","Small(EHILKMNPQV)","Large(FRWY)"))
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    #            0   1   2   3   4   5   6   7   8   9   10  11  12  13  14 15   16  17 18   19
    index = range(len(data))
    for row in index:
        sequence=data.values[row][3]
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
        FEATURES=np.vstack((FEATURES,feature))
    return FEATURES


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

def create_dipeptide_composition(data): 
    """
    Returns an array with the dipeptide composition of each example in the dataset
    """
    
    data=data
    
    aminoacids=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
    Dipeptides=[]
    for item in aminoacids:
        pep1=item
        for item in aminoacids:
            pep2=str(pep1)+str(item)
            Dipeptides.append(pep2)
    Dipeptidestuple=tuple(Dipeptides)
    
    FEATURES=np.array(Dipeptidestuple)
    
    index = range(len(data))
    for row in index:
        print("Creating Dipeptide composition. Row:%s"  % row)
        sequence=data.values[row][3]
        feature=np.array(dipeptide_composition(sequence))
        FEATURES=np.vstack((FEATURES,feature))
    return FEATURES



def create_PhobiusJobIds(data): 
    data=data
    AllPhobiusJobIds=[]
    index = range(len(data))
    for row in index:
        try:
            print (row)
            if row % 500 == 0:
                time.sleep(180)
            seq=data.values[row][3]
            PhobiusId= (P_API.run(email="danielvarzim@hotmail.com",
                       format="short",
                       stype="protein",
                       sequence=seq))
            PhobiusId=str(PhobiusId)
            PhobiusId=PhobiusId.strip("b'")
            print(PhobiusId)
            AllPhobiusJobIds.append(PhobiusId)        
        except:
            i=True
            while i:
                try:
                    print (row)
                    if row % 500 == 0:
                        time.sleep(180)
                    seq=data.values[row][3]
                    PhobiusId= (P_API.run(email="danielvarzim@hotmail.com",
                                    format="short",
                                    stype="protein",
                                    sequence=seq))
                    PhobiusId=str(PhobiusId)
                    PhobiusId=PhobiusId.strip("b'")
                    print(PhobiusId)
                    AllPhobiusJobIds.append(PhobiusId)
                    i=False
                except:
                    print("Error on row %s" % row )    
            
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/PhobiusJobIds.txt", 'wb') as f:
        pickle.dump(AllPhobiusJobIds, f)
    
def create_num_alphahelices_signalpeptide(data,PhobiusJobIds): #utilizacao desta feature em: "TransportTP- A two-phase classification approach for membrane transporter prediction and characterization"
    data=data
    FEATURES=np.array(("Number of alpha helices","Is a signal peptide present"))
    index = range(len(data))
    for row in index:
        print(row)
        try:
            res=P_API.result(job_id=PhobiusJobIds[row],result_type="out")
             
            AlphaHelices=str(res["TM"])
            numberAlphaHelices=int(AlphaHelices.strip("b'"))
             
            SignalPeptide=str(res["SP"])
            if SignalPeptide.strip("b'") == "Y" or "found":
                HaveSignalPeptide=1
            elif SignalPeptide.strip("b'") == "0":
                HaveSignalPeptide=0
            else:
                HaveSignalPeptide=0
                print("Error on row :%s" % row)
                 
            feature=np.array((numberAlphaHelices,HaveSignalPeptide))
         
            FEATURES=np.vstack((FEATURES,feature))
        except:

            try:
                res=P_API.result(job_id=PhobiusJobIds[row],result_type="out")
                AlphaHelices=str(res["TM"])
                numberAlphaHelices=int(AlphaHelices.strip("b'"))
                        
                SignalPeptide=str(res["SP"])
                if SignalPeptide.strip("b'") == "Y" or "found":
                    HaveSignalPeptide=1
                elif SignalPeptide.strip("b'") == "0":
                    HaveSignalPeptide=0
                else:
                    pass
            except:
                print("Error on row %s: result error" % row)
            
                try:
                    res=P_API.result(job_id=PhobiusJobIds[row],result_type="out")
                    AlphaHelices=str(res["TM"])
                    numberAlphaHelices=int(AlphaHelices.strip("b'"))
                            
                    SignalPeptide=str(res["SP"])
                    if SignalPeptide.strip("b'") == "Y" or "found":
                        HaveSignalPeptide=1
                    elif SignalPeptide.strip("b'") == "0":
                        HaveSignalPeptide=0
                    else:
                        pass
                            
                except:
                    numberAlphaHelices=0
                    SignalPeptide=0
                    print("Error on row %s: values =0" % row)
                         
            feature=np.array((numberAlphaHelices,HaveSignalPeptide))
                 
            FEATURES=np.vstack((FEATURES,feature))            
                   
    return FEATURES

def create_BOMPJobIds(data):
    data=data
    AllBOMPJobIds=[]
    index = range(len(data))
    for row in index:
        try:
            print (row)
            if row % 500 == 0:
                time.sleep(180)
            seq=data.values[row][3]
            seq=">SeqName\n"+seq
            BOMPId= (B_API.run(seqs=seq))
            
            print(BOMPId)
            AllBOMPJobIds.append(BOMPId)        
        except:
            i=True
            while i:
                try:
                    print (row)
                    if row % 500 == 0:
                        time.sleep(180)
                    seq=data.values[row][3]
                    seq=">SeqName\n"+seq
                    BOMPId= (B_API.run(seqs=seq))
            
                    print(BOMPId)
                    AllBOMPJobIds.append(BOMPId) 
                    i=False
                except:
                    print("Error on row %s" % row )    
            
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/BOMPJobIds.txt", 'wb') as f:
        pickle.dump(AllBOMPJobIds, f)


def create_betabarrels(data,BOMPJobIds): 
    data=data
    FEATURES=np.array(("Number of Beta-barrels"))
    index = range(len(data))
    for row in index:
        print(row)
        print(BOMPJobIds[row])
        try:
            res=B_API.result(job_id=BOMPJobIds[row])              
            feature=np.array((res))
        except:
            try:
                res=B_API.result(job_id=BOMPJobIds[row])
                feature=np.array((res))
            except:
                res=0
                feature=np.array((res))
                print("Error on row %s: result error" % row)

        FEATURES=np.vstack((FEATURES,feature))            
                   
    return FEATURES


def create_LocTree3JobIds(data):
    Entrez.email = "danielvarzim@hotmail.com"
    data=data
    AllLocTree3JobIds=[]
    AllLocTree3ReqIds=[]
    index = range(len(data))
    for row in index:
        try:
            print (row)
            if row % 500 == 0:
                time.sleep(180)
            fulldom=data.values[row][6]
            print(fulldom)
            if fulldom =="Eukaryota":
                dom ="euka"
            elif fulldom =="Bacteria":
                dom="bact"
            elif fulldom =="Archaea":
                dom="arch"
            else:
                dom="arch"
            print(dom)
            seq=data.values[row][3]
            seq=">SeqName\n"+seq
            reqid,LocTree3Id=(L_API.run(domain=dom, email="danielvarzim@hotmail.com",sequence=seq))
            print(reqid)
            print(LocTree3Id)
            print("----------------------------------")
            AllLocTree3JobIds.append(LocTree3Id)  
            AllLocTree3ReqIds.append(reqid)      
        except:
            i=True
            while i:
                try:
                    print (row)
                    if row % 500 == 0:
                        time.sleep(180)
                    fulldom=data.values[row][6]
                    if fulldom =="Eukaryota":
                        dom ="euka"
                    elif fulldom =="Bacteria":
                        dom="bact"
                    elif fulldom =="Archaea":
                        dom="arch"
                    else:
                        dom="arch"
                    print(dom)
                    seq=data.values[row][3]
                    seq=">SeqName\n"+seq
                    reqid,LocTree3Id=(L_API.run(domain=dom, email="danielvarzim@hotmail.com",sequence=seq))
                    print(reqid)
                    print(LocTree3Id)
                    print("----------------------------------")
                    AllLocTree3JobIds.append(LocTree3Id)  
                    AllLocTree3ReqIds.append(reqid)  
                    i=False
                except:
                    print("Error on row %s" % row )    
            
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Loctree3JobIds.txt", 'wb') as f:
        pickle.dump(AllLocTree3JobIds, f)
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Loctree3ReqIds.txt", 'wb') as g:
        pickle.dump(AllLocTree3ReqIds, g)


def create_location_prediction(data, Loctree3JobIds, Loctree3ReqIds):
    data=data
    FEATURES=np.array(("Protein Subcellular Location"))
    index = range(len(data))
    #namelist=[]
    for row in index:
        print(row)
        print(Loctree3JobIds[row])
        try:
            #L_API.LoadingFinished(Loctree3ReqIds[row])
            res=L_API.result(result_url=Loctree3JobIds[row])
            res=res.lower()
            print(res)
            #namelist.append(res)    
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
            elif res =="golgi apparatus":
                resfinal=8
            elif res =="golgi apparatus membrane":
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
            elif res =="inner membrane":
                resfinal=24    
            else:
                print("Error in the location name")          
            print(resfinal)
            
            feature=np.array((resfinal))
        except:

            resfinal=0
            feature=np.array((resfinal))
            print("Error on row %s: result error" % row)
        
        FEATURES=np.vstack((FEATURES,feature))
    
    #===========================================================================
    # names=list(set(namelist))
    # file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/Loctree3Names.txt", 'w+')
    # for name in names:
    #     file.write(name+"\n")
    # file.close()
    #===========================================================================
    
    return FEATURES


def create_transporter_related_pfam_domains(data, Transport_Pfam_domains):#utilizacao desta feature em: "TransportTP- A two-phase classification approach for membrane transporter prediction and characterization"
    Entrez.email = "danielvarzim@hotmail.com"
    AllPfams=[]
    AllTransporterPfams=[]   
    FEATURES=np.array(("Number of Transporter related Pfam domains"))
    index= range(len(data))#(27576)
    for row in index:
        ProteinPfams=[]
        if row % 500 == 0:
            time.sleep(180)
        numPfamDomains=0
        Uniprot_Acc = data.values[row][2]
        try:
            try:
                handle = Entrez.efetch(db="nucleotide", id=Uniprot_Acc, rettype="gb", retmode="text")
                record = SeqIO.read(handle, format="genbank")
                dic=record.annotations
                db_source = (dic.get("db_source"))
                Pfams= re.findall("(Pfam:.+?),",db_source)
                for Pfam in Pfams:
                    PF=Pfam[5:]
                    AllPfams.append(PF)
                    if PF in Transport_Pfam_domains:
                        numPfamDomains+=1
                        ProteinPfams.append(PF)
                                                    
                #print(ProteinPfams)        
                #print(Pfams)
                #print(numPfamDomains)
                
            except:
                
                url="http://www.uniprot.org/uniprot/"+str(Uniprot_Acc)
                source_code=requests.get(url)
                plain_text=source_code.text
                soup=BeautifulSoup(plain_text)
                links=""
                for link in soup.findAll("a"):
                    links+=(str(link.get("href"))+",")
                PfamLinks= re.findall("(http://pfam.xfam.org/family/.{7})",links)
                UniquePfamLinks=list(set(PfamLinks))  
                #print(UniquePfamLinks)
                for PfamLink in UniquePfamLinks:
                    PfamId=PfamLink[-7:]
                    AllPfams.append(PfamId)
                    if PfamId in Transport_Pfam_domains:
                        numPfamDomains+=1
                        ProteinPfams.append(PfamId) 
                        
                #print(ProteinPfams)        
                #print(numPfamDomains)
                
        except:
            print("Error in row %s" %(row))
        
        AllTransporterPfams.append(ProteinPfams)           
        feature=np.array((numPfamDomains))
        FEATURES=np.vstack((FEATURES,feature))
    
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Protein_Pfams.txt", 'wb') as f:
        pickle.dump(AllTransporterPfams, f)
    
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Protein_All_Pfams.txt", 'wb') as g:
        pickle.dump(AllPfams, g)
        
    print(AllTransporterPfams)  
    print(FEATURES)
    return FEATURES


def create_transporter_related_go_terms(data, Transport_GO_terms):#utilização desta feature em: "TransportTP- A two-phase classification approach for membrane transporter prediction and characterization"
    """
    https://www.ebi.ac.uk/QuickGO/GSearch?q=transport&what=GO&limit=1000
    site to get the transport related GO terms, filter: transport (1000 GO terms) saved in data/Transport_GO_terms.xlsx
    some proteins can´t be fetched by Entrez--> web crawler that scans Uniprot
    """
    Entrez.email = "danielvarzim@hotmail.com"
    FEATURES=np.array(("Number of Transporter related GO terms"))
    index= range(len(data))#(27576)
    for row in index:
        print (row)
        if row % 500 == 0:
            time.sleep(180)
        numGoTerms=0
        Uniprot_Asc = data.values[row][2]
        try:
            try:
                handle = Entrez.efetch(db="nucleotide", id=Uniprot_Asc, rettype="gb", retmode="text")
                record = SeqIO.read(handle, format="genbank")
                dic=record.annotations
                #print(dic.keys())
                db_source = (dic.get("db_source"))
                GOs= re.findall("(GO:.+?),",db_source)
                for go in GOs:
                    if go in Transport_GO_terms:
                        numGoTerms+=1
                #print(GOs)
                #print(numGoTerms)
            except:
                url="http://www.uniprot.org/uniprot/"+str(Uniprot_Asc)
                source_code=requests.get(url)
                plain_text=source_code.text
                soup=BeautifulSoup(plain_text)
                links=""
                for link in soup.findAll("a"):
                    links+=(str(link.get("href"))+",")
                GOLinks= re.findall("(http://www.ebi.ac.uk/QuickGO/GTerm.id=GO:.{7})",links)#GTerm?id alterado para GTerm.id devido a conflito com as expressões regulares
                UniqueGOLinks=list(set(GOLinks))  
                #print(UniqueGOLinks)
                for GOLink in UniqueGOLinks:
                    GOId=GOLink[-10:]
                    if GOId in Transport_GO_terms:
                        numGoTerms+=1
                #print(numGoTerms)
        except:
            print("Error in row %s" %(row))
        
        feature=np.array((numGoTerms))
        FEATURES=np.vstack((FEATURES,feature))
    print(FEATURES)
    return FEATURES

def create_single_pfam_features(Transport_Pfam_domains,data,ProteinPfams):
    t=True
    for Pfam in Transport_Pfam_domains:
        if t:
            FEATURES=np.array((str(Pfam)))
            t=False
        else:
            FEATURES=np.hstack((FEATURES,str(Pfam)))       
    
    indexy= range(len(data))#(27576)
    indexx= range(len(FEATURES))
    tam=len(FEATURES)
    for row in indexy:
        print ("ROW: %s" % row)
        FEATURES=np.vstack((FEATURES,np.zeros(tam)))
        for column in indexx:
            #print("Column: %s" % column)
            feature=FEATURES[0][column]
            #print(feature)
            for item in ProteinPfams[row]:
                #print("Item: %s" % item)
                if item == feature:
                    FEATURES[row][column]=1


    print("Features:\n %s" % FEATURES)
    return (FEATURES)



def Go_terms():
    workbook= xlrd.open_workbook("../Data/Transport_GO_terms.xlsx")
    worksheet= workbook.sheet_by_name("Transport_GO_terms")
    Transport_GO_terms=[]
    for row in range(worksheet.nrows):
        EntryText=str(worksheet.row(row))
        Entry=EntryText[11:(len(EntryText)-2)]
        Transport_GO_terms.append(Entry)
        
    return Transport_GO_terms

def Pfam_domains():

    workbook= xlrd.open_workbook("../Data/Pfam_domains.xlsx")
    worksheet= workbook.sheet_by_name("Pfam_domains")
    Pfam_domains=[]
    for row in range(worksheet.nrows):
        EntryText=str(worksheet.row(row))
        Entry=EntryText[7:(len(EntryText)-2)]
        Pfam_domains.append(Entry)
            
    return(Pfam_domains)

#===============================================================================
# def create_features_dataset(features):#Verificar esta função quando adicionar features de arrays de diferente numero de colunas
#     """
#     Inserts all the features created into a dataset
#     Returns the final dataset
#     """
#     i=True
#     for feature in features:
#         if i:
#             finalDataset=np.array((feature))
#             i=False
#         else:
#             finalDataset=np.hstack((finalDataset,feature))
#     return finalDataset
#===============================================================================




if __name__ == '__main__': 
    NumDataset=input("Qual o nome do dataset?")
    directory="../Data/Datasets/Dataset"+str(NumDataset)+"/Features"
    if not os.path.exists(directory):
        os.makedirs(directory)
         
    data=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/DataDoms.csv",sep=",")
     
        
###########################################################################    
    amino_composition=create_aminoacid_composition(data)
    amino_occurrence=create_aminoacid_occurrence(data)
                        
    CSV_amino_composition=pandas.DataFrame(amino_composition[1:], index=None, columns=amino_composition[0])
    CSV_amino_composition.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_composition.csv", sep=',')    
    CSV_amino_occurrence=pandas.DataFrame(amino_occurrence[1:], index=None, columns=amino_occurrence[0])
    CSV_amino_occurrence.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_occurrence.csv", sep=',')    
                       
    feature1 = pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_composition.csv")
    ft1=np.array(feature1)
    ft1=np.delete(ft1,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_composition.csv",ft1,delimiter=",")
                   
    feature2=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_occurrence.csv")
    ft2=np.array(feature2)
    ft2=np.delete(ft2,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_occurrence.csv",ft2,delimiter=",")
                 
################################################################################ 
    
    amino_physico_chemical_composition=create_aminoacid_physico_chemical_composition(data)
    amino_physico_chemical_occurrence=create_aminoacid_physico_chemical_occurrence(data)
                         
    CSV_amino_physico_chemical_composition=pandas.DataFrame(amino_physico_chemical_composition[1:], index=None, columns=amino_physico_chemical_composition[0])
    CSV_amino_physico_chemical_composition.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_physico_chemical_composition.csv", sep=',')    
    CSV_amino_physico_chemical_occurrence=pandas.DataFrame(amino_physico_chemical_occurrence[1:], index=None, columns=amino_physico_chemical_occurrence[0])
    CSV_amino_physico_chemical_occurrence.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_physico_chemical_occurrence.csv", sep=',')    
                        
    feature9 = pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_physico_chemical_composition.csv")
    ft9=np.array(feature9)
    ft9=np.delete(ft9,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_physico_chemical_composition.csv",ft9,delimiter=",")
                   
    feature10=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_aminoacid_physico_chemical_occurrence.csv")
    ft10=np.array(feature10)
    ft10=np.delete(ft10,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_physico_chemical_occurrence.csv",ft10,delimiter=",")
           
##################################################################################
         
    dip_composition=create_dipeptide_composition(data)
              
    CSV_dip_composition=pandas.DataFrame(dip_composition[1:],index=None,columns=dip_composition[0])
    CSV_dip_composition.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_dipeptide_composition.csv", sep=',')
             
    feature11=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_dipeptide_composition.csv")
    ft11=np.array(feature11)
    ft11=np.delete(ft11,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/dipeptide_composition.csv",ft11,delimiter=",")
             
   
###############################################################################      
    Transport_Pfam_domains=Pfam_domains()
                      
    PfamDomains=create_transporter_related_pfam_domains(data, Transport_Pfam_domains)
                    
    CSV_Pfam_domains=pandas.DataFrame(PfamDomains[1:], index=None, columns=PfamDomains[0])
    CSV_Pfam_domains.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Pfam_domains.csv",sep=",")
                   
    feature4=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Pfam_domains.csv")
    ft4=np.array(feature4)
    ft4=np.delete(ft4,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Transporter_Pfam_domains.csv",ft4,delimiter=",")
                    
################################################################################    
    Transport_Pfam_domains=Pfam_domains()
             
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Protein_Pfams.txt", 'rb') as f:
        ProteinPfams = pickle.load(f)
                   
    single_PfamDomains=create_single_pfam_features(Transport_Pfam_domains, data, ProteinPfams)
                   
    CSV_singlePfam_domains=pandas.DataFrame(single_PfamDomains[1:], index=None, columns=single_PfamDomains[0])
    CSV_singlePfam_domains.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_single_Pfam_domains.csv",sep=",")
                  
    feature5=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_single_Pfam_domains.csv")
    ft5=np.array(feature5)
    ft5=np.delete(ft5,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/single_Transporter_Pfam_domains.csv",ft5,delimiter=",")
                   
#################################################################################        
                  
    create_PhobiusJobIds(data)
                   
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/PhobiusJobIds.txt", 'rb') as f:
        PhobiusJobIds= pickle.load(f)
                    
    Alphahelices_Signalpeptide=create_num_alphahelices_signalpeptide(data,PhobiusJobIds)
                    
    CSV_Alphahelices_Signalpeptide=pandas.DataFrame(Alphahelices_Signalpeptide[1:],index=None, columns=Alphahelices_Signalpeptide[0])
    CSV_Alphahelices_Signalpeptide.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Alphahelices_Signalpeptide.csv",sep=",")
                   
    feature6=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Alphahelices_Signalpeptide.csv")
    ft6=np.array(feature6)
    ft6=np.delete(ft6,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Alphahelices_Signalpeptide.csv", ft6,delimiter=",")
                  
##################################################################################

    create_BOMPJobIds(data)
            
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/BOMPJobIds.txt", 'rb') as f:
        BOMPJobIds= pickle.load(f)
                  
    BetaBarrels=create_betabarrels(data, BOMPJobIds)
                  
    CSV_BetaBarrels=pandas.DataFrame(BetaBarrels[1:],index=None, columns=BetaBarrels[0])
    CSV_BetaBarrels.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_BetaBarrels.csv",sep=",")
                  
    feature7=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_BetaBarrels.csv")
    ft7=np.array(feature7)
    ft7=np.delete(ft7,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/BetaBarrels.csv", ft7, delimiter=",")
              
################################################################################
    create_LocTree3JobIds(data)
             
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Loctree3JobIds.txt", 'rb') as f:
        Loctree3JobIds= pickle.load(f)
    with open("../Data/Datasets/Dataset"+str(NumDataset)+"/Loctree3ReqIds.txt", 'rb') as g:
        Loctree3ReqIds= pickle.load(g)
               
    LocationPrediction=create_location_prediction(data, Loctree3JobIds,Loctree3ReqIds)
              
    CSV_LocationPred=pandas.DataFrame(LocationPrediction[1:],index=None, columns=LocationPrediction[0])
    CSV_LocationPred.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Location_Prediction.csv",sep=",")
               
    feature8=pandas.read_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Feature_Location_Prediction.csv")
    ft8=np.array(feature8)
    ft8=np.delete(ft8,0,1)
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Location_Prediction.csv", ft8, delimiter=",")
                
      
 
 
    aminoacid_composition=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_composition.csv",delimiter=",")
    aminoacid_composition.astype(np.float64)
    aminoacid_occurrence=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_occurrence.csv",delimiter=",")
    aminoacid_occurrence.astype(np.int64)
    aminoacid_physico_chemical_composition=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_physico_chemical_composition.csv",delimiter=",")
    aminoacid_physico_chemical_composition.astype(np.float64)
    aminoacid_physico_chemical_occurrence=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/aminoacid_physico_chemical_occurrence.csv",delimiter=",")
    aminoacid_physico_chemical_occurrence.astype(np.int64)
    Dipeptide_composition=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/dipeptide_composition.csv",delimiter=",")
    Dipeptide_composition.astype(np.float64)
    Transporter_Pfam_domains=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Transporter_Pfam_domains.csv", delimiter=",")
    Transporter_Pfam_domains.astype(np.int64)
    single_Transporter_Pfam_domains=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/single_Transporter_Pfam_domains.csv", delimiter=",")
    single_Transporter_Pfam_domains.astype(np.int64)
    Ahelices_Signalpeptide=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Alphahelices_Signalpeptide.csv", delimiter=",")
    Ahelices_Signalpeptide.astype(np.int64)
    Bbarrels=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/BetaBarrels.csv",delimiter=",")
    Bbarrels.astype(np.int64)
    LocPred=np.genfromtxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Location_Prediction.csv",delimiter=",")
    LocPred.astype(np.int64)
         
          
    features=np.hstack((aminoacid_composition,aminoacid_occurrence))
    features=np.hstack((features,aminoacid_physico_chemical_composition))
    features=np.hstack((features,aminoacid_physico_chemical_occurrence))
    features=np.hstack((features,Dipeptide_composition))
    features=np.hstack((features,Transporter_Pfam_domains[:,np.newaxis]))
    features=np.hstack((features,single_Transporter_Pfam_domains))
    features=np.hstack((features,Ahelices_Signalpeptide))
    features=np.hstack((features,Bbarrels[:,np.newaxis]))
    features=np.hstack((features,LocPred[:,np.newaxis]))
    features.astype(np.float64)
          
     
    np.savetxt("../Data/Datasets/Dataset"+str(NumDataset)+"/Features/Features_Dataset.csv",features,delimiter=",")    
