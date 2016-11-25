# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import random
import xlrd
from Bio import ExPASy
from Bio import SwissProt
import urllib
from Bio import SeqIO
import os

def Chose_rand():
    """
    Returns a list of 13788 randomly chosen negative cases from the 466621 total negative cases 
    """
    total_list=list(range(1,467681))
    select=13788
    random_selected= random.sample(total_list,select)
    return (random_selected)
    
def Get_Uniprot_IDs(Neg_cases):
    """
    Reads the excel file containg all the negative cases and returns the uniprot ids from the
    randomly chosen 
    """
    workbook= xlrd.open_workbook("../Data/Negative_cases/Uniprot_total_neg_cases.xlsx")
    worksheet= workbook.sheet_by_name("Uniprot_negative_cases")
    UniprotEntries=[]
    for i in Neg_cases:
        EntryText=str(worksheet.cell(i,0))
        Entry=EntryText[6:(len(EntryText)-1)]
        UniprotEntries.append(Entry)
    return UniprotEntries

 
def Save_Fastas(UniprotIDs):
    """
    Uses ExPASy go access the protein fastas using the acession numbers
    Saves a fasta file with all the 13788 randomly chosen
    negative cases
    Incomplete
    """
 
    records = []
     
     
    for accession in UniprotIDs:
        handle = ExPASy.get_sprot_raw(accession)
        record = SwissProt.read(handle)
        records.append(record)
    return records

def Save_Fastas2(UniprotIDs):
    """
    Enters the Swissprot database url using the acession numbers of the selected
    proteins, and writes the fastas in a file
    """
    file=open("../Data/Negative_cases/negative_cases.fasta","w")
    for ID in UniprotIDs:
        data=urllib.request.urlopen("http://www.uniprot.org/uniprot/%s.fasta" %ID)
        f=data.readlines()
        for lines in f:
            file.write(str(lines))
        #help(data)
    file.close()

def Save_Fastas3(Neg_cases,NumDataset):
    """
    Reads a file containing all the negative cases and writes the fastas
    of those in Neg_cases
    """
    handle = open("../Data/Negative_cases/Uniprot_total_neg_cases_pfam_dom.fasta", "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    
    directory="../Data/Datasets/Dataset"+str(NumDataset)
    if not os.path.exists(directory):
        os.makedirs(directory)
    file=open("../Data/Datasets/Dataset"+str(NumDataset)+"/negative_cases.fasta","w")
    for i in Neg_cases:
        file.write(str(records[i].format("fasta")))
    file.close
        
       
if __name__ == '__main__':
    NumDataset=input("Qual o número do dataset?")
    Neg_cases=Chose_rand()
    #UniprotIDs=Get_Uniprot_IDs(Neg_cases)
    Save_Fastas3(Neg_cases,NumDataset)
    
    
    