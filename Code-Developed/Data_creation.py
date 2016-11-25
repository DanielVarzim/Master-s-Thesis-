# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

import numpy as np
import os
from Bio import Entrez
from Bio import SeqIO
import pandas
import requests
import re
import time

def reverse(text):
    """
    Returns the text backwards
    ex: input = text; output = txet
    """
    return (text[::-1])

def get_fasta_id(record):
    """
    Returns the Id of the fasta record
    """
    return record.id

def get_uniprot_accession(record):
    """
    Returns the Uniprot Acession number of the record
    """
    y=record.id[0:3]
    if y == "gnl": #positive case
        IDs=record.id.split("|")
        ID=IDs[2]
    else: #negative case
        IDs=record.id.split("|")
        ID=IDs[1]
    return ID

def get_sequence(record):
    """
    Returns the record´s sequence
    """
    return str(record.seq)

def check_if_transporter(record):
    """
    Verifies if the record is a transporter protein (returns 1) or a negative case (returns 0)
    """
    y=record.id[0:3]
    if y == "gnl": #positive case
        return "1"
    else: return "0" #negative case

def get_TCDB_id(record):
    y=record.id[0:3]
    if y == "gnl": #positive case
        ID=str(record.id)
        rev=""
        i=len(ID)-1
        while ID[i] != "|":
            rev+= ID[i]
            i-=1
        return str(reverse(rev))
    else: 
        return "0"
    

def get_domain(record):
    try:
        Entrez.email = "danielvarzim@hotmail.com"
        Uniprot_Asc = get_uniprot_accession(record)
        handle = Entrez.efetch(db="nucleotide", id=Uniprot_Asc, rettype="gb", retmode="text")
        record = SeqIO.read(handle, format="genbank")
        dic=record.annotations
        taxonomy=dic.get("taxonomy")
        print(taxonomy[0])
        dom=taxonomy[0]
        return dom
         
    except:
        try:
            session=requests.Session()
            Uniprot_Asc = get_uniprot_accession(record)
            print(Uniprot_Asc)
            url="http://www.uniprot.org/uniprot/"+str(Uniprot_Asc)
            resp=session.get(url)
            #print(resp.text)
            Organism=str(re.findall("setCustomVar\', 2, \"Organism\", \".*?]",resp.text))
            #print(Organism)
            Org=re.findall('\d+',Organism)
            #print(Org[1])
                
            url2="http://www.uniprot.org/taxonomy/"+Org[1]
            session2=requests.Session()
            resp2=session2.get(url2)   
            job_content=resp2.text
            #print(job_content)
            Bact=re.findall("property\=\"rdfs:subClassOf\">Bacteria<", job_content)
            Eukaryota=re.findall("property\=\"rdfs:subClassOf\">Eukaryota<", job_content)
            Archaea=re.findall("property\=\"rdfs:subClassOf\">Archaea<", job_content)
            if Eukaryota:
                print(Eukaryota)
                dom="Eukaryota"
            elif Bact:
                print(Bact)
                dom="Bacteria"
            elif (Archaea):
                print(Archaea)
                dom="Archaea"
            else:
                print("none of the following: Bacteria, Eukaryota, Archaea")
                dom="Unknown"
                
            print(dom)
            return dom

        except:
            time.sleep(60)
            try:
                Uniprot_Asc = get_uniprot_accession(record)
                print(Uniprot_Asc)
                session=requests.Session()
                url="http://www.uniprot.org/uniprot/"+str(Uniprot_Asc)
                resp=session.get(url)
                #print(resp.text)
                Organism=str(re.findall("setCustomVar\', 2, \"Organism\", \".*?]",resp.text))
                #print(Organism)
                Org=re.findall('\d+',Organism)
                #print(Org[1])                            
                                
                url2="http://www.uniprot.org/taxonomy/"+Org[1]
                session2=requests.Session()
                resp2=session2.get(url2)   
                job_content=resp2.text
                #print(job_content)
                Bact=re.findall("property\=\"rdfs:subClassOf\">Bacteria<", job_content)
                Eukaryota=re.findall("property\=\"rdfs:subClassOf\">Eukaryota<", job_content)
                Archaea=re.findall("property\=\"rdfs:subClassOf\">Archaea<", job_content)
                if Eukaryota:
                    print(Eukaryota)
                    dom="Eukaryota"
                elif Bact:
                    print(Bact)
                    dom="Bacteria"
                elif (Archaea):
                    print(Archaea)
                    dom="Archaea"
                else:
                    print("none of the following: Bacteria, Eukaryota, Archaea")
                    dom="Unknown"
                                
                print(dom)
                return dom
            except:
                print("Error")
                dom="Unknown"
                return dom
    
def build_data(NumDataset):
    """
    creates a csv file named Data.csv containing all the positive and randomly chosen negative cases
    for a total of 27576 proteins (13788 of each)
    indexes are the row indexes
    columns: Fasta ID | Uniprot Acession | Sequence | Is Tranporter? | TCDB ID | Taxonomy Domain
    """
    i=1
    dataset=np.array(('Fasta ID','Uniprot Accession', 'Sequence', 'Is transporter?','TCDB ID','Taxonomy Domain'),dtype=object)
    for filename in os.listdir("../Data/Positive_cases"):
        handle = open("../Data/Positive_cases/"+filename, "rU")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        for record in records:
            data=np.array((get_fasta_id(record),get_uniprot_accession(record),get_sequence(record),check_if_transporter(record),get_TCDB_id(record),get_domain(record)),dtype=object)
            dataset=np.vstack((dataset,data))
            print(i)
            i+=1
    
    directory="../Data/Datasets/Dataset"+str(NumDataset)
    if not os.path.exists(directory):
        os.makedirs(directory)        
    handle2 = open("../Data/Datasets/Dataset"+str(NumDataset)+"/negative_cases.fasta", "rU")
    records2 = list(SeqIO.parse(handle2, "fasta"))
    handle2.close()
    for record in records2:
        data=np.array((get_fasta_id(record),get_uniprot_accession(record),get_sequence(record),check_if_transporter(record),get_TCDB_id(record),get_domain(record)),dtype=object)
        dataset=np.vstack((dataset,data))
        print(i)
        i+=1
        
    df=pandas.DataFrame(dataset[1:], index=None, columns=dataset[0])
    df.to_csv("../Data/Datasets/Dataset"+str(NumDataset)+"/DataDoms.csv", sep=',')    
        

if __name__ == '__main__':
    NumDataset=input("Qual o número do dataset?")
    build_data(NumDataset)

    