# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

from Bio import SeqIO
import re
#import os.path

"""
The SeqRecord class offers the following information as attributes:

.seq
– The sequence itself, typically a Seq object.
.id
– The primary ID used to identify the sequence – a string. In most cases this is something like an accession number.
.name
– A “common” name/id for the sequence – a string. In some cases this will be the same as the accession number, but it could also be a clone name. I think of this as being analogous to the LOCUS id in a GenBank record.
.description
– A human readable description or expressive name for the sequence – a string.
.letter_annotations
– Holds per-letter-annotations using a (restricted) dictionary of additional information about the letters in the sequence. The keys are the name of the information, and the information is contained in the value as a Python sequence (i.e. a list, tuple or string) with the same length as the sequence itself. This is often used for quality scores (e.g. Section 19.1.6) or secondary structure information (e.g. from Stockholm/PFAM alignment files).
.annotations
– A dictionary of additional information about the sequence. The keys are the name of the information, and the information is contained in the value. This allows the addition of more “unstructured” information to the sequence.
.features
– A list of SeqFeature objects with more structured information about the features on a sequence (e.g. position of genes on a genome, or domains on a protein sequence). The structure of sequence features is described below in Section 4.3.
.dbxrefs
- A list of database cross-references as strings.
"""

def reverse(text):
    """
    Returns the text backwards
    ex: input = text; output = txet
    """
    return (text[::-1])
        
def TC_ID(records):
    """
    Iterates through every record in the records file 
    and returns a list with each record TC_ID
    """
    TC_IDs=[]
    for record in records:
        ID=str(record.id)
        rev=""
        i=len(ID)-1
        while ID[i] != "|":
            rev+= ID[i]
            i-=1
        TC_IDs.append(reverse(rev))
    
    return (TC_IDs)


#verificar os 8 primeiros caracteres no TC_ID

def Divide_types(TC_IDs,records):
    """
    Searches for the 8 first characters from the TC_IDs 
    and creates a fasta file containing all the fastas with the same 8 characters
    grouping the fastas acording to the TCDB first 4 Transport classifications 
    """
    i=0
    while i < len(records):
        ID=re.findall('.{8}',TC_IDs[i])
        #print(ID)
        #print(i)
        
        file=open("../Data/Positive_cases/"+str(ID)+".fasta","a")
        file.write(str(records[i].format("fasta")))
        file.close
        i+=1
            
            
if __name__ == '__main__':
    
    handle = open("../Data/tcdb.fasta", "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    #print(len(records))
    #records contem 13788 proteinas transportadoras
    handle.close()
    
               
    #print (records[0].id)  #first record
    #print (records[0].seq)
    #print (records[0].name)
    #print (records[0].description)
    #print (records[0].letter_annotations)
    #print (records[0].annotations)
    #print (records[0].features)
    #print (records[0].dbxrefs)
    
    TC_IDs= TC_ID(records)
    #print(TC_ID(records))
    
    #m=re.search('.{8}',TC_IDs[0])
    #m=re.findall('.{8}',TC_IDs[0])
    #print(m[0])
    
    Divide_types(TC_IDs,records)
    
    
    
    
    