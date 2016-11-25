# -*- coding: utf-8 -*-

'''
Created on 14/12/2015

@author: danie
'''

from Bio import SeqIO
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
        id=str(record.id)
        rev=""
        i=len(id)-1
        while id[i] != "|":
            rev+= id[i]
            i-=1
        TC_IDs.append(reverse(rev))
    
    return (TC_IDs)

def First_division(TC_IDs,records):
    """
    Divides the records from the TCDB into 7 different files according
    to their TC_ID
    1:Channels_Pores
    2:Electrochemical_potential_driven_Transporters
    3:Primary_active_transporters
    4:Group_Translocators
    5:Transmembrane_electron_carriers
    8:Accessory_Factors_Involved_in_Transport
    9:Incompletely_Characterized_Transport_Systems
    """
    
    i=0
    file1=open("../Data/First_division/1 Channels_Pores.fasta","w")
    file2=open("../Data/First_division/2 Electrochemical_potential_driven_Transporters.fasta","w")
    file3=open("../Data/First_division/3 Primary_active_transporters.fasta","w")
    file4=open("../Data/First_division/4 Group_Translocators.fasta","w")
    file5=open("../Data/First_division/5 Transmembrane_electron_carriers.fasta","w")
    file8=open("../Data/First_division/8 Accessory_Factors_Involved_in_Transport.fasta","w+")
    file9=open("../Data/First_division/9 Incompletely_Characterized_Transport_Systems.fasta","w+")
    for TC_ID in TC_IDs:
        if TC_ID[0]=="1":
            file1.write(str(records[i].format("fasta")))
            i+=1
        elif TC_ID[0]=="2":
            file2.write(str(records[i].format("fasta")))
        
            i+=1
        elif TC_ID[0]=="3":
            file3.write(str(records[i].format("fasta")))
            
            i+=1
        elif TC_ID[0]=="4":
            file4.write(str(records[i].format("fasta")))
            
            i+=1
        elif TC_ID[0]=="5":
            file5.write(str(records[i].format("fasta")))
            
            i+=1
        elif TC_ID[0]=="8":
            file8.write(str(records[i].format("fasta")))
            
            i+=1
        elif TC_ID[0]=="9":
            file9.write(str(records[i].format("fasta")))
        
            i+=1
        else:
            print("ERROR")
    
    file1.close()
    file2.close()
    file3.close()
    file4.close()
    file5.close()
    file8.close()
    file9.close()

def Second_division(TC_IDs,records):
    """
    Divides the records from the 7different categories into different files according
    to their TC_ID
    
     1.A: α-Type Channels
     1.B: β-Barrel Porins
     1.C: Pore-Forming Toxins (Proteins and Peptides)
     1.D: Non-Ribosomally Synthesized Channels
     1.E: Holins
     1.F: Vesicle Fusion Pores
     1.G: Viral Fusion Pores
     1.H: Paracellular Channels
     1.I: Membrane-bounded Channels
     1.J: Virion Egress Pyramidal Apertures
     1.K: Phage DNA Injection Channels
     1.L: Tunneling Nanotubes, TNTs
     1.M: Membrane Fusion-mediating Spannins
     ------------------------------------------------------
     2.A: Porters (uniporters, symporters, antiporters)
     2.B: Nonribosomally synthesized porters
     2.C: Ion-gradient-driven energizers
     2.D: Transcompartment Lipid Carrier
     ------------------------------------------------------
     3.A: P-P-bond-hydrolysis-driven transporters
     3.B: Decarboxylation-driven transporters
     3.C: Methyltransfer-driven transporters
     3.D: Oxidoreduction-driven transporters
     3.E: Light absorption-driven transporters
     ------------------------------------------------------
     4.A: Phosphotransfer-driven Group Translocators (PTS)
     4.B: Nicotinamide ribonucleoside uptake transporters
     4.C: Acyl CoA ligase-coupled transporters
     4.D: Polysaccharide Synthase/Exporters
     4.E: Vacuolar Polyphosphate Polymerase-catalyzed Group Translocators
     ------------------------------------------------------
     5.A: Transmembrane 2-electron transfer carriers
     5.B: Transmembrane 1-electron transfer carriers
     ------------------------------------------------------
     8.A: Auxiliary transport proteins
     8.B: Ribosomally synthesized protein/peptide toxins/agonists that target channels and carriers
     8.C: Non-ribosomally synthesized toxins that target channels, carriers and pumps
     ------------------------------------------------------
     9.A: Recognized transporters of unknown biochemical mechanism
     9.B: Putative transport proteins
     9.C: Functionally characterized transporters lacking identified sequences
    """
    
    i=0
    file1_A=open("../Data/Second_division/1_A alpha-Type Channels.fasta","w")
    file1_B=open("../Data/Second_division/1_B beta-Barrel Porins.fasta","w")
    file1_C=open("../Data/Second_division/1_C Pore-Forming Toxins (Proteins and Peptides).fasta","w")
    file1_D=open("../Data/Second_division/1_D Non-Ribosomally Synthesized Channels.fasta","w")
    file1_E=open("../Data/Second_division/1_E Holins.fasta","w")
    file1_F=open("../Data/Second_division/1_F Vesicle Fusion Pores.fasta","w")
    file1_G=open("../Data/Second_division/1_G Viral Fusion Pores.fasta","w")
    file1_H=open("../Data/Second_division/1_H Paracellular Channels.fasta","w")
    file1_I=open("../Data/Second_division/1_I Membrane-bounded Channels.fasta","w")
    file1_J=open("../Data/Second_division/1_J Virion Egress Pyramidal Apertures.fasta","w")
    file1_K=open("../Data/Second_division/1_K Phage DNA Injection Channels.fasta","w")
    file1_L=open("../Data/Second_division/1_L Tunneling Nanotubes, TNTs.fasta","w")
    file1_M=open("../Data/Second_division/1_M Membrane Fusion-mediating Spannins.fasta","w")
    
    file2_A=open("../Data/Second_division/2_A Porters (uniporters, symporters, antiporters).fasta","w")
    file2_B=open("../Data/Second_division/2_B Nonribosomally synthesized porters.fasta","w")
    file2_C=open("../Data/Second_division/2_C Ion-gradient-driven energizers.fasta","w")
    file2_D=open("../Data/Second_division/2_D Transcompartment Lipid Carrier.fasta","w")
    
    file3_A=open("../Data/Second_division/3_A P-P-bond-hydrolysis-driven transporters.fasta","w")
    file3_B=open("../Data/Second_division/3_B Decarboxylation-driven transporters.fasta","w")
    file3_C=open("../Data/Second_division/3_C Methyltransfer-driven transporters.fasta","w")
    file3_D=open("../Data/Second_division/3_D Oxidoreduction-driven transporters.fasta","w")
    file3_E=open("../Data/Second_division/3_E Light absorption-driven transporters.fasta","w")
    
    file4_A=open("../Data/Second_division/4_A Phosphotransfer-driven Group Translocators (PTS).fasta","w")
    file4_B=open("../Data/Second_division/4_B Nicotinamide ribonucleoside uptake transporters.fasta","w")
    file4_C=open("../Data/Second_division/4_C Acyl CoA ligase-coupled transporters.fasta","w")
    file4_D=open("../Data/Second_division/4_D Polysaccharide Synthase-Exporters.fasta","w")
    file4_E=open("../Data/Second_division/4_E Vacuolar Polyphosphate Polymerase-catalyzed Group Translocators.fasta","w")
    
    file5_A=open("../Data/Second_division/5_A Transmembrane 2-electron transfer carriers.fasta","w")
    file5_B=open("../Data/Second_division/5_B Transmembrane 1-electron transfer carriers.fasta","w")
    
    file8_A=open("../Data/Second_division/8_A Auxiliary transport proteins.fasta","w+")
    file8_B=open("../Data/Second_division/8_B Ribosomally synthesized protein-peptide toxins-agonists that target channels and carriers.fasta","w+")
    file8_C=open("../Data/Second_division/8_C Functionally characterized transporters lacking identified sequences.fasta","w+")
    
    file9_A=open("../Data/Second_division/9_A Recognized transporters of unknown biochemical mechanism.fasta","w+")
    file9_B=open("../Data/Second_division/9_B Putative transport proteins.fasta","w+")
    file9_C=open("../Data/Second_division/9_C Functionally characterized transporters lacking identified sequences.fasta","w+")
    
    for TC_ID in TC_IDs:
        if TC_ID[0:3]=="1.A":
            file1_A.write(str(records[i].format("fasta")))
            i+=1
        elif TC_ID[0:3]=="1.B":
            file1_B.write(str(records[i].format("fasta")))
        
            i+=1
        elif TC_ID[0:3]=="1.C":
            file1_C.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="1.D":
            file1_D.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="1.E":
            file1_E.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="1.F":
            file1_F.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="1.G":
            file1_G.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="1.H":
            file1_H.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="1.I":
            file1_I.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="1.J":
            file1_J.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="1.K":
            file1_K.write(str(records[i].format("fasta")))
            i+=1    
            
        elif TC_ID[0:3]=="1.L":
            file1_L.write(str(records[i].format("fasta")))
            i+=1    
            
        elif TC_ID[0:3]=="1.M":
            file1_M.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="2.A":
            file2_A.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="2.B":
            file2_B.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="2.C":
            file2_C.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="2.D":
            file2_D.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="3.A":
            file3_A.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="3.B":
            file3_B.write(str(records[i].format("fasta")))
            i+=1   
        
        elif TC_ID[0:3]=="3.C":
            file3_C.write(str(records[i].format("fasta")))
            i+=1  
            
        elif TC_ID[0:3]=="3.D":
            file3_D.write(str(records[i].format("fasta")))
            i+=1 

        elif TC_ID[0:3]=="3.E":
            file3_E.write(str(records[i].format("fasta")))
            i+=1   
              
        elif TC_ID[0:3]=="4.A":
            file4_A.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="4.B":
            file4_B.write(str(records[i].format("fasta")))
            i+=1   
        
        elif TC_ID[0:3]=="4.C":
            file4_C.write(str(records[i].format("fasta")))
            i+=1 
            
        elif TC_ID[0:3]=="4.D":
            file4_D.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="4.E":
            file4_E.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="5.A":
            file5_A.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="5.B":
            file5_B.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="8.A":
            file8_A.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="8.B":
            file8_B.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="8.C":
            file8_C.write(str(records[i].format("fasta")))
            i+=1
            
        elif TC_ID[0:3]=="9.A":
            file9_A.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="9.B":
            file9_B.write(str(records[i].format("fasta")))
            i+=1
        
        elif TC_ID[0:3]=="9.C":
            file9_C.write(str(records[i].format("fasta")))
            i+=1
        else:
            print("ERROR")
    
    file1_A.close()
    file1_B.close()
    file1_C.close()
    file1_D.close()
    file1_E.close()
    file1_G.close()
    file1_I.close()
    file1_K.close()
    file1_L.close()
    file1_M.close()
    
    file2_A.close()
    file2_B.close()
    file2_C.close()
    file2_D.close()
    
    file3_A.close()
    file3_B.close()
    file3_C.close()
    file3_D.close()
    file3_E.close()
    
    file4_A.close()
    file4_B.close()
    file4_C.close()
    file4_D.close()
    file4_E.close()
    
    file5_A.close()
    file5_B.close()
    
    file8_A.close()
    file8_B.close()
    file8_C.close()
    
    file9_A.close()
    file9_B.close()
    file9_C.close()
    
if __name__ == '__main__':
    
    handle = open("../Data/tcdb.fasta", "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    #print(len(records))
    #records contem 13788 proteinas transportadoras
    handle.close()
    
               
    print (records[0].id)  #first record
    print (records[0].seq)
    print (records[0].name)
    print (records[0].description)
    print (records[0].letter_annotations)
    print (records[0].annotations)
    print (records[0].features)
    print (records[0].dbxrefs)
    
    TC_IDs= TC_ID(records)
    print(TC_ID(records))
    
    
    
    #First_division(TC_IDs, records)
    #Second_division(TC_IDs, records)

