# -*- coding: utf-8 -*-

'''
Created on 16/08/2016

@author: danie
'''


import re
import urllib
import requests
import time
import re


def run(seqs,evalue=0.0000000001):
    """
            Submit a job with the specified parameters.

            seqs         Protein sequence to analyse in fasta format.
            useblast     Use BLAST to get additional information (will increase running time)
            evalue       Select desired E-value (13 different options)
            queryfile    Submission of a local file in fasta format
    """

    url = 'http://services.cbu.uib.no/tools/bomp/handleForm'

    #opener = urllib.request.build_opener(urllib.request.HTTPHandler(debuglevel=0))

    args = {'seqs': seqs,
            'useblast':"on",
            'evalue':evalue,
            'queryfile':"",
            'SUBMIT':"Submit Search"}
    
    data = urllib.parse.urlencode(args)
    data=str.encode(data)

    #job_content = str(opener.open(url, data=data).read())
    
    session=requests.Session()
    resp=session.post(url,data=args)   
    job_content=resp.text
    #print(job_content)
    s= re.search("viewOutput\?id=.*?\"",job_content)
    job_id=s.group(0)
    job_id=job_id[:-1]

    return job_id

def result(job_id):
    
    url="http://services.cbu.uib.no/tools/bomp/"+job_id
    session=requests.Session()
    #===========================================================================
    # resp=session.get(url)
    # #print(resp.text)
    # s= re.search("Not finished \(please press reload to check again\)", resp.text)
    # print(s.group(0))
    # load=s.group(0)
    # while load == "Not finished (please press reload to check again)":
    #     print("Still loading")
    #     time.sleep(1)
    #     resp2=session.get(url)
    #     s=re.search("Not finished \(please press reload to check again\)", resp2.text)
    #     if s == None:
    #         load="Finished"
    #         print(load)
    #     else:
    #         load=s.group(0)
    #===========================================================================
    
    resp3=session.get(url)       
    r=re.search("barrel outer membrane proteins predicted is: .{5}",resp3.text)
    textres=r.group(0)
    print(textres)
    res=re.findall(r'\d+',textres)
    finalres=int(res[0])
    print(finalres)
    
    return finalres


if __name__ == "__main__":
    seqs=""">Sequence name
MQTKLLAIMLAAPVVFSSQEASASDFFGPEKISTEINLGTLSGKTKERVYEPEEGGRKVS
QLDWKYSNAAILKGAVNWELNPWLSVGAAGWTTLNSRGGNMVDQDWMDSGTPGTWTDESR
HPDTRLNYANEFDLNVKGWFLKESDYRLAIMAGYQESRYSFNATGGTYIYSENGGFRNET
GALPDKIKVIGYKQHFKIPYVGLTGNYRYDNFEFGGAFKYSGWVRGSDNDEHYVRQTTFR
SKVINQNYYSVAVNAGYYITPEAKVYIEGVWSRLTNKKGDTSLYDRSDNTSEHNNNGAGI
ENYNFITTAGLKYTF"""     

    job_id=run(seqs)
    print(job_id)
    #print(job_id)
    #result(job_id)
    
    #result("viewOutput?id=37365009")
    

    