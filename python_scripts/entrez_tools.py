from Bio import SeqIO
from Bio import Entrez
Entrez.email = "Kiranrambalig@gmail.com"
import time
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
from typing import Literal,Optional

def einfo(db = None,tries=3):
    for attempt in range(tries):
        try:
            handle = Entrez.einfo(db=db) if db else Entrez.einfo()
            record = Entrez.read(handle)
            handle.close()
            return record  
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `einfo function ` {attempt+1} failed, retrying...")
            time.sleep(2)
    raise Exception("NCBI failed after retries")

def esearch(db,term,usehistory=None,retmax = None,tries = 3):
    for attempt in range(tries):
        try:
            handle = Entrez.esearch(
                db=db,
                term=term,
                retmax = retmax,
                usehistory = usehistory
            )
            records = Entrez.read(handle)
            webenv = records.get("WebEnv")
            query_key = records.get("QueryKey")
            handle.close()
            return records, webenv, query_key
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `esearch function` {attempt+1} failed, retrying...")
            time.sleep(2)  #meaning -> stop everything for 2sec and then continue
    raise Exception("NCBI failed after retries")

def efetch(db,retmode:Literal["text","xml"]="text",
           rettype:Literal["fasta","gb","abstract"]="fasta",
           id=None,webenv=None,query_key=None,
           retmax=None,retstart=None,tries=3,):
    for attempt in range(tries):
        try:
            fetch_handle = Entrez.efetch(
                db=db,
                id = id,
                retmode = retmode,
                rettype = rettype,
                webenv = webenv,
                query_key = query_key,
                retmax = retmax,
                retstart = retstart
            )
            return fetch_handle
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `efetch function` {attempt+1} failed, retrying...")
            time.sleep(2)
    raise Exception("NCBI failed after retries")

def esummary(db,id,tries=3):
    for attempt in range(tries):
        try:
            handle = Entrez.esummary(
                db = db,
                id = id
            )
            records = Entrez.read(handle)
            handle.close()
            return records
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `esummary function` {attempt+1} failed, retrying...")
            time.sleep(2)
    raise Exception("NCBI failed after retries")

def elink(dbfrom,db,id,tries=3):
    for attempt in range(tries):
        try:
            handle = Entrez.elink(
                dbfrom = dbfrom,
                db = db,
                id = id
            )
            records = Entrez.read(handle)
            handle.close()
            return records
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `elink function ` {attempt+1} failed, retrying...")
            time.sleep(2)
    raise Exception("NCBI failed after retries")

def efetch_batch(db,webenv,query_key,total_count,batch_size,
        retmode:Literal["text","xml"]="text",
        rettype:Literal["fasta","gb","abstract"]="fasta",parse_as = "raw"):
    
    total_results = []
    for repeat in range(0,total_count,batch_size):
        print(f"Results batch from {repeat} to {repeat + batch_size} out of {total_count}... ")
        handle = efetch(
            db=db,
            id=None,
            webenv=webenv,
            query_key=query_key,
            retmode=retmode,
            rettype=rettype,
            retstart=repeat,
            retmax=batch_size)
        if parse_as == "seqio":
            records = list(SeqIO.parse(handle,rettype))
            total_results.extend(records)
        elif parse_as == "xml":
            data = Entrez.read(handle)
            total_results.append(data)
        else:
            raw_text = handle.read()
            total_results.append(raw_text)
        handle.close()
    return total_results

def epost(db,id_list,tries=3):  
    # isinstance checks the type of the variable and is safer 
    if isinstance(id_list,list):
        formatted_list = ",".join((map(str,id_list)))
        # map says take everything from the variable (here:id_list) 
        # and convert it into the said datatype(str)
        # .join() can only join str and if the id is int it will crash
    else:
        formatted_list = id_list
    for attempt in range(tries):
        try:
            
            handle = Entrez.epost(db=db,id=formatted_list)
            records = Entrez.read(handle)
            handle.close()
            webenv = records.get("WebEnv")
            query_key = records.get("QueryKey")
            return webenv, query_key
        except (RuntimeError,IncompleteRead,HTTPError,URLError):
            print(f"Attempt `epost function` {attempt+1} failed, retrying...")
            time.sleep(2)
    raise Exception("NCBI failed after retries")   
