#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#this is the import section
import Bio
from Bio import SeqIO
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import os
import re
import csv
import pandas as pd
from Bio import Entrez


# In[ ]:


'''
the below functions pulls a list of strings that contain 2 lines;
1.>fasta metadata as the description
2.the amino acid sequence of the protein in question
it requires a URL that will be shared for all these functions from the interpro API with the fasta additional information 
selected.

There is also a bit of code for creating a taxanomic look up table from the latest taxanomic codes and their scientific names
in the NCBI taxanomic database. The look up table gets saved as a pickle file. This was done to decrease the time needed to 
pull metadata on taxonomy for the InterPro entries by downloading them ahead of time and searching for information in the look 
up table. When the entries are not in the look up tabel instead the information is pulled by the Entrez system from NCBI through
biopython. You need to download the zip files from the NCBI site to run the first time: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
Then extract them to a folder to run the code for the lookup table.

The API pull code initially was partially based on the guides provided by the InterPro API webpage but was modified.

'''




# In[ ]:


### DONT RUN UNLESS YOUR REMAKING THE LOOKUP TABLE###
#how the taxonomy look up table was constructed:

taxfilepath = 'YOUR/PATH/HERE/new_taxdump'  # Update this to your actual path for taxdump

#col names in the nodes dmp file, for parsing the dmp file
nodes_cols = ['tax_id', 'parent_tax_id', 'rank', 'embl_code', 'division_id', 'inherited_div_flag',
              'genetic_code_id', 'inherited_GC_flag', 'mitochondrial_genetic_code_id', 'inherited_MGC_flag',
              'GenBank_hidden_flag', 'hidden_subtree_root_flag', 'comments','plastid genetic code id','inherited PGC flag',
              'specified_species',' hydrogenosome genetic code id',' inherited HGC flag ']
nodes = pd.read_csv(f'{taxfilepath}/nodes.dmp', sep='\t\|\t', engine='python', names=nodes_cols, usecols=['tax_id', 'parent_tax_id', 'rank'], dtype=str)

nodes['tax_id'] = nodes['tax_id'].str.strip()
nodes['parent_tax_id'] = nodes['parent_tax_id'].str.strip()
nodes['rank'] = nodes['rank'].str.strip()
print(nodes)


names_cols = ['tax_id', 'name_txt', 'unique_name', 'name_class']
names = pd.read_csv(f'{taxfilepath}/names.dmp', sep='\t\|\t', engine='python', names=names_cols, usecols=['tax_id', 'name_txt', 'name_class'], dtype=str)
names['tax_id'] = names['tax_id'].str.strip()
names['name_txt'] = names['name_txt'].str.strip()
names['name_class'] = names['name_class'].str.strip()
print(names)

scientific_names = names[names['name_class'] == 'scientific name\t|']


taxonomy = pd.merge(nodes, scientific_names, on='tax_id')

#this is a test to make sure the table works before making the pkl
print(taxonomy)
print(taxonomy[taxonomy['tax_id'] == '9606'])  # for human
print(taxonomy[taxonomy['tax_id'] == '9606']['rank'])

taxonomy.to_pickle('YOUR/PATH/HERE/new_taxdump/ncbi_2024_taxonomy_table.pkl')


# In[ ]:


# this function just pulls fastas and metadata saving things seperatly as individual fasta files

import json
import ssl
from urllib import request, error
from time import sleep
import sys

def output_fasta_list(base_url):
    context = ssl._create_unverified_context()
    next_url = base_url
    FASTA_STRINGS = [] 

    while next_url:
        HEADER_SEPARATOR = "|"
        try:
            req = request.Request(next_url, headers={"Accept": "application/json"})
            with request.urlopen(req, context=context) as res:
                if res.status == 408:  
                    sleep(61)
                    continue
                elif res.status == 204:  
                    break
                payload = json.loads(res.read().decode())
                next_url = payload.get("next") 

            for item in payload.get("results", []):
                entries = item.get("entry_subset") or item.get("entries", [])
                taxa = item.get("taxa", [{"lineage": []}])
                taxa_header = "|".join(taxa[0].get("lineage", []))
                #this part is modified to so the fasta files includes the metadata for down stream stuff I did. sperated by |
                if entries:
                    fasta_header = (
                        ">" + item["metadata"].get("accession", "unknown") + HEADER_SEPARATOR +
                        entries[0].get("accession", "unknown") + HEADER_SEPARATOR +
                        item["metadata"].get("name", "unknown") + HEADER_SEPARATOR +
                        item["metadata"]["source_organism"].get("taxId", "unknown") + HEADER_SEPARATOR +
                        item["metadata"]["source_organism"].get("scientificName", "unknown") + HEADER_SEPARATOR +
                        taxa_header
                    )
                else:
                    fasta_header = ">" + item["metadata"].get("accession", "unknown") + HEADER_SEPARATOR + item["metadata"].get("name", "unknown")

                fasta_seq = item.get("extra_fields", {}).get("sequence", "")
                fasta_string = fasta_header + "\n" + fasta_seq
                FASTA_STRINGS.append(fasta_string)

        except error.HTTPError as e:
            sys.stderr.write("Error processing URL: " + str(next_url) + "; Error: " + str(e) + "\n")
            if e.code == 408 and attempts < 3:
                attempts += 1
                sleep(61)
                continue
            else:
                break

        sleep(1)  # Rate limit delay

    return FASTA_STRINGS  # Return after fully processing all pages


# In[ ]:


#This function exports the individual fasta files to a given path
#if you are using windows you must changes the \ in the path to / or us the r''

def fasta_file_write(output_fasta_list,path):
    fasta_list = output_fasta_list
    pathvar = path
    
    for fasta in fasta_list:
        fasta_string = fasta
        fasta_string_2 = fasta
        fasta_string_name , fasta_string_seq = fasta_string.strip('>').split('\n')
        print(fasta_string_name,fasta_string_seq,'\n')
        filename_list = fasta_string_2.strip('>').split('|')
        filename = filename_list[0]
        rec = SeqRecord(
            Seq(fasta_string_seq),
            id='',
            name='',
            description = fasta_string_name)
        print('SeqRecord object : ',rec,sep='\n')
        records = [rec]
        SeqIO.write(records,'{}/{}.fasta'.format(pathvar,filename),'fasta')
        print('done with : ','{}/{}.fasta'.format(pathvar,filename))

    


# In[ ]:


#I load the NCBI pkl file after its made. Its not a bad idea to clear your variable at this point or restart your kernal for mem
tax_LOT = pd.read_pickle("YOUR/PATH/HERE/new_taxdump/ncbi_2024_taxonomy_table.pkl")


# In[ ]:


#this is the NCBI taxanomic ID handling for getting IDs and ranks of the IDs so all entries get lineage data
#these are called below
def fetch_TAXID_data(tax_id, email):
    classification = tax_id
    Entrez.email = email  # Set the email directly here
    retry_limit = 3  # Maximum number of retries for the ENTREZ PULL i set this to be 3 for speed

    attempts = 0
    while attempts < retry_limit:
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(classification), retmode="xml")
            records = Entrez.read(handle)
            for rec in records:
                rank = rec.get("Rank")
            break  
        except Exception as e:
            print(f"HTTP error: {str(e)} - Sleeping 10s and retrying for ID {classification}")
            sleep(10)
            attempts += 1  # Increment attempts
        finally:
            if 'handle' in locals() and handle:
                handle.close()  
    if attempts == retry_limit:
        print(f"Failed after {retry_limit} retries for ID {classification}")

    return rank



def fetch_LINEAGE_data(lineage_list, email):
    Entrez.email = email  
    lineage_df = {}
    
    for classification in lineage_list[2:]:  # Assuming the first two entries are to be skipped, this is root and cell orgs and skips to bacter (id:2)
        try:
            rank_series = tax_LOT[tax_LOT['tax_id'] == classification]['rank']
            if not rank_series.empty:
                rank = rank_series.iloc[0]  
                if rank:
                    lineage_df[rank] = classification
                    print(f'Rank: {rank}, Tax ID: {classification}')
                else:
                    print(f'No rank available for Tax ID: {classification}')
            else:
                print(f'No entry found for Tax ID: {classification}, fetching from NCBI...')
                rank = fetch_TAXID_data(classification, email)
                if rank and rank != "No rank available":
                    lineage_df[rank] = classification
                else:
                    print(f'Failed to fetch data for Tax ID: {classification} from NCBI.')
        except Exception as e:
            print(f'Error processing Tax ID {classification}: {str(e)}')

    return lineage_df

# testing everything still works:
list_tax = ['1', '0', '2624677', '1644055', '2157', '1236']
email = 'your_email@example.com'  # Replace with your actual email
print(fetch_LINEAGE_data(list_tax, email))


# In[ ]:


'''
this function creates a list of dictionaries. The dictionaries are Key value pairs:
"accession" : A0X009XJT10 , "IPR" : IPR009908 , etc,....
I have set it up so that it pulls what I consider to be all the useful information contained in the requested JSON
dictionary from the interpro pull
if you want to modify it open the url in the browser to see the actual JSON contents or print them with your own code
then odify the for loop that contains the metadata section (its around "metadata = {'Accession': accession,.....")
or look at InterPro API guide

'''
def fetch_data(BASE_URL,email):
    context = ssl._create_unverified_context()
    next_url = BASE_URL
    data_list = []  # List to hold each entry's data
    email = email

    while next_url:
        try:
            req = request.Request(next_url, headers={"Accept": "application/json"})
            res = request.urlopen(req, context=context)

            if res.status == 408:
                sleep(61)
                continue
            elif res.status == 204:
                break

            payload = json.loads(res.read().decode())
            next_url = payload.get("next", None)
            #check that the json isnt empty
            for i, item in enumerate(payload["results"]):
                entries = None
                if ("entry_subset" in item):
                    entries = item["entry_subset"]
                elif ("entries" in item):
                    entries = item["entries"]
                    
            #parsing the JSON from the results of the API pull
            for item in payload["results"]:
                taxa = item['taxa']
                accession = item["metadata"]["accession"]
                ipr = "{}".format(entries[0]["accession"]) 
                name = item["metadata"]["name"]
                gene = item["metadata"]["gene"]
                protein_length = "{}".format(entries[0]["protein_length"])
                in_alphafold = "{}".format(item["metadata"]["in_alphafold"])
                source_org_ID = item["metadata"]["source_organism"]["taxId"]
                source_org_name = item["metadata"]["source_organism"]["scientificName"]
                taxa_range = len(taxa[0]["lineage"])
                sequence = item["extra_fields"]["sequence"]
                
                #this is the metadata used to construct the full dataset used in the publication for every interpro IPR
                metadata = {'Accession': accession,
                                  'IPR': ipr, 
                                  'name': name,
                                  'gene': gene,
                                  'protein_length': protein_length,
                                  'in_alphafold': in_alphafold,
                                  'Source_Org_ID': source_org_ID,
                                  'Source_Org_Name': source_org_name,
                                  'taxa_range' : taxa_range,
                                  'Sequence': sequence,
                                 }
               
                lineage_df = fetch_LINEAGE_data(taxa[0]['lineage'],email)
                metadata.update(lineage_df)
                data_list.append(metadata)
                
            

            sleep(1)  

        except HTTPError as e:
            if e.code == 408:
                sleep(61)
                continue
            else:
                print(f"Error occurred: {e}")
                break

    return data_list


# In[ ]:


# this simple function just creates a dataframe from the fetch_data list 

def create_dataframe(data_list):
    df = pd.DataFrame(data_list)
    return df

# this function creates a csv file for exporting the dataframe, you dont need to use this if you dont to
#I just used this for personal ease of use later whe merging data from different API pulls.

def construct_csv(dataframe,pathvar,filename):
    path = pathvar
    name = filename
    working_df = dataframe
    header_df = list(working_df.keys())
    
    working_df.to_csv(path_or_buf="{}/{}.csv".format(path,name),sep=",",header = header_df)
    print("exported : " + name + " to : " + path)
    


# In[ ]:


#code used just for fastas this DOES NOT construct the data set
if __name__ == "__main__":
  
  pathvar = 'YOUR/PATH/HERE'
  #this url is provided by InterPro
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR012932/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"

# fasta stuff
  email = 'your@email.com'  
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)

  print("done VKOR")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR009908/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"

# fasta stuff
  email = 'your@email.com'
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)

  print("done MauE")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR003752/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"

# fasta stuff
  email = 'your@email.com'  
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)

  print("done DsbB")


# In[ ]:


#IPRs used: 
#    MauE: IPR009908
#    VKOR: IPR012932
#    DsbB: IPR003752
#URLS used for bacteria from interpro:
#    VKOR: "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR012932/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
#    MauE: "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR009908/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
#    DsbB: "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR003752/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


# In[ ]:


#code used for main data acquisition of the complete data set, it was done in multiple parts and takes a few days for all the data
#can be decreased by using thread pools.

if __name__ == "__main__":
  
  pathvar = 'YOUR/PATH/HERE'

  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR012932/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
  email = 'your@email.com'   
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)
#exporting dataset to csv

  archae_list = fetch_data(url,email)
  archaea_df = create_dataframe(archaea_list)
  name = "bacteria_vkor_set"
  construct_csv(archaea_df,pathvar,name)

  print("done VKOR")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR009908/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
  email = 'your@email.com' 
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)
#exporting dataset to csv

  archaea_list = fetch_data(url,email)
  archaea_df = create_dataframe(archaea_list)
  name = "bacteria_maue_set"
  construct_csv(archaea_df,pathvar,name)

  print("done MauE")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR003752/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
  email = 'your@email.com'  
  fasta_list = output_fasta_list(url)
  print(len(fasta_list))
  fasta_file_write(fasta_list,pathvar)
#exporting dataset to csv

  archaea_list = fetch_data(url,email)
  archaea_df = create_dataframe(archaea_list)
  name = "bacteria_dsbb_set"
  construct_csv(archaea_df,pathvar,name)

  print("done DsbB")


# In[ ]:


#this sesction is for the control set it includes IPR
# for ribsosomal_sub IPR047873: "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR047873/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
# for RNA polymeras_sub IPR010243:  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR010243/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"
# for signal peptidase IPR036286: "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR036286/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"

# signal peptidase has far more entries and takes the longest


# In[ ]:


#code used for main data acquisition
if __name__ == "__main__":
  
  pathvar = 'YOUR/PATH/HERE'

  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR047873/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


#exporting dataset to csv

  bacteria_list = fetch_data(url,email)
  bacteria_df = create_dataframe(bacteria_list)
  name = "bacteria_RIBOSOME_set"
  construct_csv(bacteria_df,pathvar,name)

  print("done Ribosomal subunit")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR010243/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


#exporting dataset to csv

  bacteria_list = fetch_data(url,email)
  bacteria_df = create_dataframe(bacteria_list)
  name = "bacteria_RIBOSOME_set"
  construct_csv(bacteria_df,pathvar,name)

  print("done RNA poly subunit")
    
########################
  pathvar = 'YOUR/PATH/HERE'
   
  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR036286/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


#exporting dataset to csv
  
  bacteria_list = fetch_data(url,email)
  bacteria_df = create_dataframe(bacteria_list)
  name = "bacteria_RIBOSOME_set"
  construct_csv(bacteria_df,pathvar,name)

  print("done singnal peptidase")


# In[ ]:


#This section is for MauG to determine its overlap with MauE
#https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR026259/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence

if __name__ == "__main__":
  
  pathvar = 'YOUR/PATH/HERE'

  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR026259/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


#exporting dataset to csv

  bacteria_list = fetch_data(url,email)
  bacteria_df = create_dataframe(bacteria_list)
  name = "bacteria_MAUG_set"
  construct_csv(bacteria_df,pathvar,name)

  print("done with MauG")


# In[ ]:


### this below section is for the MauA data and used IPR036560
### the species is tagged at the end and the entire thing is a for loop which cuts down the results from 1.2M to like less than 100K hopefully :D
#https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR036560/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"

if __name__ == "__main__":
  
  pathvar = 'YOUR/PATH/HERE'

  url =  "https://www.ebi.ac.uk:443/interpro/api/protein/UniProt/entry/InterPro/IPR036560/taxonomy/uniprot/2/?page_size=200&extra_fields=sequence"


#exporting dataset to csv

  bacteria_list = fetch_data(url,email)
  bacteria_df = create_dataframe(bacteria_list)
  name = "bacteria_MAUA_set"
  construct_csv(bacteria_df,pathvar,name)

  print("done with MauA")


# In[ ]:


#the data was then merged used as the foundation of the analysis.

