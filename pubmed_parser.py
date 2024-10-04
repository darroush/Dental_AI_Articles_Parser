from Bio import Medline, Entrez 
import pandas as pd
from tqdm import tqdm
from helper import *


# TODO:  set E-mail
Entrez.email = "mostafadesoki86@gmail.com"

# TODO:  set time interval
time_interval = '2024/09/01:2024/09/30'

# TODO:  modify research query if needed
search_query = f"""((("Artificial Intelligence"[Mesh] OR "Artificial Intelligence"[tw] OR "Machine Learning"[tw] OR "Convolutional neural network*"[tw] OR "Deep Learning"[tw] AND ({time_interval}[pdat])) AND 
                ("Dentistry"[Mesh] OR Dentistry[tw] OR Dental[tw] AND ({time_interval}[pdat])) NOT "systematic review"[Filter]) NOT "review"[Filter])"""

alldata = []

# Fetch the articles based on the query
articles = fetch_pubmed_data(search_query, max_records=100)

# Process the retrieved articles
for article in tqdm(articles):
    if article["JT"] not in journals_blacklist:
        try:
            pmid = article["PMID"]
        except: 
            pmid = ""
        try:
            author = article["FAU"]
        except:
            author = ""
        try:
            title = article["TI"]
        except:
            title = ""
        try:
            journal = article["JT"]
        except:
            journal = ""
        try:
            year = article["DP"]
        except:
            year = ""
        try:
            identifier = article["AID"][-1][:-6]
        except:
            identifier = ""
        try:
            abstract = article["AB"]
        except:
            abstract = ""

        dic = {
            "PMID" : pmid,
            "Author" : author,
            "Title" : title,
            "Journal" : journal,
            "Year" : year,
            "Abstract" : abstract,
            "DOI" : "https://doi.org/" + identifier
        }
        alldata.append(dic)

# save the data as excel file
df = pd.DataFrame(alldata)
df.to_excel(f'parser_output.xlsx', index=False)

print('finished')




    