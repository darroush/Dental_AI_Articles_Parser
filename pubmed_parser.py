from Bio import Medline, Entrez 
import pandas as pd
from tqdm import tqdm
from helper import *


# set E-mail
Entrez.email = "mostafadesoki86@gmail.com"

# set time interval
time_interval = '2024/09/01:2024/09/30'

# set research query
search_query = f"""((("Artificial Intelligence"[Mesh] OR "Artificial Intelligence"[tw] OR "Machine Learning"[tw] OR "Convolutional neural network*"[tw] OR "Deep Learning"[tw] AND ({time_interval}[pdat])) AND 
                ("Dentistry"[Mesh] OR Dentistry[tw] OR Dental[tw] AND ({time_interval}[pdat])) NOT "systematic review"[Filter]) NOT "review"[Filter])"""

# exclude Q3 and Q4 journals
journals_blacklist = ["AMIA ... Annual Symposium proceedings. AMIA Symposium",
"Bioinformation","Compendium of continuing education in dentistry (Jamesburg, N.J. : 1995)",
"Cureus",
"Journal of the Indian Society of Pedodontics and Preventive Dentistry",
"Journal of veterinary dentistry",
"Journal of visualized experiments : JoVE",
"JPMA. The Journal of the Pakistan Medical Association",
"L' Orthodontie francaise",
"Medecine sciences : M/S",
"Special care in dentistry : official publication of the American Association of Hospital Dentists, the Academy of Dentistry for the Handicapped, and the American Society for Geriatric Dentistry",
"Studies in health technology and informatics",
"Technology and health care : official journal of the European Society for Engineering and Medicine",
"The journal of contemporary dental practice",
"The Journal of forensic odonto-stomatology",
"Zhonghua kou qiang yi xue za zhi = Zhonghua kouqiang yixue zazhi = Chinese journal of stomatology"]

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




    